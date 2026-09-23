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

#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TPDGCode.h>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iterator>
#include <memory>
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

using EffMapPtr = std::variant<std::shared_ptr<TH2>, std::shared_ptr<TH3>>;

struct BoundEfficiencyMap {
  using CoordsTuple = std::tuple<float, float, float>;

  const EffMapPtr* effMap;
  CoordsTuple coords;

  BoundEfficiencyMap(const EffMapPtr& effMap, float x, float y, float z) : effMap(&effMap), coords(x, y, z) {}
  BoundEfficiencyMap(const EffMapPtr& effMap, const CoordsTuple& coords) : effMap(&effMap), coords(coords) {}

  [[nodiscard]] std::pair<float, float> getBinEfficiencyAndError() const
  {
    return std::visit(
      [this](auto&& mapPtr) -> std::pair<float, float> {
        if (!mapPtr) {
          return {1.0f, 0.0f};
        }

        const auto& [x, y, z] = coords;

        // Extract the actual histogram type (TH2 or TH3) held by the smart pointer
        using HistoType = typename std::decay_t<decltype(mapPtr)>::element_type;

        // Compile-time branching: generates the exact correct function call
        if constexpr (std::is_same_v<HistoType, TH2>) {
          int bin = mapPtr->FindFixBin(y, z); // 2D case only
          return {mapPtr->GetBinContent(bin), mapPtr->GetBinError(bin)};
        } else {
          int bin = mapPtr->FindFixBin(x, y, z); // Full 3D case
          return {mapPtr->GetBinContent(bin), mapPtr->GetBinError(bin)};
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

        using HistoType = typename std::decay_t<decltype(mapPtr)>::element_type;

        if constexpr (std::is_same_v<HistoType, TH2>) {
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
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"}; // TO BE REMOVED

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables on phi mass selection
  struct : ConfigurableGroup {
    Configurable<std::pair<float, float>> rangeMPhiSignal{"rangeMPhiSignal", {1.0095f, 1.029f}, "Phi mass range for signal extraction"};
    Configurable<std::pair<float, float>> rangeMPhiSideband{"rangeMPhiSideband", {1.1f, 1.2f}, "Phi mass range for sideband extraction"};
  } phiConfigs;

  // Configurables for K0s selection
  struct : ConfigurableGroup {
    Configurable<bool> selectK0sInSigRegion{"selectK0sInSigRegion", true, "Select K0s candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMK0sSignal{"rangeMK0sSignal", {0.47f, 0.53f}, "K0S mass range for signal extraction"};
  } k0sConfigs;

  // Configurables for Lambda selection
  struct : ConfigurableGroup {
    Configurable<bool> selectLambdaInSigRegion{"selectLambdaInSigRegion", true, "Select Lambda candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMLambdaSignal{"rangeMLambdaSignal", {1.11f, 1.12f}, "Lambda mass range for signal extraction"};
  } lambdaConfigs;

  // Configurables for Xi selection
  struct : ConfigurableGroup {
    Configurable<bool> selectXiInSigRegion{"selectXiInSigRegion", true, "Select Xi candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMXiSignal{"rangeMXiSignal", {1.31f, 1.33f}, "Xi mass range for signal extraction"};
  } xiConfigs;

  // Configurables for Omega selection
  struct : ConfigurableGroup {
    Configurable<bool> selectOmegaInSigRegion{"selectOmegaInSigRegion", true, "Select Omega candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMOmegaSignal{"rangeMOmegaSignal", {1.66f, 1.68f}, "Omega mass range for signal extraction"};
  } omegaConfigs;

  // Configurables for Pions selection
  struct : ConfigurableGroup {
    Configurable<bool> selectPionInSigRegion{"selectPionInSigRegion", true, "Select Pion candidates in signal region"};
    Configurable<float> pidTPCMax{"pidTPCMax", 3.0f, "Maximum nSigma TPC"};
    Configurable<float> pidTOFMax{"pidTOFMax", 3.0f, "Maximum nSigma TOF"};
    // Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.5f, "Minimum pT after which TOF PID is applicable"};
  } pionConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};
  Configurable<std::vector<double>> binspTPhiExt{"binspTPhiExt", {0.0, 0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi extended for MC Gen"};

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};
  Configurable<std::vector<double>> binspTK0SExt{"binspTK0SExt", {0.0, 0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S extended for MC Gen"};

  // Configurable on Lambda pT bins
  Configurable<std::vector<double>> binspTLambda{"binspTLambda", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Lambda"};
  Configurable<std::vector<double>> binspTLambdaExt{"binspTLambdaExt", {0.0, 0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Lambda extended for MC Gen"};

  // Configurable on Xi pT bins
  Configurable<std::vector<double>> binspTXi{"binspTXi", {0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Xi"};
  Configurable<std::vector<double>> binspTXiExt{"binspTXiExt", {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Xi extended for MC Gen"};

  // Configurable on Omega pT bins
  Configurable<std::vector<double>> binspTOmega{"binspTOmega", {0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Omega"};
  Configurable<std::vector<double>> binspTOmegaExt{"binspTOmegaExt", {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Omega extended for MC Gen"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for Pions"};
  Configurable<std::vector<double>> binspTPiExt{"binspTPiExt", {0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for Pions extended for MC Gen"};

  Configurable<std::array<bool, kAssocPartSize>> activeCorrelationTypes{"activeCorrelationTypes", {true, false, false, false, false, true}, "Enable correlation types: 0 - K0S, 1 - Lambda, 2 - AntiLambda, 3 - Xi, 4 - Omega, 5 - Pion"};

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
    Configurable<bool> propagateEffError{"propagateEffError", false, "Propagate efficiency error"};
    Configurable<bool> perPeriodEfficiency{"perPeriodEfficiency", false, "Use per period efficiency maps"};
  } efficiencyConfigs;

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurables for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbEfficiencyPath{"ccdbEfficiencyPath", "Users/s/scannito/Efficiencies/h3EffMap", "Correction path to file"};

  // Configurables for minimum pt selection in McGen
  struct : ConfigurableGroup {
    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<float> v0SettingMinPt{"v0SettingMinPt", 0.1f, "V0 min pt"};
    Configurable<float> lambdaSettingMinPt{"lambdaSettingMinPt", 0.2f, "Lambda min pt"};
    Configurable<float> cascadeSettingMinPt{"cascadeSettingMinPt", 0.8f, "Cascade min pt"};
    Configurable<float> omegaSettingMinPt{"omegaSettingMinPt", 0.8f, "Omega min pt"};
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};
    Configurable<bool> bypassPtCut{"bypassPtCut", false, "Bypass the minimum pt cut at MCGen level"};
  } minPtMcGenConfigs;

  // Filter on phi selected collisions
  Filter collisionFilter = (eventSelectionType == 0 && aod::lf_selection_event::defaultSel == true) ||
                           (eventSelectionType == 1 && aod::lf_selection_event::defaultSel == true && aod::lf_selection_event::phimesonSel == true);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike>>;
  using SimCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike, aod::McCollisionLabels>>;
  using MCCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::PhiStrangeEvtSelMcGen>>;

  // Slice cache and Preslices for table slicing
  SliceCache cache;

  struct : PresliceGroup {
    Preslice<aod::PhimesonCandidatesData> phiCandDataPerCollision = aod::lf_selection_phi_candidate::collisionId;
    // PresliceUnsorted<aod::PhimesonCandidatesMcReco> phiCandMcRecoPerCollision = aod::lf_selection_phi_candidate::collisionId;
    Preslice<aod::PhimesonCandidatesMcReco> phiCandMcRecoPerCollision = aod::lf_selection_phi_candidate::collisionId;

    Preslice<aod::K0sReducedCandidatesData> k0sDataPerCollision = aod::v0::collisionId;
    Preslice<aod::K0sReducedCandidatesMcReco> k0sMcRecoPerCollision = aod::v0::collisionId;

    Preslice<aod::LambdaReducedCandidatesData> lambdaDataPerCollision = aod::v0::collisionId;
    Preslice<aod::LambdaReducedCandidatesMcReco> lambdaMcRecoPerCollision = aod::v0::collisionId;

    Preslice<aod::AntiLambdaReducedCandidatesData> antiLambdaDataPerCollision = aod::v0::collisionId;
    Preslice<aod::AntiLambdaReducedCandidatesMcReco> antiLambdaMcRecoPerCollision = aod::v0::collisionId;

    Preslice<aod::XiReducedCandidatesData> xiDataPerCollision = aod::cascade::collisionId;
    Preslice<aod::XiReducedCandidatesMcReco> xiMcRecoPerCollision = aod::cascade::collisionId;

    Preslice<aod::OmegaReducedCandidatesData> omegaDataPerCollision = aod::cascade::collisionId;
    Preslice<aod::OmegaReducedCandidatesMcReco> omegaMcRecoPerCollision = aod::cascade::collisionId;

    Preslice<aod::PionTracksData> pionTrackDataPerCollision = aod::track::collisionId;
    Preslice<aod::PionTracksMcReco> pionTrackMcRecoPerCollision = aod::track::collisionId;

    Preslice<aod::McParticles> mcPartPerMcCollision = aod::mcparticle::mcCollisionId;
  } preslices;

  // Necessary service to retrieve efficiency maps from CCDB
  Service<ccdb::BasicCCDBManager> ccdb{};

  // std::shared_ptr<TH3> effMapPhi{};
  // std::array<std::shared_ptr<TH3>, kAssocPartSize> effMapsAssoc{};

  EffMapPtr effMapPhi;
  std::array<EffMapPtr, kAssocPartSize> effMapsAssoc{};

  // Binning policy and axes for mixed event
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {10, -10.0f, 10.0f}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 15.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Multiplicity percentage binning for mixing"};

  using BinningTypeVertexCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeVertexCent binningOnVertexAndCent{{axisVertexMixing, axisCentralityMixing}, true};

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
    std::vector<MiniParticle> k0sParticles;
    std::vector<MiniParticle> lambdaParticles;
    std::vector<MiniParticle> antilambdaParticles;
    std::vector<MiniParticle> xiParticles;
    std::vector<MiniParticle> omegaParticles;
    std::vector<MiniParticle> pionParticles;
  };

  // Buffer for mixed event, organized as a vector of deques, one for each multiplicity bin, containing the past events with their particles of interest needed for mixing
  std::vector<std::deque<MiniEvent>> eventBuffer;

  int mLastRunNumber{-1};

  void init(InitContext&)
  {
    AxisSpec vertexZAxis = {100, -cutZVertex, cutZVertex, "vrtx_{Z} [cm]"}; // TO BE REMOVED
    AxisSpec yAxis = {yConfigs.nBinsY, -yConfigs.cfgYAcceptance, yConfigs.cfgYAcceptance, "#it{y}"};
    AxisSpec deltayAxis = {yConfigs.nBinsDeltaY, -1.0f, 1.0f, "#Delta#it{y}"};
    AxisSpec deltaphiAxis = {72, -o2::constants::math::PIHalf, o2::constants::math::PIHalf * 3, "#Delta#varphi"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTPhiAxis = {120, 0.0f, 12.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPhiAxis{(std::vector<double>)binspTPhi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPhiAxisExt{(std::vector<double>)binspTPhiExt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massK0SAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTK0SAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxisExt{(std::vector<double>)binspTK0SExt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massLambdaAxis = {200, 1.1f, 1.13f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTLambdaAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTLambdaAxis{(std::vector<double>)binspTLambda, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTLambdaAxisExt{(std::vector<double>)binspTLambdaExt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massXiAxis = {200, 1.3f, 1.4f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTXiAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTXiAxis{(std::vector<double>)binspTXi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTXiAxisExt{(std::vector<double>)binspTXiExt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massOmegaAxis = {200, 1.65f, 1.7f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTOmegaAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTOmegaAxis{(std::vector<double>)binspTOmega, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTOmegaAxisExt{(std::vector<double>)binspTOmegaExt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec nSigmaPiAxis = {100, -10.0f, 10.0f, "N#sigma #pi"};
    AxisSpec pTPiAxis = {50, 0.0f, 5.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxisExt{(std::vector<double>)binspTPiExt, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("phi/h3PhiData", "Invariant mass of Phi in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

    histos.add("phiK0S/h6PhiK0SData", "Invariant mass of Phi vs Invariant mass of K0Short in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, massPhiAxis, massK0SAxis});
    histos.add("phiLambda/h6PhiLambdaData", "Invariant mass of Phi vs Invariant mass of Lambda in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTLambdaAxis, deltayAxis, massPhiAxis, massLambdaAxis});
    histos.add("phiXi/h6PhiXiData", "Invariant mass of Phi vs Invariant mass of Xi in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTXiAxis, deltayAxis, massPhiAxis, massXiAxis});
    histos.add("phiOmega/h6PhiOmegaData", "Invariant mass of Phi vs Invariant mass of Omega in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTOmegaAxis, deltayAxis, massPhiAxis, massOmegaAxis});
    histos.add("phiPi/h6PhiPiTPCData", "Invariant mass of Phi vs nSigmaTPC Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, massPhiAxis, nSigmaPiAxis});
    histos.add("phiPi/h6PhiPiTOFData", "Invariant mass of Phi vs nSigmaTOF Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, massPhiAxis, nSigmaPiAxis});

    histos.add("phiK0S/h6PhiK0SDataME", "Invariant mass of Phi vs Invariant mass of K0Short in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, massPhiAxis, massK0SAxis});
    histos.add("phiLambda/h6PhiLambdaDataME", "Invariant mass of Phi vs Invariant mass of Lambda in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTLambdaAxis, deltayAxis, massPhiAxis, massLambdaAxis});
    histos.add("phiXi/h6PhiXiDataME", "Invariant mass of Phi vs Invariant mass of Xi in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTXiAxis, deltayAxis, massPhiAxis, massXiAxis});
    histos.add("phiOmega/h6PhiOmegaDataME", "Invariant mass of Phi vs Invariant mass of Omega in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTOmegaAxis, deltayAxis, massPhiAxis, massOmegaAxis});
    histos.add("phiPi/h6PhiPiTPCDataME", "Invariant mass of Phi vs nSigmaTPC Pion in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, massPhiAxis, nSigmaPiAxis});
    histos.add("phiPi/h6PhiPiTOFDataME", "Invariant mass of Phi vs nSigmaTOF Pion in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, massPhiAxis, nSigmaPiAxis});

    for (const auto& label : PhiMassRegionLabels) {
      histos.add(fmt::format("phiK0S/h5PhiK0SData{}", label).c_str(), "Deltay vs deltaphi for Phi and K0Short in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiLambda/h5PhiLambdaData{}", label).c_str(), "Deltay vs deltaphi for Phi and Lambda in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTLambdaAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiXi/h5PhiXiData{}", label).c_str(), "Deltay vs deltaphi for Phi and Xi in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTXiAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiOmega/h5PhiOmegaData{}", label).c_str(), "Deltay vs deltaphi for Phi and Omega in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTOmegaAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiPi/h5PhiPiData{}", label).c_str(), "Deltay vs deltaphi for Phi and Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, deltaphiAxis});

      histos.add(fmt::format("phiK0S/h5PhiK0SDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and K0Short in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiLambda/h5PhiLambdaDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and Lambda in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTLambdaAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiXi/h5PhiXiDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and Xi in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTXiAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiOmega/h5PhiOmegaDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and Omega in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTOmegaAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiPi/h5PhiPiDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and Pion in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, deltaphiAxis});
    }

    histos.add("event/hRecoMCMultiplicityPercent", "RecoMC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/h2RecoMCVertexZvsMult", "RecoMC Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    histos.add("event/hSplitVertexZ", "Split in z-vtx", kTH1F, {{100, -5.0f, 5.0f}});
    histos.add("event/hGenMCMultiplicityPercent", "Generated MC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/hGenMCAssocRecoMultiplicityPercent", "Generated MC associated Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/h2GenMCAssocRecoVertexZvsMult", "Generated MC associated reco Vertex Z vs multiplicity", kTH2F, {vertexZAxis, binnedmultAxis});

    histos.add("phi/h4PhiMCReco", "Phi in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPhiAxis, yAxis});
    histos.add("phi/h3PhiMCGen", "Phi in MC Gen", kTH3F, {binnedmultAxis, binnedpTPhiAxis, yAxis});
    histos.add("phi/h4PhiMCGenAssocReco", "Phi in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPhiAxis, yAxis});

    histos.add("k0s/h4K0SMCReco", "K0S in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTK0SAxis, yAxis});
    histos.add("k0s/h3K0SMCGen", "K0S in MC Gen", kTH3F, {binnedmultAxis, binnedpTK0SAxis, yAxis});
    histos.add("k0s/h4K0SMCGenAssocReco", "K0S in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTK0SAxis, yAxis});

    histos.add("lambda/h4LambdaMCReco", "Lambda in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTLambdaAxis, yAxis});
    histos.add("lambda/h3LambdaMCGen", "Lambda in MC Gen", kTH3F, {binnedmultAxis, binnedpTLambdaAxis, yAxis});
    histos.add("lambda/h4LambdaMCGenAssocReco", "Lambda in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTLambdaAxis, yAxis});

    histos.add("antilambda/h4AntiLambdaMCReco", "AntiLambda in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTLambdaAxis, yAxis});
    histos.add("antilambda/h3AntiLambdaMCGen", "AntiLambda in MC Gen", kTH3F, {binnedmultAxis, binnedpTLambdaAxis, yAxis});
    histos.add("antilambda/h4AntiLambdaMCGenAssocReco", "AntiLambda in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTLambdaAxis, yAxis});

    histos.add("xi/h4XiMCReco", "Xi in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTXiAxis, yAxis});
    histos.add("xi/h3XiMCGen", "Xi in MC Gen", kTH3F, {binnedmultAxis, binnedpTXiAxis, yAxis});
    histos.add("xi/h4XiMCGenAssocReco", "Xi in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTXiAxis, yAxis});

    histos.add("omega/h4OmegaMCReco", "Omega in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTOmegaAxis, yAxis});
    histos.add("omega/h3OmegaMCGen", "Omega in MC Gen", kTH3F, {binnedmultAxis, binnedpTOmegaAxis, yAxis});
    histos.add("omega/h4OmegaMCGenAssocReco", "Omega in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTOmegaAxis, yAxis});

    histos.add("pi/h4PiMCReco", "Pion in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPiAxis, yAxis});
    histos.add("pi/h3PiMCGen", "Pion in MC Gen", kTH3F, {binnedmultAxis, binnedpTPiAxis, yAxis});
    histos.add("pi/h4PiMCGenAssocReco", "Pion in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPiAxis, yAxis});

    histos.add("phi/h3PhiMCClosureGen", "Phi in MC Gen for MC Closure Test", kTH3F, {binnedmultAxis, binnedpTPhiAxisExt, yAxis});

    histos.add("phiK0S/h5PhiK0SClosureMCGen", "Deltay vs deltaphi for Phi and K0Short in MCGen", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTK0SAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiLambda/h5PhiLambdaClosureMCGen", "Deltay vs deltaphi for Phi and Lambda in MCGen", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTLambdaAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiXi/h5PhiXiClosureMCGen", "Deltay vs deltaphi for Phi and Xi in MCGen", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTXiAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiOmega/h5PhiOmegaClosureMCGen", "Deltay vs deltaphi for Phi and Omega in MCGen", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTOmegaAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiPi/h5PhiPiClosureMCGen", "Deltay vs deltaphi for Phi and Pion in MCGen", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTPiAxisExt, deltayAxis, deltaphiAxis});

    histos.add("phiK0S/h5PhiK0SClosureMCGenME", "Deltay vs deltaphi for Phi and K0Short in MCGen ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTK0SAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiLambda/h5PhiLambdaClosureMCGenME", "Deltay vs deltaphi for Phi and Lambda in MCGen ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTLambdaAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiXi/h5PhiXiClosureMCGenME", "Deltay vs deltaphi for Phi and Xi in MCGen ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTXiAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiOmega/h5PhiOmegaClosureMCGenME", "Deltay vs deltaphi for Phi and Omega in MCGen ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTOmegaAxisExt, deltayAxis, deltaphiAxis});
    histos.add("phiPi/h5PhiPiClosureMCGenME", "Deltay vs deltaphi for Phi and Pion in MCGen ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTPiAxisExt, deltayAxis, deltaphiAxis});

    // Load efficiency maps from CCDB
    if (efficiencyConfigs.applyEfficiency) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      /*for (int i = 0; i < ParticleOfInterestSize; ++i) {
        loadEfficiencyMapFromCCDB(static_cast<ParticleOfInterest>(i));
      }*/
      if (!efficiencyConfigs.perPeriodEfficiency) {
        loadEfficiencyMaps();
      }
    }

    eventBuffer.resize(binsMult->size() - 1);
  }

  void fetchSingleEfficiencyMapFromCCDB(EffMapPtr& effMap, std::string_view particleName, std::optional<int> runNumber = std::nullopt)
  {
    std::string path = fmt::format("{}{}", ccdbEfficiencyPath.value, particleName);

    auto fetchFromCCDB = [&]<typename HistoType>() -> std::shared_ptr<HistoType> {
      auto mapPtr = runNumber ? ccdb->getForRun<HistoType>(path, *runNumber) : ccdb->get<HistoType>(path);
      return std::shared_ptr<HistoType>(mapPtr);
    };

    if (auto map3D = fetchFromCCDB.operator()<TH3D>()) {
      effMap = map3D;
      LOG(debug) << "Efficiency map (TH3) for " << particleName << " loaded from CCDB";
      return;
    }

    if (auto map2D = fetchFromCCDB.operator()<TH2D>()) {
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
      if (activeCorrelationTypes->at(i)) {
        fetchSingleEfficiencyMapFromCCDB(effMapsAssoc[i], AssocParticleLabels[i], runNumber);
      }
    }

    if (runNumber) {
      LOG(info) << "Efficiency maps successfully updated from CCDB for run " << *runNumber;
    } else {
      LOG(info) << "Global efficiency maps successfully loaded from CCDB";
    }
  }

  // Compute weight based on efficiencies
  template <typename... BoundEffMaps>
  std::pair<float, float> computeWeightAndError(const BoundEffMaps&... boundEffMaps)
  {
    if (!efficiencyConfigs.applyEfficiency) {
      return {1.0f, 0.0f};
    }

    float effTot = 1.0f;
    float relErrSqSum = 0.0f;

    auto processMap = [&](const auto& boundMap) {
      auto [eff, err] = boundMap.getBinEfficiencyAndError(); // Unpack efficiency and error from the map

      if (efficiencyConfigs.useEffInterpolation) {
        eff = boundMap.interpolateEfficiency();
        // For simplicity, we keep the error from the bin content even when interpolating, but this can be refined if needed
      }

      effTot *= eff;

      if (eff > 0.0f) {
        float mapErr = efficiencyConfigs.propagateEffError ? err : 0.0f; // Optionally propagate error, otherwise treat as zero
        relErrSqSum += (mapErr / eff) * (mapErr / eff);
      }
    };

    (processMap(boundEffMaps), ...); // Fold expression to process all bound efficiency maps

    if (effTot <= 0.0f) {
      return {1.0f, 0.0f};
    }

    float weight = 1.0f / effTot;
    float weightErr = weight * std::sqrt(relErrSqSum); // Propagate relative error to the weight

    return {weight, weightErr};
  }

  template <AssociatedParticleType PartType>
  std::pair<float, float> computeAssocWeight(float multiplicity, const auto& phiCand, const auto& assocCand)
  {
    return computeWeightAndError(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()),
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
      return minPtMcGenConfigs.v0SettingMinPt.value;
    } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
      return minPtMcGenConfigs.lambdaSettingMinPt.value;
    } else if constexpr (PartType == kXi) {
      return minPtMcGenConfigs.cascadeSettingMinPt.value;
    } else if constexpr (PartType == kOmega) {
      return minPtMcGenConfigs.omegaSettingMinPt.value;
    } else if constexpr (PartType == kPion) {
      return minPtMcGenConfigs.cMinPionPtcut.value;
    } else {
      static_assert(PartType == kK0S || PartType == kLambda || PartType == kAntiLambda ||
                      PartType == kXi || PartType == kOmega || PartType == kPion,
                    "Unsupported particle type in getMinGenPt");
      return 0.f;
    }
  }

  template <typename THist, typename... Args>
  void customFillHist(auto histId, const std::pair<float, float>& weightPair, Args... coords)
  {
    auto hist = histos.get<THist>(histId);
    if (!hist) {
      return;
    }

    // Extract the weight and its propagated uncertainty
    auto [w, wErr] = weightPair;

    // Find the global bin number for the given physical coordinates
    int bin = hist->FindFixBin(coords...);

    // Retrieve the previous bin content and its absolute error
    double prevContent = hist->GetBinContent(bin);
    double prevErr = hist->GetBinError(bin);

    // Calculate the new content by adding the current weight
    double newContent = prevContent + w;

    // Error propagation in quadrature:
    // prevErr^2 : previous accumulated variance
    // (w * w)   : statistical fluctuation added by this specific particle count
    // (wErr * wErr) : systematic uncertainty added by the efficiency map inaccuracy
    double newErr = std::sqrt(prevErr * prevErr + (w * w) + (wErr * wErr));

    // Update the histogram bin
    hist->SetBinContent(bin, newContent);
    hist->SetBinError(bin, newErr);
  }

  template <typename... Args>
  void customFillTHn(auto histId, const std::pair<float, float>& weightPair, Args... coords)
  {
    auto hist = histos.get<THnSparse>(histId);
    if (!hist) {
      return;
    }

    // Extract the weight and its propagated uncertainty
    auto [w, wErr] = weightPair;

    // Create an array of floats for the coordinates to pass to GetBin and SetBinContent
    std::array<double, sizeof...(Args)> coordArray{static_cast<double>(coords)...};

    // Find the bin number.
    // The 'true' flag is mandatory: it allocates the bin in memory if it doesn't exist yet.
    auto bin = hist->GetBin(coordArray.data(), true);

    // Retrieve the previous content and the squared error (variance)
    double prevContent = hist->GetBinContent(bin);
    double prevErr2 = hist->GetBinError2(bin);

    // Calculate the new content
    double newContent = prevContent + w;

    // Add the new variances to the accumulated variance
    // (w * w) is the statistical term, (wErr * wErr) is the map uncertainty term
    double newErr2 = prevErr2 + (w * w) + (wErr * wErr);

    // Update the THnSparse bin
    hist->SetBinContent(bin, newContent);
    hist->SetBinError2(bin, newErr2); // Note: SetBinError2 takes the squared error directly
  }

  template <AssociatedParticleType PartType, typename TAssoc>
  bool isAssocValid(TAssoc const& assoc)
  {
    if constexpr (PartType == kK0S) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && k0sConfigs.selectK0sInSigRegion;
      const auto& [minMass, maxMass] = k0sConfigs.rangeMK0sSignal.value;
      return (!efficiencyConfigs.applyEfficiency || assoc.pt() < binspTK0S->back()) && (!applyMassCut || assoc.inMassRegion(minMass, maxMass));
    } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && lambdaConfigs.selectLambdaInSigRegion;
      const auto& [minMass, maxMass] = lambdaConfigs.rangeMLambdaSignal.value;
      return (!efficiencyConfigs.applyEfficiency || assoc.pt() < binspTLambda->back()) && (!applyMassCut || assoc.inMassRegion(minMass, maxMass));
    } else if constexpr (PartType == kXi) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && xiConfigs.selectXiInSigRegion;
      const auto& [minMass, maxMass] = xiConfigs.rangeMXiSignal.value;
      return (!efficiencyConfigs.applyEfficiency || assoc.pt() < binspTXi->back()) && (!applyMassCut || assoc.inMassRegion(minMass, maxMass));
    } else if constexpr (PartType == kOmega) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && omegaConfigs.selectOmegaInSigRegion;
      const auto& [minMass, maxMass] = omegaConfigs.rangeMOmegaSignal.value;
      return (!efficiencyConfigs.applyEfficiency || assoc.pt() < binspTOmega->back()) && (!applyMassCut || assoc.inMassRegion(minMass, maxMass));
    } else if constexpr (PartType == kPion) {
      const bool applyNSigmaCut = (analysisMode == kDeltaYvsDeltaPhi) && pionConfigs.selectPionInSigRegion;
      return (!efficiencyConfigs.applyEfficiency || assoc.pt() < binspTPi->back()) &&
             (!applyNSigmaCut || assoc.inNSigmaRegion(pionConfigs.pidTPCMax, pionConfigs.pidTOFMax));
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
  void fillPhiAssocCorrelation(TPhiCand const& phiCand, TAssoc const& assoc, float multiplicity, const std::pair<float, float>& weight)
  {
    const std::array<std::pair<float, float>, kPhiMassRegions> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    if (analysisMode == kMassvsMass) {
      if constexpr (PartType == kK0S) {
        if constexpr (IsME) {
          customFillTHn(HIST("phiK0S/h6PhiK0SDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        } else {
          customFillTHn(HIST("phiK0S/h6PhiK0SData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        }
      } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
        if constexpr (IsME) {
          customFillTHn(HIST("phiLambda/h6PhiLambdaDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        } else {
          customFillTHn(HIST("phiLambda/h6PhiLambdaData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        }
      } else if constexpr (PartType == kXi) {
        if constexpr (IsME) {
          customFillTHn(HIST("phiXi/h6PhiXiDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        } else {
          customFillTHn(HIST("phiXi/h6PhiXiData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        }
      } else if constexpr (PartType == kOmega) {
        if constexpr (IsME) {
          customFillTHn(HIST("phiOmega/h6PhiOmegaDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        } else {
          customFillTHn(HIST("phiOmega/h6PhiOmegaData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.m());
        }
      } else if constexpr (PartType == kPion) {
        if constexpr (IsME) {
          customFillTHn(HIST("phiPi/h6PhiPiTPCDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.nSigmaTPC());
          customFillTHn(HIST("phiPi/h6PhiPiTOFDataME"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.nSigmaTOF());
        } else {
          customFillTHn(HIST("phiPi/h6PhiPiTPCData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.nSigmaTPC());
          customFillTHn(HIST("phiPi/h6PhiPiTOFData"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), phiCand.m(), assoc.nSigmaTOF());
        }
      }
    } else if (analysisMode == kDeltaYvsDeltaPhi) {
      std::function<void()> fillSignal;
      std::function<void()> fillSideband;

      if constexpr (PartType == kK0S) {
        if constexpr (IsME) {
          fillSignal = [&]() { customFillTHn(HIST("phiK0S/h5PhiK0SDataMESignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiK0S/h5PhiK0SDataMESideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        } else {
          fillSignal = [&]() { customFillTHn(HIST("phiK0S/h5PhiK0SDataSignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiK0S/h5PhiK0SDataSideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        }
      } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
        if constexpr (IsME) {
          fillSignal = [&]() { customFillTHn(HIST("phiLambda/h5PhiLambdaDataMESignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiLambda/h5PhiLambdaDataMESideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        } else {
          fillSignal = [&]() { customFillTHn(HIST("phiLambda/h5PhiLambdaDataSignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiLambda/h5PhiLambdaDataSideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        }
      } else if constexpr (PartType == kXi) {
        if constexpr (IsME) {
          fillSignal = [&]() { customFillTHn(HIST("phiXi/h5PhiXiDataMESignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiXi/h5PhiXiDataMESideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        } else {
          fillSignal = [&]() { customFillTHn(HIST("phiXi/h5PhiXiDataSignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiXi/h5PhiXiDataSideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        }
      } else if constexpr (PartType == kOmega) {
        if constexpr (IsME) {
          fillSignal = [&]() { customFillTHn(HIST("phiOmega/h5PhiOmegaDataMESignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiOmega/h5PhiOmegaDataMESideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        } else {
          fillSignal = [&]() { customFillTHn(HIST("phiOmega/h5PhiOmegaDataSignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiOmega/h5PhiOmegaDataSideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        }
      } else if constexpr (PartType == kPion) {
        if constexpr (IsME) {
          fillSignal = [&]() { customFillTHn(HIST("phiPi/h5PhiPiDataMESignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiPi/h5PhiPiDataMESideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        } else {
          fillSignal = [&]() { customFillTHn(HIST("phiPi/h5PhiPiDataSignal"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
          fillSideband = [&]() { customFillTHn(HIST("phiPi/h5PhiPiDataSideband"), weight, multiplicity, phiCand.pt(), assoc.pt(), phiCand.y() - assoc.y(), getDeltaPhi(phiCand.phi(), assoc.phi())); };
        }
      }

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
    if (!isAssocValid<PartType>(assoc)) {
      return;
    }

    auto weight = computeAssocWeight<PartType>(multiplicity, phiCand, assoc);
    fillPhiAssocCorrelation<PartType, IsME>(phiCand, assoc, multiplicity, weight);
  }

  template <AssociatedParticleType PartType, typename TPhiCand, typename TAssocCands>
  void processOneAssocSpeciesSE(TPhiCand const& phiCand, TAssocCands const& assocCands, float multiplicity)
  {
    if (!activeCorrelationTypes->at(PartType)) {
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
      if (efficiencyConfigs.applyEfficiency && efficiencyConfigs.applyPhiEfficiency && phiCand.pt() >= binspTPhi->back()) {
        continue;
      }

      auto weightPhi = computeWeightAndError(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()));
      customFillHist<TH3>(HIST("phi/h3PhiData"), weightPhi, multiplicity, phiCand.pt(), phiCand.m());

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
    auto tuplePhiAssoc = std::make_tuple(phiCandidates, assocCands);
    Pair<TCollisions, TPhiCands, TAssocCands, BinningTypeVertexCent> pairPhiAssoc{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiAssoc, &cache};

    for (const auto& [c1, phiCands, c2, assocRed] : pairPhiAssoc) {
      float multiplicity = c1.centFT0M();

      if (efficiencyConfigs.applyEfficiency && efficiencyConfigs.perPeriodEfficiency) {
        auto bc = c1.template bc_as<aod::BCs>();
        loadEfficiencyMaps(bc.runNumber());
      }

      for (const auto& [phiCand, assoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, assocRed))) {
        if (efficiencyConfigs.applyEfficiency && efficiencyConfigs.applyPhiEfficiency && phiCand.pt() >= binspTPhi->back()) {
          continue;
        }

        processOneAssocPair<PartType, true>(phiCand, assoc, multiplicity);
      }
    }
  }

  void processPhiK0SMEDataLike(SelCollisions const& collisions,
                               aod::PhimesonCandidatesData const& phiCandidates,
                               aod::K0sReducedCandidatesData const& k0sReduced,
                               aod::BCs const&)
  {
    processPhiAssocME<kK0S>(collisions, phiCandidates, k0sReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiK0SMEDataLike, "Process function for Phi-K0S 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiK0SMEMCWithPDG(SimCollisions const& collisions,
                                aod::PhimesonCandidatesMcReco const& phiCandidates,
                                aod::K0sReducedCandidatesMcReco const& k0sReduced,
                                aod::BCs const&)
  {
    processPhiAssocME<kK0S>(collisions, phiCandidates, k0sReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiK0SMEMCWithPDG, "Process function for Phi-K0S 2D Correlations in MC with PDG ME", true);

  void processPhiLambdaMEDataLike(SelCollisions const& collisions,
                                  aod::PhimesonCandidatesData const& phiCandidates,
                                  aod::LambdaReducedCandidatesData const& lambdaReduced,
                                  aod::BCs const&)
  {
    processPhiAssocME<kLambda>(collisions, phiCandidates, lambdaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiLambdaMEDataLike, "Process function for Phi-Lambda 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiLambdaMEMCWithPDG(SimCollisions const& collisions,
                                   aod::PhimesonCandidatesMcReco const& phiCandidates,
                                   aod::LambdaReducedCandidatesMcReco const& lambdaReduced,
                                   aod::BCs const&)
  {
    processPhiAssocME<kLambda>(collisions, phiCandidates, lambdaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiLambdaMEMCWithPDG, "Process function for Phi-Lambda 2D Correlations in MC with PDG ME", true);

  void processPhiAntiLambdaMEDataLike(SelCollisions const& collisions,
                                      aod::PhimesonCandidatesData const& phiCandidates,
                                      aod::AntiLambdaReducedCandidatesData const& antilambdaReduced,
                                      aod::BCs const&)
  {
    processPhiAssocME<kAntiLambda>(collisions, phiCandidates, antilambdaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiAntiLambdaMEDataLike, "Process function for Phi-AntiLambda 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiAntiLambdaMEMCWithPDG(SimCollisions const& collisions,
                                       aod::PhimesonCandidatesMcReco const& phiCandidates,
                                       aod::AntiLambdaReducedCandidatesMcReco const& antilambdaReduced,
                                       aod::BCs const&)
  {
    processPhiAssocME<kAntiLambda>(collisions, phiCandidates, antilambdaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiAntiLambdaMEMCWithPDG, "Process function for Phi-AntiLambda 2D Correlations in MC with PDG ME", true);

  void processPhiXiMEDataLike(SelCollisions const& collisions,
                              aod::PhimesonCandidatesData const& phiCandidates,
                              aod::XiReducedCandidatesData const& xiReduced,
                              aod::BCs const&)
  {
    processPhiAssocME<kXi>(collisions, phiCandidates, xiReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiXiMEDataLike, "Process function for Phi-Xi 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiXiMEMCWithPDG(SimCollisions const& collisions,
                               aod::PhimesonCandidatesMcReco const& phiCandidates,
                               aod::XiReducedCandidatesMcReco const& xiReduced,
                               aod::BCs const&)
  {
    processPhiAssocME<kXi>(collisions, phiCandidates, xiReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiXiMEMCWithPDG, "Process function for Phi-Xi 2D Correlations in MC with PDG ME", true);

  void processPhiOmegaMEDataLike(SelCollisions const& collisions,
                                 aod::PhimesonCandidatesData const& phiCandidates,
                                 aod::OmegaReducedCandidatesData const& omegaReduced,
                                 aod::BCs const&)
  {
    processPhiAssocME<kOmega>(collisions, phiCandidates, omegaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiOmegaMEDataLike, "Process function for Phi-Omega 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiOmegaMEMCWithPDG(SimCollisions const& collisions,
                                  aod::PhimesonCandidatesMcReco const& phiCandidates,
                                  aod::OmegaReducedCandidatesMcReco const& omegaReduced,
                                  aod::BCs const&)
  {
    processPhiAssocME<kOmega>(collisions, phiCandidates, omegaReduced);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiOmegaMEMCWithPDG, "Process function for Phi-Omega 2D Correlations in MC with PDG ME", true);

  void processPhiPionMEDataLike(SelCollisions const& collisions,
                                aod::PhimesonCandidatesData const& phiCandidates,
                                aod::PionTracksData const& pionTracks,
                                aod::BCs const&)
  {
    processPhiAssocME<kPion>(collisions, phiCandidates, pionTracks);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiPionMEDataLike, "Process function for Phi-Pion 2D Correlations in Data or MC w/o PDG ME", true);

  void processPhiPionMEMCWithPDG(SimCollisions const& collisions,
                                 aod::PhimesonCandidatesMcReco const& phiCandidates,
                                 aod::PionTracksMcReco const& pionTracks,
                                 aod::BCs const&)
  {
    processPhiAssocME<kPion>(collisions, phiCandidates, pionTracks);
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiPionMEMCWithPDG, "Process function for Phi-Pion 2D Correlations in MC with PDG ME", true);

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
            if (!activeCorrelationTypes->at(TAssocInputs::Type)) {
              return;
            }

            const auto assocThisColl = assocCands.candidates->sliceBy(getPreslice<TAssocInputs::Type>(), collision.globalIndex());

            for (const auto& assoc : assocThisColl) {
              if (!isAssocValid<TAssocInputs::Type>(assoc)) {
                continue;
              }

              histos.fill(histoKey, collision.posZ(), mcCollision.centFT0M(), assoc.pt(), assoc.y());
            }
          };

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
        auto inYAcceptance = [&]() {
          return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
        };

        auto fillGenHistos = [&](auto h3Key, auto h4Key) {
          histos.fill(h3Key, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          if (hasAssoc) {
            histos.fill(h4Key, zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          }
        };

        auto fillGenAssocSpecies = [&]<AssociatedParticleType PartType>(auto h3Key, auto h4Key) {
          if (!activeCorrelationTypes->at(PartType) || !isGenSpeciesValid<PartType>(mcParticle, false)) {
            return;
          }

          fillGenHistos(h3Key, h4Key);
        };

        if (!inYAcceptance()) {
          continue;
        }

        switch (std::abs(mcParticle.pdgCode())) {
          case o2::constants::physics::Pdg::kPhi:
            if (eventSelectionType == 0 && mcParticle.pt() >= minPtMcGenConfigs.minPhiPt) {
              fillGenHistos(HIST("phi/h3PhiMCGen"), HIST("phi/h4PhiMCGenAssocReco"));
            }
            break;
          /*case PDG_t::kK0Short:
            if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.v0SettingMinPt)
              fillGenHistos(HIST("k0s/h3K0SMCGen"), HIST("k0s/h4K0SMCGenAssocReco"));
            break;
          case PDG_t::kXiMinus:
            if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.cascadeSettingMinPt)
              fillGenHistos(HIST("xi/h3XiMCGen"), HIST("xi/h4XiMCGenAssocReco"));
            break;
          case PDG_t::kPiPlus:
            if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.cMinPionPtcut)
              fillGenHistos(HIST("pi/h3PiMCGen"), HIST("pi/h4PiMCGenAssocReco"));
            break;*/
          case getPdgCode<kK0S>():
            fillGenAssocSpecies.template operator()<kK0S>(HIST("k0s/h3K0SMCGen"), HIST("k0s/h4K0SMCGenAssocReco"));
            break;
          case getPdgCode<kLambda>():
            fillGenAssocSpecies.template operator()<kLambda>(HIST("lambda/h3LambdaMCGen"), HIST("lambda/h4LambdaMCGenAssocReco"));
            break;
          case getPdgCode<kAntiLambda>():
            fillGenAssocSpecies.template operator()<kAntiLambda>(HIST("antilambda/h3AntiLambdaMCGen"), HIST("antilambda/h4AntiLambdaMCGenAssocReco"));
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

  /*void processMCGenClosureSE(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    float multiplicity = mcCollision.centFT0M();

    std::vector<int64_t> phiIndices;
    std::vector<int64_t> k0sIndices;
    std::vector<int64_t> pionIndices;
    std::vector<std::vector<int64_t>> assocIndices;

    auto inYAcceptance = [&](const auto& mcParticle) {
      return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
    };

    for (const auto& mcParticle : mcParticles) {
      if (!inYAcceptance(mcParticle))
        continue;

      switch (std::abs(mcParticle.pdgCode())) {
        case o2::constants::physics::Pdg::kPhi:
          if (eventSelectionType == 0 && mcParticle.pt() >= minPtMcGenConfigs.minPhiPt)
            phiIndices.push_back(mcParticle.globalIndex());
          break;
        case PDG_t::kK0Short:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.v0SettingMinPt)
            k0sIndices.push_back(mcParticle.globalIndex());
          break;
        case PDG_t::kPiPlus:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.cMinPionPtcut)
            pionIndices.push_back(mcParticle.globalIndex());
          break;
        default:
          break;
      }
    }

    assocIndices.push_back(k0sIndices);
    assocIndices.push_back(pionIndices);

    for (std::size_t iTrigg{0}; iTrigg < phiIndices.size(); ++iTrigg) {
      auto& phiParticle = mcParticles.rawIteratorAt(phiIndices[iTrigg]);

      static_for<0, assocIndices.size() - 1>([&](auto i_idx) {
        constexpr unsigned int Idx = i_idx.value;

        for (std::size_t iAssoc{0}; iAssoc < assocIndices[Idx].size(); ++iAssoc) {
          auto& assocParticle = mcParticles.rawIteratorAt(assocIndices[Idx][iAssoc]);

          histos.fill(HIST("mcGenClosure/h5Phi") + HIST(AssocParticleLabels[Idx]) + HIST("ClosureGenSE"), multiplicity, phiParticle.pt(), assocParticle.pt(), phiParticle.y() - assocParticle.y(), getDeltaPhi(phiParticle.phi(), assocParticle.phi()));
        }
      });
    }
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processMCGenClosureSE, "Process function for MC Gen Closure Test in SE", false);

  void processMCGenClosureME(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    float multiplicity = mcCollision.centFT0M();

    std::vector<MiniParticle> phiParticles;
    std::vector<MiniParticle> k0sParticles;
    std::vector<MiniParticle> pionParticles;

    auto inYAcceptance = [&](const auto& mcParticle) {
      return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
    };

    for (const auto& mcParticle : mcParticles) {
      if (!inYAcceptance(mcParticle))
        continue;

      switch (std::abs(mcParticle.pdgCode())) {
        case o2::constants::physics::Pdg::kPhi:
          if (eventSelectionType == 0 && mcParticle.pt() >= minPtMcGenConfigs.minPhiPt)
            phiParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;
        case PDG_t::kK0Short:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.v0SettingMinPt)
            k0sParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;
        case PDG_t::kPiPlus:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.cMinPionPtcut)
            pionParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;
        default:
          break;
      }
    }

    if (phiParticles.empty() && k0sParticles.empty() && pionParticles.empty())
      return;

    int multBin = getCentBin(multiplicity);

    // Loop over past events in the same multiplicity bin and fill histograms with all combinations of current phi particles and past K0S and pion particles
    for (const auto& pastEvent : eventBuffer[multBin]) {
      for (const auto& phiParticle : phiParticles) {
        for (const auto& k0sParticle : pastEvent.k0sParticles) {
          histos.fill(HIST("mcGenClosure/h5PhiK0SClosureGenME"), multiplicity, phiParticle.pt, k0sParticle.pt, phiParticle.y - k0sParticle.y, getDeltaPhi(phiParticle.phi, k0sParticle.phi));
        }
        for (const auto& pionParticle : pastEvent.pionParticles) {
          histos.fill(HIST("mcGenClosure/h5PhiPiClosureGenME"), multiplicity, phiParticle.pt, pionParticle.pt, phiParticle.y - pionParticle.y, getDeltaPhi(phiParticle.phi, pionParticle.phi));
        }
      }
    }

    // Add current event to buffer
    MiniEvent currentEvent;
    currentEvent.multiplicity = multiplicity;
    currentEvent.phiParticles = std::move(phiParticles);
    currentEvent.k0sParticles = std::move(k0sParticles);
    currentEvent.pionParticles = std::move(pionParticles);

    eventBuffer[multBin].push_front(std::move(currentEvent));
    if (eventBuffer[multBin].size() > static_cast<std::size_t>(cfgNoMixedEvents.value))
      eventBuffer[multBin].pop_back();
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processMCGenClosureME, "Process function for MC Gen Closure Test in ME", false);*/

  void processMCGenClosure(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    float multiplicity = mcCollision.centFT0M();

    std::vector<MiniParticle> phiParticles;
    std::vector<MiniParticle> k0sParticles;
    std::vector<MiniParticle> lambdaParticles;
    std::vector<MiniParticle> antilambdaParticles;
    std::vector<MiniParticle> xiParticles;
    std::vector<MiniParticle> omegaParticles;
    std::vector<MiniParticle> pionParticles;

    // Preliminary loop to fill vectors of particles of interest for the current event, applying pt and y cuts
    for (const auto& mcParticle : mcParticles) {
      auto inYAcceptance = [&]() {
        return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
      };

      auto fillPartCollection = [&]<AssociatedParticleType PartType>(auto& collection) {
        if (!activeCorrelationTypes->at(PartType) || !isGenSpeciesValid<PartType>(mcParticle, minPtMcGenConfigs.bypassPtCut.value)) {
          return;
        }

        collection.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
      };

      if (!inYAcceptance()) {
        continue;
      }

      switch (std::abs(mcParticle.pdgCode())) {
        case o2::constants::physics::Pdg::kPhi:
          if (eventSelectionType == 0 && (minPtMcGenConfigs.bypassPtCut || mcParticle.pt() >= minPtMcGenConfigs.minPhiPt)) {
            phiParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          }
          break;
        /*case PDG_t::kK0Short:
          if (mcParticle.isPhysicalPrimary() && (minPtMcGenConfigs.bypassPtCut || mcParticle.pt() >= minPtMcGenConfigs.v0SettingMinPt))
            k0sParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;
        case PDG_t::kXiMinus:
          if (mcParticle.isPhysicalPrimary() && (minPtMcGenConfigs.bypassPtCut || mcParticle.pt() >= minPtMcGenConfigs.cascadeSettingMinPt))
            xiParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;
        case PDG_t::kPiPlus:
          if (mcParticle.isPhysicalPrimary() && (minPtMcGenConfigs.bypassPtCut || mcParticle.pt() >= minPtMcGenConfigs.cMinPionPtcut))
            pionParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          break;*/
        case getPdgCode<kK0S>():
          fillPartCollection.template operator()<kK0S>(k0sParticles);
          break;
        case getPdgCode<kLambda>():
          fillPartCollection.template operator()<kLambda>(lambdaParticles);
          break;
        case getPdgCode<kAntiLambda>():
          fillPartCollection.template operator()<kAntiLambda>(antilambdaParticles);
          break;
        case getPdgCode<kXi>():
          fillPartCollection.template operator()<kXi>(xiParticles);
          break;
        case getPdgCode<kOmega>():
          fillPartCollection.template operator()<kOmega>(omegaParticles);
          break;
        case getPdgCode<kPion>():
          fillPartCollection.template operator()<kPion>(pionParticles);
          break;
        default:
          break;
      }
    }

    const bool skipK0s = !activeCorrelationTypes->at(kK0S) || k0sParticles.empty();
    const bool skipLambda = !activeCorrelationTypes->at(kLambda) || lambdaParticles.empty();
    const bool skipAntiLambda = !activeCorrelationTypes->at(kAntiLambda) || antilambdaParticles.empty();
    const bool skipXi = !activeCorrelationTypes->at(kXi) || xiParticles.empty();
    const bool skipOmega = !activeCorrelationTypes->at(kOmega) || omegaParticles.empty();
    const bool skipPion = !activeCorrelationTypes->at(kPion) || pionParticles.empty();

    if (phiParticles.empty() && skipK0s && skipLambda && skipAntiLambda && skipXi && skipOmega && skipPion) {
      return;
    }

    int multBin = getCentBin(multiplicity);
    if (multBin < 0) {
      return;
    }

    // Same Event Correlations
    std::array<std::vector<MiniParticle>*, kAssocPartSize> currentAssocParticles{&k0sParticles, &lambdaParticles, &antilambdaParticles,
                                                                                 &xiParticles, &omegaParticles, &pionParticles};

    for (const auto& phiParticle : phiParticles) {
      histos.fill(HIST("phi/h3PhiMCClosureGen"), multiplicity, phiParticle.pt, phiParticle.y);

      static_for<0, AssocParticleLabels.size() - 1>([&](auto i_idx) {
        constexpr unsigned int Idx = i_idx.value;
        if (!activeCorrelationTypes->at(Idx)) {
          return;
        }

        for (const auto& assocParticle : *(currentAssocParticles[Idx])) {
          histos.fill(HIST("phi") + HIST(AssocParticleLabels[Idx]) + HIST("/h5Phi") + HIST(AssocParticleLabels[Idx]) + HIST("ClosureMCGen"),
                      multiplicity, phiParticle.pt, assocParticle.pt,
                      phiParticle.y - assocParticle.y,
                      getDeltaPhi(phiParticle.phi, assocParticle.phi));
        }
      });
    }

    // Mixed Event Correlations
    for (const auto& pastEvent : eventBuffer[multBin]) {
      std::array<const std::vector<MiniParticle>*, kAssocPartSize> pastAssocParticles{&pastEvent.k0sParticles, &pastEvent.lambdaParticles, &pastEvent.antilambdaParticles,
                                                                                      &pastEvent.xiParticles, &pastEvent.omegaParticles, &pastEvent.pionParticles};

      // Loop over past events in the same multiplicity bin and fill histograms with all combinations of current phi particles and past associated particles
      for (const auto& phiParticle : phiParticles) {
        static_for<0, AssocParticleLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int Idx = i_idx.value;
          if (!activeCorrelationTypes->at(Idx)) {
            return;
          }

          for (const auto& assocParticle : *(pastAssocParticles[Idx])) {
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
    currentEvent.k0sParticles = std::move(k0sParticles);
    currentEvent.lambdaParticles = std::move(lambdaParticles);
    currentEvent.antilambdaParticles = std::move(antilambdaParticles);
    currentEvent.xiParticles = std::move(xiParticles);
    currentEvent.omegaParticles = std::move(omegaParticles);
    currentEvent.pionParticles = std::move(pionParticles);

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
