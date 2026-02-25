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
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <Framework/BinningPolicy.h>
#include <Framework/SliceCache.h>
#include <Framework/StaticFor.h>

#include <Math/Vector4D.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TList.h>
#include <TMCProcess.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr std::array<std::string_view, 2> phiMassRegionLabels{"Signal", "Sideband"};

enum ParticleOfInterest {
  Phi = 0,
  K0S,
  Pion,
  /*PionTPC,
  PionTPCTOF*/
  ParticleOfInterestSize
};

static constexpr std::array<std::string_view, ParticleOfInterestSize> particleOfInterestLabels{"Phi", "K0S", "Pion" /*"PionTPC", "PionTPCTOF"*/};

/*
#define LIST_OF_PARTICLES_OF_INTEREST \
  X(Phi)                          \
  X(K0S)                          \
  X(Pion)                         \
  //X(PionTPC)                      \
  //X(PionTPCTOF)

enum ParticleOfInterest {
#define X(name) name,
  LIST_OF_PARTICLES_OF_INTEREST
#undef X
  ParticleOfInterestSize
};

static constexpr std::array<std::string_view, ParticleOfInterestSize> particleOfInterestLabels{
#define X(name) #name,
  LIST_OF_PARTICLES_OF_INTEREST
#undef X
};

static constexpr auto particleOfInterestLabels = std::to_array<std::string_view>({
#define X(name) #name,
  LIST_OF_PARTICLES_OF_INTEREST
#undef X
});
*/

struct BoundEfficiencyMap {
  using CoordsTuple = std::tuple<float, float, float>;

  const TH3* effMap;
  CoordsTuple coords;

  BoundEfficiencyMap(const std::shared_ptr<TH3>& effMap, float x, float y, float z) : effMap(effMap.get()), coords(x, y, z) {}
  BoundEfficiencyMap(const std::shared_ptr<TH3>& effMap, const CoordsTuple& coords) : effMap(effMap.get()), coords(coords) {}

  float getBinEfficiency() const
  {
    if (!effMap) {
      return 1.0f;
    }

    const auto& [x, y, z] = coords;
    return effMap->GetBinContent(effMap->FindFixBin(x, y, z));
  }

  float interpolateEfficiency() const
  {
    if (!effMap) {
      return 1.0f;
    }

    const auto& [x, y, z] = coords;
    return effMap->Interpolate(x, y, z);
  }
};

struct PhiStrangenessCorrelation {
  HistogramRegistry histos{"phiStrangenessCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for selection type
  Configurable<int> eventSelectionType{"eventSelectionType", 0, "Event selection type: 0 - default selection only, 1 - default + phi meson selection"};

  // Configurable for analysis mode
  Configurable<int> analysisMode{"analysisMode", 1, "Analysis mode: 0 - old method with online normalization, 1 - new method with correlations"};

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
    Configurable<bool> selectK0sInSigRegion{"selectK0sInSigRegion", false, "Select K0s candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMK0sSignal{"rangeMK0sSignal", {0.47f, 0.53f}, "K0S mass range for signal extraction"};
  } k0sConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  // Configurables for delta y selection
  struct : ConfigurableGroup {
    Configurable<int> nBinsY{"nBinsY", 20, "Number of bins in y axis"};
    Configurable<int> nBinsDeltaY{"nBinsDeltaY", 20, "Number of bins in deltay axis"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
    Configurable<std::vector<float>> cfgDeltaYAcceptanceBins{"cfgDeltaYAcceptanceBins", {0.5f}, "Rapidity acceptance bins"};
  } yConfigs;

  // Configurables to apply efficiency online and how to
  Configurable<bool> applyEfficiency{"applyEfficiency", false, "Use efficiency for filling histograms"};
  Configurable<bool> useEffInterpolation{"useEffInterpolation", false, "If true, interpolates efficiency map, else uses bin center"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurables for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbEfficiencyPath{"ccdbEfficiencyPath", "Users/s/scannito/Efficiencies", "Correction path to file"};

  // Configurables for minimum pt selection in McGen
  struct : ConfigurableGroup {
    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<float> v0SettingMinPt{"v0SettingMinPt", 0.1f, "V0 min pt"};
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};
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
    PresliceUnsorted<aod::PhimesonCandidatesMcReco> phiCandMcRecoPerCollision = aod::lf_selection_phi_candidate::collisionId;

    Preslice<aod::K0sReducedCandidatesData> k0sDataPerCollision = aod::v0::collisionId;
    Preslice<aod::K0sReducedCandidatesMcReco> k0sMcRecoPerCollision = aod::v0::collisionId;

    Preslice<aod::PionTracksData> pionTrackDataPerCollision = aod::track::collisionId;
    Preslice<aod::PionTracksMcReco> pionTrackMcRecoPerCollision = aod::track::collisionId;

    Preslice<aod::McParticles> mcPartPerMcCollision = aod::mcparticle::mcCollisionId;
  } preslices;

  // Necessary service to retrieve efficiency maps from CCDB
  Service<ccdb::BasicCCDBManager> ccdb;

  std::array<std::shared_ptr<TH3>, ParticleOfInterestSize> effMaps{};

  // Binning policy and axes for mixed event
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {20, -10, 10}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {20, 0, 100}, "Multiplicity percentil binning for mixing"};

  using BinningTypeVertexCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeVertexCent binningOnVertexAndCent{{axisVertexMixing, axisCentralityMixing}, true};

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
    AxisSpec pTK0SAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPiAxis = {50, 0.0f, 5.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("phi/h3PhiData", "Invariant mass of Phi in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    for (const auto& label : phiMassRegionLabels) {
      histos.add(fmt::format("phiK0S/h5PhiK0SData{}", label).c_str(), "Deltay vs deltaphi for Phi and K0Short in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phiPi/h5PhiPiData{}", label).c_str(), "Deltay vs deltaphi for Phi and Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, deltaphiAxis});

      histos.add(fmt::format("phiK0S/h5PhiK0SDataME{}", label).c_str(), "Deltay vs deltaphi for Phi and K0Short in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
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

    histos.add("pi/h4PiMCReco", "Pion in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPiAxis, yAxis});
    histos.add("pi/h3PiMCGen", "Pion in MC Gen", kTH3F, {binnedmultAxis, binnedpTPiAxis, yAxis});
    histos.add("pi/h4PiMCGenAssocReco", "Pion in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPiAxis, yAxis});

    histos.add("pi/h2RecMCDCAxyPrimPi", "Dcaxy distribution vs pt for Primary Pions", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("pi/h2RecMCDCAxySecWeakDecayPi", "Dcaz distribution vs pt for Secondary Pions from Weak Decay", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("pi/h2RecMCDCAxySecMaterialPi", "Dcaxy distribution vs pt for Secondary Pions from Material", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});

    // Load efficiency maps from CCDB
    if (applyEfficiency) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      for (int i = 0; i < ParticleOfInterestSize; ++i) {
        loadEfficiencyMapFromCCDB(static_cast<ParticleOfInterest>(i));
      }
    }
  }

  void loadEfficiencyMapFromCCDB(ParticleOfInterest poi)
  {
    effMaps[poi] = std::shared_ptr<TH3>(ccdb->get<TH3F>(fmt::format("{}/h3EffMap{}", ccdbEfficiencyPath.value, particleOfInterestLabels[poi])));
    if (!effMaps[poi])
      LOG(fatal) << "Could not load efficiency map for " << particleOfInterestLabels[poi] << "!";
    LOG(info) << "Efficiency map for " << particleOfInterestLabels[poi] << " loaded from CCDB";
  }

  // Compute weight based on efficiencies
  template <typename... BoundEffMaps>
  float computeWeight(const BoundEffMaps&... boundEffMaps)
  {
    if (!applyEfficiency)
      return 1.0f;

    float totalEfficiency = ((useEffInterpolation ? boundEffMaps.interpolateEfficiency() : boundEffMaps.getBinEfficiency()) * ...);

    return totalEfficiency <= 0.0f ? 1.0f : 1.0f / totalEfficiency;
  }

  float getDeltaPhi(float phiTrigger, float phiAssociated)
  {
    return RecoDecay::constrainAngle(phiTrigger - phiAssociated, -o2::constants::math::PIHalf);
  }

  void processPhiK0SPionData(SelCollisions::iterator const& collision, aod::PhimesonCandidatesData const& phiCandidates, aod::K0sReducedCandidatesData const& k0sReduced, aod::PionTracksData const& pionTracks)
  {
    float multiplicity = collision.centFT0M();

    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    for (const auto& phiCand : phiCandidates) {
      float weightPhi = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()));

      histos.fill(HIST("phi/h3PhiData"), multiplicity, phiCand.pt(), phiCand.m(), weightPhi);

      static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
        constexpr unsigned int i = i_idx.value;

        const auto& [minMassPhi, maxMassPhi] = phiMassRegions[i];
        if (!phiCand.inMassRegion(minMassPhi, maxMassPhi))
          return;

        // Loop over all reduced K0S candidates
        for (const auto& k0s : k0sReduced) {
          if (k0sConfigs.selectK0sInSigRegion) {
            const auto& [minMassK0s, maxMassK0s] = k0sConfigs.rangeMK0sSignal.value;
            if (!k0s.inMassRegion(minMassK0s, maxMassK0s))
              continue;
          }

          float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                             BoundEfficiencyMap(effMaps[K0S], multiplicity, k0s.pt(), k0s.y()));

          histos.fill(HIST("phiK0S/h5PhiK0SData") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), k0s.pt(), phiCand.y() - k0s.y(), getDeltaPhi(phiCand.phi(), k0s.phi()), weightPhiK0S);
        }

        // Loop over all primary pion candidates
        for (const auto& pionTrack : pionTracks) {
          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                              BoundEfficiencyMap(effMaps[Pion], multiplicity, pionTrack.pt(), pionTrack.y()));

          histos.fill(HIST("phiPi/h5PhiPiData") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), pionTrack.pt(), phiCand.y() - pionTrack.y(), getDeltaPhi(phiCand.phi(), pionTrack.phi()), weightPhiPion);
        }
      });
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiK0SPionData, "Process function for Phi-K0S and Phi-Pion Deltay and Deltaphi 2D Correlations in Data", true);

  void processPhiK0SDataME(SelCollisions const& collisions, aod::PhimesonCandidatesData const& phiCandidates, aod::K0sReducedCandidatesData const& k0sReduced)
  {
    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    auto tuplePhiK0S = std::make_tuple(phiCandidates, k0sReduced);
    Pair<SelCollisions, aod::PhimesonCandidatesData, aod::K0sReducedCandidatesData, BinningTypeVertexCent> pairPhiK0S{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiK0S, &cache};

    for (const auto& [c1, phiCands, c2, k0sRed] : pairPhiK0S) {

      float multiplicity = c1.centFT0M();

      for (const auto& [phiCand, k0s] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, k0sRed))) {
        static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int i = i_idx.value;

          const auto& [minMassPhi, maxMassPhi] = phiMassRegions[i];
          if (!phiCand.inMassRegion(minMassPhi, maxMassPhi))
            return;

          if (k0sConfigs.selectK0sInSigRegion) {
            const auto& [minMassK0s, maxMassK0s] = k0sConfigs.rangeMK0sSignal.value;
            if (!k0s.inMassRegion(minMassK0s, maxMassK0s))
              return;
          }

          float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                             BoundEfficiencyMap(effMaps[K0S], multiplicity, k0s.pt(), k0s.y()));

          histos.fill(HIST("phiK0S/h5PhiK0SDataME") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), k0s.pt(), phiCand.y() - k0s.y(), getDeltaPhi(phiCand.phi(), k0s.phi()), weightPhiK0S);
        });
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiK0SDataME, "Process function for Phi-K0S and Deltay and Deltaphi 2D Correlations in Data ME", false);

  void processPhiPionDataME(SelCollisions const& collisions, aod::PhimesonCandidatesData const& phiCandidates, aod::PionTracksData const& pionTracks)
  {
    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    auto tuplePhiPion = std::make_tuple(phiCandidates, pionTracks);
    Pair<SelCollisions, aod::PhimesonCandidatesData, aod::PionTracksData, BinningTypeVertexCent> pairPhiPion{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiPion, &cache};

    for (const auto& [c1, phiCands, c2, piTracks] : pairPhiPion) {

      float multiplicity = c1.centFT0M();

      for (const auto& [phiCand, piTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, piTracks))) {
        static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int i = i_idx.value;

          const auto& [minMassPhi, maxMassPhi] = phiMassRegions[i];
          if (!phiCand.inMassRegion(minMassPhi, maxMassPhi))
            return;

          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                              BoundEfficiencyMap(effMaps[Pion], multiplicity, piTrack.pt(), piTrack.y()));

          histos.fill(HIST("phiPi/h5PhiPiDataME") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), piTrack.pt(), phiCand.y() - piTrack.y(), getDeltaPhi(phiCand.phi(), piTrack.phi()), weightPhiPion);
        });
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiPionDataME, "Process function for Phi-Pion Deltay and Deltaphi 2D Correlations in Data ME", false);

  void processParticleEfficiency(MCCollisions const& mcCollisions, SimCollisions const& collisions, aod::PhimesonCandidatesMcReco const& phiCandidates, aod::K0sReducedCandidatesMcReco const& k0sReduced, aod::PionTracksMcReco const& pionTracks, aod::McParticles const& mcParticles)
  {
    // std::vector<std::vector<int>> collsGrouped(mcCollisions.size());
    std::unordered_map<int, std::vector<int>> collsGrouped;
    collsGrouped.reserve(mcCollisions.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
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

          const auto k0sThisColl = k0sReduced.sliceBy(preslices.k0sMcRecoPerCollision, collision.globalIndex());
          const auto pionTracksThisColl = pionTracks.sliceBy(preslices.pionTrackMcRecoPerCollision, collision.globalIndex());

          for (const auto& k0s : k0sThisColl) {
            histos.fill(HIST("k0s/h4K0SMCReco"), collision.posZ(), mcCollision.centFT0M(), k0s.pt(), k0s.y());
          }

          for (const auto& pionTrack : pionTracksThisColl) {
            histos.fill(HIST("pi/h4PiMCReco"), collision.posZ(), mcCollision.centFT0M(), pionTrack.pt(), pionTrack.y());
          }

          numberAssocColls++;
        }
      }

      histos.fill(HIST("event/hGenMCMultiplicityPercent"), mcCollision.centFT0M());

      const bool hasAssoc = (numberAssocColls > 0);
      const float zVtxRef = hasAssoc ? zVtxs[0] : 0.0f;

      //////TOBECHANGED//////
      if (hasAssoc) {
        if (zVtxs.size() > 1) {
          for (size_t i = 1; i < zVtxs.size(); ++i) {
            histos.fill(HIST("event/hSplitVertexZ"), zVtxs[i] - zVtxRef);
          }
        }

        histos.fill(HIST("event/hGenMCAssocRecoMultiplicityPercent"), mcCollision.centFT0M());
        histos.fill(HIST("event/h2GenMCAssocRecoVertexZvsMult"), zVtxRef, mcCollision.centFT0M());
      }
      ///////////////////////

      const auto mcParticlesThisMcColl = mcParticles.sliceBy(preslices.mcPartPerMcCollision, mcCollision.globalIndex());

      auto inYAcceptance = [&](const auto& mcParticle) {
        return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
      };

      auto fillGenHistos = [&](auto h3Key, auto h4Key, const auto& mcParticle) {
        histos.fill(h3Key, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
        if (hasAssoc)
          histos.fill(h4Key, zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
      };

      for (const auto& mcParticle : mcParticlesThisMcColl /*| std::views::filter(inYAcceptance)*/) {
        if (!inYAcceptance(mcParticle))
          continue;

        switch (std::abs(mcParticle.pdgCode())) {
          case o2::constants::physics::Pdg::kPhi:
            if (eventSelectionType == 0 && mcParticle.pt() >= minPtMcGenConfigs.minPhiPt)
              fillGenHistos(HIST("phi/h3PhiMCGen"), HIST("phi/h4PhiMCGenAssocReco"), mcParticle);
            break;
          case PDG_t::kK0Short:
            if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.v0SettingMinPt)
              fillGenHistos(HIST("k0s/h3K0SMCGen"), HIST("k0s/h4K0SMCGenAssocReco"), mcParticle);
            break;
          case PDG_t::kPiPlus:
            if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= minPtMcGenConfigs.cMinPionPtcut)
              fillGenHistos(HIST("pi/h3PiMCGen"), HIST("pi/h4PiMCGenAssocReco"), mcParticle);
            break;
          default:
            break;
        }
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processParticleEfficiency, "Process function for Efficiency Computation for Particles of Interest", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiStrangenessCorrelation>(cfgc)};
}
