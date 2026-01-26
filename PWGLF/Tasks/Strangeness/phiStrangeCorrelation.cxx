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
  Configurable<int> selectionType{"selectionType", 1, "Selection type: 0 - default selection only, 1 - default + phi meson selection"};

  // Configurable for analysis mode
  Configurable<int> analysisMode{"analysisMode", 1, "Analysis mode: 0 - old method with online normalization, 1 - new method with correlations"};

  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"}; // TO BE REMOVED

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for tracks selection
  struct : ConfigurableGroup {
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<std::vector<float>> cMaxDCArToPVPhi{"cMaxDCArToPVPhi", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Phi"};

    Configurable<bool> cfgIsTOFChecked{"cfgIsTOFChecked", true, "Is TOF checked in PID for pions"};
    Configurable<std::vector<float>> cMaxDCArToPVPion{"cMaxDCArToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Pions"};
    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};
    Configurable<std::vector<float>> cMaxDCAzToPVPion{"cMaxDCAzToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAz cut to PV for Pions"};

    Configurable<float> nSigmaCutTPCKa{"nSigmaCutTPCKa", 2.0f, "Value of the TPC Nsigma cut for Kaons"};
    Configurable<float> nSigmaCutCombinedKa{"nSigmaCutCombinedKa", 2.0f, "Value of the TPC and TOF Nsigma cut for Kaons"};

    Configurable<float> nSigmaCutTPCPrimPion{"nSigmaCutTPCPrimPion", 2.0f, "Value of the TPC Nsigma cut for primary Pions"};
    Configurable<float> nSigmaCutTPCSecPion{"nSigmaCutTPCSecPion", 4.0f, "Value of the TPC Nsigma cut for secondary Pions"};
    Configurable<float> nSigmaCutCombinedPi{"nSigmaCutCombinedPi", 2.0f, "Value of the TPC and TOF Nsigma cut for Pions"};
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
    Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};
    Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};

    Configurable<bool> forceTOF{"forceTOF", false, "force the TOF signal for the PID"};
    Configurable<float> tofPIDThreshold{"tofPIDThreshold", 0.5, "minimum pT after which TOF PID is applicable"};
    Configurable<std::vector<int>> trkPIDspecies{"trkPIDspecies", std::vector<int>{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton}, "Trk sel: Particles species for PID, proton, pion, kaon"};
    Configurable<std::vector<float>> pidTPCMax{"pidTPCMax", std::vector<float>{2.0f, 2.0f, 2.0f}, "maximum nSigma TPC"};
    Configurable<std::vector<float>> pidTOFMax{"pidTOFMax", std::vector<float>{2.0f, 2.0f, 2.0f}, "maximum nSigma TOF"};
  } trackConfigs;

  // Configurables on phi selection
  struct : ConfigurableGroup {
    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<std::pair<float, float>> rangeMPhiSignal{"rangeMPhiSignal", {1.0095f, 1.029f}, "Phi mass range for signal extraction"};
    Configurable<std::pair<float, float>> rangeMPhiSideband{"rangeMPhiSideband", {1.1f, 1.2f}, "Phi mass range for sideband extraction"};
  } phiConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};

  // Configurables for V0 selection
  struct : ConfigurableGroup {
    Configurable<float> v0SettingCosPA{"v0SettingCosPA", 0.98f, "V0 CosPA"};
    Configurable<float> v0SettingRadius{"v0SettingRadius", 0.5f, "v0radius"};
    Configurable<float> v0SettingDCAV0Dau{"v0SettingDCAV0Dau", 1.0f, "DCA V0 Daughters"};
    Configurable<float> v0SettingDCAPosToPV{"v0SettingDCAPosToPV", 0.1f, "DCA Pos To PV"};
    Configurable<float> v0SettingDCANegToPV{"v0SettingDCANegToPV", 0.1f, "DCA Neg To PV"};
    Configurable<float> v0SettingMinPt{"v0SettingMinPt", 0.1f, "V0 min pt"};

    Configurable<bool> cfgFurtherV0Selection{"cfgFurtherV0Selection", false, "Further V0 selection"};
    Configurable<float> ctauK0s{"ctauK0s", 20.0f, "C tau K0s(cm)"};
    Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2f, "parameter Armenteros Cut"};
    Configurable<float> v0rejK0s{"v0rejK0s", 0.005f, "V0 rej K0s"};

    Configurable<float> lowMK0S{"lowMK0S", 0.48f, "Lower limit on K0Short mass"};
    Configurable<float> upMK0S{"upMK0S", 0.52f, "Upper limit on K0Short mass"};
  } v0Configs;

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

  // Constants
  double massPi = o2::constants::physics::MassPiPlus;
  double massK0S = o2::constants::physics::MassK0Short;
  double massLambda = o2::constants::physics::MassLambda0;

  // Filter on phi selected collisions
  Filter collisionFilter = aod::lf_selection_phi_collision::phimesonSel == true;

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter v0PreFilter = (nabs(aod::v0data::dcapostopv) > v0Configs.v0SettingDCAPosToPV && nabs(aod::v0data::dcanegtopv) > v0Configs.v0SettingDCANegToPV && aod::v0data::dcaV0daughters < v0Configs.v0SettingDCAV0Dau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhimesonSelectionData>>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  // using SimCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhimesonSelectionMcReco, aod::McCollisionLabels>>;
  using MCCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::PhimesonSelectionMcGen>>;

  // Defining the type of the V0s and corresponding daughter tracks for data and MC
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Join<FullV0s, aod::McV0Labels>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  // using FilteredTracks = soa::Filtered<FullTracks>;
  // using FilteredMCTracks = soa::Filtered<FullMCTracks>;

  // Preslice for manual slicing
  struct : PresliceGroup {
    // Preslice<SimCollisions> collPerMCCollision = aod::mccollisionlabel::mcCollisionId;
    Preslice<FullMCV0s> v0PerCollision = aod::v0::collisionId;
    Preslice<FullMCTracks> trackPerCollision = aod::track::collisionId;
    // Preslice<aod::PhimesonCandidatesData> phiCandDataPerCollision = aod::lf_selection_phi_candidate::collisionId;
    // PresliceUnsorted<SimCollisions> collPerMCCollision = aod::mccollisionlabel::mcCollisionId;
    PresliceUnsorted<aod::PhimesonCandidatesMcReco> phiCandPerCollision = aod::lf_selection_phi_candidate::collisionId;

    // Preslice<aod::McParticles> mcPartPerMCCollision = aod::mcparticle::mcCollisionId;
  } preslices;

  // Slice cache for mixed event
  SliceCache cache;

  // Necessary service to retrieve efficiency maps from CCDB
  Service<ccdb::BasicCCDBManager> ccdb;

  // Efficiency maps
  /*std::shared_ptr<TH3> effMapPhi = nullptr;
  std::shared_ptr<TH3> effMapK0S = nullptr;
  std::shared_ptr<TH3> effMapPionTPC = nullptr;
  std::shared_ptr<TH3> effMapPionTPCTOF = nullptr;*/

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

    // histos.add("phiK0S/h5PhiK0SDataNewProc", "2D Invariant mass of Phi and K0Short in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, massPhiAxis});
    // histos.add("phiPi/h5PhiPiTPCDataNewProc", "Phi Invariant mass vs Pion nSigma TPC in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
    // histos.add("phiPi/h5PhiPiTOFDataNewProc", "Phi Invariant mass vs Pion nSigma TOF in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});

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

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC)
      return false;

    if (std::abs(track.eta()) > trackConfigs.etaMax)
      return false;
    return true;
  }

  // V0 selection
  template <bool isMC, typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& collision)
  {
    using V0DauTrackType = std::conditional_t<isMC, V0DauMCTracks, V0DauTracks>;

    const auto& posDaughterTrack = v0.template posTrack_as<V0DauTrackType>();
    const auto& negDaughterTrack = v0.template negTrack_as<V0DauTrackType>();

    // const auto& posDaughterTrack = v0.template posTrack_as<V0DauTracks>();
    // const auto& negDaughterTrack = v0.template negTrack_as<V0DauTracks>();

    if (!selectionTrackStrangeness(posDaughterTrack) || !selectionTrackStrangeness(negDaughterTrack))
      return false;

    if constexpr (!isMC) {
      if (std::abs(posDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion)
        return false;
      if (std::abs(negDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCSecPion)
        return false;
    }

    if (v0.v0cosPA() < v0Configs.v0SettingCosPA)
      return false;
    if (v0.v0radius() < v0Configs.v0SettingRadius)
      return false;
    if (v0.pt() < v0Configs.v0SettingMinPt)
      return false;

    if (v0Configs.cfgFurtherV0Selection) {
      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0S > v0Configs.ctauK0s)
        return false;
      if (v0.qtarm() < (v0Configs.paramArmenterosCut * std::abs(v0.alpha())))
        return false;
      if (std::abs(v0.mLambda() - massLambda) < v0Configs.v0rejK0s)
        return false;
    }

    if (std::abs(v0.yK0Short()) > yConfigs.cfgYAcceptance)
      return false;

    return true;
  }

  // PID selection for Pions
  template <typename T>
  bool pidSelectionPion(const T& track)
  {
    for (size_t speciesIndex = 0; speciesIndex < trackConfigs.trkPIDspecies->size(); ++speciesIndex) {
      auto const& pid = trackConfigs.trkPIDspecies->at(speciesIndex);
      auto nSigmaTPC = aod::pidutils::tpcNSigma(pid, track);

      if (trackConfigs.forceTOF && !track.hasTOF()) {
        return false;
      }

      if (speciesIndex == 0) { // First species logic
        if (std::abs(nSigmaTPC) >= trackConfigs.pidTPCMax->at(speciesIndex)) {
          return false; // TPC check failed
        }
        if (trackConfigs.forceTOF || (track.pt() >= trackConfigs.tofPIDThreshold && track.hasTOF())) {
          auto nSigmaTOF = aod::pidutils::tofNSigma(pid, track);
          if (std::abs(nSigmaTOF) >= trackConfigs.pidTOFMax->at(speciesIndex)) {
            return false; // TOF check failed
          }
        }
      } else {                                                                // Other species logic
        if (std::abs(nSigmaTPC) < trackConfigs.pidTPCMax->at(speciesIndex)) { // Check TPC nSigma  first
          if (track.hasTOF()) {
            auto nSigmaTOF = aod::pidutils::tofNSigma(pid, track);
            if (std::abs(nSigmaTOF) < trackConfigs.pidTOFMax->at(speciesIndex)) {
              return false; // Reject if both TPC and TOF are within thresholds
            }
          } else {
            return false; // Reject if only TPC is within threshold and TOF is unavailable
          }
        }
      }
    }

    return true;
  }

  // Track selection for Pions
  template <typename T>
  bool selectionPion(const T& track)
  {
    if (!track.isGlobalTrackWoDCA())
      return false;

    if (track.itsNCls() < trackConfigs.minITSnCls)
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinPionPtcut)
      return false;

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPion->at(0) + (trackConfigs.cMaxDCArToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPion->at(2))))
      return false;
    if (trackConfigs.cfgIsDCAzParameterized) {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVPion->at(0) + (trackConfigs.cMaxDCAzToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCAzToPVPion->at(2))))
        return false;
    } else {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
        return false;
    }

    if (trackConfigs.cfgIsTOFChecked && track.pt() >= trackConfigs.tofPIDThreshold && !track.hasTOF())
      return false;

    if (analysisMode == 1 && !pidSelectionPion(track))
      return false;

    /*
    if (analysisMode == 1) {
      if (track.pt() < trackConfigs.tofPIDThreshold && std::abs(track.tpcNSigmaPi()) >= trackConfigs.nSigmaCutTPCPrimPion)
        return false;
      if (trackConfigs.cfgIsTOFChecked && track.pt() >= trackConfigs.tofPIDThreshold && (std::pow(track.tofNSigmaPi(), 2) + std::pow(track.tpcNSigmaPi(), 2)) >= std::pow(trackConfigs.nSigmaCutCombinedPi, 2))
        return false;
    }
    */

    if (std::abs(track.rapidity(massPi)) > yConfigs.cfgYAcceptance)
      return false;

    return true;
  }

  void processPhiK0SPionData(SelCollisions::iterator const& collision, aod::PhimesonCandidatesData const& phiCandidates, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    float multiplicity = collision.centFT0M();

    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    for (const auto& phiCand : phiCandidates) {
      float weightPhi = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()));

      // float weightPhi = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()));

      histos.fill(HIST("phi/h3PhiData"), multiplicity, phiCand.pt(), phiCand.m(), weightPhi);

      static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
        constexpr unsigned int i = i_idx.value;

        const auto& [minMass, maxMass] = phiMassRegions[i];
        if (!phiCand.inMassRegion(minMass, maxMass))
          return;

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          // Cut on V0 dynamic columns
          if (!selectionV0<false>(v0, collision))
            continue;

          float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                             BoundEfficiencyMap(effMaps[K0S], multiplicity, v0.pt(), v0.yK0Short()));

          /*float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()),
                                             BoundEfficiencyMap(effMapK0S, multiplicity, v0.pt(), v0.yK0Short()));*/

          histos.fill(HIST("phiK0S/h5PhiK0SData") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), v0.pt(), phiCand.y() - v0.yK0Short(), getDeltaPhi(phiCand.phi(), v0.phi()), weightPhiK0S);
        }

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion(track))
            continue;

          // auto Pion = track.pt() < trackConfigs.tofPIDThreshold ? PionTPC : PionTPCTOF;

          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                              BoundEfficiencyMap(effMaps[Pion], multiplicity, track.pt(), track.rapidity(massPi)));

          /*auto effMapPion = track.pt() < trackConfigs.tofPIDThreshold ? effMapPionTPC : effMapPionTPCTOF;

          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()),
                                              BoundEfficiencyMap(effMapPion, multiplicity, track.pt(), track.rapidity(massPi)));*/

          histos.fill(HIST("phiPi/h5PhiPiData") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), track.pt(), phiCand.y() - track.rapidity(massPi), getDeltaPhi(phiCand.phi(), track.phi()), weightPhiPion);
        }
      });
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiK0SPionData, "Process function for Phi-K0S and Phi-Pion Deltay and Deltaphi 2D Correlations in Data", true);

  void processPhiK0SDataME(SelCollisions const& collisions, aod::PhimesonCandidatesData const& phiCandidates, FullV0s const& V0s, V0DauTracks const&)
  {
    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    auto tuplePhiV0 = std::make_tuple(phiCandidates, V0s);
    Pair<SelCollisions, aod::PhimesonCandidatesData, FullV0s, BinningTypeVertexCent> pairPhiK0S{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiV0, &cache};

    for (const auto& [c1, phiCands, c2, v0s] : pairPhiK0S) {

      float multiplicity = c1.centFT0M();

      for (const auto& [phiCand, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, v0s))) {
        static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int i = i_idx.value;

          const auto& [minMass, maxMass] = phiMassRegions[i];
          if (!phiCand.inMassRegion(minMass, maxMass))
            return;

          if (!selectionV0<false>(v0, c2))
            return;

          float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                             BoundEfficiencyMap(effMaps[K0S], multiplicity, v0.pt(), v0.yK0Short()));

          histos.fill(HIST("phiK0S/h5PhiK0SDataME") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), v0.pt(), phiCand.y() - v0.yK0Short(), getDeltaPhi(phiCand.phi(), v0.phi()), weightPhiK0S);
        });
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiK0SDataME, "Process function for Phi-K0S and Deltay and Deltaphi 2D Correlations in Data ME", false);

  void processPhiPionDataME(SelCollisions const& collisions, aod::PhimesonCandidatesData const& phiCandidates, FullTracks const& fullTracks)
  {
    const std::array<std::pair<float, float>, 2> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    auto tuplePhiPion = std::make_tuple(phiCandidates, fullTracks);
    Pair<SelCollisions, aod::PhimesonCandidatesData, FullTracks, BinningTypeVertexCent> pairPhiPion{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiPion, &cache};

    for (const auto& [c1, phiCands, c2, tracks] : pairPhiPion) {

      float multiplicity = c1.centFT0M();

      for (const auto& [phiCand, track] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, tracks))) {
        static_for<0, phiMassRegionLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int i = i_idx.value;

          const auto& [minMass, maxMass] = phiMassRegions[i];
          if (!phiCand.inMassRegion(minMass, maxMass))
            return;

          if (!selectionPion(track))
            return;

          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMaps[Phi], multiplicity, phiCand.pt(), phiCand.y()),
                                              BoundEfficiencyMap(effMaps[Pion], multiplicity, track.pt(), track.rapidity(massPi)));

          histos.fill(HIST("phiPi/h5PhiPiDataME") + HIST(phiMassRegionLabels[i]), multiplicity, phiCand.pt(), track.pt(), phiCand.y() - track.rapidity(massPi), getDeltaPhi(phiCand.phi(), track.phi()), weightPhiPion);
        });
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiPionDataME, "Process function for Phi-Pion Deltay and Deltaphi 2D Correlations in Data ME", false);

  void processParticleEfficiency(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, aod::McParticles const& mcParticles, aod::PhimesonCandidatesMcReco const& phiCandidatesMcReco)
  {
    uint16_t numberAssocColls{0};
    std::vector<float> zVtxs;

    // const auto collsThisMCColl = collisions.sliceBy(preslices.collPerMCCollision, mcCollision.globalIndex());

    for (const auto& collision : collisions) {
      histos.fill(HIST("event/hRecoMCMultiplicityPercent"), mcCollision.centFT0M());
      histos.fill(HIST("event/h2RecoMCVertexZvsMult"), collision.posZ(), mcCollision.centFT0M());

      zVtxs.push_back(collision.posZ());

      if (selectionType == 0) {
        const auto phiCandidatesThisColl = phiCandidatesMcReco.sliceBy(preslices.phiCandPerCollision, collision.globalIndex());
        for (const auto& phiCand : phiCandidatesThisColl) {
          histos.fill(HIST("phi/h4PhiMCReco"), collision.posZ(), mcCollision.centFT0M(), phiCand.pt(), phiCand.y());
        }
      }

      const auto v0sThisColl = V0s.sliceBy(preslices.v0PerCollision, collision.globalIndex());
      const auto fullMCTracksThisColl = fullMCTracks.sliceBy(preslices.trackPerCollision, collision.globalIndex());

      for (const auto& v0 : v0sThisColl) {
        if (!selectionV0<true>(v0, collision))
          continue;

        if (!v0.has_mcParticle())
          continue;

        const auto& v0McParticle = mcParticles.rawIteratorAt(v0.mcParticleId());
        if (std::abs(v0McParticle.pdgCode()) != PDG_t::kK0Short || !v0McParticle.isPhysicalPrimary())
          continue;

        histos.fill(HIST("k0s/h4K0SMCReco"), collision.posZ(), mcCollision.centFT0M(), v0McParticle.pt(), v0McParticle.y());
      }

      for (const auto& track : fullMCTracksThisColl) {
        if (!selectionPion(track))
          continue;

        if (!track.has_mcParticle())
          continue;

        const auto& trackMcParticle = mcParticles.rawIteratorAt(track.mcParticleId());
        if (std::abs(trackMcParticle.pdgCode()) != PDG_t::kPiPlus)
          continue;

        if (trackMcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("pi/h2RecMCDCAxyPrimPi"), track.pt(), track.dcaXY());
        } else {
          if (trackMcParticle.getProcess() == TMCProcess::kPDecay) { // Selection of secondary pions from weak decay
            histos.fill(HIST("pi/h2RecMCDCAxySecWeakDecayPi"), track.pt(), track.dcaXY());
          } else { // Selection of secondary pions from material interactions
            histos.fill(HIST("pi/h2RecMCDCAxySecMaterialPi"), track.pt(), track.dcaXY());
          }
          continue;
        }

        histos.fill(HIST("pi/h4PiMCReco"), collision.posZ(), mcCollision.centFT0M(), trackMcParticle.pt(), trackMcParticle.y());
      }

      numberAssocColls++;
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

    auto inYAcceptance = [&](const auto& mcParticle) {
      return std::abs(mcParticle.y()) <= yConfigs.cfgYAcceptance;
    };

    auto fillGenHistos = [&](auto h3Key, auto h4Key, const auto& mcParticle) {
      histos.fill(h3Key, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
      if (hasAssoc)
        histos.fill(h4Key, zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
    };

    for (const auto& mcParticle : mcParticles /*| std::views::filter(inYAcceptance)*/) {
      if (!inYAcceptance(mcParticle))
        continue;

      switch (std::abs(mcParticle.pdgCode())) {
        case o2::constants::physics::Pdg::kPhi:
          if (selectionType == 0 && mcParticle.pt() >= phiConfigs.minPhiPt)
            fillGenHistos(HIST("phi/h3PhiMCGen"), HIST("phi/h4PhiMCGenAssocReco"), mcParticle);
          break;
        case PDG_t::kK0Short:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= v0Configs.v0SettingMinPt)
            fillGenHistos(HIST("k0s/h3K0SMCGen"), HIST("k0s/h4K0SMCGenAssocReco"), mcParticle);
          break;
        case PDG_t::kPiPlus:
          if (mcParticle.isPhysicalPrimary() && mcParticle.pt() >= trackConfigs.cMinPionPtcut)
            fillGenHistos(HIST("pi/h3PiMCGen"), HIST("pi/h4PiMCGenAssocReco"), mcParticle);
          break;
        default:
          break;
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processParticleEfficiency, "Process function for Efficiency Computation for Particles of Interest", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiStrangenessCorrelation>(cfgc)};
}
