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

/// \file doubleOmegaTrackIndexSkimCreator.cxx
/// \brief Table producer for V0s and Cascade tables from tracks selected by a track-level BDT.
///        Designed to improve the efficiency in the double Omega search, the BDT is supposed to
///        discriminate between primary tracks and daughters of Lambda and Omega decays,
///        regardless of the fact that those particles are from a double Omega or not.
///        A table producer for the AO2Ds required for the training of the BDT is also provided.
///
/// \author Giorgio Alberto Lucia (giorgio.alberto.lucia@cern.ch)

#include "PWGLF/DataModel/LFDoubleOmegaTrackTables.h"
#include "PWGLF/Utils/doubleOmegaBuilderHelper.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <GPUROOTCartesianFwd.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::pwglf::doubleomega;

using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using TracksData = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using TracksMC = soa::Join<TracksData, aod::McTrackLabels>;

struct doubleOmegaTrackIndexSkimCreator {
  Produces<aod::DOmegaTrackTrain> trackTable;
  Produces<aod::DOmegaV0s> v0Table;
  Produces<aod::DOmegaCascades> cascadeTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;

  static constexpr int kDoubleOmegaPdg = 1060020020;

  enum Role : int8_t {
    kBackground = 0,
    kOmegaProton,
    kOmegaPion,
    kOmegaKaon,
    kDirectLambdaProton,
    kDirectLambdaPion,
    kDirectKaon,
    kOtherLambdaProton,
    kOtherLambdaPion,
    kOtherOmegaKaon,
    kNRoles
  };
  static constexpr std::array<const char*, kNRoles> RoleLabels{
    "Background", "p #leftarrow #Lambda #leftarrow #Omega #leftarrow #Omega#Omega", "#pi #leftarrow #Lambda #leftarrow #Omega #leftarrow #Omega#Omega",
    "K #leftarrow #Omega #leftarrow #Omega#Omega", "p #leftarrow #Lambda #leftarrow #Omega#Omega", "#pi #leftarrow #Lambda #leftarrow #Omega#Omega",
    "K #leftarrow #Omega#Omega", "p #leftarrow other #Lambda", "#pi #leftarrow other #Lambda", "K #leftarrow other #Omega"};

  // Species of the chain starts in hGeneratedAncestors (charge conjugates included)
  enum Species : int8_t {
    kSpeciesDoubleOmega = 0,
    kSpeciesOmega,
    kSpeciesXi,
    kSpeciesLambda,
    kSpeciesPion,
    kSpeciesKaon,
    kSpeciesProton,
    kNSpecies
  };
  static constexpr std::array<const char*, kNSpecies> SpeciesLabels{
    "#Omega#Omega", "#Omega", "#Xi^{-}, #Xi^{0}", "#Lambda", "#pi^{#pm}", "K^{#pm}", "p"};
  static constexpr int kMaxChainDepth = 10;

  // Classes of the track BDT used for the candidates (candOpts.classIndices, classThresholds)
  enum TrackClass : uint8_t {
    kClassProton = 0,
    kClassPion,
    kClassKaon,
    kNTrackClasses
  };

  Preslice<TracksData> tracksDataPerCollision = aod::track::collisionId;
  Preslice<TracksMC> tracksMCPerCollision = aod::track::collisionId;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Material correction for the propagation to the primary vertex"};
  Configurable<bool> cfgApplyEventSelection{"cfgApplyEventSelection", true, "Apply the standard collision event selection"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum |z| of the primary vertex (cm)"};

  // Track preselection, the baseline of the selection: it should be looser than any selection studied with it.
  // No track-quality cut (e.g. TPC crossed rows): the quality variables are stored and left to the selections.
  Configurable<bool> requireTPC{"requireTPC", true, "Require TPC (PID of all the daughters)"};
  Configurable<float> etaMax{"etaMax", 1.f, "Maximum |eta|"};
  Configurable<float> genEtaMax{"genEtaMax", 1.f, "Maximum |eta| of the generated daughters (efficiency denominator)"};
  ConfigurableAxis ptAxis{"ptAxis", {100, 0.f, 10.f}, "pT binning of the efficiency and purity histograms (GeV/c)"};
  Configurable<float> ancestorRapidityMax{"ancestorRapidityMax", 0.5f, "Maximum |y| of the generated chain starts (hGeneratedAncestors)"};
  Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT (GeV/c)"};

  // Downsampling (roles 1-6 are always stored)
  Configurable<float> fractionBackground{"fractionBackground", 0.01f, "Fraction of the background tracks (role 0) stored"};
  Configurable<float> fractionOtherDaughters{"fractionOtherDaughters", 0.1f, "Fraction of the daughters of other Lambdas and Omegas (roles 7-9) stored"};
  Configurable<int> randomSeed{"randomSeed", 0, "Seed of the downsampling (0: random)"};

  // V0 and cascade candidates (processCandidates), from the preselected tracks above selected by the track BDT
  struct : ConfigurableGroup {
    std::string prefix = "candOpts";
    Configurable<std::string> onnxFileName{"onnxFileName", "model.onnx", "ONNX file of the track BDT (local path, or file name on CCDB)"};
    Configurable<bool> loadModelFromCCDB{"loadModelFromCCDB", false, "Load the ONNX model from CCDB"};
    Configurable<std::string> modelPathCCDB{"modelPathCCDB", "Users/g/galucia/DoubleOmega/TrackBDT", "CCDB path of the model"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "Timestamp of the model on CCDB"};
    Configurable<std::vector<std::string>> inputFeatures{"inputFeatures", std::vector<std::string>{"fPt", "fTpcInnerParam", "fTpcNSigmaPr", "fTpcNSigmaPi", "fTpcNSigmaKa", "fTofNSigmaPr", "fTofNSigmaPi", "fTofNSigmaKa", "fDcaXY", "fDcaZ", "fTpcNClsCrossedRows", "fTpcNClsShared", "fTpcChi2NCl", "fIsPVContributor"}, "Model inputs, in order (model_inputs of config.json of the training; names of PWGLF/Utils/doubleOmegaBuilderHelper.h)"};
    Configurable<int> nClasses{"nClasses", 4, "Number of model outputs (classes)"};
    Configurable<std::vector<int>> classIndices{"classIndices", std::vector<int>{1, 2, 3}, "Output index of the proton, pion and kaon classes (class_names of config.json)"};
    Configurable<std::vector<float>> classThresholds{"classThresholds", std::vector<float>{0.5f, 0.5f, 0.5f}, "Minimum probability of the proton, pion and kaon classes"};
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.05f, "Lambda mass window of the V0s, from the momenta at the innermost update (GeV/c^2; negative: no selection)"};
    Configurable<float> omegaMassWindow{"omegaMassWindow", 0.1f, "Omega mass window of the cascades, from the momenta at the innermost update and the Lambda mass (GeV/c^2; negative: no selection)"};
  } candOpts;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::analysis::MlResponse<float> mlResponse; // track BDT
  std::vector<int> mlFeatureIndices;          // model inputs in the feature array of computeTrackFeatures
  std::array<int, kNTrackClasses> mlClassIndex{};
  std::array<float, kNTrackClasses> mlThreshold{};

  int mRunNumber = 0;
  std::mt19937 mRandom;
  std::uniform_real_distribution<float> mUniform{0.f, 1.f};
  std::unordered_map<int64_t, Role> mRoles; // MC particle index -> role, per dataframe
  std::unordered_set<int64_t> mAcceptedMcCollisions;
  int64_t mDataframe = 0;

  enum Stage : uint8_t {
    kPreselected = 0,
    kStored,
    kNStages
  };

  void init(InitContext const&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    mRandom.seed(randomSeed.value != 0 ? static_cast<unsigned>(randomSeed.value) : std::random_device{}());

    histos.add("hEvents", ";;Collisions", HistType::kTH1F, {{2, -0.5f, 1.5f}});
    histos.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(1, "All");
    histos.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(2, "Selected");
    const AxisSpec roleAxis{kNRoles, -0.5f, static_cast<float>(kNRoles) - 0.5f, "Role"};
    auto hTracks = histos.add<TH3>("hTracks", ";Role;Stage;#it{p}_{T} (GeV/#it{c})", HistType::kTH3F, {roleAxis, {kNStages, -0.5f, static_cast<float>(kNStages) - 0.5f}, ptAxis});
    hTracks->GetYaxis()->SetBinLabel(kPreselected + 1, "Preselected");
    hTracks->GetYaxis()->SetBinLabel(kStored + 1, "Stored");
    auto hGenerated = histos.add<TH2>("hGenerated", ";Role;#it{p}_{T}^{gen} (GeV/#it{c})", HistType::kTH2F, {roleAxis, ptAxis});
    auto hGeneratedAncestorPt = histos.add<TH3>("hGeneratedAncestorPt", ";Role;#it{p}_{T}^{gen} (GeV/#it{c});#it{p}_{T}^{gen} of the first particle of the decay chain (GeV/#it{c})",
                                                HistType::kTH3F, {roleAxis, ptAxis, ptAxis});
    auto hGeneratedAncestors = histos.add<TH2>("hGeneratedAncestors", ";Species;#it{p}_{T}^{gen} (GeV/#it{c})", HistType::kTH2F,
                                               {{kNSpecies, -0.5f, static_cast<float>(kNSpecies) - 0.5f, "Species"}, ptAxis});
    for (int species = 0; species < kNSpecies; ++species) {
      hGeneratedAncestors->GetXaxis()->SetBinLabel(species + 1, SpeciesLabels[species]);
    }
    for (auto* axis : {hTracks->GetXaxis(), hGenerated->GetXaxis(), hGeneratedAncestorPt->GetXaxis()}) {
      for (int role = 0; role < kNRoles; ++role) {
        axis->SetBinLabel(role + 1, RoleLabels[role]);
      }
    }

    if (doprocessCandidates) {
      initTrackMl();
      histos.add("Candidates/hEvents", ";;Collisions", HistType::kTH1F, {{2, -0.5f, 1.5f}});
      histos.get<TH1>(HIST("Candidates/hEvents"))->GetXaxis()->SetBinLabel(1, "All");
      histos.get<TH1>(HIST("Candidates/hEvents"))->GetXaxis()->SetBinLabel(2, "Selected");
      auto hPools = histos.add<TH2>("Candidates/hPools", ";Class;#it{p}_{T} (GeV/#it{c})", HistType::kTH2F, {{2 * kNTrackClasses, -0.5f, static_cast<float>(2 * kNTrackClasses) - 0.5f}, ptAxis});
      for (const auto& [bin, label] : {std::pair{1, "p"}, std::pair{2, "#bar{p}"}, std::pair{3, "#pi^{+}"}, std::pair{4, "#pi^{-}"}, std::pair{5, "K^{+}"}, std::pair{6, "K^{-}"}}) {
        hPools->GetXaxis()->SetBinLabel(bin, label);
      }
      histos.add("Candidates/hMassLambda", ";#it{p}_{T} (GeV/#it{c});Sign;#it{M}(p#pi) (GeV/#it{c}^{2})", HistType::kTH3F, {ptAxis, {2, -2.f, 2.f}, {200, 1.05f, 1.25f}});
      histos.add("Candidates/hMassOmega", ";#it{p}_{T} (GeV/#it{c});Sign;#it{M}(#LambdaK) (GeV/#it{c}^{2})", HistType::kTH3F, {ptAxis, {2, -2.f, 2.f}, {200, 1.57f, 1.77f}});
    }
  }

  void initTrackMl()
  {
    mlFeatureIndices = trackFeatureIndices(candOpts.inputFeatures.value);
    const auto& indices = candOpts.classIndices.value;
    const auto& thresholds = candOpts.classThresholds.value;
    if (indices.size() != kNTrackClasses || thresholds.size() != kNTrackClasses) {
      LOG(fatal) << "candOpts.classIndices and candOpts.classThresholds need one entry per class (proton, pion, kaon)";
    }
    for (int iClass = 0; iClass < kNTrackClasses; ++iClass) {
      if (indices[iClass] < 0 || indices[iClass] >= candOpts.nClasses) {
        LOG(fatal) << "candOpts.classIndices: index " << indices[iClass] << " out of the " << candOpts.nClasses << " model outputs";
      }
      mlClassIndex[iClass] = indices[iClass];
      mlThreshold[iClass] = thresholds[iClass];
    }

    // One model for all pT; the class thresholds are applied here, not by MlResponse.
    const int nClasses = candOpts.nClasses;
    std::vector<double> noCuts(nClasses, 0.);
    std::vector<std::string> labelsClasses;
    labelsClasses.reserve(nClasses);
    for (int iClass = 0; iClass < nClasses; ++iClass) {
      labelsClasses.push_back(Form("class %d", iClass));
    }
    mlResponse.configure(std::vector<double>{0., 1.e6}, LabeledArray<double>{noCuts.data(), 1, static_cast<uint32_t>(nClasses), {"all pT"}, labelsClasses},
                         std::vector<int>(nClasses, o2::cuts_ml::CutNot), static_cast<uint8_t>(nClasses));
    if (candOpts.loadModelFromCCDB) {
      ccdbApi.init(ccdburl);
      mlResponse.setModelPathsCCDB(std::vector<std::string>{candOpts.onnxFileName.value}, ccdbApi, std::vector<std::string>{candOpts.modelPathCCDB.value}, candOpts.timestampCCDB);
    } else {
      mlResponse.setModelPathsLocal(std::vector<std::string>{candOpts.onnxFileName.value});
    }
    mlResponse.init();
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Could not retrieve GRPMagField for timestamp " << bc.timestamp();
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", bc.timestamp()));
      if (!lut) {
        LOG(fatal) << "Could not retrieve material LUT for timestamp " << bc.timestamp();
      }
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    mRunNumber = bc.runNumber();
  }

  // isCandidates: counted in Candidates/hEvents (processCandidates) instead of hEvents (processMC)
  template <bool isCandidates = false, class C>
  bool acceptCollision(C const& collision)
  {
    initCCDB(collision.template bc_as<aod::BCsWithTimestamps>());
    const auto fillEvents = [&](int bin) {
      if constexpr (isCandidates) {
        histos.fill(HIST("Candidates/hEvents"), bin);
      } else {
        histos.fill(HIST("hEvents"), bin);
      }
    };
    fillEvents(0);
    if (cfgApplyEventSelection &&
        (!collision.sel8() || std::abs(collision.posZ()) > zVtxMax)) {
      return false;
    }
    fillEvents(1);
    return true;
  }

  template <class T>
  bool preselectTrack(T const& track) const
  {
    return track.sign() != 0 &&
           (!requireTPC.value || track.hasTPC()) &&
           std::abs(track.eta()) <= etaMax.value &&
           track.pt() >= ptMin.value;
  }

  bool keep(Role role)
  {
    switch (role) {
      case kBackground:
        return mUniform(mRandom) < fractionBackground;
      case kOtherLambdaProton:
      case kOtherLambdaPion:
      case kOtherOmegaKaon:
        return mUniform(mRandom) < fractionOtherDaughters;
      default:
        return true;
    }
  }

  // ---------------------------------------------------------------------------------------------
  // MC roles

  // Daughter of a given PDG code, -1 if none
  template <class McParticle>
  static int64_t findDaughter(McParticle const& particle, int pdgCode)
  {
    for (const auto& daughter : particle.template daughters_as<aod::McParticles>()) {
      if (daughter.pdgCode() == pdgCode) {
        return daughter.globalIndex();
      }
    }
    return -1;
  }

  // First particle of the decay chain of an MC particle: the double Omega if it is an ancestor, otherwise the
  // first physical primary going up the mothers (the particle itself if it is a physical primary)
  static int64_t chainStart(aod::McParticles const& mcParticles, int64_t id)
  {
    for (int depth = 0; depth < kMaxChainDepth; ++depth) {
      const auto particle = mcParticles.rawIteratorAt(id);
      if (!particle.has_mothers()) {
        break;
      }
      const int64_t motherId = particle.mothersIds()[0];
      if (motherId < 0 || motherId >= mcParticles.size()) {
        break;
      }
      if (std::abs(mcParticles.rawIteratorAt(motherId).pdgCode()) == kDoubleOmegaPdg) {
        return motherId;
      }
      if (particle.isPhysicalPrimary()) {
        break;
      }
      id = motherId;
    }
    return id;
  }

  // Species of a chain start in hGeneratedAncestors, -1 if none
  static int species(int pdgCode)
  {
    switch (std::abs(pdgCode)) {
      case kDoubleOmegaPdg:
        return kSpeciesDoubleOmega;
      case PDG_t::kOmegaMinus:
        return kSpeciesOmega;
      case PDG_t::kXiMinus:
      case 3322: // Xi0
        return kSpeciesXi;
      case PDG_t::kLambda0:
        return kSpeciesLambda;
      case PDG_t::kPiPlus:
        return kSpeciesPion;
      case PDG_t::kKPlus:
        return kSpeciesKaon;
      case PDG_t::kProton:
        return kSpeciesProton;
      default:
        return -1;
    }
  }

  // Proton and pion of a Lambda (sign: of the Lambda PDG code)
  void assignLambdaDaughters(aod::McParticles const& mcParticles, int64_t lambdaId, int sign, Role protonRole, Role pionRole)
  {
    if (lambdaId < 0) {
      return;
    }
    const auto lambda = mcParticles.rawIteratorAt(lambdaId);
    if (const auto id = findDaughter(lambda, sign * PDG_t::kProton); id >= 0) {
      mRoles.emplace(id, protonRole);
    }
    if (const auto id = findDaughter(lambda, -sign * PDG_t::kPiPlus); id >= 0) {
      mRoles.emplace(id, pionRole);
    }
  }

  // Roles of the MC particles of the dataframe. The double-Omega chains are assigned first, so that their
  // Lambdas and Omega are not counted again as "other" ones (emplace does not overwrite).
  void assignRoles(aod::McParticles const& mcParticles)
  {
    mRoles.clear();
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.pdgCode()) != kDoubleOmegaPdg) {
        continue;
      }
      const int sign = particle.pdgCode() > 0 ? 1 : -1;
      const int64_t omegaId = findDaughter(particle, sign * PDG_t::kOmegaMinus);
      const int64_t directLambdaId = findDaughter(particle, sign * PDG_t::kLambda0);
      if (const auto id = findDaughter(particle, -sign * PDG_t::kKPlus); id >= 0) {
        mRoles.emplace(id, kDirectKaon);
      }
      assignLambdaDaughters(mcParticles, directLambdaId, sign, kDirectLambdaProton, kDirectLambdaPion);
      if (omegaId >= 0) {
        const auto omega = mcParticles.rawIteratorAt(omegaId);
        if (const auto id = findDaughter(omega, -sign * PDG_t::kKPlus); id >= 0) {
          mRoles.emplace(id, kOmegaKaon);
        }
        assignLambdaDaughters(mcParticles, findDaughter(omega, sign * PDG_t::kLambda0), sign, kOmegaProton, kOmegaPion);
      }
    }
    for (const auto& particle : mcParticles) {
      const int pdgCode = particle.pdgCode();
      const int sign = pdgCode > 0 ? 1 : -1;
      if (std::abs(pdgCode) == PDG_t::kLambda0) {
        assignLambdaDaughters(mcParticles, particle.globalIndex(), sign, kOtherLambdaProton, kOtherLambdaPion);
      } else if (std::abs(pdgCode) == PDG_t::kOmegaMinus) {
        if (const auto id = findDaughter(particle, -sign * PDG_t::kKPlus); id >= 0) {
          mRoles.emplace(id, kOtherOmegaKaon);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------------------------
  // Filling

  template <class C, class T>
  void fillTrack(C const& collision, T const& track, Role role, int pdgCode, int motherPdgCode, bool isPhysicalPrimary, float ptMC, float productionRadius, int64_t mcKey,
                 int ancestorPdgCode, float ancestorPtMC)
  {
    // DCA to the primary vertex with the uncertainties of the track and of the vertex
    auto trackParCov = getTrackParCov(track);
    const o2::dataformats::VertexBase primaryVertex{o2::math_utils::Point3D<float>{collision.posX(), collision.posY(), collision.posZ()},
                                                    std::array<float, 6>{collision.covXX(), collision.covXY(), collision.covYY(),
                                                                         collision.covXZ(), collision.covYZ(), collision.covZZ()}};
    o2::dataformats::DCA dca;
    float dcaXY = -999.f, dcaZ = -999.f, sigmaDcaXY = -999.f, sigmaDcaZ = -999.f;
    if (o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCov, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dca)) {
      dcaXY = dca.getY();
      dcaZ = dca.getZ();
      sigmaDcaXY = std::sqrt(std::max(0.f, static_cast<float>(dca.getSigmaY2())));
      sigmaDcaZ = std::sqrt(std::max(0.f, static_cast<float>(dca.getSigmaZ2())));
    }
    const uint8_t detectorMap = static_cast<uint8_t>(track.hasITS()) | (static_cast<uint8_t>(track.hasTPC()) << 1) |
                                (static_cast<uint8_t>(track.hasTRD()) << 2) | (static_cast<uint8_t>(track.hasTOF()) << 3);
    trackTable(track.pt(), track.eta(), track.phi(), static_cast<int8_t>(track.sign()),
               track.tpcInnerParam(), track.tpcSignal(),
               track.tpcNSigmaPr(), track.tpcNSigmaPi(), track.tpcNSigmaKa(),
               track.hasTOF() ? track.tofNSigmaPr() : -999.f, track.hasTOF() ? track.tofNSigmaPi() : -999.f, track.hasTOF() ? track.tofNSigmaKa() : -999.f,
               dcaXY, dcaZ, sigmaDcaXY, sigmaDcaZ,
               track.itsClusterMap(), track.itsClusterSizes(), track.itsChi2NCl(),
               static_cast<int16_t>(track.tpcNClsFound()), static_cast<int16_t>(track.tpcNClsCrossedRows()),
               track.tpcNClsShared(), track.tpcChi2NCl(),
               detectorMap, track.isPVContributor(),
               static_cast<uint16_t>(collision.numContrib()), collision.posZ(),
               role, pdgCode, motherPdgCode, isPhysicalPrimary, ptMC, productionRadius, mcKey,
               ancestorPdgCode, ancestorPtMC);
  }

  void processMC(CollisionsMC const& collisions, TracksMC const& tracks, aod::McParticles const& mcParticles,
                 aod::McCollisions const&, aod::BCsWithTimestamps const&)
  {
    ++mDataframe;
    assignRoles(mcParticles);
    mAcceptedMcCollisions.clear();
    for (const auto& collision : collisions) {
      if (!acceptCollision(collision)) {
        continue;
      }
      if (collision.has_mcCollision()) {
        mAcceptedMcCollisions.insert(collision.mcCollisionId());
      }
      for (const auto& track : tracks.sliceBy(tracksMCPerCollision, collision.globalIndex())) {
        Role role = kBackground;
        int pdgCode = 0, motherPdgCode = 0;
        bool isPhysicalPrimary = false;
        float ptMC = -999.f, productionRadius = -999.f, ancestorPtMC = -999.f;
        int64_t mcKey = -1;
        int ancestorPdgCode = 0;
        if (track.has_mcParticle()) {
          const auto mcParticle = track.template mcParticle_as<aod::McParticles>();
          if (const auto found = mRoles.find(mcParticle.globalIndex()); found != mRoles.end()) {
            role = found->second;
          }
          pdgCode = mcParticle.pdgCode();
          for (const auto& mother : mcParticle.template mothers_as<aod::McParticles>()) {
            motherPdgCode = mother.pdgCode();
            break;
          }
          isPhysicalPrimary = mcParticle.isPhysicalPrimary();
          ptMC = mcParticle.pt();
          productionRadius = std::hypot(mcParticle.vx(), mcParticle.vy());
          mcKey = (mDataframe << 32) | mcParticle.globalIndex();
          const auto ancestor = mcParticles.rawIteratorAt(chainStart(mcParticles, mcParticle.globalIndex()));
          ancestorPdgCode = ancestor.pdgCode();
          ancestorPtMC = ancestor.pt();
        }
        if (!preselectTrack(track)) {
          continue;
        }
        histos.fill(HIST("hTracks"), role, kPreselected, track.pt());
        if (keep(role)) {
          histos.fill(HIST("hTracks"), role, kStored, track.pt());
          fillTrack(collision, track, role, pdgCode, motherPdgCode, isPhysicalPrimary, ptMC, productionRadius, mcKey,
                    ancestorPdgCode, ancestorPtMC);
        }
      }
    }

    // Generated daughters (efficiency denominator)
    for (const auto& [id, role] : mRoles) {
      const auto particle = mcParticles.rawIteratorAt(id);
      if (!mAcceptedMcCollisions.contains(particle.mcCollisionId()) || std::abs(particle.eta()) > genEtaMax) {
        continue;
      }
      histos.fill(HIST("hGenerated"), role, particle.pt());
      histos.fill(HIST("hGeneratedAncestorPt"), role, particle.pt(), mcParticles.rawIteratorAt(chainStart(mcParticles, id)).pt());
    }

    // Generated chain starts: the MC pT shape of the double Omegas and of the primary hadrons
    for (const auto& particle : mcParticles) {
      const int index = species(particle.pdgCode());
      if (index < 0 || (index != kSpeciesDoubleOmega && !particle.isPhysicalPrimary()) ||
          !mAcceptedMcCollisions.contains(particle.mcCollisionId()) || std::abs(particle.y()) > ancestorRapidityMax) {
        continue;
      }
      histos.fill(HIST("hGeneratedAncestors"), index, particle.pt());
    }
  }
  PROCESS_SWITCH(doubleOmegaTrackIndexSkimCreator, processMC, "Store the tracks of MC with their role (signal: double-Omega daughters, background downsampled)", true);

  // ---------------------------------------------------------------------------------------------
  // V0 and cascade candidates

  struct PoolTrack {
    int64_t id = -1;
    std::array<float, 3> momentum{};
  };
  // Tracks of a collision selected by the BDT, per class and charge ([0] negative, [1] positive)
  using Pools = std::array<std::array<std::vector<PoolTrack>, 2>, kNTrackClasses>;

  struct StoredV0 {
    int64_t index = -1; // row of DOmegaV0s
    int64_t pionId = -1;
    int8_t sign = 0;
    std::array<float, 3> momentum{};
  };

  template <class C, class T>
  Pools fillPools(C const& collision, T const& collisionTracks)
  {
    Pools pools;
    std::array<float, kNTrackFeatures> features{};
    std::vector<float> inputs(mlFeatureIndices.size());
    for (const auto& track : collisionTracks) {
      if (!preselectTrack(track)) {
        continue;
      }
      computeTrackFeatures(collision, track, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), features);
      for (std::size_t iInput = 0; iInput < mlFeatureIndices.size(); ++iInput) {
        inputs[iInput] = features[mlFeatureIndices[iInput]];
      }
      const auto scores = mlResponse.getModelOutput(inputs, 0);
      const int iCharge = track.sign() > 0 ? 1 : 0;
      for (int iClass = 0; iClass < kNTrackClasses; ++iClass) {
        if (scores[mlClassIndex[iClass]] <= mlThreshold[iClass]) {
          continue;
        }
        pools[iClass][iCharge].push_back({track.globalIndex(), {track.px(), track.py(), track.pz()}});
        histos.fill(HIST("Candidates/hPools"), 2 * iClass + 1 - iCharge, track.pt());
      }
    }
    return pools;
  }

  template <class C>
  void buildCandidates(C const& collision, Pools const& pools)
  {
    constexpr auto MassLambda = o2::constants::physics::MassLambda0;
    constexpr auto MassOmega = o2::constants::physics::MassOmegaMinus;
    std::vector<StoredV0> v0s;

    // Lambda: p pi-, anti-Lambda: pbar pi+
    for (const int8_t sign : {int8_t{1}, int8_t{-1}}) {
      const auto& protons = pools[kClassProton][sign > 0 ? 1 : 0];
      const auto& pions = pools[kClassPion][sign > 0 ? 0 : 1];
      for (const auto& proton : protons) {
        for (const auto& pion : pions) {
          const auto momentum = RecoDecay::pVec(proton.momentum, pion.momentum);
          const float mass = RecoDecay::m(std::array{proton.momentum, pion.momentum}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
          if (candOpts.lambdaMassWindow >= 0.f && std::abs(mass - MassLambda) > candOpts.lambdaMassWindow) {
            continue;
          }
          histos.fill(HIST("Candidates/hMassLambda"), RecoDecay::pt(momentum), sign, mass);
          const int64_t posId = sign > 0 ? proton.id : pion.id;
          const int64_t negId = sign > 0 ? pion.id : proton.id;
          v0Table(collision.globalIndex(), posId, negId, sign);
          v0s.push_back({v0Table.lastIndex(), pion.id, sign, momentum});
        }
      }
    }

    // Omega-: Lambda K-, Omega+: anti-Lambda K+; the bachelor has the charge of the pion of the V0
    for (const auto& v0 : v0s) {
      const int8_t omegaSign = -v0.sign;
      for (const auto& kaon : pools[kClassKaon][omegaSign > 0 ? 1 : 0]) {
        if (kaon.id == v0.pionId) {
          continue;
        }
        const float mass = RecoDecay::m(std::array{v0.momentum, kaon.momentum}, std::array{MassLambda, o2::constants::physics::MassKaonCharged});
        if (candOpts.omegaMassWindow >= 0.f && std::abs(mass - MassOmega) > candOpts.omegaMassWindow) {
          continue;
        }
        histos.fill(HIST("Candidates/hMassOmega"), RecoDecay::pt(RecoDecay::pVec(v0.momentum, kaon.momentum)), omegaSign, mass);
        cascadeTable(collision.globalIndex(), v0.index, kaon.id, omegaSign);
      }
    }
  }

  void processCandidates(Collisions const& collisions, TracksData const& tracks, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      if (!acceptCollision<true>(collision)) {
        continue;
      }
      const auto pools = fillPools(collision, tracks.sliceBy(tracksDataPerCollision, collision.globalIndex()));
      buildCandidates(collision, pools);
    }
  }
  PROCESS_SWITCH(doubleOmegaTrackIndexSkimCreator, processCandidates, "V0 and cascade candidates (DOmegaV0s, DOmegaCascades) from the tracks selected by the track BDT, data or MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<doubleOmegaTrackIndexSkimCreator>(cfgc)};
}
