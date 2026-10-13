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

#include "PWGLF/DataModel/LFDoubleOmegaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <GPU/GPUROOTCartesianFwd.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <chrono>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsTable = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>;
using Collisions = CollisionsTable::iterator;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults, aod::McCollisionLabels>;
using FullCascades = aod::CascDataExt;
// Innermost-update tracks: DCAs to the PV are recomputed in this task by propagation,
// so neither the propagated tracks nor TracksDCA are needed.
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFKa>;
using TracksFullMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFKa, aod::McTrackLabels>;

struct DoubleOmegaCandidate {
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float x = -999.f;
  float y = -999.f;
  float z = -999.f;
  float pvX = -999.f;
  float pvY = -999.f;
  float pvZ = -999.f;
  float cosPAOmega = -999.f;
  float cosPADirectLambda = -999.f;
  float cosPADoubleOmega = -999.f;
  float cosPAOmegaToPV = -999.f;
  float cosPADirectLambdaToPV = -999.f;
  float dcaXYOmegaToPV = -999.f;
  float dcaXYDirectLambdaToPV = -999.f;
  float dcaXYDirectKaonToPV = -999.f;
  float dcaXYOmegaKaonToPV = -999.f;
  float dcaXYOmegaProtonToPV = -999.f;
  float dcaXYOmegaPionToPV = -999.f;
  float dcaXYDirectProtonToPV = -999.f;
  float dcaXYDirectPionToPV = -999.f;
  float tpcNSigmaDirectKaon = -999.f;
  float tofNSigmaDirectKaon = -999.f;
  uint32_t daughterDetectorMap = 0;
  float mass = -999.f;
  float massOmega = -999.f;
  float massXi = -999.f;
  int8_t sign = 0;
};

struct DoubleOmegaMCInfo {
  int64_t motherId = -1;
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float decayLength = -999.f;
  int pdgCode = 0;
  // Generated vertices: double-Omega production point (primary vertex), double-Omega
  // decay point, and Omega decay point.
  std::array<float, 3> genProductionVertex{-999.f, -999.f, -999.f};
  std::array<float, 3> genDecayVertex{-999.f, -999.f, -999.f};
  std::array<float, 3> genOmegaDecayVertex{-999.f, -999.f, -999.f};
  int8_t sign = 0;
};

struct LambdaCandidate {
  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float mass = -999.f;
};

struct BuiltLambdaCandidate {
  LambdaCandidate candidate;
  o2::track::TrackParCov parentTrack;
  std::array<float, 3> decayVertex{};
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
  uint32_t daughterDetectorMap = 0; // Proton, pion; three detector bits each.
  float dcaXYProtonToPV = -999.f;
  float dcaXYPionToPV = -999.f;
};

struct BuiltOmegaCandidate {
  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  float massOmega = -999.f;
  float massXi = -999.f;
  int8_t sign = 0;
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
  int64_t bachelorId = -1;
  o2::track::TrackParCov parentTrack;
  uint32_t daughterDetectorMap = 0; // Kaon, proton, pion; three detector bits each.
  float dcaXYKaonToPV = -999.f;
  float dcaXYProtonToPV = -999.f;
  float dcaXYPionToPV = -999.f;
};

struct BuiltKaonCandidate {
  o2::track::TrackParCov track;
  int64_t trackId = -1;
  int8_t sign = 0;
  float tpcNSigma = -999.f; // -999 when TPC is unavailable.
  float tofNSigma = -999.f; // -999 when TOF is unavailable.
  uint32_t detectorMap = 0;
};

// Omega, Lambda and kaon candidates of one collision, built once and paired afterwards.
// Lambdas and kaons are stored per charge of the Omega they can be paired with (see chargeIndex).
struct CollisionCandidates {
  std::array<float, 3> primaryVertex{};
  std::vector<BuiltOmegaCandidate> omegas;
  std::array<std::vector<BuiltLambdaCandidate>, 2> lambdas;
  std::array<std::vector<BuiltKaonCandidate>, 2> kaons;
};

// Same-event double-Omega candidate used as a seed of the candidate-level event mixing. Seeds are
// kept in a pool across dataframes, so they store copies of their daughters, with the collision and
// track indices made unique across dataframes (see uniqueIndex).
struct MixingSeed {
  int64_t collisionId = -1;
  std::array<float, 3> primaryVertex{};
  float bz = 0.f;
  BuiltOmegaCandidate omega;
  BuiltLambdaCandidate lambda;
  BuiltKaonCandidate kaon;
};

struct doubleOmegaTreeCreator {
  Produces<aod::DoubleOmegaTable> doubleOmegaTable;
  Produces<aod::DoubleOmegaTableMC> doubleOmegaTableMC;
  Produces<aod::DoubleOmegaMixedTable> doubleOmegaMixedTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  static constexpr int kDoubleOmegaPdg = 1060020020;
  enum FindabilityStep : uint8_t {
    kAllGenerated,
    kOmegaLambdaKaon,
    kTwoLambdasTwoKaons,
    kFinalState,
    kFindable,
    kFindableSelectedTracks,
    kFindableCascadeAndV0,
    kFindableCascadeAndV0SelectedTracks,
    kNFindabilitySteps
  };
  enum MixingType : uint8_t {
    kSameEvent = 0,
    kMixedOmega,  // Lambda and kaon from the same collision, Omega from another one
    kMixedLambda, // Omega and kaon from the same collision, Lambda from another one
    kMixedKaon,   // Omega and Lambda from the same collision, kaon from another one
    kMixedAll,    // Omega, Lambda and kaon from three different collisions
    kNMixingTypes
  };
  std::vector<int64_t> reconstructedDoubleOmegaIds;

  Preslice<TracksFullMC> tracksPerCollision = aod::track::collisionId;
  Preslice<TracksFull> tracksPerCollisionData = aod::track::collisionId;
  Preslice<aod::Cascades> rawCascadesPerCollision = aod::cascade::collisionId;
  Preslice<aod::V0s> rawV0sPerCollision = aod::v0::collisionId;

  int mRunNumber = 0;
  float mBz = 0.f;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3Body;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<bool> cfgApplyEventSelection{"cfgApplyEventSelection", true, "Apply the standard collision event selection"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Material correction for the raw V0/cascade fits"};

  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};
  ConfigurableAxis massOmegaAxis{"massOmegaAxis", {400, o2::constants::physics::MassOmegaMinus - 0.05f, o2::constants::physics::MassOmegaMinus + 0.05f}, "binning for the Omega invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {400, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};
  ConfigurableAxis massLambdaAxis{"massLambdaAxis", {400, o2::constants::physics::MassLambda0 - 0.05f, o2::constants::physics::MassLambda0 + 0.05f}, "binning for the Lambda invariant-mass"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  Configurable<float> cascPtMin{"cascPtMin", 1.f, "minimum (anti)cascade pT (GeV/c)"};
  Configurable<float> cascPtMax{"cascPtMax", 5.f, "maximum (anti)cascade pT (GeV/c)"};

  Configurable<float> minNCrossedRows{"minNCrossedRows", 100, "Minimum number of crossed TPC rows"};
  Configurable<float> minNITSClus{"minNITSClus", 0., "Minimum number of ITS clusters"};
  Configurable<float> minNTPCClus{"minNTPCClus", 80, "Minimum number of TPC clusters"};
  Configurable<float> maxNSharedTPCClus{"maxNSharedTPCClus", 5, "Maximum number of shared TPC clusters"};

  Configurable<float> nSigmaTPCCut{"nSigmaTPCCut", 3.f, "Number of sigmas for the TPC PID"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.05f, "Minimum DCA of a cascade bachelor to the primary vertex"};
  Configurable<float> dcaKaonToPV{"dcaKaonToPV", 0.05f, "Minimum absolute transverse DCA of the direct kaon to the primary vertex"};
  Configurable<float> dcaOmegaToPV{"dcaOmegaToPV", 0.f, "Minimum absolute transverse DCA of the Omega to the primary vertex"};
  Configurable<float> dcaDirectLambdaToPV{"dcaDirectLambdaToPV", 0.f, "Minimum absolute transverse DCA of the direct Lambda to the primary vertex"};
  Configurable<float> dcaV0DauToPV{"dcaV0DauToPV", 0.05f, "Minimum DCA of Lambda daughters to the primary vertex"};
  Configurable<float> dcaV0Bach{"dcaV0Bach", 1.f, "Maximum DCA between the V0 and cascade bachelor"};
  Configurable<float> dcaLambdaDaughters{"dcaLambdaDaughters", 1.f, "Maximum DCA between Lambda daughters"};
  Configurable<float> mXiWindow{"mXiWindow", 0.02f, "Xi mass window used by the cascade compatibility mode"};
  Configurable<float> mOmegaWindow{"mOmegaWindow", 0.01f, "Omega mass window"};
  Configurable<float> mLambdaWindow{"mLambdaWindow", 0.01f, "Lambda mass window"};
  Configurable<float> maxDoubleOmegaMass{"maxDoubleOmegaMass", 3.6f, "Maximum double-Omega invariant mass (GeV/c^2)"};
  Configurable<float> minCosPAOmega{"minCosPAOmega", -1.f, "Minimum Omega cosPA relative to the double-Omega decay vertex"};
  Configurable<float> minCosPADirectLambda{"minCosPADirectLambda", -1.f, "Minimum direct-Lambda cosPA relative to the double-Omega decay vertex"};
  Configurable<float> minCosPADoubleOmega{"minCosPADoubleOmega", -1.f, "Minimum double-Omega cosPA relative to the primary vertex"};
  Configurable<float> minDoubleOmegaDecayRadius{"minDoubleOmegaDecayRadius", 1.f, "Minimum double-Omega transverse decay radius in cm"};

  struct : ConfigurableGroup {
    std::string prefix = "mixingOpts";
    Configurable<int> nMixingCollisions{"nMixingCollisions", 5, "Number of neighbouring collisions in the same bin each collision is mixed with"};
    Configurable<int> mixingTypes{"mixingTypes", 0b1111, "Bitmask of the mixed backgrounds: 1 Omega, 2 Lambda, 4 kaon from another collision, 8 all from different collisions"};
    Configurable<float> maxDeltaPVz{"maxDeltaPVz", 1.f, "Maximum difference of the primary-vertex z among the mixed collisions (cm)"};
    ConfigurableAxis binsVtxZ{"binsVtxZ", {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - primary-vertex z (cm)"};
    ConfigurableAxis binsMultiplicity{"binsMultiplicity", {VARIABLE_WIDTH, 0.f, 20.f, 50.f, 100.f, 200.f, 500.f, 1.e5f}, "Mixing bins - FT0M multiplicity"};
    // Candidate-level mixing (processDataMixingCandidates)
    Configurable<int> nMixingCandidates{"nMixingCandidates", 10, "Number of same-event candidates from other collisions each candidate is mixed with, per (R, phi, charge) bin"};
    Configurable<int> poolDepth{"poolDepth", 50, "Number of same-event candidates kept in the mixing pool per (R, phi, charge) bin, across dataframes"};
    Configurable<int> monitoringInterval{"monitoringInterval", 50, "Log the mixing pool size, memory and timing every this many dataframes (0: never)"};
    ConfigurableAxis binsRadius{"binsRadius", {VARIABLE_WIDTH, 0.f, 0.25f, 0.5f, 1.f, 2.f, 4.f, 10.f, 40.f}, "Mixing bins - double-Omega decay radius (cm)"};
    ConfigurableAxis binsPhi{"binsPhi", {VARIABLE_WIDTH, -o2::constants::math::PI, -2.f * o2::constants::math::PI / 3.f, -o2::constants::math::PI / 3.f, 0.f, o2::constants::math::PI / 3.f, 2.f * o2::constants::math::PI / 3.f, o2::constants::math::PI}, "Mixing bins - double-Omega momentum azimuth (rad)"};
  } mixingOpts;

  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <class Track>
  static uint32_t detectorMap(Track const& track)
  {
    return static_cast<uint32_t>(track.hasITS()) |
           (static_cast<uint32_t>(track.hasTPC()) << 1) |
           (static_cast<uint32_t>(track.hasTOF()) << 2);
  }

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.itsNCls() < minNITSClus ||
        track.tpcNClsFound() < minNTPCClus ||
        track.tpcNClsCrossedRows() < minNCrossedRows ||
        track.tpcNClsCrossedRows() < 0.8f * track.tpcNClsFindable() ||
        track.tpcNClsShared() > maxNSharedTPCClus) {
      LOG(debug) << "Track failed selection: itsNCls=" << static_cast<int>(track.itsNCls()) << ", tpcNClsFound=" << track.tpcNClsFound() << ", tpcNClsCrossedRows=" << static_cast<int>(track.tpcNClsCrossedRows()) << ", tpcNClsFindable=" << static_cast<int>(track.tpcNClsFindable()) << ", tpcNClsShared=" << static_cast<int>(track.tpcNClsShared());
      return false;
    }
    return true;
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
    mBz = o2::base::Propagator::Instance()->getNominalBz();
    fitter.setBz(mBz);
    fitter3Body.setBz(mBz);

    if (static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>("GLO/Param/MatLUT", bc.timestamp()));
      if (!lut) {
        LOG(fatal) << "Could not retrieve material LUT for timestamp " << bc.timestamp();
      }
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
    fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value));
    fitter3Body.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value));

    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDoubleOmega,fOmegaXi");
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  void printConfig() const
  {
    LOG(info) << "doubleOmegaTreeCreator runtime configuration (boolean values: 0=false, 1=true):";
    const auto printValue = [](auto const& config) {
      LOG(info) << config.name << " = " << config.value;
    };
    const auto printAxis = [](ConfigurableAxis const& axis) {
      std::ostringstream values;
      values << '[';
      for (size_t i = 0; i < axis.value.size(); ++i) {
        if (i > 0) {
          values << ", ";
        }
        values << axis.value[i];
      }
      values << ']';
      LOG(info) << axis.name << " = " << values.str();
    };

    printValue(cfgSkimmedProcessing);
    printValue(cfgApplyEventSelection);
    printValue(ccdburl);
    printValue(cfgMaterialCorrection);
    printAxis(zVtxAxis);
    printAxis(massOmegaAxis);
    printAxis(massXiAxis);
    printAxis(massLambdaAxis);
    printAxis(momAxis);
    printValue(zVtxMax);
    printValue(etaMax);
    printValue(cascPtMin);
    printValue(cascPtMax);
    printValue(minNCrossedRows);
    printValue(minNITSClus);
    printValue(minNTPCClus);
    printValue(maxNSharedTPCClus);
    printValue(nSigmaTPCCut);
    printValue(dcaBachToPV);
    printValue(dcaKaonToPV);
    printValue(dcaOmegaToPV);
    printValue(dcaDirectLambdaToPV);
    printValue(dcaV0DauToPV);
    printValue(dcaV0Bach);
    printValue(dcaLambdaDaughters);
    printValue(mXiWindow);
    printValue(mOmegaWindow);
    printValue(mLambdaWindow);
    printValue(maxDoubleOmegaMass);
    printValue(minCosPAOmega);
    printValue(minCosPADirectLambda);
    printValue(minCosPADoubleOmega);
    printValue(minDoubleOmegaDecayRadius);
    printValue(mixingOpts.nMixingCollisions);
    printValue(mixingOpts.mixingTypes);
    printValue(mixingOpts.maxDeltaPVz);
    printAxis(mixingOpts.binsVtxZ);
    printAxis(mixingOpts.binsMultiplicity);
    printValue(mixingOpts.nMixingCandidates);
    printValue(mixingOpts.poolDepth);
    printValue(mixingOpts.monitoringInterval);
    printAxis(mixingOpts.binsRadius);
    printAxis(mixingOpts.binsPhi);
    printValue(doprocessData);
    printValue(doprocessDataFromCascades);
    printValue(doprocessMC);
    printValue(doprocessFindableTracks);
    printValue(doprocessDataMixing);
    printValue(doprocessDataMixingCandidates);
  }

  void init(InitContext const&)
  {
    printConfig();

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMaxDZIni(4.);
    fitter.setMinParamChange(1.e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setUseAbsDCA(true);

    fitter3Body.setPropagateToPCA(true);
    fitter3Body.setMaxR(200.);
    fitter3Body.setMaxDZIni(4.);
    fitter3Body.setMinParamChange(1.e-3);
    fitter3Body.setMinRelChi2Change(0.9);
    fitter3Body.setUseAbsDCA(true);

    zorroSummary.setObject(zorro.getZorroSummary());

    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/massXi", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add<TH2>("QA/massOmega", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + K) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massOmegaAxis});
    histos.add<TH2>("QA/massLambda", ";#it{p}_{T} (GeV/#it{c});#it{M}(p + #pi) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massLambdaAxis});
    histos.add("MC/generatedAndFindable", ";Double-#Omega;Entries", HistType::kTH1F, {{kNFindabilitySteps, -0.5f, static_cast<float>(kNFindabilitySteps) - 0.5f}});
    auto generatedAndFindable = histos.get<TH1>(HIST("MC/generatedAndFindable"));
    generatedAndFindable->GetXaxis()->SetBinLabel(kAllGenerated + 1, "All generated");
    generatedAndFindable->GetXaxis()->SetBinLabel(kOmegaLambdaKaon + 1, "#Omega + #Lambda + K");
    generatedAndFindable->GetXaxis()->SetBinLabel(kTwoLambdasTwoKaons + 1, "2#Lambda + 2K");
    generatedAndFindable->GetXaxis()->SetBinLabel(kFinalState + 1, "2p + 2#pi^{-} + 2K^{-} (and c.c.)");
    generatedAndFindable->GetXaxis()->SetBinLabel(kFindable + 1, "Findable");
    generatedAndFindable->GetXaxis()->SetBinLabel(kFindableSelectedTracks + 1, "Findable, selected tracks");
    generatedAndFindable->GetXaxis()->SetBinLabel(kFindableCascadeAndV0 + 1, "Findable cascade + V0");
    generatedAndFindable->GetXaxis()->SetBinLabel(kFindableCascadeAndV0SelectedTracks + 1, "Findable cascade + V0, selected tracks");

    if (doprocessDataMixingCandidates) {
      histos.add("Mixing/hSeedRadiusPhi", ";Double-#Omega decay radius (cm);Double-#Omega #varphi (rad);Same-event candidates", HistType::kTH2F, {mixingOpts.binsRadius, mixingOpts.binsPhi});
      // Resource monitoring of the candidate-level mixing, one entry per dataframe.
      histos.add("Mixing/hProcessingTime", ";Processing time per dataframe (ms);Dataframes", HistType::kTH1F, {{500, 0.f, 5000.f}});
      histos.add("Mixing/hPoolSeeds", ";Seeds in the mixing pool;Dataframes", HistType::kTH1F, {{500, 0.f, 50000.f}});
      histos.add("Mixing/hPoolMemory", ";Estimated mixing-pool memory (MB);Dataframes", HistType::kTH1F, {{500, 0.f, 50.f}});
      histos.add("Mixing/hProcessingTimeVsPoolSeeds", ";Seeds in the mixing pool;Processing time per dataframe (ms)", HistType::kTH2F, {{100, 0.f, 50000.f}, {100, 0.f, 5000.f}});
    }
    if (doprocessDataMixing || doprocessDataMixingCandidates) {
      histos.add("Mixing/hMixedCollisions", ";;Entries", HistType::kTH1F, {{4, -0.5f, 3.5f}});
      auto mixedCollisions = histos.get<TH1>(HIST("Mixing/hMixedCollisions"));
      mixedCollisions->GetXaxis()->SetBinLabel(1, "Pairs");
      mixedCollisions->GetXaxis()->SetBinLabel(2, "Pairs, #Delta#it{z}_{PV}");
      mixedCollisions->GetXaxis()->SetBinLabel(3, "Triples");
      mixedCollisions->GetXaxis()->SetBinLabel(4, "Triples, #Delta#it{z}_{PV}");
      histos.add("Mixing/hMixedCandidates", ";;Entries", HistType::kTH1F, {{kNMixingTypes - 1, 0.5f, kNMixingTypes - 0.5f}});
      auto mixedCandidates = histos.get<TH1>(HIST("Mixing/hMixedCandidates"));
      mixedCandidates->GetXaxis()->SetBinLabel(kMixedOmega, "#Omega mixed");
      mixedCandidates->GetXaxis()->SetBinLabel(kMixedLambda, "#Lambda mixed");
      mixedCandidates->GetXaxis()->SetBinLabel(kMixedKaon, "K mixed");
      mixedCandidates->GetXaxis()->SetBinLabel(kMixedAll, "All mixed");
    }
  }

  template <class Track>
  bool isSelectedKaon(Track const& track, int8_t charge)
  {
    return track.sign() == charge &&
           selectTrack(track) &&
           std::abs(track.tpcNSigmaKa()) <= nSigmaTPCCut;
  }

  static float invariantMass2Body(std::array<float, 3> const& momentum1, float mass1,
                                  std::array<float, 3> const& momentum2, float mass2)
  {
    const float momentum1Squared = momentum1[0] * momentum1[0] + momentum1[1] * momentum1[1] + momentum1[2] * momentum1[2];
    const float momentum2Squared = momentum2[0] * momentum2[0] + momentum2[1] * momentum2[1] + momentum2[2] * momentum2[2];
    const float energy = std::sqrt(momentum1Squared + mass1 * mass1) + std::sqrt(momentum2Squared + mass2 * mass2);
    const std::array<float, 3> totalMomentum{
      momentum1[0] + momentum2[0], momentum1[1] + momentum2[1], momentum1[2] + momentum2[2]};
    const float massSquared = energy * energy - totalMomentum[0] * totalMomentum[0] - totalMomentum[1] * totalMomentum[1] - totalMomentum[2] * totalMomentum[2];
    return std::sqrt(std::max(0.f, massSquared));
  }

  // Transverse DCA of a track to the collision primary vertex, computed by propagating
  // a copy of the track parameters. Returns false if the propagation fails.
  template <class C, class Track>
  bool dcaXYToPV(C const& collision, Track const& track, float& dcaXY) const
  {
    auto trackParCov = getTrackParCov(track);
    const o2::math_utils::Point3D<float> primaryVertex{collision.posX(), collision.posY(), collision.posZ()};
    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
    std::array<float, 2> dca{};
    if (!o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCov, 2.f, matCorr, &dca)) {
      return false;
    }
    dcaXY = dca[0];
    return true;
  }

  template <class C, class Track>
  bool buildLambda(C const& collision, Track const& posTrack, Track const& negTrack, int8_t charge, BuiltLambdaCandidate& builtLambda)
  {
    if (posTrack.sign() <= 0 || negTrack.sign() >= 0 ||
        !selectTrack(posTrack) || !selectTrack(negTrack)) {
      return false;
    }

    const bool isMatter = charge < 0;
    const float protonNSigma = isMatter ? posTrack.tpcNSigmaPr() : negTrack.tpcNSigmaPr();
    const float pionNSigma = isMatter ? negTrack.tpcNSigmaPi() : posTrack.tpcNSigmaPi();
    if (std::abs(protonNSigma) > nSigmaTPCCut || std::abs(pionNSigma) > nSigmaTPCCut) {
      return false;
    }
    float dcaXYPos = -999.f;
    float dcaXYNeg = -999.f;
    if (!dcaXYToPV(collision, posTrack, dcaXYPos) || !dcaXYToPV(collision, negTrack, dcaXYNeg)) {
      LOG(debug) << "Propagation to the PV failed for Lambda daughters " << posTrack.globalIndex() << ", " << negTrack.globalIndex();
      return false;
    }
    if (std::abs(dcaXYPos) < dcaV0DauToPV || std::abs(dcaXYNeg) < dcaV0DauToPV) {
      LOG(debug) << "Lambda dau don't pass selections: " << "nsigmaPR: " << protonNSigma << ", pionNSigma: " << pionNSigma << ", dcaXYPos: " << dcaXYPos << ", dcaXYNeg: " << dcaXYNeg;
      return false;
    }

    auto posTrackParCov = getTrackParCov(posTrack);
    auto negTrackParCov = getTrackParCov(negTrack);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(posTrackParCov, negTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting Lambda daughter tracks " << posTrack.globalIndex() << ", " << negTrack.globalIndex();
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    std::array<float, 3> posMomentum{};
    std::array<float, 3> negMomentum{};
    fitter.getTrack(0).getPxPyPzGlo(posMomentum);
    fitter.getTrack(1).getPxPyPzGlo(negMomentum);
    const std::array<float, 3> lambdaMomentum{
      posMomentum[0] + negMomentum[0],
      posMomentum[1] + negMomentum[1],
      posMomentum[2] + negMomentum[2]};
    const auto& fittedVertex = fitter.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const float mass = isMatter ? invariantMass2Body(posMomentum, o2::constants::physics::MassProton,
                                                     negMomentum, o2::constants::physics::MassPionCharged)
                                : invariantMass2Body(posMomentum, o2::constants::physics::MassPionCharged,
                                                     negMomentum, o2::constants::physics::MassProton);
    const float dcaDaughters = std::sqrt(std::abs(fitter.getChi2AtPCACandidate()));
    const float eta = etaFromMomentum(lambdaMomentum[0], lambdaMomentum[1], lambdaMomentum[2]);
    if (std::abs(eta) > etaMax ||
        dcaDaughters > dcaLambdaDaughters ||
        std::abs(mass - o2::constants::physics::MassLambda0) > mLambdaWindow) {
      return false;
    }

    builtLambda.candidate = {
      .px = lambdaMomentum[0],
      .py = lambdaMomentum[1],
      .pz = lambdaMomentum[2],
      .mass = mass};
    builtLambda.parentTrack = fitter.createParentTrackParCov(0);
    builtLambda.decayVertex = decayVertex;
    builtLambda.posTrackId = posTrack.globalIndex();
    builtLambda.negTrackId = negTrack.globalIndex();
    // Keep particle-role ordering for matter and antimatter alike.
    builtLambda.daughterDetectorMap = isMatter ? detectorMap(posTrack) | (detectorMap(negTrack) << 3)
                                               : detectorMap(negTrack) | (detectorMap(posTrack) << 3);
    builtLambda.dcaXYProtonToPV = isMatter ? dcaXYPos : dcaXYNeg;
    builtLambda.dcaXYPionToPV = isMatter ? dcaXYNeg : dcaXYPos;
    return true;
  }

  template <class C, class Track>
  bool buildOmega(C const& collision, Track const& posTrack, Track const& negTrack, Track const& bachelor, BuiltOmegaCandidate& builtOmega)
  {
    if (bachelor.sign() == 0 ||
        !selectTrack(bachelor) ||
        std::abs(bachelor.tpcNSigmaKa()) > nSigmaTPCCut) {
      return false;
    }
    float dcaXYBachelor = -999.f;
    if (!dcaXYToPV(collision, bachelor, dcaXYBachelor) || std::abs(dcaXYBachelor) < dcaBachToPV) {
      return false;
    }

    BuiltLambdaCandidate lambda;
    if (!buildLambda(collision, posTrack, negTrack, bachelor.sign(), lambda)) {
      return false;
    }

    auto v0TrackParCov = lambda.parentTrack;
    auto bachelorTrackParCov = getTrackParCov(bachelor);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(v0TrackParCov, bachelorTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting Omega daughter tracks " << posTrack.globalIndex() << ", " << negTrack.globalIndex() << ", " << bachelor.globalIndex();
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    std::array<float, 3> lambdaMomentum{};
    std::array<float, 3> bachelorMomentum{};
    fitter.getTrack(0).getPxPyPzGlo(lambdaMomentum);
    fitter.getTrack(1).getPxPyPzGlo(bachelorMomentum);
    const std::array<float, 3> omegaMomentum{
      lambdaMomentum[0] + bachelorMomentum[0],
      lambdaMomentum[1] + bachelorMomentum[1],
      lambdaMomentum[2] + bachelorMomentum[2]};
    const auto& fittedVertex = fitter.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const float massOmega = invariantMass2Body(lambdaMomentum, o2::constants::physics::MassLambda0,
                                               bachelorMomentum, o2::constants::physics::MassKaonCharged);
    const float massXi = invariantMass2Body(lambdaMomentum, o2::constants::physics::MassLambda0,
                                            bachelorMomentum, o2::constants::physics::MassPionCharged);
    const float dcaDaughters = std::sqrt(std::abs(fitter.getChi2AtPCACandidate()));
    const float eta = etaFromMomentum(omegaMomentum[0], omegaMomentum[1], omegaMomentum[2]);
    if (dcaDaughters > dcaV0Bach ||
        std::abs(eta) > etaMax ||
        std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > mOmegaWindow) {
      return false;
    }

    builtOmega = {
      omegaMomentum[0],
      omegaMomentum[1],
      omegaMomentum[2],
      decayVertex[0],
      decayVertex[1],
      decayVertex[2],
      massOmega,
      massXi,
      static_cast<int8_t>(bachelor.sign()),
      posTrack.globalIndex(),
      negTrack.globalIndex(),
      bachelor.globalIndex(),
      fitter.createParentTrackParCov(),
      detectorMap(bachelor) | (lambda.daughterDetectorMap << 3),
      dcaXYBachelor,
      lambda.dcaXYProtonToPV,
      lambda.dcaXYPionToPV};
    return true;
  }

  static float etaFromMomentum(float px, float py, float pz)
  {
    const float pt = std::hypot(px, py);
    return pt > 0.f ? std::asinh(pz / pt) : 0.f;
  }

  static int chargeIndex(int8_t charge)
  {
    return charge < 0 ? 0 : 1;
  }

  template <class Track>
  static BuiltKaonCandidate buildKaon(Track const& track)
  {
    return {getTrackParCov(track),
            track.globalIndex(),
            static_cast<int8_t>(track.sign()),
            track.hasTPC() ? track.tpcNSigmaKa() : -999.f,
            track.hasTOF() ? track.tofNSigmaKa() : -999.f,
            detectorMap(track)};
  }

  static bool sharesTracks(BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda, BuiltKaonCandidate const& kaon)
  {
    const std::array<int64_t, 6> trackIds{omega.posTrackId, omega.negTrackId, omega.bachelorId, lambda.posTrackId, lambda.negTrackId, kaon.trackId};
    for (std::size_t i = 0; i < trackIds.size(); ++i) {
      for (std::size_t j = i + 1; j < trackIds.size(); ++j) {
        if (trackIds[i] == trackIds[j]) {
          return true;
        }
      }
    }
    return false;
  }

  bool buildDoubleOmega(std::array<float, 3> const& primaryVertex, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& directLambda, BuiltKaonCandidate const& kaon, DoubleOmegaCandidate& cand)
  {
    auto omegaTrackParCov = omega.parentTrack;
    auto lambdaTrackParCov = directLambda.parentTrack;
    auto kaonTrackParCov = kaon.track;
    int nCandidates = 0;
    try {
      nCandidates = fitter3Body.process(omegaTrackParCov, lambdaTrackParCov, kaonTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting the double-Omega candidate";
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }

    auto omegaTrackAtVertex = fitter3Body.getTrack(0);
    auto lambdaTrackAtVertex = fitter3Body.getTrack(1);
    auto kaonTrackAtVertex = fitter3Body.getTrack(2);
    std::array<float, 3> omegaMomentum{};
    std::array<float, 3> lambdaMomentum{};
    std::array<float, 3> kaonMomentum{};
    omegaTrackAtVertex.getPxPyPzGlo(omegaMomentum);
    lambdaTrackAtVertex.getPxPyPzGlo(lambdaMomentum);
    kaonTrackAtVertex.getPxPyPzGlo(kaonMomentum);

    const std::array<float, 3> totalMomentum{
      omegaMomentum[0] + lambdaMomentum[0] + kaonMomentum[0],
      omegaMomentum[1] + lambdaMomentum[1] + kaonMomentum[1],
      omegaMomentum[2] + lambdaMomentum[2] + kaonMomentum[2]};
    const auto& fittedVertex = fitter3Body.getPCACandidate();
    const std::array<float, 3> decayVertex{
      static_cast<float>(fittedVertex[0]), static_cast<float>(fittedVertex[1]), static_cast<float>(fittedVertex[2])};
    const o2::math_utils::Point3D<float> primaryVertexPoint{primaryVertex[0], primaryVertex[1], primaryVertex[2]};
    const std::array<float, 3> omegaDecayVertex{omega.x, omega.y, omega.z};

    std::array<float, 2> dcaOmega{};
    std::array<float, 2> dcaDirectLambda{};
    std::array<float, 2> dcaDirectKaon{};
    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value);
    if (!o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, omegaTrackAtVertex, 2.f, matCorr, &dcaOmega) ||
        !o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, lambdaTrackAtVertex, 2.f, matCorr, &dcaDirectLambda) ||
        !o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertexPoint, kaonTrackAtVertex, 2.f, matCorr, &dcaDirectKaon)) {
      return false;
    }

    const auto momentumSquared = [](std::array<float, 3> const& momentum) {
      return momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2];
    };
    const float energy = std::sqrt(momentumSquared(omegaMomentum) + o2::constants::physics::MassOmegaMinus * o2::constants::physics::MassOmegaMinus) +
                         std::sqrt(momentumSquared(lambdaMomentum) + o2::constants::physics::MassLambda0 * o2::constants::physics::MassLambda0) +
                         std::sqrt(momentumSquared(kaonMomentum) + o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged);
    const float massSquared = energy * energy - totalMomentum[0] * totalMomentum[0] - totalMomentum[1] * totalMomentum[1] - totalMomentum[2] * totalMomentum[2];

    cand.pt = std::hypot(totalMomentum[0], totalMomentum[1]);
    cand.eta = etaFromMomentum(totalMomentum[0], totalMomentum[1], totalMomentum[2]);
    cand.phi = std::atan2(totalMomentum[1], totalMomentum[0]);
    cand.x = decayVertex[0];
    cand.y = decayVertex[1];
    cand.z = decayVertex[2];
    cand.pvX = primaryVertex[0];
    cand.pvY = primaryVertex[1];
    cand.pvZ = primaryVertex[2];
    cand.cosPAOmega = RecoDecay::cpa(decayVertex, omegaDecayVertex, omegaMomentum);
    cand.cosPADirectLambda = RecoDecay::cpa(decayVertex, directLambda.decayVertex, lambdaMomentum);
    cand.cosPADoubleOmega = RecoDecay::cpa(primaryVertex, decayVertex, totalMomentum);
    // PV pointing angles use the daughter momenta at their own decay vertices.
    cand.cosPAOmegaToPV = RecoDecay::cpa(primaryVertex, omegaDecayVertex, std::array<float, 3>{omega.px, omega.py, omega.pz});
    cand.cosPADirectLambdaToPV = RecoDecay::cpa(primaryVertex, directLambda.decayVertex, std::array<float, 3>{directLambda.candidate.px, directLambda.candidate.py, directLambda.candidate.pz});
    cand.dcaXYOmegaToPV = dcaOmega[0];
    cand.dcaXYDirectLambdaToPV = dcaDirectLambda[0];
    cand.dcaXYDirectKaonToPV = dcaDirectKaon[0];
    cand.dcaXYOmegaKaonToPV = omega.dcaXYKaonToPV;
    cand.dcaXYOmegaProtonToPV = omega.dcaXYProtonToPV;
    cand.dcaXYOmegaPionToPV = omega.dcaXYPionToPV;
    cand.dcaXYDirectProtonToPV = directLambda.dcaXYProtonToPV;
    cand.dcaXYDirectPionToPV = directLambda.dcaXYPionToPV;
    cand.tpcNSigmaDirectKaon = kaon.tpcNSigma;
    cand.tofNSigmaDirectKaon = kaon.tofNSigma;
    cand.daughterDetectorMap = omega.daughterDetectorMap |
                               (directLambda.daughterDetectorMap << 9) |
                               (kaon.detectorMap << 15);
    cand.mass = std::sqrt(std::max(0.f, massSquared));
    if (cand.cosPAOmega < minCosPAOmega ||
        cand.cosPADirectLambda < minCosPADirectLambda ||
        cand.cosPADoubleOmega < minCosPADoubleOmega ||
        std::hypot(decayVertex[0], decayVertex[1]) < minDoubleOmegaDecayRadius ||
        std::abs(cand.dcaXYOmegaToPV) < dcaOmegaToPV ||
        std::abs(cand.dcaXYDirectLambdaToPV) < dcaDirectLambdaToPV ||
        std::abs(cand.dcaXYDirectKaonToPV) < dcaKaonToPV ||
        cand.mass > maxDoubleOmegaMass) {
      return false;
    }
    cand.massOmega = omega.massOmega;
    cand.massXi = omega.massXi;
    cand.sign = omega.sign;
    return true;
  }

  template <class Track>
  int64_t getLambdaMCLabel(Track const& posTrack, Track const& negTrack, int8_t sign)
  {
    if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
      return -1;
    }

    auto mcPosTrack = posTrack.template mcParticle_as<aod::McParticles>();
    auto mcNegTrack = negTrack.template mcParticle_as<aod::McParticles>();
    const int expectedPosPdg = sign < 0 ? PDG_t::kProton : PDG_t::kPiPlus;
    const int expectedNegPdg = sign < 0 ? PDG_t::kPiMinus : PDG_t::kProtonBar;
    const int expectedLambdaPdg = sign < 0 ? PDG_t::kLambda0 : PDG_t::kLambda0Bar;
    if (mcPosTrack.pdgCode() != expectedPosPdg || mcNegTrack.pdgCode() != expectedNegPdg) {
      return -1;
    }

    for (const auto& posMother : mcPosTrack.template mothers_as<aod::McParticles>()) {
      if (posMother.pdgCode() != expectedLambdaPdg) {
        continue;
      }
      for (const auto& negMother : mcNegTrack.template mothers_as<aod::McParticles>()) {
        if (posMother.globalIndex() == negMother.globalIndex()) {
          return posMother.globalIndex();
        }
      }
    }
    return -1;
  }

  template <class Track, class McParticles>
  int64_t getOmegaMCLabel(Track const& posTrack, Track const& negTrack, Track const& bachelorTrack,
                          int8_t sign, McParticles const& mcParticles)
  {
    const int64_t lambdaLabel = getLambdaMCLabel(posTrack, negTrack, sign);
    if (lambdaLabel < 0 || !bachelorTrack.has_mcParticle()) {
      return -1;
    }

    auto mcLambda = mcParticles.rawIteratorAt(lambdaLabel);
    auto mcBachelor = bachelorTrack.template mcParticle_as<aod::McParticles>();
    const int expectedBachelorPdg = sign < 0 ? PDG_t::kKMinus : PDG_t::kKPlus;
    const int expectedOmegaPdg = sign < 0 ? PDG_t::kOmegaMinus : PDG_t::kOmegaPlusBar;
    if (mcBachelor.pdgCode() != expectedBachelorPdg) {
      return -1;
    }

    for (const auto& lambdaMother : mcLambda.template mothers_as<aod::McParticles>()) {
      if (lambdaMother.pdgCode() != expectedOmegaPdg) {
        continue;
      }
      for (const auto& bachelorMother : mcBachelor.template mothers_as<aod::McParticles>()) {
        if (lambdaMother.globalIndex() == bachelorMother.globalIndex()) {
          return lambdaMother.globalIndex();
        }
      }
    }
    return -1;
  }

  template <class Track, class McParticles>
  bool getMCInfo(DoubleOmegaMCInfo& mcInfo,
                 Track const& omegaPosTrack, Track const& omegaNegTrack, Track const& omegaBachelorTrack,
                 Track const& lambdaPosTrack, Track const& lambdaNegTrack, Track const& directKaonTrack,
                 int8_t sign, McParticles const& mcParticles)
  {
    const int64_t omegaLabel = getOmegaMCLabel(omegaPosTrack, omegaNegTrack, omegaBachelorTrack, sign, mcParticles);
    const int64_t lambdaLabel = getLambdaMCLabel(lambdaPosTrack, lambdaNegTrack, sign);
    if (omegaLabel < 0 || lambdaLabel < 0 || !directKaonTrack.has_mcParticle()) {
      return false;
    }

    auto mcOmega = mcParticles.rawIteratorAt(omegaLabel);
    auto mcLambda = mcParticles.rawIteratorAt(lambdaLabel);
    auto mcKaon = directKaonTrack.template mcParticle_as<aod::McParticles>();
    const int expectedKaonPdg = sign < 0 ? PDG_t::kKMinus : PDG_t::kKPlus;
    const int expectedDoubleOmegaPdg = sign < 0 ? kDoubleOmegaPdg : -kDoubleOmegaPdg;
    if (mcKaon.pdgCode() != expectedKaonPdg) {
      return false;
    }

    for (const auto& omegaMother : mcOmega.template mothers_as<aod::McParticles>()) {
      if (omegaMother.pdgCode() != expectedDoubleOmegaPdg) {
        continue;
      }
      for (const auto& lambdaMother : mcLambda.template mothers_as<aod::McParticles>()) {
        if (omegaMother.globalIndex() != lambdaMother.globalIndex()) {
          continue;
        }
        for (const auto& kaonMother : mcKaon.template mothers_as<aod::McParticles>()) {
          if (omegaMother.globalIndex() == kaonMother.globalIndex()) {
            mcInfo.motherId = omegaMother.globalIndex();
            mcInfo.pt = omegaMother.pt();
            mcInfo.eta = omegaMother.eta();
            mcInfo.phi = omegaMother.phi();
            mcInfo.decayLength = std::hypot(mcOmega.vx() - omegaMother.vx(),
                                            mcOmega.vy() - omegaMother.vy(),
                                            mcOmega.vz() - omegaMother.vz());
            mcInfo.pdgCode = omegaMother.pdgCode();
            mcInfo.sign = sign;
            return true;
          }
        }
      }
    }
    return false;
  }

  template <class McParticle>
  bool getGeneratedMCInfo(DoubleOmegaMCInfo& info, McParticle const& particle)
  {
    if (std::abs(particle.pdgCode()) != kDoubleOmegaPdg) {
      return false;
    }

    const bool isMatter = particle.pdgCode() > 0;
    const int expectedOmegaPdg = isMatter ? PDG_t::kOmegaMinus : PDG_t::kOmegaPlusBar;
    const int expectedLambdaPdg = isMatter ? PDG_t::kLambda0 : PDG_t::kLambda0Bar;
    const int expectedKaonPdg = isMatter ? PDG_t::kKMinus : PDG_t::kKPlus;
    bool foundOmega = false;
    bool foundLambda = false;
    bool foundKaon = false;
    std::array<float, 3> decayVertex{0.f, 0.f, 0.f};
    std::array<float, 3> omegaDecayVertex{-999.f, -999.f, -999.f};

    for (const auto& daughter : particle.template daughters_as<aod::McParticles>()) {
      if (daughter.pdgCode() == expectedOmegaPdg) {
        foundOmega = true;
        decayVertex = {daughter.vx(), daughter.vy(), daughter.vz()};
        // The Omega decay point is the production point of any of its daughters.
        for (const auto& omegaDaughter : daughter.template daughters_as<aod::McParticles>()) {
          omegaDecayVertex = {omegaDaughter.vx(), omegaDaughter.vy(), omegaDaughter.vz()};
          break;
        }
      } else if (daughter.pdgCode() == expectedLambdaPdg) {
        foundLambda = true;
      } else if (daughter.pdgCode() == expectedKaonPdg) {
        foundKaon = true;
      }
    }
    if (!foundOmega || !foundLambda || !foundKaon) {
      return false;
    }

    info.motherId = particle.globalIndex();
    info.pt = particle.pt();
    info.eta = particle.eta();
    info.phi = particle.phi();
    info.decayLength = std::hypot(decayVertex[0] - particle.vx(),
                                  decayVertex[1] - particle.vy(),
                                  decayVertex[2] - particle.vz());
    info.pdgCode = particle.pdgCode();
    info.genProductionVertex = {particle.vx(), particle.vy(), particle.vz()};
    info.genDecayVertex = decayVertex;
    info.genOmegaDecayVertex = omegaDecayVertex;
    info.sign = isMatter ? -1 : 1;
    return true;
  }

  template <class McParticle>
  bool getLambdaFinalStateIds(McParticle const& lambda, int sign, std::array<int64_t, 2>& finalStateIds)
  {
    bool foundProton = false;
    bool foundPion = false;
    for (const auto& daughter : lambda.template daughters_as<aod::McParticles>()) {
      if (daughter.pdgCode() == sign * 2212 && !foundProton) {
        finalStateIds[0] = daughter.globalIndex();
        foundProton = true;
      } else if (daughter.pdgCode() == -sign * 211 && !foundPion) {
        finalStateIds[1] = daughter.globalIndex();
        foundPion = true;
      }
    }
    return foundProton && foundPion;
  }

  void writeDataCandidate(DoubleOmegaCandidate const& cand)
  {
    doubleOmegaTable(cand.pt * cand.sign,
                     cand.eta,
                     cand.phi,
                     cand.x,
                     cand.y,
                     cand.z,
                     cand.pvX,
                     cand.pvY,
                     cand.pvZ,
                     cand.cosPAOmega,
                     cand.cosPADirectLambda,
                     cand.cosPADoubleOmega,
                     cand.cosPAOmegaToPV,
                     cand.cosPADirectLambdaToPV,
                     cand.dcaXYOmegaToPV,
                     cand.dcaXYDirectLambdaToPV,
                     cand.dcaXYDirectKaonToPV,
                     cand.dcaXYOmegaKaonToPV,
                     cand.dcaXYOmegaProtonToPV,
                     cand.dcaXYOmegaPionToPV,
                     cand.dcaXYDirectProtonToPV,
                     cand.dcaXYDirectPionToPV,
                     cand.tpcNSigmaDirectKaon,
                     cand.tofNSigmaDirectKaon,
                     cand.daughterDetectorMap,
                     cand.mass,
                     cand.massOmega,
                     cand.massXi);
  }

  void writeMixedCandidate(DoubleOmegaCandidate const& cand, MixingType mixingType)
  {
    doubleOmegaMixedTable(cand.pt * cand.sign,
                          cand.eta,
                          cand.phi,
                          cand.x,
                          cand.y,
                          cand.z,
                          cand.pvX,
                          cand.pvY,
                          cand.pvZ,
                          cand.cosPAOmega,
                          cand.cosPADirectLambda,
                          cand.cosPADoubleOmega,
                          cand.cosPAOmegaToPV,
                          cand.cosPADirectLambdaToPV,
                          cand.dcaXYOmegaToPV,
                          cand.dcaXYDirectLambdaToPV,
                          cand.dcaXYDirectKaonToPV,
                          cand.dcaXYOmegaKaonToPV,
                          cand.dcaXYOmegaProtonToPV,
                          cand.dcaXYOmegaPionToPV,
                          cand.dcaXYDirectProtonToPV,
                          cand.dcaXYDirectPionToPV,
                          cand.tpcNSigmaDirectKaon,
                          cand.tofNSigmaDirectKaon,
                          cand.daughterDetectorMap,
                          cand.mass,
                          cand.massOmega,
                          cand.massXi,
                          mixingType);
  }

  void writeMCCandidate(DoubleOmegaCandidate const& cand, DoubleOmegaMCInfo const& mcInfo, bool isReco)
  {
    doubleOmegaTableMC(cand.pt * cand.sign,
                       cand.eta,
                       cand.phi,
                       cand.x,
                       cand.y,
                       cand.z,
                       cand.pvX,
                       cand.pvY,
                       cand.pvZ,
                       cand.cosPAOmega,
                       cand.cosPADirectLambda,
                       cand.cosPADoubleOmega,
                       cand.cosPAOmegaToPV,
                       cand.cosPADirectLambdaToPV,
                       cand.dcaXYOmegaToPV,
                       cand.dcaXYDirectLambdaToPV,
                       cand.dcaXYDirectKaonToPV,
                       cand.dcaXYOmegaKaonToPV,
                       cand.dcaXYOmegaProtonToPV,
                       cand.dcaXYOmegaPionToPV,
                       cand.dcaXYDirectProtonToPV,
                       cand.dcaXYDirectPionToPV,
                       cand.tpcNSigmaDirectKaon,
                       cand.tofNSigmaDirectKaon,
                       cand.daughterDetectorMap,
                       cand.mass,
                       cand.massOmega,
                       cand.massXi,
                       mcInfo.pt * mcInfo.sign,
                       mcInfo.eta,
                       mcInfo.phi,
                       mcInfo.decayLength,
                       mcInfo.pdgCode,
                       mcInfo.motherId,
                       isReco);
  }

  template <class C>
  void logVertexComparison(C const& collision, DoubleOmegaCandidate const& cand,
                           BuiltOmegaCandidate const& omega, DoubleOmegaMCInfo const& mcInfo)
  {
    auto logVertex = [](const char* label, std::array<float, 3> const& reco, std::array<float, 3> const& gen) {
      LOG(debug) << label
                 << ": reco = (" << reco[0] << ", " << reco[1] << ", " << reco[2] << ") cm"
                 << " | gen = (" << gen[0] << ", " << gen[1] << ", " << gen[2] << ") cm"
                 << " | reco - gen = (" << reco[0] - gen[0] << ", " << reco[1] - gen[1] << ", " << reco[2] - gen[2] << ") cm";
    };
    LOG(debug) << "---------------------------------------- vertices, mother " << mcInfo.motherId
               << ", pdgCode = " << mcInfo.pdgCode << ", gen pt = " << mcInfo.pt << " GeV/c";
    logVertex("PV          ", {cand.pvX, cand.pvY, cand.pvZ}, mcInfo.genProductionVertex);
    logVertex("doubleOmega ", {cand.x, cand.y, cand.z}, mcInfo.genDecayVertex);
    logVertex("Omega       ", {omega.x, omega.y, omega.z}, mcInfo.genOmegaDecayVertex);
    LOG(debug) << "decay length: reco = "
               << std::hypot(cand.x - cand.pvX, cand.y - cand.pvY, cand.z - cand.pvZ)
               << " cm | gen = " << mcInfo.decayLength << " cm";
    LOG(debug) << "collision " << collision.globalIndex() << " (direct kaon), numContrib = " << collision.numContrib();
  }

  // Builds the Omega, Lambda and kaon candidates of a collision. Lambdas and kaons are built only if the
  // collision has at least one Omega, unless they are needed for the event mixing.
  template <class C, class T, class Cascades, class V0s>
  CollisionCandidates buildCollisionCandidates(C const& collision, T const& tracks, Cascades const& cascades, V0s const& v0s, bool requireOmega)
  {
    CollisionCandidates candidates;
    candidates.primaryVertex = {collision.posX(), collision.posY(), collision.posZ()};
    for (const auto& rawCascade : cascades) {
      const auto omegaV0 = rawCascade.template v0_as<aod::V0s>();
      BuiltOmegaCandidate omega;
      if (!buildOmega(collision, omegaV0.template posTrack_as<T>(), omegaV0.template negTrack_as<T>(), rawCascade.template bachelor_as<T>(), omega)) {
        continue;
      }
      const float omegaPt = std::hypot(omega.px, omega.py);
      histos.fill(HIST("QA/massXi"), omegaPt, omega.massXi);
      histos.fill(HIST("QA/massOmega"), omegaPt, omega.massOmega);
      candidates.omegas.push_back(omega);
    }
    if (requireOmega && candidates.omegas.empty()) {
      return candidates;
    }

    for (const int8_t charge : {int8_t{-1}, int8_t{1}}) {
      for (const auto& rawV0 : v0s) {
        BuiltLambdaCandidate lambda;
        if (!buildLambda(collision, rawV0.template posTrack_as<T>(), rawV0.template negTrack_as<T>(), charge, lambda)) {
          continue;
        }
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);
        candidates.lambdas[chargeIndex(charge)].push_back(lambda);
      }
      for (const auto& track : tracks) {
        if (isSelectedKaon(track, charge)) {
          candidates.kaons[chargeIndex(charge)].push_back(buildKaon(track));
        }
      }
    }
    return candidates;
  }

  // Pairs the Omegas, Lambdas and kaons of (possibly different) collisions. The double-Omega candidate
  // is built with respect to the primary vertex of the reference collision.
  template <class F>
  void combine(CollisionCandidates const& reference, CollisionCandidates const& omegaSource,
               CollisionCandidates const& lambdaSource, CollisionCandidates const& kaonSource, F&& onCandidate)
  {
    for (const auto& omega : omegaSource.omegas) {
      for (const auto& lambda : lambdaSource.lambdas[chargeIndex(omega.sign)]) {
        for (const auto& kaon : kaonSource.kaons[chargeIndex(omega.sign)]) {
          // Tracks can also be shared by V0s and cascades associated to different collisions.
          if (sharesTracks(omega, lambda, kaon)) {
            continue;
          }
          DoubleOmegaCandidate cand;
          if (!buildDoubleOmega(reference.primaryVertex, omega, lambda, kaon, cand)) {
            continue;
          }
          onCandidate(cand, omega, lambda, kaon);
        }
      }
    }
  }

  template <class C, class T, class Cascades>
  void fillFromCascades(C const& collision, T const&, Cascades const& cascades)
  {
    for (const auto& omegaRow : cascades) {
      BuiltOmegaCandidate omega;
      if (!buildOmega(collision, omegaRow.template posTrack_as<T>(), omegaRow.template negTrack_as<T>(), omegaRow.template bachelor_as<T>(), omega)) {
        continue;
      }
      const float omegaPt = std::hypot(omega.px, omega.py);
      histos.fill(HIST("QA/massXi"), omegaPt, omega.massXi);
      histos.fill(HIST("QA/massOmega"), omegaPt, omega.massOmega);
      const std::array<int64_t, 3> omegaTrackIds{omega.posTrackId, omega.negTrackId, omega.bachelorId};

      for (const auto& lambdaKaonSource : cascades) {
        if (lambdaKaonSource.globalIndex() == omegaRow.globalIndex() ||
            lambdaKaonSource.sign() != omega.sign) {
          continue;
        }

        const std::array<int64_t, 3> sourceTrackIds{
          lambdaKaonSource.posTrackId(), lambdaKaonSource.negTrackId(), lambdaKaonSource.bachelorId()};
        bool sharesTrack = false;
        for (const auto& omegaTrackId : omegaTrackIds) {
          if (std::find(sourceTrackIds.begin(), sourceTrackIds.end(), omegaTrackId) != sourceTrackIds.end()) {
            sharesTrack = true;
            break;
          }
        }
        if (sharesTrack) {
          continue;
        }

        auto kaon = lambdaKaonSource.template bachelor_as<T>();
        if (!isSelectedKaon(kaon, omega.sign)) {
          continue;
        }
        BuiltLambdaCandidate lambda;
        if (!buildLambda(collision, lambdaKaonSource.template posTrack_as<T>(), lambdaKaonSource.template negTrack_as<T>(), omega.sign, lambda)) {
          continue;
        }
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);
        DoubleOmegaCandidate cand;
        if (buildDoubleOmega({collision.posX(), collision.posY(), collision.posZ()}, omega, lambda, buildKaon(kaon), cand)) {
          writeDataCandidate(cand);
        }
      }
    }
  }

  // Always take the vertex the direct kaon is associated with. To be improved
  template <class Kaon>
  int64_t selectFindableCollision(Kaon const& kaon)
  {
    return kaon.collisionId();
  }

  template <class C, class T>
  bool fillFromSelectedTracks(C const& collisions, T const& tracks,
                              std::array<int, 6> const& daughterTrackIds,
                              DoubleOmegaMCInfo const& mcInfo)
  {
    // Daughter slots: Omega kaon, Omega proton/pion, direct proton/pion, direct kaon.
    // Swap proton and pion slots for antimatter to pass positive/negative tracks.
    const bool isMatter = mcInfo.pdgCode > 0;
    const int8_t charge = isMatter ? -1 : 1;
    const auto bachelor = tracks.rawIteratorAt(daughterTrackIds[0]);
    if (bachelor.sign() != charge) {
      return false;
    }
    const auto omegaPosTrack = tracks.rawIteratorAt(daughterTrackIds[isMatter ? 1 : 2]);
    const auto omegaNegTrack = tracks.rawIteratorAt(daughterTrackIds[isMatter ? 2 : 1]);
    const auto lambdaPosTrack = tracks.rawIteratorAt(daughterTrackIds[isMatter ? 3 : 4]);
    const auto lambdaNegTrack = tracks.rawIteratorAt(daughterTrackIds[isMatter ? 4 : 3]);
    const auto kaon = tracks.rawIteratorAt(daughterTrackIds[5]);
    if (!isSelectedKaon(kaon, charge)) {
      return false;
    }
    const int64_t collisionId = selectFindableCollision(kaon);
    if (collisionId < 0) {
      return false;
    }
    const auto collision = collisions.rawIteratorAt(collisionId);
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    BuiltOmegaCandidate omega;
    if (!buildOmega(collision, omegaPosTrack, omegaNegTrack, bachelor, omega)) {
      return false;
    }
    const float omegaPt = std::hypot(omega.px, omega.py);
    histos.fill(HIST("QA/massXi"), omegaPt, omega.massXi);
    histos.fill(HIST("QA/massOmega"), omegaPt, omega.massOmega);
    BuiltLambdaCandidate lambda;
    if (!buildLambda(collision, lambdaPosTrack, lambdaNegTrack, charge, lambda)) {
      return false;
    }
    histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);
    DoubleOmegaCandidate cand;
    if (!buildDoubleOmega({collision.posX(), collision.posY(), collision.posZ()}, omega, lambda, buildKaon(kaon), cand)) {
      return false;
    }
    logVertexComparison(collision, cand, omega, mcInfo);
    writeMCCandidate(cand, mcInfo, true);
    return true;
  }

  template <class C>
  bool acceptCollision(C const& collision, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (cfgApplyEventSelection &&
        (!collision.sel8() ||
         std::abs(collision.posZ()) > zVtxMax ||
         !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
         !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
      return false;
    }
    if (cfgSkimmedProcessing) {
      zorro.isSelected(bc.globalBC());
    }
    histos.fill(HIST("QA/zVtx"), collision.posZ());
    return true;
  }

  void mix(CollisionCandidates const& reference, CollisionCandidates const& omegaSource,
           CollisionCandidates const& lambdaSource, CollisionCandidates const& kaonSource, MixingType mixingType)
  {
    if (!(mixingOpts.mixingTypes & (1 << (mixingType - 1)))) {
      return;
    }
    combine(reference, omegaSource, lambdaSource, kaonSource, [&](DoubleOmegaCandidate const& cand, auto const&...) {
      histos.fill(HIST("Mixing/hMixedCandidates"), mixingType);
      writeMixedCandidate(cand, mixingType);
    });
  }

  // Two candidates from one collision and one from the other, in both directions. The collision
  // providing two candidates is the reference one.
  void mixPair(CollisionCandidates const& candidates0, CollisionCandidates const& candidates1)
  {
    mix(candidates0, candidates1, candidates0, candidates0, kMixedOmega);
    mix(candidates1, candidates0, candidates1, candidates1, kMixedOmega);
    mix(candidates0, candidates0, candidates1, candidates0, kMixedLambda);
    mix(candidates1, candidates1, candidates0, candidates1, kMixedLambda);
    mix(candidates0, candidates0, candidates0, candidates1, kMixedKaon);
    mix(candidates1, candidates1, candidates1, candidates0, kMixedKaon);
  }

  // All the assignments of the Omega, Lambda and kaon to three different collisions. The collision
  // providing the Omega is the reference one.
  void mixTriple(CollisionCandidates const& candidates0, CollisionCandidates const& candidates1, CollisionCandidates const& candidates2)
  {
    const std::array<const CollisionCandidates*, 3> sources{&candidates0, &candidates1, &candidates2};
    std::array<int, 3> roles{0, 1, 2}; // Sources of the Omega, Lambda and kaon.
    do {
      const auto& omegaSource = *sources[roles[0]];
      mix(omegaSource, omegaSource, *sources[roles[1]], *sources[roles[2]], kMixedAll);
    } while (std::next_permutation(roles.begin(), roles.end()));
  }

  void processData(Collisions const& collision,
                   TracksFull const& tracks,
                   aod::V0s const& v0s,
                   aod::Cascades const& cascades,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (!acceptCollision(collision, bcs)) {
      return;
    }
    const auto candidates = buildCollisionCandidates(collision, tracks, cascades, v0s, true);
    combine(candidates, candidates, candidates, candidates,
            [&](DoubleOmegaCandidate const& cand, auto const&...) { writeDataCandidate(cand); });
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processData, "Reconstruct Omega + Lambda + kaon from raw V0/cascade indices", true);

  void processDataFromCascades(Collisions const& collision,
                               TracksFull const& tracks,
                               FullCascades const& cascades,
                               aod::BCsWithTimestamps const& bcs)
  {
    if (acceptCollision(collision, bcs)) {
      fillFromCascades(collision, tracks, cascades);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processDataFromCascades, "Reconstruct Omega + Lambda + kaon from two cascade rows", false);

  void processMC(CollisionsMC const& collisions,
                 TracksFullMC const& tracks,
                 aod::V0s const& v0s,
                 aod::Cascades const& cascades,
                 aod::McParticles const& mcParticles,
                 aod::BCsWithTimestamps const& bcs)
  {
    reconstructedDoubleOmegaIds.clear();

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!acceptCollision(collision, bcs)) {
        continue;
      }
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(rawCascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(rawV0sPerCollision, collision.globalIndex());
      v0sThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&v0s);
      const auto candidates = buildCollisionCandidates(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, true);
      combine(candidates, candidates, candidates, candidates,
              [&](DoubleOmegaCandidate const& cand, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda, BuiltKaonCandidate const& kaon) {
                DoubleOmegaMCInfo mcInfo;
                if (getMCInfo(mcInfo,
                              tracks.rawIteratorAt(omega.posTrackId), tracks.rawIteratorAt(omega.negTrackId), tracks.rawIteratorAt(omega.bachelorId),
                              tracks.rawIteratorAt(lambda.posTrackId), tracks.rawIteratorAt(lambda.negTrackId), tracks.rawIteratorAt(kaon.trackId),
                              omega.sign, mcParticles)) {
                  reconstructedDoubleOmegaIds.push_back(mcInfo.motherId);
                  writeMCCandidate(cand, mcInfo, true);
                }
              });
    }
    for (const auto& mcParticle : mcParticles) {
      DoubleOmegaMCInfo mcInfo;
      if (!getGeneratedMCInfo(mcInfo, mcParticle)) {
        continue;
      }

      if (std::find(reconstructedDoubleOmegaIds.begin(), reconstructedDoubleOmegaIds.end(), mcInfo.motherId) != reconstructedDoubleOmegaIds.end()) {
        continue;
      }
      writeMCCandidate(DoubleOmegaCandidate{}, mcInfo, false);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processMC, "Reconstruct and truth-match Omega + Lambda + kaon in MC", false);

  void processFindableTracks(soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
                             TracksFullMC const& tracks,
                             aod::Cascades const& cascades,
                             aod::V0s const& v0s,
                             aod::McParticles const& mcParticles,
                             aod::McCollisions const&,
                             aod::BCsWithTimestamps const&)
  {
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.pdgCode()) != kDoubleOmegaPdg) {
        continue;
      }
      histos.fill(HIST("MC/generatedAndFindable"), kAllGenerated);

      const int sign = mcParticle.pdgCode() > 0 ? 1 : -1;
      int64_t omegaId = -1;
      int64_t directLambdaId = -1;
      int64_t directKaonId = -1;

      // First stage: double-Omega -> Omega + Lambda + kaon.
      for (const auto& daughter : mcParticle.daughters_as<aod::McParticles>()) {
        if (daughter.pdgCode() == sign * 3334 && omegaId < 0) {
          omegaId = daughter.globalIndex();
        } else if (daughter.pdgCode() == sign * 3122 && directLambdaId < 0) {
          directLambdaId = daughter.globalIndex();
        } else if (daughter.pdgCode() == -sign * 321 && directKaonId < 0) {
          directKaonId = daughter.globalIndex();
        }
      }
      if (omegaId < 0 || directLambdaId < 0 || directKaonId < 0) {
        continue;
      }
      histos.fill(HIST("MC/generatedAndFindable"), kOmegaLambdaKaon);

      // Second stage: Omega -> Lambda + kaon, giving two Lambdas and two kaons.
      const auto omega = mcParticles.rawIteratorAt(omegaId);
      int64_t omegaLambdaId = -1;
      int64_t omegaKaonId = -1;
      for (const auto& daughter : omega.daughters_as<aod::McParticles>()) {
        if (daughter.pdgCode() == sign * 3122 && omegaLambdaId < 0) {
          omegaLambdaId = daughter.globalIndex();
        } else if (daughter.pdgCode() == -sign * 321 && omegaKaonId < 0) {
          omegaKaonId = daughter.globalIndex();
        }
      }
      if (omegaLambdaId < 0 || omegaKaonId < 0) {
        continue;
      }
      histos.fill(HIST("MC/generatedAndFindable"), kTwoLambdasTwoKaons);

      // Third stage: both Lambdas -> proton + pion.
      const auto omegaLambda = mcParticles.rawIteratorAt(omegaLambdaId);
      const auto directLambda = mcParticles.rawIteratorAt(directLambdaId);
      std::array<int64_t, 2> omegaLambdaFinalStateIds{};
      std::array<int64_t, 2> directLambdaFinalStateIds{};
      if (!getLambdaFinalStateIds(omegaLambda, sign, omegaLambdaFinalStateIds)) {
        continue;
      }
      if (!getLambdaFinalStateIds(directLambda, sign, directLambdaFinalStateIds)) {
        continue;
      }
      histos.fill(HIST("MC/generatedAndFindable"), kFinalState);

      const std::array<int64_t, 6> finalStateIds{
        omegaKaonId,
        omegaLambdaFinalStateIds[0],
        omegaLambdaFinalStateIds[1],
        directLambdaFinalStateIds[0],
        directLambdaFinalStateIds[1],
        directKaonId};
      std::array<bool, 6> reconstructedDaughters{};
      std::array<bool, 6> selectedReconstructedDaughters{};
      std::array<int, 6> daughterTrackIds{};
      daughterTrackIds.fill(-1);
      // Match by generated daughter, independently of reconstructed collision association.
      std::array<int, 6> selectedTrackIds{};
      selectedTrackIds.fill(-1);
      for (const auto& track : tracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto daughterId = std::find(finalStateIds.begin(), finalStateIds.end(), track.mcParticleId());
        if (daughterId != finalStateIds.end()) {
          const auto daughterIndex = std::distance(finalStateIds.begin(), daughterId);
          reconstructedDaughters[daughterIndex] = true;
          if (daughterTrackIds[daughterIndex] < 0) {
            daughterTrackIds[daughterIndex] = track.globalIndex();
          }
          if (selectTrack(track)) {
            selectedReconstructedDaughters[daughterIndex] = true;
            if (selectedTrackIds[daughterIndex] < 0) {
              selectedTrackIds[daughterIndex] = track.globalIndex();
            }
          }
        }
      }
      auto printDaughterTracks = [&](const char* tag) {
        LOG(info) << "---------------------------------------- " << tag;
        for (std::size_t iDaughter = 0; iDaughter < daughterTrackIds.size(); ++iDaughter) {
          LOG(info) << "+++++++";
          const auto daughterTrack = tracks.rawIteratorAt(daughterTrackIds[iDaughter]);
          const auto daughterMCParticle = daughterTrack.mcParticle_as<aod::McParticles>();
          LOG(info) << "Dau" << iDaughter + 1
                    << ": hasITS = " << daughterTrack.hasITS()
                    << ", hasTPC = " << daughterTrack.hasTPC()
                    << ", hasTOF = " << daughterTrack.hasTOF()
                    << ", hasTRD = " << daughterTrack.hasTRD()
                    << ", isITSAfterburner = " << daughterTrack.isITSAfterburner();
          LOG(info) << "Dau" << iDaughter + 1
                    << ": eta = " << daughterTrack.eta()
                    << ", pt = " << daughterTrack.pt()
                    << ", CollisionID = " << daughterTrack.collisionId()
                    << ", production point = (" << daughterMCParticle.vx()
                    << ", " << daughterMCParticle.vy()
                    << ", " << daughterMCParticle.vz() << ") cm";
        }
      };

      bool allDaughtersReconstructed = true;
      bool allDaughtersSelected = true;
      for (std::size_t iDaughter = 0; iDaughter < reconstructedDaughters.size(); ++iDaughter) {
        allDaughtersReconstructed &= reconstructedDaughters[iDaughter];
        allDaughtersSelected &= selectedReconstructedDaughters[iDaughter];
      }
      if (allDaughtersReconstructed) {
        histos.fill(HIST("MC/generatedAndFindable"), kFindable);
        bool omegaInCascade = false;
        int64_t cascadeCollisionId = -1;
        for (const auto& cascade : cascades) {
          const auto cascadeV0 = cascade.v0();
          const std::array<int64_t, 3> cascadeTrackIds{
            cascadeV0.posTrackId(), cascadeV0.negTrackId(), cascade.bachelorId()};
          bool allOmegaTracksInCascade = true;
          for (std::size_t iDaughter = 0; iDaughter < 3; ++iDaughter) {
            if (std::find(cascadeTrackIds.begin(), cascadeTrackIds.end(), daughterTrackIds[iDaughter]) == cascadeTrackIds.end()) {
              allOmegaTracksInCascade = false;
              break;
            }
          }
          if (allOmegaTracksInCascade) {
            omegaInCascade = true;
            cascadeCollisionId = cascade.collisionId();
            break;
          }
        }

        bool directLambdaInV0 = false;
        int64_t directLambdaCollisionId = -1;
        for (const auto& v0 : v0s) {
          const std::array<int64_t, 2> v0TrackIds{v0.posTrackId(), v0.negTrackId()};
          bool allDirectLambdaTracksInV0 = true;
          for (std::size_t iDaughter = 3; iDaughter < 5; ++iDaughter) {
            if (std::find(v0TrackIds.begin(), v0TrackIds.end(), daughterTrackIds[iDaughter]) == v0TrackIds.end()) {
              allDirectLambdaTracksInV0 = false;
              break;
            }
          }
          if (allDirectLambdaTracksInV0) {
            directLambdaInV0 = true;
            directLambdaCollisionId = v0.collisionId();
            break;
          }
        }

        if (omegaInCascade && directLambdaInV0) {
          histos.fill(HIST("MC/generatedAndFindable"), kFindableCascadeAndV0);
          if (allDaughtersSelected) {
            histos.fill(HIST("MC/generatedAndFindable"), kFindableCascadeAndV0SelectedTracks);
          }
          const auto directKaonTrack = tracks.rawIteratorAt(daughterTrackIds[5]);
          LOG(info) << "Findable collision IDs: cascade = " << cascadeCollisionId
                    << ", direct Lambda = " << directLambdaCollisionId
                    << ", direct kaon = " << directKaonTrack.collisionId();
          if (cascadeCollisionId == directLambdaCollisionId && cascadeCollisionId == directKaonTrack.collisionId()) {
            LOG(info) << "All findable daughter tracks are in the same collision.";
          }
          printDaughterTracks("findable cascade and V0");
        }
      }
      if (allDaughtersSelected) {
        histos.fill(HIST("MC/generatedAndFindable"), kFindableSelectedTracks);
        printDaughterTracks("findable selected tracks");
        DoubleOmegaMCInfo mcInfo;
        if (!getGeneratedMCInfo(mcInfo, mcParticle)) {
          continue;
        }
        const bool isReco = fillFromSelectedTracks(collisions, tracks, selectedTrackIds, mcInfo);
        if (!isReco) {
          writeMCCandidate(DoubleOmegaCandidate{}, mcInfo, false);
        }
      }
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processFindableTracks, "Count findable double-Omega decays and reconstruct selected daughter tracks", false);

  // Mixed-event background within a dataframe. The collisions are binned in primary-vertex z and
  // multiplicity, the candidates of each collision are built once and then paired.
  void processDataMixing(CollisionsTable const& collisions,
                         TracksFull const& tracks,
                         aod::V0s const& v0s,
                         aod::Cascades const& cascades,
                         aod::BCsWithTimestamps const& bcs)
  {
    std::unordered_map<int64_t, CollisionCandidates> candidatesPerCollision;
    for (const auto& collision : collisions) {
      if (!acceptCollision(collision, bcs)) {
        continue;
      }
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionData, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(rawCascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(rawV0sPerCollision, collision.globalIndex());
      v0sThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&v0s);
      candidatesPerCollision.emplace(collision.globalIndex(),
                                     buildCollisionCandidates(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, false));
    }

    MixingBinning binning{{mixingOpts.binsVtxZ, mixingOpts.binsMultiplicity}, true};
    for (const auto& [collision0, collision1] : selfPairCombinations(binning, mixingOpts.nMixingCollisions, -1, collisions)) {
      const auto candidates0 = candidatesPerCollision.find(collision0.globalIndex());
      const auto candidates1 = candidatesPerCollision.find(collision1.globalIndex());
      if (candidates0 == candidatesPerCollision.end() || candidates1 == candidatesPerCollision.end()) {
        continue;
      }
      histos.fill(HIST("Mixing/hMixedCollisions"), 0);
      if (std::abs(collision0.posZ() - collision1.posZ()) > mixingOpts.maxDeltaPVz) {
        continue;
      }
      histos.fill(HIST("Mixing/hMixedCollisions"), 1);
      mixPair(candidates0->second, candidates1->second);
    }

    if (!(mixingOpts.mixingTypes & (1 << (kMixedAll - 1)))) {
      return;
    }
    for (const auto& [collision0, collision1, collision2] : selfTripleCombinations(binning, mixingOpts.nMixingCollisions, -1, collisions)) {
      const auto candidates0 = candidatesPerCollision.find(collision0.globalIndex());
      const auto candidates1 = candidatesPerCollision.find(collision1.globalIndex());
      const auto candidates2 = candidatesPerCollision.find(collision2.globalIndex());
      if (candidates0 == candidatesPerCollision.end() || candidates1 == candidatesPerCollision.end() || candidates2 == candidatesPerCollision.end()) {
        continue;
      }
      histos.fill(HIST("Mixing/hMixedCollisions"), 2);
      const auto [minPVz, maxPVz] = std::minmax({collision0.posZ(), collision1.posZ(), collision2.posZ()});
      if (maxPVz - minPVz > mixingOpts.maxDeltaPVz) {
        continue;
      }
      histos.fill(HIST("Mixing/hMixedCollisions"), 3);
      mixTriple(candidates0->second, candidates1->second, candidates2->second);
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processDataMixing, "Mixed-event Omega + Lambda + kaon background within a dataframe", false);

  // Bin index of value on a ConfigurableAxis (variable or fixed width), -1 outside the axis.
  static int findBin(std::vector<double> const& axis, double value)
  {
    if (axis.size() < 3) {
      return -1;
    }
    if (axis[0] == VARIABLE_WIDTH) {
      const auto edge = std::upper_bound(axis.begin() + 1, axis.end(), value);
      if (edge == axis.begin() + 1 || edge == axis.end()) {
        return -1;
      }
      return static_cast<int>(std::distance(axis.begin() + 2, edge));
    }
    const int nBins = static_cast<int>(axis[0]);
    if (value < axis[1] || value >= axis[2]) {
      return -1;
    }
    return static_cast<int>((value - axis[1]) / (axis[2] - axis[1]) * nBins);
  }

  static int nBins(std::vector<double> const& axis)
  {
    return axis[0] == VARIABLE_WIDTH ? static_cast<int>(axis.size()) - 2 : static_cast<int>(axis[0]);
  }

  // Daughter track indices identifying a combination. Seeds sharing a daughter would produce the same
  // combination more than once: each combination is written only once.
  using MixedCombination = std::array<int64_t, 6>;

  // Mixing pool of same-event candidates per (radius, phi, charge) bin, kept across dataframes.
  std::unordered_map<int, std::deque<MixingSeed>> mSeedPool;
  std::size_t mPoolSeeds = 0;
  int64_t mDataframeCounter = 0;
  double mProcessingTimeSum = 0.;

  // Indices are only unique within a dataframe: the dataframe counter goes in the upper bits.
  int64_t uniqueIndex(int64_t index) const
  {
    return index < 0 ? index : (mDataframeCounter << 32) | index;
  }

  MixingSeed makeSeed(int64_t collisionId, std::array<float, 3> const& primaryVertex, BuiltOmegaCandidate const& omega,
                      BuiltLambdaCandidate const& lambda, BuiltKaonCandidate const& kaon) const
  {
    MixingSeed seed{uniqueIndex(collisionId), primaryVertex, mBz, omega, lambda, kaon};
    for (auto* id : {&seed.omega.posTrackId, &seed.omega.negTrackId, &seed.omega.bachelorId,
                     &seed.lambda.posTrackId, &seed.lambda.negTrackId, &seed.kaon.trackId}) {
      *id = uniqueIndex(*id);
    }
    return seed;
  }

  // Builds and writes one mixed candidate, with the primary vertex of the reference seed.
  void mixSeedDaughters(MixingSeed const& reference, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda,
                        BuiltKaonCandidate const& kaon, MixingType mixingType, std::set<MixedCombination>& written)
  {
    if (!(mixingOpts.mixingTypes & (1 << (mixingType - 1)))) {
      return;
    }
    const MixedCombination combination{omega.posTrackId, omega.negTrackId, omega.bachelorId, lambda.posTrackId, lambda.negTrackId, kaon.trackId};
    if (!written.insert(combination).second || sharesTracks(omega, lambda, kaon)) {
      return;
    }
    DoubleOmegaCandidate cand;
    if (!buildDoubleOmega(reference.primaryVertex, omega, lambda, kaon, cand)) {
      return;
    }
    histos.fill(HIST("Mixing/hMixedCandidates"), mixingType);
    writeMixedCandidate(cand, mixingType);
  }

  // One daughter of each seed replaced by the corresponding daughter of the other seed. The seed
  // keeping two of its daughters provides the primary vertex.
  void mixSeedPair(MixingSeed const& seed0, MixingSeed const& seed1, std::set<MixedCombination>& written)
  {
    mixSeedDaughters(seed0, seed1.omega, seed0.lambda, seed0.kaon, kMixedOmega, written);
    mixSeedDaughters(seed1, seed0.omega, seed1.lambda, seed1.kaon, kMixedOmega, written);
    mixSeedDaughters(seed0, seed0.omega, seed1.lambda, seed0.kaon, kMixedLambda, written);
    mixSeedDaughters(seed1, seed1.omega, seed0.lambda, seed1.kaon, kMixedLambda, written);
    mixSeedDaughters(seed0, seed0.omega, seed0.lambda, seed1.kaon, kMixedKaon, written);
    mixSeedDaughters(seed1, seed1.omega, seed1.lambda, seed0.kaon, kMixedKaon, written);
  }

  // Two daughters of a seed replaced by those of two other seeds: Omega, Lambda and kaon from three
  // different collisions, in all the assignments. The seed providing the Omega is the reference.
  void mixSeedTriple(MixingSeed const& seed0, MixingSeed const& seed1, MixingSeed const& seed2, std::set<MixedCombination>& written)
  {
    const std::array<const MixingSeed*, 3> seeds{&seed0, &seed1, &seed2};
    std::array<int, 3> roles{0, 1, 2}; // Seeds providing the Omega, Lambda and kaon.
    do {
      const auto& omegaSeed = *seeds[roles[0]];
      mixSeedDaughters(omegaSeed, omegaSeed.omega, seeds[roles[1]]->lambda, seeds[roles[2]]->kaon, kMixedAll, written);
    } while (std::next_permutation(roles.begin(), roles.end()));
  }

  // Seeds can be mixed if they come from different collisions with close primary vertices, recorded
  // with the same magnetic field.
  bool canMix(MixingSeed const& seed0, MixingSeed const& seed1) const
  {
    return seed0.collisionId != seed1.collisionId &&
           std::abs(seed0.bz - seed1.bz) < 1.e-3f &&
           std::abs(seed0.primaryVertex[2] - seed1.primaryVertex[2]) <= mixingOpts.maxDeltaPVz.value;
  }

  // Mixes a new seed with the most recent compatible seeds of its pool bin: pairs with each of them,
  // triples with each pair of them coming from different collisions.
  void mixWithPool(MixingSeed const& seed, std::deque<MixingSeed> const& pool, std::set<MixedCombination>& written)
  {
    const auto nPartners = static_cast<std::size_t>(mixingOpts.nMixingCandidates.value);
    std::vector<const MixingSeed*> partners;
    for (auto partner = pool.rbegin(); partner != pool.rend() && partners.size() < nPartners; ++partner) {
      histos.fill(HIST("Mixing/hMixedCollisions"), 0);
      if (!canMix(seed, *partner)) {
        continue;
      }
      histos.fill(HIST("Mixing/hMixedCollisions"), 1);
      partners.push_back(&*partner);
      mixSeedPair(seed, *partner, written);
    }
    if (!(mixingOpts.mixingTypes & (1 << (kMixedAll - 1)))) {
      return;
    }
    for (std::size_t j = 0; j < partners.size(); ++j) {
      for (std::size_t k = j + 1; k < partners.size(); ++k) {
        histos.fill(HIST("Mixing/hMixedCollisions"), 2);
        if (!canMix(*partners[j], *partners[k])) {
          continue;
        }
        histos.fill(HIST("Mixing/hMixedCollisions"), 3);
        mixSeedTriple(seed, *partners[j], *partners[k], written);
      }
    }
  }

  void monitorMixing(std::chrono::steady_clock::time_point const& start)
  {
    const double processingTime = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    mProcessingTimeSum += processingTime;
    // Lower bound: deque and hash-map bookkeeping are not included.
    const double poolMemory = static_cast<double>(mPoolSeeds * sizeof(MixingSeed)) / (1024. * 1024.);
    histos.fill(HIST("Mixing/hProcessingTime"), processingTime);
    histos.fill(HIST("Mixing/hPoolSeeds"), mPoolSeeds);
    histos.fill(HIST("Mixing/hPoolMemory"), poolMemory);
    histos.fill(HIST("Mixing/hProcessingTimeVsPoolSeeds"), mPoolSeeds, processingTime);
    if (mixingOpts.monitoringInterval > 0 && mDataframeCounter % mixingOpts.monitoringInterval == 0) {
      LOGP(info, "Candidate mixing, dataframe {}: {} seeds in {} pool bins ({:.2f} MB), {:.1f} ms in this dataframe, {:.1f} ms on average",
           mDataframeCounter, mPoolSeeds, mSeedPool.size(), poolMemory, processingTime, mProcessingTimeSum / mDataframeCounter);
    }
  }

  // Candidate-level mixed-event background, as in decay3bodybuilder. The same-event candidates
  // (seeds) are binned in decay radius, momentum azimuth and charge, so seeds in the same bin have
  // their decay vertex in the same region of the detector and point in the same direction. Mixed
  // candidates are built by swapping one (or two) daughters between seeds of different collisions,
  // which keeps the vertex topology without moving any track. Seeds are pooled across dataframes,
  // since a single dataframe rarely contains seeds from more than one collision.
  void processDataMixingCandidates(CollisionsTable const& collisions,
                                   TracksFull const& tracks,
                                   aod::V0s const& v0s,
                                   aod::Cascades const& cascades,
                                   aod::BCsWithTimestamps const& bcs)
  {
    const auto start = std::chrono::steady_clock::now();
    ++mDataframeCounter;
    const int nPhiBins = nBins(mixingOpts.binsPhi.value);
    const auto poolDepth = static_cast<std::size_t>(mixingOpts.poolDepth.value);
    std::set<MixedCombination> written; // Every combination has a daughter from this dataframe.

    for (const auto& collision : collisions) {
      if (!acceptCollision(collision, bcs)) {
        continue;
      }
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionData, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(rawCascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(rawV0sPerCollision, collision.globalIndex());
      v0sThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&v0s);
      const auto candidates = buildCollisionCandidates(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, true);

      // The seeds of this collision are mixed with the pool first and added to it afterwards, so
      // they are never mixed with each other.
      std::vector<std::pair<int, MixingSeed>> newSeeds;
      combine(candidates, candidates, candidates, candidates,
              [&](DoubleOmegaCandidate const& cand, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& lambda, BuiltKaonCandidate const& kaon) {
                const float radius = std::hypot(cand.x, cand.y);
                const int radiusBin = findBin(mixingOpts.binsRadius.value, radius);
                const int phiBin = findBin(mixingOpts.binsPhi.value, cand.phi);
                if (radiusBin < 0 || phiBin < 0) {
                  return;
                }
                histos.fill(HIST("Mixing/hSeedRadiusPhi"), radius, cand.phi);
                const int key = (radiusBin * nPhiBins + phiBin) * 2 + chargeIndex(cand.sign);
                newSeeds.emplace_back(key, makeSeed(collision.globalIndex(), candidates.primaryVertex, omega, lambda, kaon));
              });
      for (const auto& [key, seed] : newSeeds) {
        mixWithPool(seed, mSeedPool[key], written);
      }
      for (auto& [key, seed] : newSeeds) {
        auto& pool = mSeedPool[key];
        pool.push_back(std::move(seed));
        ++mPoolSeeds;
        if (pool.size() > poolDepth) {
          pool.pop_front();
          --mPoolSeeds;
        }
      }
    }
    monitorMixing(start);
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processDataMixingCandidates, "Candidate-level mixed-event background, seeds binned in decay radius and azimuth and pooled across dataframes", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<doubleOmegaTreeCreator>(cfgc)};
}
