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
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <GPU/GPUROOTCartesianFwd.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsTable = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>;
using Collisions = CollisionsTable::iterator;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults, aod::McCollisionLabels>;
using FullCascades = aod::CascDataExt;
using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;
using TracksFullIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::McTrackLabels>;
using TracksFullMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::McTrackLabels>;

struct DoubleOmegaCandidate {
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float x = -999.f;
  float y = -999.f;
  float z = -999.f;
  float cosPAOmega = -999.f;
  float cosPADirectLambda = -999.f;
  float cosPADoubleOmega = -999.f;
  float dcaXYOmegaToPV = -999.f;
  float dcaZOmegaToPV = -999.f;
  float dcaXYDirectLambdaToPV = -999.f;
  float dcaZDirectLambdaToPV = -999.f;
  float dcaXYDirectKaonToPV = -999.f;
  float dcaZDirectKaonToPV = -999.f;
  float mass = -999.f;
  float massOmega = -999.f;
  float massXi = -999.f;
};

struct DoubleOmegaMCInfo {
  int64_t motherId = -1;
  float pt = -999.f;
  float eta = -999.f;
  float phi = -999.f;
  float decayLength = -999.f;
  int pdgCode = 0;
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
  int64_t v0Id = -1;
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
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
  int64_t cascadeId = -1;
  int64_t posTrackId = -1;
  int64_t negTrackId = -1;
  int64_t bachelorId = -1;
  o2::track::TrackParCov parentTrack;
};

struct doubleOmegaTreeCreator {
  Produces<aod::DoubleOmegaTable> doubleOmegaTable;
  Produces<aod::DoubleOmegaTableMC> doubleOmegaTableMC;
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
    kNFindabilitySteps
  };
  std::vector<int64_t> reconstructedDoubleOmegaIds;

  Preslice<TracksFullIUMC> tracksIUPerCollision = aod::track::collisionId;
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

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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

    LOG(info) << "Retrieved GRP for timestamp " << bc.timestamp() << " with magnetic field " << mBz << " kG";
    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDoubleOmega,fOmegaXi");
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
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

  template <class T, class V0>
  bool buildLambda(T const&, V0 const& v0, int8_t charge, BuiltLambdaCandidate& builtLambda)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    if (!selectTrack(posTrack) || !selectTrack(negTrack)) {
      return false;
    }

    const bool isMatter = charge < 0;
    const float protonNSigma = isMatter ? posTrack.tpcNSigmaPr() : negTrack.tpcNSigmaPr();
    const float pionNSigma = isMatter ? negTrack.tpcNSigmaPi() : posTrack.tpcNSigmaPi();
    if (std::abs(protonNSigma) > nSigmaTPCCut || std::abs(pionNSigma) > nSigmaTPCCut ||
        std::abs(posTrack.dcaXY()) < dcaV0DauToPV || std::abs(negTrack.dcaXY()) < dcaV0DauToPV) {
      return false;
    }

    auto posTrackParCov = getTrackParCov(posTrack);
    auto negTrackParCov = getTrackParCov(negTrack);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(posTrackParCov, negTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting raw V0 " << v0.globalIndex();
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
    builtLambda.v0Id = v0.globalIndex();
    builtLambda.posTrackId = v0.posTrackId();
    builtLambda.negTrackId = v0.negTrackId();
    return true;
  }

  template <class T, class Casc>
  bool buildOmega(T const& tracks, Casc const& cascade, BuiltOmegaCandidate& builtOmega)
  {
    auto bachelor = cascade.template bachelor_as<T>();
    if (bachelor.sign() == 0 ||
        !selectTrack(bachelor) ||
        std::abs(bachelor.tpcNSigmaKa()) > nSigmaTPCCut ||
        std::abs(bachelor.dcaXY()) < dcaBachToPV) {
      return false;
    }

    auto v0 = cascade.template v0_as<aod::V0s>();
    BuiltLambdaCandidate lambda;
    if (!buildLambda(tracks, v0, bachelor.sign(), lambda)) {
      return false;
    }

    auto v0TrackParCov = lambda.parentTrack;
    auto bachelorTrackParCov = getTrackParCov(bachelor);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(v0TrackParCov, bachelorTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while fitting raw cascade " << cascade.globalIndex();
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
      cascade.globalIndex(),
      v0.posTrackId(),
      v0.negTrackId(),
      cascade.bachelorId(),
      fitter.createParentTrackParCov()};
    return true;
  }

  template <class T, class Casc>
  bool buildOmegaFromBuilder(T const& tracks, Casc const& cascade, BuiltOmegaCandidate& builtOmega)
  {
    auto bachelor = cascade.template bachelor_as<T>();
    if (bachelor.sign() == 0 ||
        !selectTrack(bachelor) ||
        std::abs(bachelor.tpcNSigmaKa()) > nSigmaTPCCut ||
        std::abs(bachelor.dcaXY()) < dcaBachToPV) {
      return false;
    }

    BuiltLambdaCandidate lambda;
    if (!buildLambda(tracks, cascade, bachelor.sign(), lambda)) {
      return false;
    }

    auto v0TrackParCov = lambda.parentTrack;
    auto bachelorTrackParCov = getTrackParCov(bachelor);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(v0TrackParCov, bachelorTrackParCov);
    } catch (...) {
      LOG(error) << "Exception while refitting cascade " << cascade.globalIndex();
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
      cascade.globalIndex(),
      cascade.posTrackId(),
      cascade.negTrackId(),
      cascade.bachelorId(),
      fitter.createParentTrackParCov()};
    return true;
  }

  static float etaFromMomentum(float px, float py, float pz)
  {
    const float pt = std::hypot(px, py);
    return pt > 0.f ? std::asinh(pz / pt) : 0.f;
  }

  template <class C, class Kaon>
  bool buildDoubleOmega(C const& collision, BuiltOmegaCandidate const& omega, BuiltLambdaCandidate const& directLambda, Kaon const& kaon, DoubleOmegaCandidate& cand)
  {
    auto omegaTrackParCov = omega.parentTrack;
    auto lambdaTrackParCov = directLambda.parentTrack;
    auto kaonTrackParCov = getTrackParCov(kaon);
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
    const std::array<float, 3> primaryVertex{collision.posX(), collision.posY(), collision.posZ()};
    const o2::math_utils::Point3D<float> primaryVertexPoint{collision.posX(), collision.posY(), collision.posZ()};
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
    cand.cosPAOmega = RecoDecay::cpa(decayVertex, omegaDecayVertex, omegaMomentum);
    cand.cosPADirectLambda = RecoDecay::cpa(decayVertex, directLambda.decayVertex, lambdaMomentum);
    cand.cosPADoubleOmega = RecoDecay::cpa(primaryVertex, decayVertex, totalMomentum);
    cand.dcaXYOmegaToPV = dcaOmega[0];
    cand.dcaZOmegaToPV = dcaOmega[1];
    cand.dcaXYDirectLambdaToPV = dcaDirectLambda[0];
    cand.dcaZDirectLambdaToPV = dcaDirectLambda[1];
    cand.dcaXYDirectKaonToPV = dcaDirectKaon[0];
    cand.dcaZDirectKaonToPV = dcaDirectKaon[1];
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
    const int expectedPosPdg = sign < 0 ? 2212 : 211;
    const int expectedNegPdg = sign < 0 ? -211 : -2212;
    const int expectedLambdaPdg = sign < 0 ? 3122 : -3122;
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
    const int expectedBachelorPdg = sign < 0 ? -321 : 321;
    const int expectedOmegaPdg = sign < 0 ? 3334 : -3334;
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
    const int expectedKaonPdg = sign < 0 ? -321 : 321;
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
    const int expectedOmegaPdg = isMatter ? 3334 : -3334;
    const int expectedLambdaPdg = isMatter ? 3122 : -3122;
    const int expectedKaonPdg = isMatter ? -321 : 321;
    bool foundOmega = false;
    bool foundLambda = false;
    bool foundKaon = false;
    std::array<float, 3> decayVertex{0.f, 0.f, 0.f};

    for (const auto& daughter : particle.template daughters_as<aod::McParticles>()) {
      if (daughter.pdgCode() == expectedOmegaPdg) {
        foundOmega = true;
        decayVertex = {daughter.vx(), daughter.vy(), daughter.vz()};
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
    doubleOmegaTable(cand.pt,
                     cand.eta,
                     cand.phi,
                     cand.x,
                     cand.y,
                     cand.z,
                     cand.cosPAOmega,
                     cand.cosPADirectLambda,
                     cand.cosPADoubleOmega,
                     cand.dcaXYOmegaToPV,
                     cand.dcaZOmegaToPV,
                     cand.dcaXYDirectLambdaToPV,
                     cand.dcaZDirectLambdaToPV,
                     cand.dcaXYDirectKaonToPV,
                     cand.dcaZDirectKaonToPV,
                     cand.mass,
                     cand.massOmega,
                     cand.massXi);
  }

  void writeMCCandidate(DoubleOmegaCandidate const& cand, DoubleOmegaMCInfo const& mcInfo, bool isReco)
  {
    doubleOmegaTableMC(cand.pt,
                       cand.eta,
                       cand.phi,
                       cand.x,
                       cand.y,
                       cand.z,
                       cand.cosPAOmega,
                       cand.cosPADirectLambda,
                       cand.cosPADoubleOmega,
                       cand.dcaXYOmegaToPV,
                       cand.dcaZOmegaToPV,
                       cand.dcaXYDirectLambdaToPV,
                       cand.dcaZDirectLambdaToPV,
                       cand.dcaXYDirectKaonToPV,
                       cand.dcaZDirectKaonToPV,
                       cand.mass,
                       cand.massOmega,
                       cand.massXi,
                       mcInfo.pt,
                       mcInfo.eta,
                       mcInfo.phi,
                       mcInfo.decayLength,
                       mcInfo.pdgCode,
                       isReco);
  }

  template <bool isMC, class C, class T, class Cascades, class V0s>
  void fillFromRawTables(C const& collision, T const& tracks, Cascades const& cascades, V0s const& v0s,
                         aod::McParticles const* mcParticles = nullptr)
  {
    for (const auto& rawCascade : cascades) {
      BuiltOmegaCandidate omega;
      if (!buildOmega(tracks, rawCascade, omega)) {
        continue;
      }
      const float omegaPt = std::hypot(omega.px, omega.py);
      histos.fill(HIST("QA/massXi"), omegaPt, omega.massXi);
      histos.fill(HIST("QA/massOmega"), omegaPt, omega.massOmega);
      const std::array<int64_t, 3> omegaTrackIds{omega.posTrackId, omega.negTrackId, omega.bachelorId};

      for (const auto& rawV0 : v0s) {
        BuiltLambdaCandidate lambda;
        if (!buildLambda(tracks, rawV0, omega.sign, lambda)) {
          continue;
        }
        if (std::find(omegaTrackIds.begin(), omegaTrackIds.end(), lambda.posTrackId) != omegaTrackIds.end() ||
            std::find(omegaTrackIds.begin(), omegaTrackIds.end(), lambda.negTrackId) != omegaTrackIds.end()) {
          continue;
        }
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);

        for (const auto& kaon : tracks) {
          if (!isSelectedKaon(kaon, omega.sign) ||
              kaon.globalIndex() == lambda.posTrackId ||
              kaon.globalIndex() == lambda.negTrackId ||
              std::find(omegaTrackIds.begin(), omegaTrackIds.end(), kaon.globalIndex()) != omegaTrackIds.end()) {
            continue;
          }
          DoubleOmegaCandidate cand;
          if (!buildDoubleOmega(collision, omega, lambda, kaon, cand)) {
            continue;
          }
          if constexpr (isMC) {
            if (mcParticles == nullptr) {
              continue;
            }
            auto omegaV0 = rawCascade.template v0_as<aod::V0s>();
            auto omegaPosTrack = omegaV0.template posTrack_as<T>();
            auto omegaNegTrack = omegaV0.template negTrack_as<T>();
            auto omegaBachelorTrack = rawCascade.template bachelor_as<T>();
            auto lambdaPosTrack = rawV0.template posTrack_as<T>();
            auto lambdaNegTrack = rawV0.template negTrack_as<T>();
            DoubleOmegaMCInfo mcInfo;
            if (getMCInfo(mcInfo,
                          omegaPosTrack, omegaNegTrack, omegaBachelorTrack,
                          lambdaPosTrack, lambdaNegTrack, kaon,
                          omega.sign, *mcParticles)) {
              reconstructedDoubleOmegaIds.push_back(mcInfo.motherId);
              writeMCCandidate(cand, mcInfo, true);
            }
          } else {
            writeDataCandidate(cand);
          }
        }
      }
    }
  }

  template <class C, class T, class Cascades>
  void fillFromCascades(C const& collision, T const& tracks, Cascades const& cascades)
  {
    for (const auto& omegaRow : cascades) {
      BuiltOmegaCandidate omega;
      if (!buildOmegaFromBuilder(tracks, omegaRow, omega)) {
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
        if (!buildLambda(tracks, lambdaKaonSource, omega.sign, lambda)) {
          continue;
        }
        histos.fill(HIST("QA/massLambda"), std::hypot(lambda.candidate.px, lambda.candidate.py), lambda.candidate.mass);
        DoubleOmegaCandidate cand;
        if (buildDoubleOmega(collision, omega, lambda, kaon, cand)) {
          writeDataCandidate(cand);
        }
      }
    }
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

  void processData(Collisions const& collision,
                   TracksFullIU const& tracks,
                   aod::V0s const& v0s,
                   aod::Cascades const& cascades,
                   aod::BCsWithTimestamps const& bcs)
  {
    if (acceptCollision(collision, bcs)) {
      fillFromRawTables<false>(collision, tracks, cascades, v0s);
    }
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
                 TracksFullIUMC const& tracks,
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
      auto tracksThisCollision = tracks.sliceBy(tracksIUPerCollision, collision.globalIndex());
      auto cascadesThisCollision = cascades.sliceBy(rawCascadesPerCollision, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(rawV0sPerCollision, collision.globalIndex());
      v0sThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&tracks);
      cascadesThisCollision.bindExternalIndices(&v0s);
      fillFromRawTables<true>(collision, tracksThisCollision, cascadesThisCollision, v0sThisCollision, &mcParticles);
    }
    LOG(info) << "Found " << reconstructedDoubleOmegaIds.size() << " reconstructed double-Omega candidates in MC";
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

  void processFindableTracks(TracksFullMC const& tracks,
                             aod::Cascades const& cascades,
                             aod::V0s const& v0s,
                             aod::McParticles const& mcParticles)
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
      std::array<int64_t, 6> daughterTrackIds{};
      daughterTrackIds.fill(-1);
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
          }
        }
      }
      bool allDaughtersReconstructed = true;
      bool allDaughtersSelected = true;
      for (std::size_t iDaughter = 0; iDaughter < reconstructedDaughters.size(); ++iDaughter) {
        allDaughtersReconstructed &= reconstructedDaughters[iDaughter];
        allDaughtersSelected &= selectedReconstructedDaughters[iDaughter];
      }
      if (allDaughtersReconstructed) {
        histos.fill(HIST("MC/generatedAndFindable"), kFindable);
        bool omegaInCascade = false;
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
            break;
          }
        }

        bool directLambdaInV0 = false;
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
            break;
          }
        }

        if (omegaInCascade && directLambdaInV0) {
          histos.fill(HIST("MC/generatedAndFindable"), kFindableCascadeAndV0);
        }
      }
      if (allDaughtersSelected) {
        histos.fill(HIST("MC/generatedAndFindable"), kFindableSelectedTracks);
        LOG(debug) << "----------------------------------------";
        for (std::size_t iDaughter = 0; iDaughter < daughterTrackIds.size(); ++iDaughter) {
          LOG(debug) << "+++++++";
          const auto daughterTrack = tracks.rawIteratorAt(daughterTrackIds[iDaughter]);
          const auto daughterMCParticle = daughterTrack.mcParticle_as<aod::McParticles>();
          LOG(debug) << "Dau" << iDaughter + 1
                     << ": hasITS = " << daughterTrack.hasITS()
                     << ", hasTPC = " << daughterTrack.hasTPC()
                     << ", hasTOF = " << daughterTrack.hasTOF()
                     << ", hasTRD = " << daughterTrack.hasTRD()
                     << ", isITSAfterburner = " << daughterTrack.isITSAfterburner();
          LOG(debug) << "Dau" << iDaughter + 1
                     << ": eta = " << daughterTrack.eta()
                     << ", pt = " << daughterTrack.pt()
                     << ", CollisionID = " << daughterTrack.collisionId()
                     << ", production point = (" << daughterMCParticle.vx()
                     << ", " << daughterMCParticle.vy()
                     << ", " << daughterMCParticle.vz() << ") cm";
        }
      }
    }
  }
  PROCESS_SWITCH(doubleOmegaTreeCreator, processFindableTracks, "Count generated and track-findable double-Omega decays", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<doubleOmegaTreeCreator>(cfgc)};
}
