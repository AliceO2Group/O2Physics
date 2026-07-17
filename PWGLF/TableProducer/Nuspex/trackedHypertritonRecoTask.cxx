// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file trackedHypertritonRecoTask.cxx
/// \brief Build two- and three-body hypertriton analysis tables from strangeness-tracked candidates

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "PWGLF/DataModel/LFHypernucleiTables.h"
#include "PWGLF/DataModel/LFPIDTOFGenericTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/Utils/decay3bodyBuilderHelper.h"
#include "PWGLF/Utils/pidTOFGeneric.h"

#include "Common/Core/MetadataHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH1.h>
#include <TH2.h>

#include <KFParticle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{
using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::EvTimeTOFFT0>;
using Tracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                         aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe,
                         aod::TOFSignal, aod::TOFEvTime, aod::EvTimeTOFFT0ForTrack>;

constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
const std::vector<std::string> particleName{"He3"};
o2::common::core::MetadataHelper metadataInfo;

enum ZorroTrigger : size_t {
  kHe = 0,
  kTracked3Body,
  kNZorroTriggers
};
} // namespace

struct TrackedHypertritonRecoTask {
  Produces<aod::DataHypCands> dataHypCands;
  Produces<aod::Vtx3BodyDatas> vtx3BodyDatas;
  Produces<aod::Vtx3BodyCovs> vtx3BodyCovs;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> requireSel8{"requireSel8", true, "Require the standard sel8 event selection"};
  Configurable<float> maxAbsZvtx{"maxAbsZvtx", 10.f, "Maximum absolute primary-vertex z; negative disables the cut"};
  Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Enable Zorro accounting and trigger selection for skimmed data"};
  Configurable<std::string> zorroCCDBPath{"zorroCCDBPath", "EventFiltering/Zorro/", "Base path of the Zorro CCDB objects"};
  Configurable<int> zorroBCTolerance{"zorroBCTolerance", 100, "Zorro BC matching tolerance"};

  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the magnetic-field object"};
  Configurable<double> bzInput{"bz", -999., "Magnetic field in kG; values below -990 use CCDB"};

  struct : ConfigurableGroup {
    std::string prefix = "twoBody";
    Configurable<bool> useSelections{"useSelections", false, "Apply the two-body candidate selections"};
    Configurable<float> maxEtaDaughters{"maxEtaDaughters", 1.e10f, "Maximum absolute daughter eta"};
    Configurable<float> minTPCNClsHe{"minTPCNClsHe", -1.f, "Minimum He TPC clusters"};
    Configurable<float> minTPCNClsPi{"minTPCNClsPi", -1.f, "Minimum pion TPC clusters"};
    Configurable<float> minTPCCrossedRowsHe{"minTPCCrossedRowsHe", -1.f, "Minimum He TPC crossed rows"};
    Configurable<float> minTPCCrossedRowsPi{"minTPCCrossedRowsPi", -1.f, "Minimum pion TPC crossed rows"};
    Configurable<float> minTPCInnerParamHe{"minTPCInnerParamHe", -1.f, "Minimum He TPC rigidity"};
    Configurable<float> minPt{"minPt", -1.f, "Minimum candidate transverse momentum"};
    Configurable<float> massWindow{"massWindow", 1.e10f, "Half-width of the hypertriton mass window"};
    Configurable<float> maxDcaDaughters{"maxDcaDaughters", 1.e10f, "Maximum DCA between daughters"};
    Configurable<float> minCosPA{"minCosPA", -2.f, "Minimum cosine of the pointing angle"};
    Configurable<float> minDcaHeToPV{"minDcaHeToPV", -1.f, "Minimum absolute He DCA to PV"};
    Configurable<float> minDcaPiToPV{"minDcaPiToPV", -1.f, "Minimum absolute pion DCA to PV"};
    Configurable<float> minNSigmaHe{"minNSigmaHe", -1.e10f, "Minimum He TPC n-sigma"};
    Configurable<bool> compensatePIDinTracking{"compensatePIDinTracking", true, "Correct the TPC inner parameter for charge-two PID in tracking"};
    Configurable<LabeledArray<double>> betheBlochParams{"betheBlochParams", {betheBlochDefault[0], 1, 6, particleName, betheBlochParNames}, "TPC Bethe-Bloch parameters for He3"};
  } twoBody;

  struct : ConfigurableGroup {
    std::string prefix = "threeBody";
    Configurable<bool> useKFParticle{"useKFParticle", false, "Use KFParticle to build the three-body candidate"};
    Configurable<bool> setTopologicalConstraint{"setTopologicalConstraint", false, "Apply the KFParticle topological constraint"};
    Configurable<bool> useSelections{"useSelections", false, "Apply the standard three-body candidate selections"};
    Configurable<bool> useChi2Selection{"useChi2Selection", false, "Apply the candidate chi2 selection"};
    Configurable<bool> useTPCforPion{"useTPCforPion", false, "Require TPC information for the pion"};
    Configurable<bool> acceptTPCOnly{"acceptTPCOnly", false, "Accept TPC-only daughter tracks"};
    Configurable<bool> askOnlyITSMatch{"askOnlyITSMatch", true, "Require an ITS match when rejecting TPC-only tracks"};
    Configurable<bool> calculateCovariance{"calculateCovariance", true, "Calculate candidate and daughter covariance matrices"};
    Configurable<float> maxEtaDaughters{"maxEtaDaughters", 0.9f, "Maximum absolute daughter eta"};
    Configurable<int> minTPCNClProton{"minTPCNClProton", 90, "Minimum proton TPC clusters"};
    Configurable<int> minTPCNClPion{"minTPCNClPion", 70, "Minimum pion TPC clusters"};
    Configurable<int> minTPCNClDeuteron{"minTPCNClDeuteron", 100, "Minimum deuteron TPC clusters"};
    Configurable<float> minDCAProtonToPV{"minDCAProtonToPV", 0.1f, "Minimum proton DCA to PV"};
    Configurable<float> minDCAPionToPV{"minDCAPionToPV", 0.1f, "Minimum pion DCA to PV"};
    Configurable<float> minDCADeuteronToPV{"minDCADeuteronToPV", 0.1f, "Minimum deuteron DCA to PV"};
    Configurable<float> minPtProton{"minPtProton", 0.3f, "Minimum proton pT"};
    Configurable<float> minPtPion{"minPtPion", 0.1f, "Minimum pion pT"};
    Configurable<float> minPtDeuteron{"minPtDeuteron", 0.6f, "Minimum deuteron pT"};
    Configurable<float> maxPtProton{"maxPtProton", 5.f, "Maximum proton pT"};
    Configurable<float> maxPtPion{"maxPtPion", 1.2f, "Maximum pion pT"};
    Configurable<float> maxPtDeuteron{"maxPtDeuteron", 10.f, "Maximum deuteron pT"};
    Configurable<float> maxTPCNSigma{"maxTPCNSigma", 5.f, "Maximum absolute daughter TPC n-sigma"};
    Configurable<float> minTOFNSigmaDeuteron{"minTOFNSigmaDeuteron", -5.f, "Minimum deuteron TOF n-sigma"};
    Configurable<float> maxTOFNSigmaDeuteron{"maxTOFNSigmaDeuteron", 5.f, "Maximum deuteron TOF n-sigma"};
    Configurable<float> minPDeuteronUseTOF{"minPDeuteronUseTOF", 1.f, "Minimum deuteron momentum at which TOF PID is required"};
    Configurable<float> maxDCADaughtersToSVAverage{"maxDCADaughtersToSVAverage", 0.5f, "Maximum average daughter DCA to SV"};
    Configurable<float> maxRapidity{"maxRapidity", 1.f, "Maximum absolute candidate rapidity"};
    Configurable<float> minPt{"minPt", 2.f, "Minimum candidate pT"};
    Configurable<float> maxPt{"maxPt", 5.f, "Maximum candidate pT"};
    Configurable<float> minMass{"minMass", 2.96f, "Minimum candidate mass"};
    Configurable<float> maxMass{"maxMass", 3.04f, "Maximum candidate mass"};
    Configurable<float> minCtau{"minCtau", 0.f, "Minimum candidate c tau"};
    Configurable<float> maxCtau{"maxCtau", 100.f, "Maximum candidate c tau"};
    Configurable<float> minCosPA{"minCosPA", 0.9f, "Minimum candidate cosine of the pointing angle"};
    Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum candidate chi2"};
  } threeBody;

  o2::vertexing::DCAFitterN<2> fitter2Body;
  o2::pwglf::decay3bodyBuilderHelper builder3Body;
  o2::pid::tof::TOFResoParamsV3 tofResponse;
  o2::aod::pidtofgeneric::TOFCalibConfig tofCalibConfig;
  o2::aod::pidtofgeneric::TofPidNewCollision<Tracks::iterator> deuteronTOFPID;

  int runNumber = 0;
  float bz = 0.f;
  std::array<float, 6> bbParamsHe{};
  std::vector<bool> goodCollision;
  std::vector<std::array<bool, kNZorroTriggers>> zorroDecision;

  void init(InitContext& context)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(zorroCCDBPath);
    zorro.setBCtolerance(zorroBCTolerance);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter2Body.setPropagateToPCA(true);
    fitter2Body.setMaxR(200.f);
    fitter2Body.setMinParamChange(1.e-3f);
    fitter2Body.setMinRelChi2Change(0.9f);
    fitter2Body.setMaxDZIni(1.e9f);
    fitter2Body.setMaxChi2(1.e9f);
    fitter2Body.setUseAbsDCA(true);

    deuteronTOFPID.SetPidType(o2::track::PID::Deuteron);
    tofCalibConfig.metadataInfo = metadataInfo;
    tofCalibConfig.inheritFromBaseTask(context);
    tofCalibConfig.initSetup(tofResponse, ccdb);

    auto& selections = builder3Body.decay3bodyselections;
    selections.maxEtaDaughters = threeBody.maxEtaDaughters;
    selections.minTPCNClProton = threeBody.minTPCNClProton;
    selections.minTPCNClPion = threeBody.minTPCNClPion;
    selections.minTPCNClDeuteron = threeBody.minTPCNClDeuteron;
    selections.minDCAProtonToPV = threeBody.minDCAProtonToPV;
    selections.minDCAPionToPV = threeBody.minDCAPionToPV;
    selections.minDCADeuteronToPV = threeBody.minDCADeuteronToPV;
    selections.minDCAProtonToPVprop = threeBody.minDCAProtonToPV;
    selections.minDCAPionToPVprop = threeBody.minDCAPionToPV;
    selections.minDCADeuteronToPVprop = threeBody.minDCADeuteronToPV;
    selections.minPtProton = threeBody.minPtProton;
    selections.minPtPion = threeBody.minPtPion;
    selections.minPtDeuteron = threeBody.minPtDeuteron;
    selections.maxPtProton = threeBody.maxPtProton;
    selections.maxPtPion = threeBody.maxPtPion;
    selections.maxPtDeuteron = threeBody.maxPtDeuteron;
    selections.maxTPCnSigma = threeBody.maxTPCNSigma;
    selections.minTOFnSigmaDeuteron = threeBody.minTOFNSigmaDeuteron;
    selections.maxTOFnSigmaDeuteron = threeBody.maxTOFNSigmaDeuteron;
    selections.minPDeuteronUseTOF = threeBody.minPDeuteronUseTOF;
    selections.maxDCADauToSVaverage = threeBody.maxDCADaughtersToSVAverage;
    selections.maxRapidity = threeBody.maxRapidity;
    selections.minPt = threeBody.minPt;
    selections.maxPt = threeBody.maxPt;
    selections.minMass = threeBody.minMass;
    selections.maxMass = threeBody.maxMass;
    selections.minCtau = threeBody.minCtau;
    selections.maxCtau = threeBody.maxCtau;
    selections.minCosPA = threeBody.minCosPA;
    selections.maxChi2 = threeBody.maxChi2;

    registry.add("events", "Event selection;selection;events", HistType::kTH1D, {{2, -0.5, 1.5}});
    auto events = registry.get<TH1>(HIST("events"));
    events->GetXaxis()->SetBinLabel(1, "all");
    events->GetXaxis()->SetBinLabel(2, "selected");
    registry.add("zorroEvents", "Zorro accounting;trigger;selection", HistType::kTH2D,
                 {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    auto zorroEvents = registry.get<TH2>(HIST("zorroEvents"));
    zorroEvents->GetXaxis()->SetBinLabel(1, "fHe");
    zorroEvents->GetXaxis()->SetBinLabel(2, "fTracked3Body");
    zorroEvents->GetYaxis()->SetBinLabel(1, "before sel8");
    zorroEvents->GetYaxis()->SetBinLabel(2, "after sel8");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    if (skimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe,fTracked3Body");
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Missing magnetic-field object " << grpmagPath << " for timestamp " << bc.timestamp();
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    bz = bzInput < -990. ? o2::base::Propagator::Instance()->getNominalBz() : bzInput;
    o2::base::Propagator::Instance()->setNominalBz(bz);
    fitter2Body.setBz(bz);
    builder3Body.fitterV0.setBz(bz);
    builder3Body.fitter3body.setBz(bz);
#ifdef HomogeneousField
    KFParticle::SetField(bz);
#endif

    for (size_t i = 0; i < bbParamsHe.size(); ++i) {
      bbParamsHe[i] = twoBody.betheBlochParams->get("He3", betheBlochParNames[i].c_str());
    }
    tofCalibConfig.processSetup(tofResponse, ccdb, bc);
    runNumber = bc.runNumber();
  }

  template <typename TCollision>
  void selectCollisions(TCollision const& collisions)
  {
    goodCollision.assign(collisions.size(), false);
    zorroDecision.assign(collisions.size(), {});
    for (const auto& collision : collisions) {
      registry.fill(HIST("events"), 0.);
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (skimmedProcessing) {
        zorro.isSelected(bc.globalBC(), zorroBCTolerance);
        auto decisions = zorro.getTriggerOfInterestResults();
        for (size_t i = 0; i < std::min(decisions.size(), static_cast<size_t>(kNZorroTriggers)); ++i) {
          zorroDecision[collision.globalIndex()][i] = decisions[i];
          if (decisions[i]) {
            registry.fill(HIST("zorroEvents"), static_cast<double>(i), 0.);
          }
        }
      }

      if (requireSel8 && !collision.sel8()) {
        continue;
      }
      if (maxAbsZvtx >= 0.f && std::abs(collision.posZ()) > maxAbsZvtx) {
        continue;
      }
      goodCollision[collision.globalIndex()] = true;
      registry.fill(HIST("events"), 1.);
      if (skimmedProcessing) {
        for (size_t i = 0; i < kNZorroTriggers; ++i) {
          if (zorroDecision[collision.globalIndex()][i]) {
            registry.fill(HIST("zorroEvents"), static_cast<double>(i), 1.);
          }
        }
      }
    }
  }

  template <typename TTrack>
  float nSigmaHe3(TTrack const& track) const
  {
    const bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    const float rigidity = heliumPID && twoBody.compensatePIDinTracking.value ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();
    const float expected = o2::common::BetheBlochAleph(rigidity * 2.f / static_cast<float>(constants::physics::MassHelium3),
                                                       bbParamsHe[0], bbParamsHe[1], bbParamsHe[2], bbParamsHe[3], bbParamsHe[4]);
    return (track.tpcSignal() - expected) / (expected * bbParamsHe[5]);
  }

  template <typename TTrack, typename TCollision>
  void buildTwoBody(TTrack const& heTrack, TTrack const& piTrack, TCollision const& collision, float trackedClSize)
  {
    auto heTrackCov = getTrackParCov(heTrack);
    auto piTrackCov = getTrackParCov(piTrack);
    int nCandidates = 0;
    try {
      nCandidates = fitter2Body.process(heTrackCov, piTrackCov);
    } catch (...) {
      LOG(error) << "Exception while fitting a tracked two-body candidate";
      return;
    }
    if (nCandidates == 0) {
      return;
    }

    std::array<float, 3> heMomentum{};
    std::array<float, 3> piMomentum{};
    fitter2Body.getTrack(0).getPxPyPzGlo(heMomentum);
    fitter2Body.getTrack(1).getPxPyPzGlo(piMomentum);
    for (auto& component : heMomentum) {
      component *= 2.f;
    }
    std::array<float, 3> momentum{heMomentum[0] + piMomentum[0], heMomentum[1] + piMomentum[1], heMomentum[2] + piMomentum[2]};
    if (twoBody.useSelections && std::hypot(momentum[0], momentum[1]) < twoBody.minPt) {
      return;
    }

    const float heMomentum2 = RecoDecay::sumOfSquares(heMomentum[0], heMomentum[1], heMomentum[2]);
    const float piMomentum2 = RecoDecay::sumOfSquares(piMomentum[0], piMomentum[1], piMomentum[2]);
    const float candidateMomentum2 = RecoDecay::sumOfSquares(momentum[0], momentum[1], momentum[2]);
    const float heEnergy = std::sqrt(heMomentum2 + constants::physics::MassHelium3 * constants::physics::MassHelium3);
    const float piEnergy = std::sqrt(piMomentum2 + constants::physics::MassPionCharged * constants::physics::MassPionCharged);
    const float mass = std::sqrt((heEnergy + piEnergy) * (heEnergy + piEnergy) - candidateMomentum2);
    if (twoBody.useSelections && std::abs(mass - constants::physics::MassHyperTriton) > twoBody.massWindow) {
      return;
    }

    const auto& secondaryVertex = fitter2Body.getPCACandidate();
    const std::array<float, 3> primaryVertex{collision.posX(), collision.posY(), collision.posZ()};
    const std::array<float, 3> decayVertex{static_cast<float>(secondaryVertex[0]), static_cast<float>(secondaryVertex[1]), static_cast<float>(secondaryVertex[2])};
    const float dcaDaughters = std::sqrt(fitter2Body.getChi2AtPCACandidate());
    if (twoBody.useSelections && (dcaDaughters > twoBody.maxDcaDaughters || RecoDecay::cpa(primaryVertex, decayVertex, momentum) < twoBody.minCosPA)) {
      return;
    }

    std::array<float, 2> dcaInfo{};
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, heTrackCov, 2.f, fitter2Body.getMatCorrType(), &dcaInfo);
    const float dcaHe = dcaInfo[0];
    o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, piTrackCov, 2.f, fitter2Body.getMatCorrType(), &dcaInfo);
    const float dcaPi = dcaInfo[0];
    if (twoBody.useSelections && (std::abs(dcaHe) < twoBody.minDcaHeToPV || std::abs(dcaPi) < twoBody.minDcaPiToPV)) {
      return;
    }

    const bool heliumPID = heTrack.pidForTracking() == o2::track::PID::Helium3 || heTrack.pidForTracking() == o2::track::PID::Alpha;
    const float tpcMomentumHe = heliumPID && twoBody.compensatePIDinTracking ? heTrack.tpcInnerParam() / 2.f : heTrack.tpcInnerParam();
    float tofMass = 0.f;
    if (heTrack.hasTOF()) {
      float beta = o2::pid::tof::Beta::GetBeta(heTrack);
      beta = std::clamp(beta, 1.e-4f, 1.f - 1.e-6f);
      tofMass = 2.f * tpcMomentumHe * std::sqrt(1.f / (beta * beta) - 1.f);
    }
    uint8_t flags = static_cast<uint8_t>((heTrack.pidForTracking() & 0xf) << 4);
    flags |= static_cast<uint8_t>(piTrack.pidForTracking() & 0xf);

    dataHypCands(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                 collision.posX(), collision.posY(), collision.posZ(),
                 runNumber, heTrack.sign() > 0,
                 std::hypot(heMomentum[0], heMomentum[1]), std::atan2(heMomentum[1], heMomentum[0]), RecoDecay::eta(heMomentum),
                 std::hypot(piMomentum[0], piMomentum[1]), std::atan2(piMomentum[1], piMomentum[0]), RecoDecay::eta(piMomentum),
                 secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                 dcaDaughters, dcaHe, dcaPi,
                 nSigmaHe3(heTrack), heTrack.tpcNClsFound(), piTrack.tpcNClsFound(),
                 static_cast<int16_t>(heTrack.tpcNClsFindable()) - heTrack.tpcNClsFindableMinusPID(),
                 static_cast<int16_t>(piTrack.tpcNClsFindable()) - piTrack.tpcNClsFindableMinusPID(),
                 heTrack.tpcNClsCrossedRows(), piTrack.tpcNClsCrossedRows(),
                 tpcMomentumHe, piTrack.tpcInnerParam(), heTrack.tpcSignal(), piTrack.tpcSignal(),
                 heTrack.tpcChi2NCl(), heTrack.itsChi2NCl(), piTrack.itsChi2NCl(), tofMass,
                 heTrack.itsClusterSizes(), piTrack.itsClusterSizes(), flags, static_cast<int>(trackedClSize));
  }

  template <typename TCollision, typename TTrack>
  double deuteronTOFNSigma(TCollision const& collision, TTrack const& track)
  {
    if (!track.has_collision() || !track.hasTOF()) {
      return -999.;
    }
    auto originalCollision = track.template collision_as<Collisions>();
    return deuteronTOFPID.GetTOFNSigma(tofResponse, track, originalCollision, collision);
  }

  void fillThreeBodyTables()
  {
    const auto& candidate = builder3Body.decay3body;
    vtx3BodyDatas(candidate.sign,
                  candidate.mass, candidate.massV0,
                  candidate.position[0], candidate.position[1], candidate.position[2],
                  candidate.momentum[0], candidate.momentum[1], candidate.momentum[2],
                  candidate.chi2, candidate.trackedClSize,
                  candidate.momProton[0], candidate.momProton[1], candidate.momProton[2],
                  candidate.momPion[0], candidate.momPion[1], candidate.momPion[2],
                  candidate.momDeuteron[0], candidate.momDeuteron[1], candidate.momDeuteron[2],
                  candidate.trackDCAxyToPV[0], candidate.trackDCAxyToPV[1], candidate.trackDCAxyToPV[2],
                  candidate.trackDCAToPV[0], candidate.trackDCAToPV[1], candidate.trackDCAToPV[2],
                  candidate.trackDCAxyToPVprop[0], candidate.trackDCAxyToPVprop[1], candidate.trackDCAxyToPVprop[2],
                  candidate.trackDCAToPVprop[0], candidate.trackDCAToPVprop[1], candidate.trackDCAToPVprop[2],
                  candidate.daughterDCAtoSV[0], candidate.daughterDCAtoSV[1], candidate.daughterDCAtoSV[2],
                  candidate.daughterDCAtoSVaverage, candidate.cosPA, candidate.ctau,
                  candidate.tpcNsigma[0], candidate.tpcNsigma[1], candidate.tpcNsigma[2], candidate.tpcNsigma[3],
                  candidate.tofNsigmaDeuteron,
                  candidate.averageITSClSize[0], candidate.averageITSClSize[1], candidate.averageITSClSize[2],
                  candidate.tpcNCl[0], candidate.tpcNCl[1], candidate.tpcNCl[2], candidate.pidForTrackingDeuteron);
    vtx3BodyCovs(candidate.covProton, candidate.covPion, candidate.covDeuteron, candidate.covariance);
  }

  void processData(Collisions const& collisions,
                   aod::V0s const& /*v0s*/,
                   aod::Decay3Bodys const& /*decay3Bodys*/,
                   aod::TrackedV0s const& trackedV0s,
                   aod::Tracked3Bodys const& tracked3Bodys,
                   Tracks const& /*tracks*/,
                   aod::BCsWithTimestamps const&)
  {
    selectCollisions(collisions);

    for (const auto& trackedV0 : trackedV0s) {
      const auto v0 = trackedV0.v0_as<aod::V0s>();
      if (v0.collisionId() < 0 || !goodCollision[v0.collisionId()] || (skimmedProcessing && !zorroDecision[v0.collisionId()][kHe])) {
        continue;
      }
      const auto collision = v0.collision_as<Collisions>();
      const auto positiveTrack = v0.posTrack_as<Tracks>();
      const auto negativeTrack = v0.negTrack_as<Tracks>();
      const float nSigmaPositive = nSigmaHe3(positiveTrack);
      const float nSigmaNegative = nSigmaHe3(negativeTrack);
      const bool positiveTrackedAsHe = positiveTrack.pidForTracking() == o2::track::PID::Helium3 || positiveTrack.pidForTracking() == o2::track::PID::Alpha;
      const bool negativeTrackedAsHe = negativeTrack.pidForTracking() == o2::track::PID::Helium3 || negativeTrack.pidForTracking() == o2::track::PID::Alpha;
      const bool positiveIsHe = positiveTrackedAsHe != negativeTrackedAsHe ? positiveTrackedAsHe : nSigmaPositive > nSigmaNegative;
      const float selectedNSigmaHe = positiveIsHe ? nSigmaPositive : nSigmaNegative;
      if (twoBody.useSelections && selectedNSigmaHe < twoBody.minNSigmaHe) {
        continue;
      }
      const auto heTrack = positiveIsHe ? positiveTrack : negativeTrack;
      const auto piTrack = positiveIsHe ? negativeTrack : positiveTrack;
      if (twoBody.useSelections &&
          (std::abs(heTrack.eta()) > twoBody.maxEtaDaughters || std::abs(piTrack.eta()) > twoBody.maxEtaDaughters ||
           heTrack.tpcNClsFound() < twoBody.minTPCNClsHe || piTrack.tpcNClsFound() < twoBody.minTPCNClsPi ||
           heTrack.tpcNClsCrossedRows() < twoBody.minTPCCrossedRowsHe || piTrack.tpcNClsCrossedRows() < twoBody.minTPCCrossedRowsPi)) {
        continue;
      }
      const bool heliumPID = heTrack.pidForTracking() == o2::track::PID::Helium3 || heTrack.pidForTracking() == o2::track::PID::Alpha;
      const float tpcMomentumHe = heliumPID && twoBody.compensatePIDinTracking ? heTrack.tpcInnerParam() / 2.f : heTrack.tpcInnerParam();
      if (twoBody.useSelections && tpcMomentumHe < twoBody.minTPCInnerParamHe) {
        continue;
      }
      buildTwoBody(heTrack, piTrack, collision, trackedV0.itsClsSize());
    }

    for (const auto& tracked3Body : tracked3Bodys) {
      const auto decay3Body = tracked3Body.decay3Body_as<aod::Decay3Bodys>();
      if (decay3Body.collisionId() < 0 || !goodCollision[decay3Body.collisionId()] || (skimmedProcessing && !zorroDecision[decay3Body.collisionId()][kTracked3Body])) {
        continue;
      }
      const auto collision = decay3Body.collision_as<Collisions>();
      const auto trackPositive = decay3Body.track0_as<Tracks>();
      const auto trackNegative = decay3Body.track1_as<Tracks>();
      const auto trackDeuteron = decay3Body.track2_as<Tracks>();
      const auto trackProton = trackDeuteron.sign() > 0 ? trackPositive : trackNegative;
      const auto trackPion = trackDeuteron.sign() > 0 ? trackNegative : trackPositive;
      if (builder3Body.buildDecay3BodyCandidate(collision, trackProton, trackPion, trackDeuteron,
                                                decay3Body.globalIndex(), deuteronTOFNSigma(collision, trackDeuteron), tracked3Body.itsClsSize(),
                                                threeBody.useKFParticle, threeBody.setTopologicalConstraint,
                                                threeBody.useSelections, threeBody.useChi2Selection, threeBody.useTPCforPion,
                                                threeBody.acceptTPCOnly, threeBody.askOnlyITSMatch, threeBody.calculateCovariance)) {
        fillThreeBodyTables();
      }
    }
  }
  PROCESS_SWITCH(TrackedHypertritonRecoTask, processData, "Process tracked hypertriton candidates in data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<TrackedHypertritonRecoTask>(cfgc)};
}
