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
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>

#include <KFParticle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{
using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::EvTimeTOFFT0>;
using Tracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                         aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe,
                         aod::TOFSignal, aod::TOFEvTime, aod::EvTimeTOFFT0ForTrack>;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::EvTimeTOFFT0, aod::McCollisionLabels>;
using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                           aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe,
                           aod::TOFSignal, aod::TOFEvTime, aod::EvTimeTOFFT0ForTrack, aod::McTrackLabels>;

constexpr std::array<double, 6> betheBlochDefault{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32};
constexpr double UseCCDBMagneticFieldThreshold = -990.;
const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
const std::vector<std::string> particleName{"He3"};

enum ZorroTrigger : std::size_t {
  kHe = 0,
  kTracked3Body,
  kNZorroTriggers
};
} // namespace

struct TrackedHypertritonRecoTask {
  o2::common::core::MetadataHelper metadataInfo{};

  Produces<aod::DataHypCands> dataHypCands;
  Produces<aod::MCHypCands> mcHypCands;
  Produces<aod::Vtx3BodyDatas> vtx3BodyDatas;
  Produces<aod::Vtx3BodyCovs> vtx3BodyCovs;
  Produces<aod::Vtx3BodyTrackedInfo> vtx3BodyTrackedInfo;
  Produces<aod::McVtx3BodyDatas> mcVtx3BodyDatas;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> requireSel8{"requireSel8", true, "Require the standard sel8 event selection"};
  Configurable<float> maxAbsZvtx{"maxAbsZvtx", 10.f, "Maximum absolute primary-vertex z; negative disables the cut"};
  Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Enable Zorro accounting and trigger selection for skimmed data"};
  Configurable<std::string> zorroCCDBPath{"zorroCCDBPath", "EventFiltering/Zorro/", "Base path of the Zorro CCDB objects"};
  Configurable<int> zorroBCTolerance{"zorroBCTolerance", 100, "Zorro BC matching tolerance"};

  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"}; // o2-linter: disable=name/configurable (Keep the standard CCDB option name.)
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the magnetic-field object"};
  Configurable<double> bzInput{"bz", -999., "Magnetic field in kG; values below -990 use CCDB"}; // o2-linter: disable=name/configurable (Keep the established magnetic-field option name.)

  struct : ConfigurableGroup {
    std::string prefix = "mc";
    Configurable<bool> storeBackground{"storeBackground", true, "Store reconstructed MC candidates not matched to a hypertriton"};
  } mc;

  struct : ConfigurableGroup {
    std::string prefix = "twoBody";
    Configurable<bool> useKFParticle{"useKFParticle", false, "Use KFParticle to build the two-body candidate"};
    Configurable<bool> useSelections{"useSelections", false, "Apply the two-body candidate selections"};
    Configurable<bool> kfSetTopologicalConstraint{"kfSetTopologicalConstraint", false, "Set topological vertex constraint in case of KFParticle reconstruction"};
    Configurable<float> maxEtaDaughters{"maxEtaDaughters", 1.e10f, "Maximum absolute daughter eta"};
    Configurable<float> minTPCNClsHe{"minTPCNClsHe", -1.f, "Minimum He TPC clusters"};
    Configurable<float> minTPCNClsPi{"minTPCNClsPi", -1.f, "Minimum pion TPC clusters"};
    Configurable<float> minTPCCrossedRowsHe{"minTPCCrossedRowsHe", -1.f, "Minimum He TPC crossed rows"};
    Configurable<float> minTPCCrossedRowsPi{"minTPCCrossedRowsPi", -1.f, "Minimum pion TPC crossed rows"};
    Configurable<float> minTPCInnerParamHe{"minTPCInnerParamHe", -1.f, "Minimum He TPC rigidity"};
    Configurable<float> minPt{"minPt", -1.f, "Minimum candidate transverse momentum"};
    Configurable<float> massWindow{"massWindow", 1.e10f, "Half-width of the hypertriton mass window"};
    Configurable<float> maxChi2{"maxChi2", 1.e10f, "KFParticle: Maximum SV chi2, DCA fitter: Maximum DCA between daughters"};
    Configurable<float> minCosPA{"minCosPA", -2.f, "Minimum cosine of the pointing angle"};
    Configurable<float> minDcaHeToPV{"minDcaHeToPV", -1.f, "Minimum absolute He DCA to PV"};
    Configurable<float> minDcaPiToPV{"minDcaPiToPV", -1.f, "Minimum absolute pion DCA to PV"};
    Configurable<float> minNSigmaHe{"minNSigmaHe", -1.e10f, "Minimum He TPC n-sigma"};
    Configurable<bool> compensatePIDinTracking{"compensatePIDinTracking", true, "Correct the TPC inner parameter for charge-two PID in tracking"};
    Configurable<LabeledArray<double>> betheBlochParams{"betheBlochParams", {betheBlochDefault.data(), 1, 6, particleName, betheBlochParNames}, "TPC Bethe-Bloch parameters for He3"};
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
    Configurable<float> maxITSDCAxytrackToSV{"maxITSDCAxytrackToSV", 10.0, "Maximum distance of ITS matched track to SV in xy"};
    Configurable<float> maxITSDCAztrackToSV{"maxITSDCAztrackToSV", 10.0, "Maximum distance of ITS matched track to SV in z"};
  } threeBody;

  o2::vertexing::DCAFitterN<2> fitter2Body;
  o2::pwglf::decay3bodyBuilderHelper builder3Body;
  o2::pid::tof::TOFResoParamsV3 tofResponse;
  o2::aod::pidtofgeneric::TOFCalibConfig tofCalibConfig;
  o2::aod::pidtofgeneric::TofPidNewCollision<Tracks::iterator> deuteronTOFPID;
  o2::aod::pidtofgeneric::TofPidNewCollision<TracksMC::iterator> deuteronTOFPIDMC;

  int runNumber = 0;
  float bz = 0.f;
  std::array<float, 6> bbParamsHe{};
  std::vector<bool> goodCollision;
  std::vector<std::array<bool, kNZorroTriggers>> zorroDecision;
  std::vector<int> recoCollisionForMC;
  std::vector<bool> survivedMCEventSelection;

  struct v0Candidate {
    // daughter properties
    std::array<float, 3> momHelium{};
    std::array<float, 3> momPion{};
    std::array<float, 3> posHelium{};
    std::array<float, 3> posPion{};
    // vertex properties
    float mass;
    float chi2;
    float cosPA;
    std::array<float, 3> decayVertex{};
    std::array<float, 3> momentum{};
  };

  v0Candidate v0;

  struct TwoBodyMCInfo {
    float genPt = -1.f;
    float genPhi = -1.f;
    float genEta = -1.f;
    float genPtHe3 = -1.f;
    std::array<float, 3> genDecayVertex{-1.f, -1.f, -1.f};
    bool isSignal = false;
    bool isRecoMCCollision = false;
    bool survivedEventSelection = false;
    uint8_t fakeHeITSLayerMap = 0;
    int motherLabel = -1;
    int statusCode = 0;
  };

  struct ThreeBodyMCInfo {
    std::array<float, 3> genMomentum{-1.f, -1.f, -1.f};
    std::array<float, 3> genDecayVertex{-1.f, -1.f, -1.f};
    float genCt = -1.f;
    float genPhi = -1.f;
    float genEta = -1.f;
    float genRapidity = -1.f;
    float genMomentumProton = -1.f;
    float genMomentumPion = -1.f;
    float genMomentumDeuteron = -1.f;
    float genPtProton = -1.f;
    float genPtPion = -1.f;
    float genPtDeuteron = -1.f;
    bool isReco = true;
    int motherLabel = -1;
    int motherPdgCode = 0;
    int protonPdgCode = -1;
    int pionPdgCode = -1;
    int deuteronPdgCode = -1;
    bool isDeuteronPrimary = false;
    bool survivedEventSelection = false;
  };

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
    deuteronTOFPIDMC.SetPidType(o2::track::PID::Deuteron);
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
    bz = bzInput < UseCCDBMagneticFieldThreshold ? o2::base::Propagator::Instance()->getNominalBz() : bzInput;
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
  void selectCollisions(TCollision const& collisions, bool applyZorro)
  {
    goodCollision.assign(collisions.size(), false);
    zorroDecision.assign(collisions.size(), {});
    for (const auto& collision : collisions) {
      registry.fill(HIST("events"), 0.);
      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyZorro) {
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
      if (applyZorro) {
        for (size_t i = 0; i < kNZorroTriggers; ++i) {
          if (zorroDecision[collision.globalIndex()][i]) {
            registry.fill(HIST("zorroEvents"), static_cast<double>(i), 1.);
          }
        }
      }
    }
  }

  template <typename TCollisions>
  void prepareMCEventInformation(TCollisions const& collisions, aod::McCollisions const& mcCollisions)
  {
    recoCollisionForMC.assign(mcCollisions.size(), -1);
    survivedMCEventSelection.assign(mcCollisions.size(), false);
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const int mcCollisionId = collision.mcCollisionId();
      if (mcCollisionId < 0 || mcCollisionId >= static_cast<int>(mcCollisions.size())) {
        continue;
      }
      if (recoCollisionForMC[mcCollisionId] < 0) {
        recoCollisionForMC[mcCollisionId] = collision.globalIndex();
      }
      if (goodCollision[collision.globalIndex()]) {
        recoCollisionForMC[mcCollisionId] = collision.globalIndex();
        survivedMCEventSelection[mcCollisionId] = true;
      }
    }
  }

  template <typename TParticleA, typename TParticleB>
  int findCommonMother(TParticleA const& particleA, TParticleB const& particleB) const
  {
    if (!particleA.has_mothers() || !particleB.has_mothers()) {
      return -1;
    }
    for (const auto& motherA : particleA.template mothers_as<aod::McParticles>()) {
      for (const auto& motherB : particleB.template mothers_as<aod::McParticles>()) {
        if (motherA.globalIndex() == motherB.globalIndex()) {
          return motherA.globalIndex();
        }
      }
    }
    return -1;
  }

  template <typename TParticleA, typename TParticleB, typename TParticleC>
  int findCommonMother(TParticleA const& particleA, TParticleB const& particleB, TParticleC const& particleC) const
  {
    if (!particleA.has_mothers() || !particleB.has_mothers() || !particleC.has_mothers()) {
      return -1;
    }
    for (const auto& motherA : particleA.template mothers_as<aod::McParticles>()) {
      for (const auto& motherB : particleB.template mothers_as<aod::McParticles>()) {
        if (motherA.globalIndex() != motherB.globalIndex()) {
          continue;
        }
        for (const auto& motherC : particleC.template mothers_as<aod::McParticles>()) {
          if (motherA.globalIndex() == motherC.globalIndex()) {
            return motherA.globalIndex();
          }
        }
      }
    }
    return -1;
  }

  template <typename TTrack, typename TCollision>
  TwoBodyMCInfo getTwoBodyMCInfo(TTrack const& heTrack, TTrack const& piTrack, TCollision const& collision, aod::McParticles const& mcParticles) const
  {
    TwoBodyMCInfo info;
    if (collision.has_mcCollision()) {
      info.survivedEventSelection = survivedMCEventSelection[collision.mcCollisionId()];
    }
    if (!heTrack.has_mcParticle() || !piTrack.has_mcParticle()) {
      return info;
    }

    const auto mcHe = heTrack.template mcParticle_as<aod::McParticles>();
    const auto mcPi = piTrack.template mcParticle_as<aod::McParticles>();
    info.fakeHeITSLayerMap = heTrack.mcMask() & 0x7f;
    const int motherLabel = findCommonMother(mcHe, mcPi);
    if (motherLabel < 0) {
      return info;
    }
    const auto mother = mcParticles.rawIteratorAt(motherLabel);
    const int sign = mother.pdgCode() > 0 ? 1 : -1;
    if (std::abs(mother.pdgCode()) != constants::physics::Pdg::kHyperTriton ||
        mcHe.pdgCode() != sign * constants::physics::Pdg::kHelium3 || mcPi.pdgCode() != -sign * PDG_t::kPiPlus) {
      return info;
    }

    info.isSignal = true;
    info.motherLabel = motherLabel;
    info.genPt = sign * mother.pt();
    info.genPhi = mother.phi();
    info.genEta = mother.eta();
    info.genPtHe3 = mcHe.pt();
    info.genDecayVertex = {mcHe.vx() - mother.vx(), mcHe.vy() - mother.vy(), mcHe.vz() - mother.vz()};
    info.statusCode = mcHe.getProcess();
    if (mother.mcCollisionId() >= 0 && mother.mcCollisionId() < static_cast<int>(recoCollisionForMC.size())) {
      info.isRecoMCCollision = recoCollisionForMC[mother.mcCollisionId()] >= 0;
      info.survivedEventSelection = survivedMCEventSelection[mother.mcCollisionId()];
    }
    return info;
  }

  template <typename TTrack, typename TCollision>
  ThreeBodyMCInfo getThreeBodyMCInfo(TTrack const& trackProton, TTrack const& trackPion, TTrack const& trackDeuteron,
                                     TCollision const& collision, aod::McParticles const& mcParticles) const
  {
    ThreeBodyMCInfo info;
    if (collision.has_mcCollision()) {
      info.survivedEventSelection = survivedMCEventSelection[collision.mcCollisionId()];
    }
    if (!trackProton.has_mcParticle() || !trackPion.has_mcParticle() || !trackDeuteron.has_mcParticle()) {
      return info;
    }

    const auto mcProton = trackProton.template mcParticle_as<aod::McParticles>();
    const auto mcPion = trackPion.template mcParticle_as<aod::McParticles>();
    const auto mcDeuteron = trackDeuteron.template mcParticle_as<aod::McParticles>();
    info.protonPdgCode = mcProton.pdgCode();
    info.pionPdgCode = mcPion.pdgCode();
    info.deuteronPdgCode = mcDeuteron.pdgCode();
    info.isDeuteronPrimary = mcDeuteron.isPhysicalPrimary();
    info.genMomentumProton = mcProton.p();
    info.genMomentumPion = mcPion.p();
    info.genMomentumDeuteron = mcDeuteron.p();
    info.genPtProton = mcProton.pt();
    info.genPtPion = mcPion.pt();
    info.genPtDeuteron = mcDeuteron.pt();

    const int motherLabel = findCommonMother(mcProton, mcPion, mcDeuteron);
    if (motherLabel < 0) {
      return info;
    }
    const auto mother = mcParticles.rawIteratorAt(motherLabel);
    const int sign = mother.pdgCode() > 0 ? 1 : -1;
    if (std::abs(mother.pdgCode()) != constants::physics::Pdg::kHyperTriton ||
        mcProton.pdgCode() != sign * PDG_t::kProton || mcPion.pdgCode() != -sign * PDG_t::kPiPlus ||
        mcDeuteron.pdgCode() != sign * constants::physics::Pdg::kDeuteron) {
      return info;
    }

    info.motherLabel = motherLabel;
    info.motherPdgCode = mother.pdgCode();
    info.genMomentum = {mother.px(), mother.py(), mother.pz()};
    info.genDecayVertex = {mcProton.vx(), mcProton.vy(), mcProton.vz()};
    info.genCt = RecoDecay::sqrtSumOfSquares(mcProton.vx() - mother.vx(), mcProton.vy() - mother.vy(), mcProton.vz() - mother.vz()) * constants::physics::MassHyperTriton / (mother.p() + 1.e-10f);
    info.genPhi = mother.phi();
    info.genEta = mother.eta();
    info.genRapidity = mother.y();
    if (mother.mcCollisionId() >= 0 && mother.mcCollisionId() < static_cast<int>(survivedMCEventSelection.size())) {
      info.survivedEventSelection = survivedMCEventSelection[mother.mcCollisionId()];
    }
    return info;
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

  template <typename TCollision, typename TTrack, typename TTrackParCov>
  void fit2BodyWithKF(TCollision const& collision,
                      TTrack const& trackHelium,
                      TTrack const& trackPion,
                      TTrackParCov const& trackHeliumCov,
                      TTrackParCov const& trackPionCov)
  {
    // initialise KF primary vertex
    KFVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle kfpv(kfpVertex);

    // create KFParticle objects
    KFParticle kfpHelium, kfpPion;
    // helium
    std::array<float, 3> xyz, pxpypz;
    float xyzpxpypz[6];
    trackHeliumCov.getPxPyPzGlo(pxpypz);
    trackHeliumCov.getXYZGlo(xyz);
    for (int i = 0; i < 3; ++i) {
      xyzpxpypz[i] = xyz[i];
      xyzpxpypz[i + 3] = pxpypz[i] * 2;
    }
    std::array<float, 21> cv{};
    trackHeliumCov.getCovXYZPxPyPzGlo(cv);
    KFParticle kfHelium;
    kfHelium.Create(xyzpxpypz, cv.data(), trackHelium.sign() * 2, constants::physics::MassHelium3);
    // pion
    kfpPion = createKFParticleFromTrackParCov(trackPionCov, trackPion.sign(), constants::physics::MassPionCharged);

    // construct V0 vertex
    KFParticle KFV0;
    int nDaughtersV0 = 2;
    const KFParticle* DaughtersV0[2] = {&kfpHelium, &kfpPion};
    KFV0.SetConstructMethod(2);
    try {
      KFV0.Construct(DaughtersV0, nDaughtersV0);
    } catch (std::runtime_error& e) {
      LOG(debug) << "Failed to create V0 vertex." << e.what();
      return;
    }

    // topological constraint
    if (twoBody.kfSetTopologicalConstraint) {
      KFV0.SetProductionVertex(kfpv);
      KFV0.TransportToDecayVertex();
    }

    // get vertex position and momentum
    v0.decayVertex[0] = KFV0.GetX();
    v0.decayVertex[1] = KFV0.GetY();
    v0.decayVertex[2] = KFV0.GetZ();
    v0.momentum[0] = KFV0.GetPx();
    v0.momentum[1] = KFV0.GetPy();
    v0.momentum[2] = KFV0.GetPz();

    // transport all daughter tracks to hypertriton vertex
    // float position[3];
    // for (int i; i < 3; i++) {
    //   position[i] = v0.decayVertex[i];
    // }
    kfpHelium.TransportToPoint(v0.decayVertex.data());
    kfpPion.TransportToPoint(v0.decayVertex.data());

    // daughter positions
    v0.posHelium[0] = kfpHelium.GetX();
    v0.posHelium[1] = kfpHelium.GetY();
    v0.posHelium[2] = kfpHelium.GetZ();
    v0.posPion[0] = kfpPion.GetX();
    v0.posPion[1] = kfpPion.GetY();
    v0.posPion[2] = kfpPion.GetZ();

    // daughter momenta
    v0.momHelium[0] = kfpHelium.GetPx();
    v0.momHelium[1] = kfpHelium.GetPy();
    v0.momHelium[2] = kfpHelium.GetPz();
    v0.momPion[0] = kfpPion.GetPx();
    v0.momPion[1] = kfpPion.GetPy();
    v0.momPion[2] = kfpPion.GetPz();

    // candidate mass
    float mass, massErr;
    KFV0.GetMass(mass, massErr);
    v0.mass = mass;

    // vertex chi2
    v0.chi2 = KFV0.GetChi2() / KFV0.GetNDF();
  }

  template <typename TTrackParCov>
  void fit2bodyWithDCAFitter(TTrackParCov const& trackHeliumCov,
                             TTrackParCov const& trackPionCov)
  {
    int nCandidates = 0;
    try {
      nCandidates = fitter2Body.process(trackHeliumCov, trackPionCov);
    } catch (...) {
      LOG(error) << "Exception while fitting a tracked two-body candidate";
      return;
    }
    if (nCandidates == 0) {
      return;
    }

    // get daughter momenta
    fitter2Body.getTrack(0).getPxPyPzGlo(v0.momHelium);
    fitter2Body.getTrack(1).getPxPyPzGlo(v0.momPion);

    for (std::size_t i = 0; i < v0.momHelium.size(); ++i) {
      v0.momHelium[i] *= 2.f;
    }

    // compute candidate momentum
    v0.momentum = {v0.momHelium[0] + v0.momPion[0], v0.momHelium[1] + v0.momPion[1], v0.momHelium[2] + v0.momPion[2]};

    // compute candidate mass
    const float heMomentum2 = RecoDecay::sumOfSquares(v0.momHelium[0], v0.momHelium[1], v0.momHelium[2]);
    const float piMomentum2 = RecoDecay::sumOfSquares(v0.momPion[0], v0.momPion[1], v0.momPion[2]);
    const float candidateMomentum2 = RecoDecay::sumOfSquares(v0.momentum[0], v0.momentum[1], v0.momentum[2]);
    const float heEnergy = std::sqrt(heMomentum2 + constants::physics::MassHelium3 * constants::physics::MassHelium3);
    const float piEnergy = std::sqrt(piMomentum2 + constants::physics::MassPionCharged * constants::physics::MassPionCharged);
    v0.mass = std::sqrt((heEnergy + piEnergy) * (heEnergy + piEnergy) - candidateMomentum2);

    // get SV position
    const auto& secondaryVertex = fitter2Body.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      v0.decayVertex[i] = secondaryVertex[i];
    }
    v0.chi2 = std::sqrt(fitter2Body.getChi2AtPCACandidate());
  }

  template <typename TTrack, typename TCollision, typename TFillCandidate>
  void buildTwoBody(TTrack const& heTrack, TTrack const& piTrack, TCollision const& collision, float trackedClSize, TFillCandidate const& fillCandidate)
  {
    const std::array<float, 3> primaryVertex{collision.posX(), collision.posY(), collision.posZ()};

    auto heTrackCov = getTrackParCov(heTrack);
    auto piTrackCov = getTrackParCov(piTrack);

    if (twoBody.useKFParticle) {
      fit2BodyWithKF(collision, heTrack, piTrack, heTrackCov, piTrackCov);
    } else {
      fit2bodyWithDCAFitter(heTrackCov, piTrackCov);
    }

    v0.cosPA = RecoDecay::cpa(primaryVertex, v0.decayVertex, v0.momentum);
    if (twoBody.useSelections &&
        (std::hypot(v0.momentum[0], v0.momentum[1]) < twoBody.minPt ||
         std::abs(v0.mass - constants::physics::MassHyperTriton) > twoBody.massWindow ||
         v0.chi2 > twoBody.maxChi2 ||
         v0.cosPA < twoBody.minCosPA)) {
      return;
    }

    // Do propagation with Propagator including material interactions in all cases (KFParticle propagation would not include material)
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
    auto flags = static_cast<uint8_t>((heTrack.pidForTracking() & 0xf) << 4);
    flags |= static_cast<uint8_t>(piTrack.pidForTracking() & 0xf);

    fillCandidate(collision.centFT0A(), collision.centFT0C(), collision.centFT0M(),
                  collision.posX(), collision.posY(), collision.posZ(),
                  runNumber, heTrack.sign() > 0,
                  std::hypot(v0.momHelium[0], v0.momHelium[1]), std::atan2(v0.momHelium[1], v0.momHelium[0]), RecoDecay::eta(v0.momHelium),
                  std::hypot(v0.momPion[0], v0.momPion[1]), std::atan2(v0.momPion[1], v0.momPion[0]), RecoDecay::eta(v0.momPion),
                  v0.decayVertex[0], v0.decayVertex[1], v0.decayVertex[2],
                  v0.chi2, dcaHe, dcaPi,
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

  template <typename TTrack, typename TCollision>
  double deuteronTOFNSigmaMC(TCollision const& collision, TTrack const& track)
  {
    if (!track.has_collision() || !track.hasTOF()) {
      return -999.;
    }
    auto originalCollision = track.template collision_as<CollisionsMC>();
    return deuteronTOFPIDMC.GetTOFNSigma(tofResponse, track, originalCollision, collision);
  }

  void fillThreeBodyTables()
  {
    const auto& candidate = builder3Body.decay3body;
    vtx3BodyDatas(static_cast<float>(candidate.sign),
                  candidate.mass, candidate.massV0,
                  candidate.position[0], candidate.position[1], candidate.position[2],
                  candidate.momentum[0], candidate.momentum[1], candidate.momentum[2],
                  candidate.chi2, candidate.trackedClSize,
                  candidate.momProton[0], candidate.momProton[1], candidate.momProton[2],
                  candidate.momPion[0], candidate.momPion[1], candidate.momPion[2],
                  candidate.momDeuteron[0], candidate.momDeuteron[1], candidate.momDeuteron[2],
                  candidate.xProton, candidate.xPion, candidate.xDeuteron,
                  candidate.trackDCAxyToPV[0], candidate.trackDCAxyToPV[1], candidate.trackDCAxyToPV[2],
                  candidate.trackDCAToPV[0], candidate.trackDCAToPV[1], candidate.trackDCAToPV[2],
                  candidate.trackDCAxyToPVprop[0], candidate.trackDCAxyToPVprop[1], candidate.trackDCAxyToPVprop[2],
                  candidate.trackDCAToPVprop[0], candidate.trackDCAToPVprop[1], candidate.trackDCAToPVprop[2],
                  candidate.daughterDCAtoSV[0], candidate.daughterDCAtoSV[1], candidate.daughterDCAtoSV[2],
                  candidate.daughterDCAtoSVaverage, candidate.cosPA, candidate.ctau,
                  candidate.tpcNsigma[0], candidate.tpcNsigma[1], candidate.tpcNsigma[2], candidate.tpcNsigma[3],
                  static_cast<float>(candidate.tofNsigmaDeuteron),
                  candidate.averageITSClSize[0], candidate.averageITSClSize[1], candidate.averageITSClSize[2],
                  static_cast<int>(candidate.tpcNCl[0]), static_cast<int>(candidate.tpcNCl[1]), static_cast<int>(candidate.tpcNCl[2]),
                  static_cast<uint32_t>(candidate.pidForTrackingDeuteron));
    vtx3BodyCovs(candidate.covProton, candidate.covPion, candidate.covDeuteron, candidate.covariance);
    vtx3BodyTrackedInfo(candidate.itsTrackDCAToSV[0], candidate.itsTrackDCAToSV[1]);
  }

  void fillThreeBodyMCTable(ThreeBodyMCInfo const& info)
  {
    const auto& candidate = builder3Body.decay3body;
    mcVtx3BodyDatas(static_cast<float>(candidate.sign),
                    candidate.mass, candidate.massV0,
                    candidate.position[0], candidate.position[1], candidate.position[2],
                    candidate.momentum[0], candidate.momentum[1], candidate.momentum[2],
                    candidate.chi2, candidate.trackedClSize,
                    candidate.momProton[0], candidate.momProton[1], candidate.momProton[2],
                    candidate.momPion[0], candidate.momPion[1], candidate.momPion[2],
                    candidate.momDeuteron[0], candidate.momDeuteron[1], candidate.momDeuteron[2],
                    candidate.xProton, candidate.xPion, candidate.xDeuteron,
                    candidate.trackDCAxyToPV[0], candidate.trackDCAxyToPV[1], candidate.trackDCAxyToPV[2],
                    candidate.trackDCAToPV[0], candidate.trackDCAToPV[1], candidate.trackDCAToPV[2],
                    candidate.trackDCAxyToPVprop[0], candidate.trackDCAxyToPVprop[1], candidate.trackDCAxyToPVprop[2],
                    candidate.trackDCAToPVprop[0], candidate.trackDCAToPVprop[1], candidate.trackDCAToPVprop[2],
                    candidate.daughterDCAtoSV[0], candidate.daughterDCAtoSV[1], candidate.daughterDCAtoSV[2],
                    candidate.daughterDCAtoSVaverage, candidate.cosPA, candidate.ctau,
                    candidate.tpcNsigma[0], candidate.tpcNsigma[1], candidate.tpcNsigma[2], candidate.tpcNsigma[3],
                    static_cast<float>(candidate.tofNsigmaDeuteron),
                    candidate.averageITSClSize[0], candidate.averageITSClSize[1], candidate.averageITSClSize[2],
                    static_cast<int>(candidate.tpcNCl[0]), static_cast<int>(candidate.tpcNCl[1]), static_cast<int>(candidate.tpcNCl[2]),
                    static_cast<uint32_t>(candidate.pidForTrackingDeuteron),
                    info.genMomentum[0], info.genMomentum[1], info.genMomentum[2],
                    info.genDecayVertex[0], info.genDecayVertex[1], info.genDecayVertex[2],
                    info.genCt, info.genPhi, info.genEta, info.genRapidity,
                    info.genMomentumProton, info.genMomentumPion, info.genMomentumDeuteron,
                    info.genPtProton, info.genPtPion, info.genPtDeuteron,
                    static_cast<int>(info.isReco), info.motherLabel, info.motherPdgCode,
                    info.protonPdgCode, info.pionPdgCode, info.deuteronPdgCode,
                    info.isDeuteronPrimary, static_cast<int>(info.survivedEventSelection));
  }

  void fillGeneratedThreeBodyMCTable(ThreeBodyMCInfo const& info)
  {
    mcVtx3BodyDatas(-1.f,
                    -1.f, -1.f,       // mass, massV0
                    -1.f, -1.f, -1.f, // position
                    -1.f, -1.f, -1.f, // momentum
                    -1.f, -1.f,       // chi2, trackedClSize
                    -1.f, -1.f, -1.f, // proton momentum
                    -1.f, -1.f, -1.f, // pion momentum
                    -1.f, -1.f, -1.f, // deuteron momentum
                    -1.f, -1.f, -1.f, // daughter x at inner update
                    -1.f, -1.f, -1.f, // track DCAxy to PV
                    -1.f, -1.f, -1.f, // track DCA to PV
                    -1.f, -1.f, -1.f, // propagated track DCAxy to PV
                    -1.f, -1.f, -1.f, // propagated track DCA to PV
                    -1.f, -1.f, -1.f, // daughter DCA to SV
                    -1.f, -1.f, -1.f, // average daughter DCA, cosPA, ctau
                    -1.f, -1.f, -1.f, -1.f,
                    -1.f,
                    -1.f, -1.f, -1.f,
                    -1, -1, -1, std::numeric_limits<uint32_t>::max(),
                    info.genMomentum[0], info.genMomentum[1], info.genMomentum[2],
                    info.genDecayVertex[0], info.genDecayVertex[1], info.genDecayVertex[2],
                    info.genCt, info.genPhi, info.genEta, info.genRapidity,
                    info.genMomentumProton, info.genMomentumPion, info.genMomentumDeuteron,
                    info.genPtProton, info.genPtPion, info.genPtDeuteron,
                    0, info.motherLabel, info.motherPdgCode,
                    info.protonPdgCode, info.pionPdgCode, info.deuteronPdgCode,
                    info.isDeuteronPrimary, static_cast<int>(info.survivedEventSelection));
  }

  void processData(Collisions const& collisions,
                   aod::V0s const& /*v0s*/,
                   aod::Decay3Bodys const& /*decay3Bodys*/,
                   aod::TrackedV0s const& trackedV0s,
                   aod::Tracked3Bodys const& tracked3Bodys,
                   Tracks const& /*tracks*/,
                   aod::BCsWithTimestamps const&)
  {
    selectCollisions(collisions, skimmedProcessing);

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
      buildTwoBody(heTrack, piTrack, collision, trackedV0.itsClsSize(), [&](auto... values) { dataHypCands(values...); });
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

      // get DCA of ITS track to SV
      const auto itsTrack = tracked3Body.itsTrack_as<Tracks>();
      auto itsTrackParCov = getTrackParCov(itsTrack);
      std::array<float, 2> dcaInfoItsTrack{};
      o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, itsTrackParCov, 2.f, fitter2Body.getMatCorrType(), &dcaInfoItsTrack);
      builder3Body.decay3body.itsTrackDCAToSV[0] = dcaInfoItsTrack[0];
      builder3Body.decay3body.itsTrackDCAToSV[1] = dcaInfoItsTrack[1];
      if (threeBody.useSelections && (builder3Body.decay3body.itsTrackDCAToSV[0] > threeBody.maxITSDCAxytrackToSV || builder3Body.decay3body.itsTrackDCAToSV[1] > threeBody.maxITSDCAztrackToSV)) {
        continue;
      }

      if (builder3Body.buildDecay3BodyCandidate(collision, trackProton, trackPion, trackDeuteron,
                                                decay3Body.globalIndex(), deuteronTOFNSigma(collision, trackDeuteron), tracked3Body.itsClsSize(),
                                                threeBody.useKFParticle, threeBody.setTopologicalConstraint,
                                                threeBody.useSelections, threeBody.useChi2Selection, threeBody.useTPCforPion,
                                                threeBody.acceptTPCOnly, threeBody.askOnlyITSMatch, threeBody.calculateCovariance)) {
        fillThreeBodyTables();
      }
    }
  }

  void processMC(CollisionsMC const& collisions,
                 aod::McCollisions const& mcCollisions,
                 aod::V0s const& /*v0s*/,
                 aod::Decay3Bodys const& /*decay3Bodys*/,
                 aod::TrackedV0s const& trackedV0s,
                 aod::Tracked3Bodys const& tracked3Bodys,
                 TracksMC const& /*tracks*/,
                 aod::BCsWithTimestamps const&,
                 aod::McParticles const& mcParticles)
  {
    selectCollisions(collisions, false);
    prepareMCEventInformation(collisions, mcCollisions);
    std::vector<bool> reconstructedTwoBody(mcParticles.size(), false);
    std::vector<bool> reconstructedThreeBody(mcParticles.size(), false);

    for (const auto& trackedV0 : trackedV0s) {
      const auto v0 = trackedV0.v0_as<aod::V0s>();
      if (v0.collisionId() < 0 || !goodCollision[v0.collisionId()]) {
        continue;
      }
      const auto collision = v0.collision_as<CollisionsMC>();
      const auto positiveTrack = v0.posTrack_as<TracksMC>();
      const auto negativeTrack = v0.negTrack_as<TracksMC>();
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

      const auto mcInfo = getTwoBodyMCInfo(heTrack, piTrack, collision, mcParticles);
      if (!mcInfo.isSignal && !mc.storeBackground) {
        continue;
      }
      buildTwoBody(heTrack, piTrack, collision, trackedV0.itsClsSize(), [&](auto... values) {
        mcHypCands(values...,
                   mcInfo.genPt, mcInfo.genPhi, mcInfo.genEta, mcInfo.genPtHe3,
                   mcInfo.genDecayVertex[0], mcInfo.genDecayVertex[1], mcInfo.genDecayVertex[2],
                   true, mcInfo.fakeHeITSLayerMap, mcInfo.isSignal,
                   mcInfo.isRecoMCCollision, mcInfo.survivedEventSelection, true, mcInfo.statusCode);
        if (mcInfo.isSignal && mcInfo.motherLabel >= 0) {
          reconstructedTwoBody[mcInfo.motherLabel] = true;
        }
      });
    }

    for (const auto& tracked3Body : tracked3Bodys) {
      const auto decay3Body = tracked3Body.decay3Body_as<aod::Decay3Bodys>();
      if (decay3Body.collisionId() < 0 || !goodCollision[decay3Body.collisionId()]) {
        continue;
      }
      const auto collision = decay3Body.collision_as<CollisionsMC>();
      const auto trackPositive = decay3Body.track0_as<TracksMC>();
      const auto trackNegative = decay3Body.track1_as<TracksMC>();
      const auto trackDeuteron = decay3Body.track2_as<TracksMC>();
      const auto trackProton = trackDeuteron.sign() > 0 ? trackPositive : trackNegative;
      const auto trackPion = trackDeuteron.sign() > 0 ? trackNegative : trackPositive;
      if (!builder3Body.buildDecay3BodyCandidate(collision, trackProton, trackPion, trackDeuteron,
                                                 decay3Body.globalIndex(), deuteronTOFNSigmaMC(collision, trackDeuteron), tracked3Body.itsClsSize(),
                                                 threeBody.useKFParticle, threeBody.setTopologicalConstraint,
                                                 threeBody.useSelections, threeBody.useChi2Selection, threeBody.useTPCforPion,
                                                 threeBody.acceptTPCOnly, threeBody.askOnlyITSMatch, threeBody.calculateCovariance)) {
        continue;
      }
      const auto mcInfo = getThreeBodyMCInfo(trackProton, trackPion, trackDeuteron, collision, mcParticles);
      if (mcInfo.motherLabel < 0 && !mc.storeBackground) {
        continue;
      }
      fillThreeBodyMCTable(mcInfo);
      if (mcInfo.motherLabel >= 0) {
        reconstructedThreeBody[mcInfo.motherLabel] = true;
      }
    }

    for (const auto& mother : mcParticles) {
      if (std::abs(mother.pdgCode()) != constants::physics::Pdg::kHyperTriton) {
        continue;
      }
      const int sign = mother.pdgCode() > 0 ? 1 : -1;
      bool hasHelium3 = false;
      bool hasPionTwoBody = false;
      bool hasProton = false;
      bool hasPionThreeBody = false;
      bool hasDeuteron = false;
      std::array<float, 3> heliumMomentum{};
      std::array<float, 3> twoBodyDecayVertex{};
      ThreeBodyMCInfo threeBodyInfo;
      threeBodyInfo.isReco = false;
      threeBodyInfo.motherLabel = mother.globalIndex();
      threeBodyInfo.motherPdgCode = mother.pdgCode();
      threeBodyInfo.genMomentum = {mother.px(), mother.py(), mother.pz()};
      threeBodyInfo.genPhi = mother.phi();
      threeBodyInfo.genEta = mother.eta();
      threeBodyInfo.genRapidity = mother.y();
      if (mother.mcCollisionId() >= 0 && mother.mcCollisionId() < static_cast<int>(survivedMCEventSelection.size())) {
        threeBodyInfo.survivedEventSelection = survivedMCEventSelection[mother.mcCollisionId()];
      }
      int twoBodyStatusCode = 0;

      for (const auto& daughter : mother.daughters_as<aod::McParticles>()) {
        if (daughter.pdgCode() == sign * constants::physics::Pdg::kHelium3) {
          hasHelium3 = true;
          heliumMomentum = {daughter.px(), daughter.py(), daughter.pz()};
          twoBodyDecayVertex = {daughter.vx(), daughter.vy(), daughter.vz()};
          twoBodyStatusCode = daughter.getProcess();
        } else if (daughter.pdgCode() == -sign * PDG_t::kPiPlus) {
          hasPionTwoBody = true;
          hasPionThreeBody = true;
          threeBodyInfo.pionPdgCode = daughter.pdgCode();
          threeBodyInfo.genMomentumPion = daughter.p();
          threeBodyInfo.genPtPion = daughter.pt();
        } else if (daughter.pdgCode() == sign * PDG_t::kProton) {
          hasProton = true;
          threeBodyInfo.protonPdgCode = daughter.pdgCode();
          threeBodyInfo.genMomentumProton = daughter.p();
          threeBodyInfo.genPtProton = daughter.pt();
          threeBodyInfo.genDecayVertex = {daughter.vx(), daughter.vy(), daughter.vz()};
        } else if (daughter.pdgCode() == sign * constants::physics::Pdg::kDeuteron) {
          hasDeuteron = true;
          threeBodyInfo.deuteronPdgCode = daughter.pdgCode();
          threeBodyInfo.genMomentumDeuteron = daughter.p();
          threeBodyInfo.genPtDeuteron = daughter.pt();
          threeBodyInfo.isDeuteronPrimary = daughter.isPhysicalPrimary();
        }
      }

      if (hasHelium3 && hasPionTwoBody && !reconstructedTwoBody[mother.globalIndex()]) {
        float centralityFT0A = -1.f;
        float centralityFT0C = -1.f;
        float centralityFT0M = -1.f;
        float primaryVertexX = -1.f;
        float primaryVertexY = -1.f;
        float primaryVertexZ = -1.f;
        bool isRecoMCCollision = false;
        bool survivedEventSelection = false;
        if (mother.mcCollisionId() >= 0 && mother.mcCollisionId() < static_cast<int>(recoCollisionForMC.size())) {
          const int recoCollisionId = recoCollisionForMC[mother.mcCollisionId()];
          isRecoMCCollision = recoCollisionId >= 0;
          survivedEventSelection = survivedMCEventSelection[mother.mcCollisionId()];
          if (isRecoMCCollision) {
            const auto collision = collisions.rawIteratorAt(recoCollisionId);
            centralityFT0A = collision.centFT0A();
            centralityFT0C = collision.centFT0C();
            centralityFT0M = collision.centFT0M();
            primaryVertexX = collision.posX();
            primaryVertexY = collision.posY();
            primaryVertexZ = collision.posZ();
          }
        }
        mcHypCands(centralityFT0A, centralityFT0C, centralityFT0M,
                   primaryVertexX, primaryVertexY, primaryVertexZ,
                   runNumber, mother.pdgCode() > 0,
                   -1.f, -1.f, -1.f,
                   -1.f, -1.f, -1.f,
                   -1.f, -1.f, -1.f,
                   -1.f, -1.f, -1.f,
                   -1.f, -1, -1, -1, -1, -1, -1,
                   -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f,
                   0u, 0u, 0u, -1,
                   sign * mother.pt(), mother.phi(), mother.eta(), std::hypot(heliumMomentum[0], heliumMomentum[1]),
                   twoBodyDecayVertex[0] - mother.vx(), twoBodyDecayVertex[1] - mother.vy(), twoBodyDecayVertex[2] - mother.vz(),
                   false, 0u, true, isRecoMCCollision, survivedEventSelection, true, twoBodyStatusCode);
      }

      if (hasProton && hasPionThreeBody && hasDeuteron && !reconstructedThreeBody[mother.globalIndex()]) {
        threeBodyInfo.genCt = RecoDecay::sqrtSumOfSquares(threeBodyInfo.genDecayVertex[0] - mother.vx(),
                                                          threeBodyInfo.genDecayVertex[1] - mother.vy(),
                                                          threeBodyInfo.genDecayVertex[2] - mother.vz()) *
                              constants::physics::MassHyperTriton / (mother.p() + 1.e-10f);
        fillGeneratedThreeBodyMCTable(threeBodyInfo);
      }
    }
  }
  PROCESS_SWITCH(TrackedHypertritonRecoTask, processData, "Process tracked hypertriton candidates in data", true);
  PROCESS_SWITCH(TrackedHypertritonRecoTask, processMC, "Process tracked hypertriton candidates in Monte Carlo", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto metadataInfo = o2::common::core::MetadataHelper{};
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<TrackedHypertritonRecoTask>(cfgc, std::move(metadataInfo))};
}
