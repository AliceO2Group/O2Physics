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

/// \file nucleiThreeBodyBuilder.cxx
/// \brief Generic He3-hadron-hadron secondary-vertex table producer
/// \author ALICE Collaboration

#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/CCDB/EventSelectionParams.h"
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
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{
using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionsPbPb = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using Tracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                         aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                         aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                         aod::pidTOFFullHe>;
using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                           aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                           aod::pidTOFFullHe, aod::McTrackLabels>;

constexpr int He3Pdg = o2::constants::physics::Pdg::kHelium3;
constexpr float InvalidFloat = -999.f;
constexpr int InvalidLabel = -1;
constexpr int NProngs = 3;
constexpr float He3Charge = 2.f;
constexpr float UseCCDBBzThreshold = -990.f;
constexpr int He3PoolId = 1;
constexpr int Daughter1PoolId = 2;
constexpr int Daughter2PoolId = 3;
constexpr std::array<double, 6> DefaultHe3BetheBloch{-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09};
const std::vector<std::string> BetheBlochNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
const std::vector<std::string> He3Name{"He3"};

enum MCMatchStatus : uint8_t {
  MissingDaughterLabel = 0,
  NoCommonImmediateMother,
  CommonImmediateMother,
  ConfiguredDecay
};

bool isSupportedHadronPdg(int pdg)
{
  const int absPdg = std::abs(pdg);
  return absPdg == PDG_t::kPiPlus ||
         absPdg == PDG_t::kKPlus || absPdg == PDG_t::kProton;
}
} // namespace

struct NucleiThreeBodyBuilder {
  Produces<aod::NucleiThreeBodyDatas> outputData;
  Produces<aod::McNucleiThreeBodyDatas> outputMC;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  struct : ConfigurableGroup {
    std::string prefix{"event"};
    Configurable<bool> requireSel8{"requireSel8", true, "Require sel8 (main pp default)"};
    Configurable<float> maxAbsZvtx{"maxAbsZvtx", 10.f, "Maximum absolute PV z; negative disables"};
    Configurable<bool> requireNoITSROFrameBorder{"requireNoITSROFrameBorder", false, "Require kNoITSROFrameBorder"};
    Configurable<bool> requireNoTimeFrameBorder{"requireNoTimeFrameBorder", false, "Require kNoTimeFrameBorder"};
    Configurable<bool> requireNoSameBunchPileup{"requireNoSameBunchPileup", false, "Require kNoSameBunchPileup"};
    Configurable<bool> requireGoodZvtxFT0vsPV{"requireGoodZvtxFT0vsPV", false, "Require kIsGoodZvtxFT0vsPV"};
  } event;

  struct : ConfigurableGroup {
    std::string prefix{"pbpb"};
    Configurable<float> minCentFT0C{"minCentFT0C", 0.f, "Minimum FT0C centrality"};
    Configurable<float> maxCentFT0C{"maxCentFT0C", 100.f, "Maximum FT0C centrality"};
  } pbpb;

  struct : ConfigurableGroup {
    std::string prefix{"he3"};
    Configurable<float> minTPCInnerParam{"minTPCInnerParam", 0.8f, "Minimum He3 TPC rigidity"};
    Configurable<float> maxAbsTPCNSigma{"maxAbsTPCNSigma", 5.f, "Maximum absolute custom He3 TPC n-sigma"};
    Configurable<bool> requirePIDForTracking{"requirePIDForTracking", false, "Require Helium3 or Alpha PID during tracking"};
    Configurable<bool> compensatePIDinTracking{"compensatePIDinTracking", true, "Divide TPC inner parameter by two for charge-two tracking PID"};
    Configurable<LabeledArray<double>> betheBlochParams{"betheBlochParams", {DefaultHe3BetheBloch.data(), 1, 6, He3Name, BetheBlochNames}, "TPC Bethe-Bloch parameters for He3"};
  } he3;

  struct : ConfigurableGroup {
    std::string prefix{"track"};
    Configurable<float> maxAbsEta{"maxAbsEta", 0.9f, "Maximum absolute eta for every daughter"};
    Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 70, "Minimum TPC crossed rows for every daughter"};
    Configurable<float> minAbsDCAxyToPV{"minAbsDCAxyToPV", 0.f, "Minimum absolute DCAxy to the assigned/candidate PV for every daughter; zero disables"};
    Configurable<float> minAbsDCAzToPV{"minAbsDCAzToPV", 0.f, "Minimum absolute DCAz to the assigned/candidate PV for every daughter; zero disables"};
  } trackSelections;

  struct : ConfigurableGroup {
    std::string prefix{"decay"};
    Configurable<int> daughter1Pdg{"daughter1Pdg", PDG_t::kProton, "Signed pi/K/p PDG relative to a positive He3; used for charge pairing, optional PID selection, and MC matching"};
    Configurable<int> daughter2Pdg{"daughter2Pdg", -PDG_t::kPiPlus, "Signed pi/K/p PDG relative to a positive He3; used for charge pairing, optional PID selection, and MC matching"};
    Configurable<int> motherPdg{"motherPdg", 1010020040, "Mother PDG for MC matching; zero accepts either common mother"};
  } decay;

  struct : ConfigurableGroup {
    std::string prefix{"vertex"};
    Configurable<float> minCosPA{"minCosPA", 0.9f, "Minimum cosine of pointing angle; below -1 disables"};
  } vertexSelections;

  struct : ConfigurableGroup {
    std::string prefix{"candidate"};
    Configurable<bool> applySelections{"applySelections", false, "Apply the additional candidate selections"};
    Configurable<float> maxChi2{"maxChi2", 100.f, "Maximum DCA-fitter chi2"};
    Configurable<float> maxDCADaughters{"maxDCADaughters", 1.f, "Maximum square root of fitter chi2"};
    Configurable<float> maxAbsTPCNSigmaDaughter1{"maxAbsTPCNSigmaDaughter1", 5.f, "Maximum daughter 1 TPC n-sigma for its configured hypothesis"};
    Configurable<float> maxAbsTOFNSigmaDaughter1{"maxAbsTOFNSigmaDaughter1", 5.f, "Maximum daughter 1 TOF n-sigma; applied only with TOF"};
    Configurable<float> maxAbsTPCNSigmaDaughter2{"maxAbsTPCNSigmaDaughter2", 5.f, "Maximum daughter 2 TPC n-sigma for its configured hypothesis"};
    Configurable<float> maxAbsTOFNSigmaDaughter2{"maxAbsTOFNSigmaDaughter2", 5.f, "Maximum daughter 2 TOF n-sigma; applied only with TOF"};
    Configurable<float> minDecayLength{"minDecayLength", 0.f, "Minimum decay length"};
    Configurable<float> maxDecayLength{"maxDecayLength", 100.f, "Maximum decay length"};
  } candidateSelections;

  struct : ConfigurableGroup {
    std::string prefix{"pool"};
    Configurable<float> timeMargin{"timeMargin", 800.f, "Track-collision compatibility margin in ns"};
    Configurable<bool> skipAmbiguousTracks{"skipAmbiguousTracks", false, "Reject unassigned ambiguous tracks"};
  } pool;

  struct : ConfigurableGroup {
    std::string prefix{"zorro"};
    Configurable<bool> enabled{"enabled", false, "Enable fHe trigger selection and accounting for skimmed pp data"};
    Configurable<std::string> ccdbPath{"ccdbPath", "EventFiltering/Zorro/", "Base path of Zorro CCDB objects"};
    Configurable<int> bcTolerance{"bcTolerance", 100, "Zorro BC matching tolerance"};
  } zorroOptions;

  struct : ConfigurableGroup {
    std::string prefix{"mc"};
    Configurable<bool> requireConfiguredDecayForGen{"requireConfiguredDecayForGen", true, "Fill generated quantities only for the configured daughters and optional mother PDG"};
  } mc;

  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "CCDB URL"}; // o2-linter: disable=name/configurable (Established CCDB option name.)
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of magnetic-field object"};
  Configurable<float> bzInput{"bz", -999.f, "Magnetic field in kG; below -990 uses CCDB"}; // o2-linter: disable=name/configurable (Established magnetic-field option name.)

  o2::vertexing::DCAFitterN<3> fitter;
  svPoolCreator3Body tripletPool{He3PoolId, Daughter1PoolId, Daughter2PoolId};
  int runNumber = -1;
  float bz = 0.f;
  std::array<float, 6> he3BetheBloch{};
  std::vector<bool> selectedCollisions;
  std::vector<bool> collisionsWithHe3;
  std::vector<bool> identifiedHe3Tracks;
  std::vector<bool> preselectedTracks;

  struct RecoCandidate {
    std::array<int8_t, 3> signs{};
    float pt = 0.f;
    float eta = 0.f;
    float phi = 0.f;
    std::array<float, 3> daughterMomentum{};
    float chi2 = 0.f;
    float dcaDaughters = 0.f;
    float cosPA = 0.f;
    float decayLength = 0.f;
    std::array<float, 3> dcaXYToPV{};
    std::array<float, 3> dcaZToPV{};
    float tpcNSigmaHe3 = InvalidFloat;
    float tofNSigmaHe3 = InvalidFloat;
    std::array<std::array<float, 3>, 2> tpcNSigmaHadron{};
    std::array<std::array<float, 3>, 2> tofNSigmaHadron{};
    std::array<uint16_t, 3> tpcNCls{};
    std::array<uint16_t, 3> tpcCrossedRows{};
    std::array<uint8_t, 3> itsNCls{};
  };

  struct MCInfo {
    std::array<int, 3> pdgs{};
    uint8_t matchStatus = MissingDaughterLabel;
    int motherPdg = 0;
    float motherPt = InvalidFloat;
    float motherEta = InvalidFloat;
    float motherPhi = InvalidFloat;
    float decayLength = InvalidFloat;
  };

  void init(InitContext const&)
  {
    const int enabledProcesses = static_cast<int>(doprocessData) + static_cast<int>(doprocessPbPb) + static_cast<int>(doprocessMC);
    if (enabledProcesses != 1) {
      LOG(fatal) << "Enable exactly one of processData, processPbPb, and processMC";
    }
    if (zorroOptions.enabled && !doprocessData) {
      LOG(fatal) << "Zorro is supported only by the pp data process path";
    }
    if (!isSupportedHadronPdg(decay.daughter1Pdg) || !isSupportedHadronPdg(decay.daughter2Pdg)) {
      LOG(fatal) << "Signed daughter PDG codes must identify a pion, kaon, or proton";
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(zorroOptions.ccdbPath);
    zorro.setBCtolerance(zorroOptions.bcTolerance);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.f);
    fitter.setMinParamChange(1.e-3f);
    fitter.setMinRelChi2Change(0.9f);
    fitter.setMaxDZIni(1.e9f);
    fitter.setMaxChi2(1.e9f);
    fitter.setUseAbsDCA(true);

    tripletPool.setTimeMargin(pool.timeMargin);
    if (pool.skipAmbiguousTracks) {
      tripletPool.setSkipAmbiTracks();
    }
    for (size_t i = 0; i < he3BetheBloch.size(); ++i) {
      he3BetheBloch[i] = he3.betheBlochParams->get("He3", BetheBlochNames[i].c_str());
    }

    registry.add("events", "Event and He3 preselection;step;events", HistType::kTH1D, {{5, -0.5, 4.5}});
    auto events = registry.get<TH1>(HIST("events"));
    events->GetXaxis()->SetBinLabel(1, "all");
    events->GetXaxis()->SetBinLabel(2, "Zorro fHe");
    events->GetXaxis()->SetBinLabel(3, "event selected");
    events->GetXaxis()->SetBinLabel(4, "identified He3");
    events->GetXaxis()->SetBinLabel(5, "candidate filled");
    registry.add("zorroEvents", "Zorro accounting;selection;events", HistType::kTH1D, {{2, -0.5, 1.5}});
    auto zorroEvents = registry.get<TH1>(HIST("zorroEvents"));
    zorroEvents->GetXaxis()->SetBinLabel(1, "fHe before event selection");
    zorroEvents->GetXaxis()->SetBinLabel(2, "fHe after event selection");
    registry.add("triplets", "Triplet building;step;count", HistType::kTH1D, {{3, -0.5, 2.5}});
    auto triplets = registry.get<TH1>(HIST("triplets"));
    triplets->GetXaxis()->SetBinLabel(1, "pool combinations");
    triplets->GetXaxis()->SetBinLabel(2, "fit succeeded");
    triplets->GetXaxis()->SetBinLabel(3, "stored");
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    if (zorroOptions.enabled) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
    auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    if (!grpmag) {
      LOG(fatal) << "Missing magnetic-field object " << grpmagPath << " for timestamp " << bc.timestamp();
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);
    bz = bzInput < UseCCDBBzThreshold ? o2::base::Propagator::Instance()->getNominalBz() : bzInput;
    o2::base::Propagator::Instance()->setNominalBz(bz);
    fitter.setBz(bz);
    runNumber = bc.runNumber();
  }

  template <typename TTrack>
  float correctedHe3Rigidity(TTrack const& track)
  {
    const bool chargeTwoPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return chargeTwoPID && he3.compensatePIDinTracking ? track.tpcInnerParam() / He3Charge : track.tpcInnerParam();
  }

  template <typename TTrack>
  float tpcNSigmaHe3(TTrack const& track)
  {
    const float rigidity = correctedHe3Rigidity(track);
    const float expected = o2::common::BetheBlochAleph(static_cast<float>(rigidity * He3Charge / constants::physics::MassHelium3),
                                                       he3BetheBloch[0], he3BetheBloch[1], he3BetheBloch[2],
                                                       he3BetheBloch[3], he3BetheBloch[4]);
    return (track.tpcSignal() - expected) / (expected * he3BetheBloch[5]);
  }

  template <typename TTrack>
  bool isIdentifiedHe3(TTrack const& track)
  {
    if (!track.hasTPC() || correctedHe3Rigidity(track) < he3.minTPCInnerParam ||
        std::abs(tpcNSigmaHe3(track)) > he3.maxAbsTPCNSigma) {
      return false;
    }
    if (he3.requirePIDForTracking &&
        track.pidForTracking() != o2::track::PID::Helium3 &&
        track.pidForTracking() != o2::track::PID::Alpha) {
      return false;
    }
    return true;
  }

  template <typename TCollision>
  bool passesEventSelection(TCollision const& collision)
  {
    if (event.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (event.maxAbsZvtx >= 0.f && std::abs(collision.posZ()) > event.maxAbsZvtx) {
      return false;
    }
    if (event.requireNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (event.requireNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (event.requireNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (event.requireGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    return true;
  }

  template <bool applyCentrality, typename TCollisions>
  void selectCollisions(TCollisions const& collisions, bool applyZorro)
  {
    selectedCollisions.assign(collisions.size(), false);
    collisionsWithHe3.assign(collisions.size(), false);
    for (const auto& collision : collisions) {
      registry.fill(HIST("events"), 0.);
      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (applyZorro) {
        if (!zorro.isSelected(bc.globalBC(), zorroOptions.bcTolerance)) {
          continue;
        }
        registry.fill(HIST("events"), 1.);
        registry.fill(HIST("zorroEvents"), 0.);
      }
      if (!passesEventSelection(collision)) {
        continue;
      }
      if constexpr (applyCentrality) {
        if (collision.centFT0C() < pbpb.minCentFT0C || collision.centFT0C() > pbpb.maxCentFT0C) {
          continue;
        }
      }
      selectedCollisions[collision.globalIndex()] = true;
      registry.fill(HIST("events"), 2.);
      if (applyZorro) {
        registry.fill(HIST("zorroEvents"), 1.);
      }
    }
  }

  template <typename TCollision, typename TTrack>
  std::array<float, 2> dcaToPV(TCollision const& collision, TTrack const& track)
  {
    auto trackPar = getTrackParCov(track);
    std::array<float, 2> dca{};
    o2::base::Propagator::Instance()->propagateToDCABxByBz(
      {collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, fitter.getMatCorrType(), &dca);
    return dca;
  }

  template <typename TTrack>
  bool passesTrackSelections(TTrack const& track, std::array<float, 2> const& dca)
  {
    return std::abs(track.eta()) <= trackSelections.maxAbsEta &&
           track.tpcNClsCrossedRows() >= trackSelections.minTPCCrossedRows &&
           std::abs(dca[0]) >= trackSelections.minAbsDCAxyToPV &&
           std::abs(dca[1]) >= trackSelections.minAbsDCAzToPV;
  }

  template <typename TTrack>
  static std::array<float, 3> tpcHadronNSigma(TTrack const& track)
  {
    return {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
  }

  template <typename TTrack>
  static std::array<float, 3> tofHadronNSigma(TTrack const& track)
  {
    if (!track.hasTOF()) {
      return {InvalidFloat, InvalidFloat, InvalidFloat};
    }
    return {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
  }

  static float selectedNSigma(std::array<float, 3> const& values, int pdg)
  {
    switch (std::abs(pdg)) {
      case PDG_t::kPiPlus:
        return values[0];
      case PDG_t::kKPlus:
        return values[1];
      case PDG_t::kProton:
        return values[2];
      default:
        return 0.f;
    }
  }

  template <typename TCollision, typename TTrack>
  bool buildCandidate(TCollision const& collision, TTrack const& he3Track,
                      TTrack const& daughter1Track, TTrack const& daughter2Track,
                      RecoCandidate& out)
  {
    const std::array tracks{he3Track, daughter1Track, daughter2Track};
    for (int i = 0; i < NProngs; ++i) {
      const auto dca = dcaToPV(collision, tracks[i]);
      if (!passesTrackSelections(tracks[i], dca)) {
        return false;
      }
      out.dcaXYToPV[i] = dca[0];
      out.dcaZToPV[i] = dca[1];
    }

    auto he3TrackFit = getTrackParCov(he3Track);
    auto daughter1TrackFit = getTrackParCov(daughter1Track);
    auto daughter2TrackFit = getTrackParCov(daughter2Track);
    int nCandidates = 0;
    try {
      nCandidates = fitter.process(he3TrackFit, daughter1TrackFit, daughter2TrackFit);
    } catch (...) {
      LOG(error) << "Exception while fitting a nuclei three-body candidate";
      return false;
    }
    if (nCandidates == 0) {
      return false;
    }
    registry.fill(HIST("triplets"), 1.);

    out.signs = {static_cast<int8_t>(he3Track.sign() > 0 ? 1 : -1),
                 static_cast<int8_t>(daughter1Track.sign() > 0 ? 1 : -1),
                 static_cast<int8_t>(daughter2Track.sign() > 0 ? 1 : -1)};

    std::array<std::array<float, 3>, 3> daughterMomenta{};
    std::array<float, 3> momentum{};
    for (int i = 0; i < NProngs; ++i) {
      fitter.getTrack(i).getPxPyPzGlo(daughterMomenta[i]);
    }
    for (int component = 0; component < NProngs; ++component) {
      daughterMomenta[0][component] *= He3Charge;
      momentum[component] = daughterMomenta[0][component] +
                            daughterMomenta[1][component] +
                            daughterMomenta[2][component];
    }
    for (int i = 0; i < NProngs; ++i) {
      out.daughterMomentum[i] = RecoDecay::sqrtSumOfSquares(
        daughterMomenta[i][0], daughterMomenta[i][1], daughterMomenta[i][2]);
    }
    out.pt = RecoDecay::sqrtSumOfSquares(momentum[0], momentum[1]);
    out.eta = RecoDecay::eta(momentum);
    out.phi = RecoDecay::phi(momentum);
    const float momentumMagnitude = RecoDecay::sqrtSumOfSquares(momentum[0], momentum[1], momentum[2]);

    out.chi2 = fitter.getChi2AtPCACandidate();
    out.dcaDaughters = std::sqrt(std::max(0.f, out.chi2));
    const auto& secondaryVertex = fitter.getPCACandidate();
    const std::array<float, 3> flight{
      static_cast<float>(secondaryVertex[0] - collision.posX()),
      static_cast<float>(secondaryVertex[1] - collision.posY()),
      static_cast<float>(secondaryVertex[2] - collision.posZ())};
    out.decayLength = RecoDecay::sqrtSumOfSquares(flight[0], flight[1], flight[2]);
    out.cosPA = (flight[0] * momentum[0] + flight[1] * momentum[1] + flight[2] * momentum[2]) /
                (out.decayLength * momentumMagnitude + 1.e-10f);
    if (out.cosPA < vertexSelections.minCosPA) {
      return false;
    }

    out.tpcNSigmaHe3 = tpcNSigmaHe3(he3Track);
    out.tofNSigmaHe3 = he3Track.hasTOF() ? he3Track.tofNSigmaHe() : InvalidFloat;
    out.tpcNSigmaHadron[0] = tpcHadronNSigma(daughter1Track);
    out.tpcNSigmaHadron[1] = tpcHadronNSigma(daughter2Track);
    out.tofNSigmaHadron[0] = tofHadronNSigma(daughter1Track);
    out.tofNSigmaHadron[1] = tofHadronNSigma(daughter2Track);
    for (int i = 0; i < NProngs; ++i) {
      out.tpcNCls[i] = tracks[i].tpcNClsFound();
      out.tpcCrossedRows[i] = tracks[i].tpcNClsCrossedRows();
      out.itsNCls[i] = tracks[i].itsNCls();
    }

    if (!candidateSelections.applySelections) {
      return true;
    }
    if (out.chi2 > candidateSelections.maxChi2 ||
        out.dcaDaughters > candidateSelections.maxDCADaughters ||
        out.decayLength < candidateSelections.minDecayLength ||
        out.decayLength > candidateSelections.maxDecayLength ||
        std::abs(selectedNSigma(out.tpcNSigmaHadron[0], decay.daughter1Pdg)) > candidateSelections.maxAbsTPCNSigmaDaughter1 ||
        std::abs(selectedNSigma(out.tpcNSigmaHadron[1], decay.daughter2Pdg)) > candidateSelections.maxAbsTPCNSigmaDaughter2 ||
        (daughter1Track.hasTOF() &&
         std::abs(selectedNSigma(out.tofNSigmaHadron[0], decay.daughter1Pdg)) > candidateSelections.maxAbsTOFNSigmaDaughter1) ||
        (daughter2Track.hasTOF() &&
         std::abs(selectedNSigma(out.tofNSigmaHadron[1], decay.daughter2Pdg)) > candidateSelections.maxAbsTOFNSigmaDaughter2)) {
      return false;
    }
    return true;
  }

  void fillDataTable(RecoCandidate const& c)
  {
    outputData(c.signs[0], c.signs[1], c.signs[2],
               c.pt, c.eta, c.phi,
               c.daughterMomentum[0], c.daughterMomentum[1], c.daughterMomentum[2],
               c.chi2, c.dcaDaughters, c.cosPA, c.decayLength,
               c.dcaXYToPV[0], c.dcaZToPV[0],
               c.dcaXYToPV[1], c.dcaZToPV[1],
               c.dcaXYToPV[2], c.dcaZToPV[2],
               c.tpcNSigmaHe3, c.tofNSigmaHe3,
               c.tpcNSigmaHadron[0][0], c.tpcNSigmaHadron[0][1], c.tpcNSigmaHadron[0][2],
               c.tofNSigmaHadron[0][0], c.tofNSigmaHadron[0][1], c.tofNSigmaHadron[0][2],
               c.tpcNSigmaHadron[1][0], c.tpcNSigmaHadron[1][1], c.tpcNSigmaHadron[1][2],
               c.tofNSigmaHadron[1][0], c.tofNSigmaHadron[1][1], c.tofNSigmaHadron[1][2],
               c.tpcNCls[0], c.tpcNCls[1], c.tpcNCls[2],
               c.tpcCrossedRows[0], c.tpcCrossedRows[1], c.tpcCrossedRows[2],
               c.itsNCls[0], c.itsNCls[1], c.itsNCls[2]);
  }

  template <typename TTrack>
  MCInfo getMCInfo(TTrack const& he3Track, TTrack const& daughter1Track,
                   TTrack const& daughter2Track, aod::McParticles const& mcParticles)
  {
    MCInfo info;
    if (!he3Track.has_mcParticle() || !daughter1Track.has_mcParticle() || !daughter2Track.has_mcParticle()) {
      return info;
    }
    const auto he3Particle = he3Track.template mcParticle_as<aod::McParticles>();
    const auto daughter1Particle = daughter1Track.template mcParticle_as<aod::McParticles>();
    const auto daughter2Particle = daughter2Track.template mcParticle_as<aod::McParticles>();
    info.pdgs = {he3Particle.pdgCode(), daughter1Particle.pdgCode(), daughter2Particle.pdgCode()};
    info.matchStatus = NoCommonImmediateMother;

    int commonMotherLabel = InvalidLabel;
    for (const auto& he3Mother : he3Particle.template mothers_as<aod::McParticles>()) {
      for (const auto& daughter1Mother : daughter1Particle.template mothers_as<aod::McParticles>()) {
        if (he3Mother.globalIndex() != daughter1Mother.globalIndex()) {
          continue;
        }
        for (const auto& daughter2Mother : daughter2Particle.template mothers_as<aod::McParticles>()) {
          if (he3Mother.globalIndex() == daughter2Mother.globalIndex()) {
            commonMotherLabel = he3Mother.globalIndex();
            break;
          }
        }
      }
    }
    if (commonMotherLabel < 0) {
      return info;
    }

    const auto selectedMother = mcParticles.rawIteratorAt(commonMotherLabel);
    info.matchStatus = CommonImmediateMother;
    info.motherPdg = selectedMother.pdgCode();
    const int he3Sign = he3Track.sign() > 0 ? 1 : -1;
    const bool expectedDaughters = info.pdgs[0] == he3Sign * He3Pdg &&
                                   info.pdgs[1] == he3Sign * decay.daughter1Pdg &&
                                   info.pdgs[2] == he3Sign * decay.daughter2Pdg;
    const bool motherAccepted = decay.motherPdg == 0 ||
                                info.motherPdg == he3Sign * std::abs(decay.motherPdg.value);
    const bool configuredDecay = expectedDaughters && motherAccepted;
    if (configuredDecay) {
      info.matchStatus = ConfiguredDecay;
    }
    if ((!mc.requireConfiguredDecayForGen || configuredDecay) && motherAccepted) {
      const std::array<float, 3> motherMomentum{selectedMother.px(), selectedMother.py(), selectedMother.pz()};
      info.motherPt = RecoDecay::sqrtSumOfSquares(motherMomentum[0], motherMomentum[1]);
      info.motherEta = RecoDecay::eta(motherMomentum);
      info.motherPhi = RecoDecay::phi(motherMomentum);
      info.decayLength = RecoDecay::sqrtSumOfSquares(he3Particle.vx() - selectedMother.vx(),
                                                     he3Particle.vy() - selectedMother.vy(),
                                                     he3Particle.vz() - selectedMother.vz());
    }
    return info;
  }

  void fillMCTable(RecoCandidate const& c, MCInfo const& m)
  {
    outputMC(c.signs[0], c.signs[1], c.signs[2],
             c.pt, c.eta, c.phi,
             c.daughterMomentum[0], c.daughterMomentum[1], c.daughterMomentum[2],
             c.chi2, c.dcaDaughters, c.cosPA, c.decayLength,
             c.dcaXYToPV[0], c.dcaZToPV[0],
             c.dcaXYToPV[1], c.dcaZToPV[1],
             c.dcaXYToPV[2], c.dcaZToPV[2],
             c.tpcNSigmaHe3, c.tofNSigmaHe3,
             c.tpcNSigmaHadron[0][0], c.tpcNSigmaHadron[0][1], c.tpcNSigmaHadron[0][2],
             c.tofNSigmaHadron[0][0], c.tofNSigmaHadron[0][1], c.tofNSigmaHadron[0][2],
             c.tpcNSigmaHadron[1][0], c.tpcNSigmaHadron[1][1], c.tpcNSigmaHadron[1][2],
             c.tofNSigmaHadron[1][0], c.tofNSigmaHadron[1][1], c.tofNSigmaHadron[1][2],
             c.tpcNCls[0], c.tpcNCls[1], c.tpcNCls[2],
             c.tpcCrossedRows[0], c.tpcCrossedRows[1], c.tpcCrossedRows[2],
             c.itsNCls[0], c.itsNCls[1], c.itsNCls[2],
             m.pdgs[0], m.pdgs[1], m.pdgs[2],
             m.matchStatus, m.motherPdg,
             m.motherPt, m.motherEta, m.motherPhi, m.decayLength);
  }

  template <bool isMC, bool isPbPb, typename TCollisions, typename TTracks>
  void build(TCollisions const& collisions, TTracks const& tracks,
             aod::AmbiguousTracks const& ambiguousTracks,
             aod::BCsWithTimestamps const& bcs,
             aod::McParticles const* mcParticles = nullptr)
  {
    selectCollisions<isPbPb>(collisions, zorroOptions.enabled && !isMC);
    identifiedHe3Tracks.assign(tracks.size(), false);
    preselectedTracks.assign(tracks.size(), false);

    // Apply the common track preselection before identifying He3 candidates.
    // Consequently, the triplet pools are never built unless a selected event
    // contains at least one He3 that also passes the minimal track selections.
    for (const auto& track : tracks) {
      if (!track.has_collision() || track.collisionId() < 0 ||
          !selectedCollisions[track.collisionId()]) {
        continue;
      }
      const auto collision = collisions.rawIteratorAt(track.collisionId());
      const auto dca = dcaToPV(collision, track);
      if (!passesTrackSelections(track, dca)) {
        continue;
      }
      preselectedTracks[track.globalIndex()] = true;
      if (isIdentifiedHe3(track)) {
        identifiedHe3Tracks[track.globalIndex()] = true;
        collisionsWithHe3[track.collisionId()] = true;
      }
    }
    for (size_t i = 0; i < collisionsWithHe3.size(); ++i) {
      if (collisionsWithHe3[i]) {
        registry.fill(HIST("events"), 3.);
      }
    }

    tripletPool.clearPools();
    tripletPool.fillBC2Coll(collisions, bcs);
    for (const auto& track : tracks) {
      if (!track.has_collision() || track.collisionId() < 0 ||
          !preselectedTracks[track.globalIndex()] ||
          !collisionsWithHe3[track.collisionId()]) {
        continue;
      }
      if (identifiedHe3Tracks[track.globalIndex()]) {
        tripletPool.appendTrackCand(track, collisions, He3PoolId, ambiguousTracks, bcs);
      } else {
        tripletPool.appendTrackCand(track, collisions, Daughter1PoolId, ambiguousTracks, bcs);
        tripletPool.appendTrackCand(track, collisions, Daughter2PoolId, ambiguousTracks, bcs);
      }
    }

    const bool daughter1LikeSign = decay.daughter1Pdg > 0;
    const bool daughter2LikeSign = decay.daughter2Pdg > 0;
    auto& triplets = tripletPool.getSVCandPool(collisions, daughter1LikeSign, daughter2LikeSign);
    for (const auto& triplet : triplets) {
      registry.fill(HIST("triplets"), 0.);
      const auto he3Track = tracks.rawIteratorAt(triplet.tr0Idx);
      const auto daughter1Track = tracks.rawIteratorAt(triplet.tr1Idx);
      const auto daughter2Track = tracks.rawIteratorAt(triplet.tr2Idx);
      for (int collisionIndex = triplet.collBracket.getMin(); collisionIndex <= triplet.collBracket.getMax(); ++collisionIndex) {
        if (!selectedCollisions[collisionIndex] || !collisionsWithHe3[collisionIndex]) {
          continue;
        }
        const auto collision = collisions.rawIteratorAt(collisionIndex);
        RecoCandidate reco;
        if (!buildCandidate(collision, he3Track, daughter1Track, daughter2Track, reco)) {
          continue;
        }
        if constexpr (isMC) {
          fillMCTable(reco, getMCInfo(he3Track, daughter1Track, daughter2Track, *mcParticles));
        } else {
          fillDataTable(reco);
        }
        registry.fill(HIST("events"), 4.);
        registry.fill(HIST("triplets"), 2.);
      }
    }
  }

  void processData(Collisions const& collisions, Tracks const& tracks,
                   aod::AmbiguousTracks const& ambiguousTracks,
                   aod::BCsWithTimestamps const& bcs)
  {
    build<false, false>(collisions, tracks, ambiguousTracks, bcs);
  }
  PROCESS_SWITCH(NucleiThreeBodyBuilder, processData, "Process pp data", true);

  void processPbPb(CollisionsPbPb const& collisions, Tracks const& tracks,
                   aod::AmbiguousTracks const& ambiguousTracks,
                   aod::BCsWithTimestamps const& bcs)
  {
    build<false, true>(collisions, tracks, ambiguousTracks, bcs);
  }
  PROCESS_SWITCH(NucleiThreeBodyBuilder, processPbPb, "Process Pb-Pb data", false);

  void processMC(CollisionsMC const& collisions, TracksMC const& tracks,
                 aod::AmbiguousTracks const& ambiguousTracks,
                 aod::BCsWithTimestamps const& bcs,
                 aod::McParticles const& mcParticles)
  {
    build<true, false>(collisions, tracks, ambiguousTracks, bcs, &mcParticles);
  }
  PROCESS_SWITCH(NucleiThreeBodyBuilder, processMC, "Process reconstructed MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<NucleiThreeBodyBuilder>(cfgc)};
}
