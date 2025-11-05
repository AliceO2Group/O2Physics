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

/// \file resonancesGfwFlow.cxx
/// \brief PID flow for resonances using the generic framework
/// \author Preet Bhanjan Pati <preet.bhanjan.pati@cern.ch>

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/GenericFramework/Core/GFWWeightsList.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "Math/Vector4D.h"
#include <TF1.h>
#include <TProfile.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

namespace
{
std::vector<std::shared_ptr<TProfile>> refV2;
std::vector<std::shared_ptr<TProfile3D>> phiV2;
std::vector<std::shared_ptr<TProfile3D>> k0V2;
std::vector<std::shared_ptr<TProfile3D>> lambdaV2;

std::vector<std::vector<std::shared_ptr<TProfile>>> refBoot;
std::vector<std::vector<std::shared_ptr<TProfile3D>>> phiBoot;
std::vector<std::vector<std::shared_ptr<TProfile3D>>> k0Boot;
std::vector<std::vector<std::shared_ptr<TProfile3D>>> lambdaBoot;
} // namespace

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::analysis::genericframework
{
GFWRegions regions;
GFWCorrConfigs configs;
} // namespace o2::analysis::genericframework

template <typename T>
auto projectMatrix(Array2D<T> const& mat, std::array<float, 6>& array1, std::array<float, 6>& array2, std::array<float, 6>& array3)
{
  for (auto j = 0; j < static_cast<int>(mat.cols); ++j) {
    array1[j] = mat(0, j);
    array2[j] = mat(1, j);
    array3[j] = mat(2, j);
  }
  return;
}
template <typename T, typename P>
auto readMatrix(Array2D<T> const& mat, P& array)
{
  for (auto i = 0; i < static_cast<int>(mat.rows); ++i) {
    for (auto j = 0; j < static_cast<int>(mat.cols); ++j) {
      array[i][j] = mat(i, j);
    }
  }

  return;
}

using namespace o2::analysis::genericframework;

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
static constexpr int LongArrayInt[3][20] = {{1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {2, 2, 2, -2, -2, -2, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {3, 3, 3, -3, -3, -3, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}};

struct ResonancesGfwFlow {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  enum OutputSpecies {
    K0 = 0,
    LAMBDA = 1,
    PHI = 2,
    ANLAMBDA = 3,
    REF = 4,
    kCount_OutputSpecies
  };
  enum Particles {
    PIONS,
    KAONS,
    PROTONS
  };
  enum ParticleCuts {
    kCosPA = 0,
    kMassMin,
    kMassMax,
    kPosTrackPt,
    kNegTrackPt,
    kDCAPosToPVMin,
    kDCANegToPVMin,
    kLifeTime,
    kRadiusMin,
    kRadiusMax,
    kRapidity
  };
  enum ParticleSwitches {
    kUseParticle = 0,
    kUseCosPA,
    kMassBins,
    kDCABetDaug,
    kUseProperLifetime,
    kUseV0Radius
  };
  enum EventCutTypes {
    kFilteredEvents = 0,
    kAfterSel8,
    kUseNoTimeFrameBorder,
    kUseNoITSROFrameBorder,
    kUseNoSameBunchPileup,
    kUseGoodZvtxFT0vsPV,
    kUseNoCollInTimeRangeStandard,
    kUseGoodITSLayersAll,
    kUseNoCollInRofStandard,
    kUseNoHighMultCollInPrevRof,
    kUseOccupancy,
    kUseMultCorrCut,
    kUseT0AV0ACut,
    kUseVertexITSTPC,
    kUseTVXinTRD
  };
  enum TrackCutTypes {
    kFilteredTracks = 0,
    kUseGlobalTracks,
    kUsePvContributor,
    kItsClustersCut,
    kHasTpcSignal,
    kTpcClustersCut,
    kTpcCrossedRowsCut
  };

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 50, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCrossRows, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 5, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Minimum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Maximum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgUseGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgFakeKaonCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaVtxz, bool, true, "Use Phi, Eta, VertexZ dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiPtCent, bool, false, "Use Phi, Pt, Centrality dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaPt, bool, false, "Use Phi, Eta, Pt dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgUseBootStrap, bool, true, "Use bootstrap for error estimation")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, true, "Use track density efficiency correction")
  O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  O2_DEFINE_CONFIGURABLE(cfgUseLsPhi, bool, true, "Use LikeSign for Phi v2")
  O2_DEFINE_CONFIGURABLE(cfgUseOnlyTPC, bool, true, "Use only TPC PID for daughter selection")
  O2_DEFINE_CONFIGURABLE(cfgUseStrictPID, bool, true, "Use strict PID cuts for TPC")
  O2_DEFINE_CONFIGURABLE(cfgUseAsymmetricPID, bool, false, "Use asymmetric PID cuts")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")

  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
  Configurable<std::vector<int>> cfgUseEventCuts{"cfgUseEventCuts", std::vector<int>{1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, "Switch for various event cuts [Filtered Events, Sel8, kNoTimeFrameBorder, kNoITSROFrameBorder, kNoSameBunchPileup, kIsGoodZvtxFT0vsPV, kNoCollInTimeRangeStandard, kIsGoodITSLayersAll, kNoCollInRofStandard, kNoHighMultCollInPrevRof, Occupancy, Multiplicity correlation, T0AV0A 3 sigma cut, kIsVertexITSTPC, kTVXinTRD]"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 3, 6, {"TPC", "TOF", "ITS"}, {"pos_pi", "pos_ka", "pos_pr", "neg_pi", "neg_ka", "neg_pr"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  Configurable<LabeledArray<float>> resonanceCuts{"resonanceCuts", {LongArrayFloat[0], 3, 11, {"K0", "Lambda", "Phi"}, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "Rapidity"}}, "Labeled array (float) for various cuts on resonances"};
  Configurable<LabeledArray<int>> resonanceSwitches{"resonanceSwitches", {LongArrayInt[0], 3, 6, {"K0", "Lambda", "Phi"}, {"UseParticle", "UseCosPA", "NMassBins", "DCABetDaug", "UseProperLifetime", "UseV0Radius"}}, "Labeled array (int) for various cuts on resonances"};

  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN08", "refP08", "refFull", "poiNphi", "poiPphi", "poifullphi", "olNphi", "olPphi", "olfullphi", "poiNk0", "poiPk0", "poifullk0", "olNk0", "olPk0", "olfullk0", "poiNlam", "poiPlam", "poifulllam", "olNlam", "olPlam", "olfulllam", "poiNantilam", "poiPantilam", "poifullantilam", "olNantilam", "olPantilam", "olfullantilam"}, {-0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8, -0.8, 0.4, -0.8}, {-0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8, -0.4, 0.8, 0.8}, {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 2, 2, 2, 32, 32, 32, 4, 4, 4, 64, 64, 64, 8, 8, 8, 128, 128, 128, 16, 16, 16, 256, 256, 256}}, "Configurations for GFW regions"};
  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"refN08 {2} refP08 {-2}", "refN08 {2 2} refP08 {-2 -2}", "poiNphi refN08 | olNphi {2} refP08 {-2}", "poiNphi refN08 | olNphi {2 2} refP08 {-2 -2}", "poiPphi refP08 | olPphi {2} refN08 {-2}", "poiPphi refP08 | olPphi {2 2} refN08 {-2 -2}", "poiNk0 refN08 | olNk0 {2} refP08 {-2}", "poiNk0 refN08 | olNk0 {2 2} refP08 {-2 -2}", "poiPk0 refP08 | olPk0 {2} refN08 {-2}", "poiPk0 refP08 | olPk0 {2 2} refN08 {-2 -2}", "poiNlam refN08 | olNlam {2} refP08 {-2}", "poiNlam refN08 | olNlam {2 2} refP08 {-2 -2}", "poiPlam refP08 | olPlam {2} refN08 {-2}", "poiPlam refP08 | olPlam {2 2} refN08 {-2 -2}", "poiNantilam refN08 | olNantilam {2} refP08 {-2}", "poiNantilam refN08 | olNantilam {2 2} refP08 {-2 -2}", "poiPantilam refP08 | olPantilam {2} refN08 {-2}", "poiPantilam refP08 | olPantilam {2 2} refN08 {-2 -2}"}, {"Ref08Gap22", "Ref08Gap24", "PhiF08Gap22", "PhiF08Gap24", "PhiB08Gap22", "PhiB08Gap24", "K0F08Gap22", "K0F08Gap24", "K0B08Gap22", "K0B08Gap24", "LamF08Gap22", "LamF08Gap24", "LamB08Gap22", "LamB08Gap24", "AnLamF08Gap22", "AnLamF08Gap24", "AnLamB08Gap22", "AnLamB08Gap24"}, {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, "Configurations for each correlation to calculate"};

  // Defining configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  SliceCache cache;
  Partition<AodTracksWithoutBayes> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<AodTracksWithoutBayes> negTracks = aod::track::signed1Pt < 0.0f;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::array<std::array<float, 11>, 3> resoCutVals;
  std::array<std::array<int, 6>, 3> resoSwitchVals;
  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;
  std::vector<int> eventCuts;

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TAxis* fPhiMassAxis;
  TAxis* fK0MassAxis;
  TAxis* fLambdaMassAxis;
  TRandom3* fRndm = new TRandom3(0);

  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  // local track density correction
  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;

  // Additional Event selection cuts - Copy from flowGenericFramework.cxx
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  void init(InitContext const&)
  {
    int64_t noLaterThan = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // Initilizing ccdb
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan);

    LOGF(info, "flowGenericFramework::init()");
    regions.SetNames(cfgRegions->GetNames());
    regions.SetEtaMin(cfgRegions->GetEtaMin());
    regions.SetEtaMax(cfgRegions->GetEtaMax());
    regions.SetpTDifs(cfgRegions->GetpTDifs());
    regions.SetBitmasks(cfgRegions->GetBitmasks());
    configs.SetCorrs(cfgCorrConfig->GetCorrs());
    configs.SetHeads(cfgCorrConfig->GetHeads());
    configs.SetpTDifs(cfgCorrConfig->GetpTDifs());
    configs.SetpTCorrMasks(cfgCorrConfig->GetpTCorrMasks());
    regions.Print();
    configs.Print();

    projectMatrix(nSigmas->getData(), tpcNsigmaCut, tofNsigmaCut, itsNsigmaCut);
    readMatrix(resonanceCuts->getData(), resoCutVals);
    readMatrix(resonanceSwitches->getData(), resoSwitchVals);
    eventCuts = cfgUseEventCuts;

    AxisSpec singleCount = {1, 0, 1};
    AxisSpec axisK0Mass = {resoSwitchVals[K0][kMassBins], resoCutVals[K0][kMassMin], resoCutVals[K0][kMassMax]};
    AxisSpec axisLambdaMass = {resoSwitchVals[LAMBDA][kMassBins], resoCutVals[LAMBDA][kMassMin], resoCutVals[LAMBDA][kMassMax]};
    AxisSpec axisPhiMass = {resoSwitchVals[PHI][kMassBins], resoCutVals[PHI][kMassMin], resoCutVals[PHI][kMassMax]};

    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});

    refBoot.resize(cfgNbootstrap);
    phiBoot.resize(cfgNbootstrap);
    k0Boot.resize(cfgNbootstrap);
    lambdaBoot.resize(cfgNbootstrap);

    // Defining histograms to store correlations
    for (auto i = 0; i < configs.GetSize(); ++i) {
      if (resoSwitchVals[PHI][kUseParticle] && configs.GetHeads()[i].starts_with("Phi")) {
        phiV2.push_back(histos.add<TProfile3D>(Form("h%spt", configs.GetHeads()[i].c_str()), "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}}));
        if (cfgUseBootStrap) {
          for (int j = 0; j < cfgNbootstrap; ++j) {
            phiBoot[j].push_back(histos.add<TProfile3D>(Form("BootStrap/h%spt_boot_%d", configs.GetHeads()[i].c_str(), j), "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}}));
          }
        } // end of bootstrap condition
      } // end of phi loop

      if (resoSwitchVals[K0][kUseParticle] && configs.GetHeads()[i].starts_with("K0")) {
        k0V2.push_back(histos.add<TProfile3D>(Form("h%spt", configs.GetHeads()[i].c_str()), "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}}));
        if (cfgUseBootStrap) {
          for (int j = 0; j < cfgNbootstrap; ++j) {
            k0Boot[j].push_back(histos.add<TProfile3D>(Form("BootStrap/h%spt_boot_%d", configs.GetHeads()[i].c_str(), j), "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}}));
          }
        } // end of bootstrap condition
      } // end of K0 loop

      if (resoSwitchVals[LAMBDA][kUseParticle] && (configs.GetHeads()[i].starts_with("Lam") || configs.GetHeads()[i].starts_with("AnLam"))) {
        lambdaV2.push_back(histos.add<TProfile3D>(Form("h%spt", configs.GetHeads()[i].c_str()), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}}));
        if (cfgUseBootStrap) {
          for (int j = 0; j < cfgNbootstrap; ++j) {
            lambdaBoot[j].push_back(histos.add<TProfile3D>(Form("BootStrap/h%spt_boot_%d", configs.GetHeads()[i].c_str(), j), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}}));
          }
        } // end of bootstrap condition
      } // end of lambda loop

      if (configs.GetHeads()[i].starts_with("Ref")) {
        refV2.push_back(histos.add<TProfile>(Form("h%s", configs.GetHeads()[i].c_str()), "", {HistType::kTProfile, {axisMultiplicity}}));
        if (cfgUseBootStrap) {
          for (int j = 0; j < cfgNbootstrap; ++j) {
            refBoot[j].push_back(histos.add<TProfile>(Form("BootStrap/h%s_boot_%d", configs.GetHeads()[i].c_str(), j), "", {HistType::kTProfile, {axisMultiplicity}}));
          }
        } // end of bootstrap condition
      } // end of ref loop

    } // end of configs loop

    if (cfgUseLsPhi) {
      histos.add("hLsPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});
    }

    if (resoSwitchVals[PHI][kUseParticle]) {
      histos.add("KaPlusTPC", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("KaMinusTPC", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("KaPlusTOF", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("KaMinusTOF", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hPhiPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hPhiEta", "", {HistType::kTH1D, {axisEta}});
      histos.add("hPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});

      histos.add("hPhiCount", "Number of Phi;; Count", {HistType::kTH1D, {{5, 0, 5}}});
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(1, "Phi candidates");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(2, "Daughter track selection");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(3, "Fake Kaon");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(4, "CosPA");
      histos.get<TH1>(HIST("hPhiCount"))->GetXaxis()->SetBinLabel(5, "Rapidity cut");
    }
    if (resoSwitchVals[K0][kUseParticle]) {
      histos.add("PiPlusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiPlusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});
      histos.add("hK0Mass_sparse", "", {HistType::kTHnSparseF, {{axisK0Mass, axisPt, axisMultiplicity}}});
      histos.add("hK0s", "", {HistType::kTH1D, {singleCount}});

      histos.add("hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{10, 0, 10}}});
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(1, "K0 candidates");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(3, "Mass cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(7, "V0radius");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(8, "CosPA");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
      histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(10, "Daughter track selection");
    }
    if (resoSwitchVals[LAMBDA][kUseParticle]) {
      histos.add("PrPlusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PrPlusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});
      histos.add("hLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, axisPt, axisMultiplicity}}});
      histos.add("PiPlusTPC_AL", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PrMinusTPC_AL", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiPlusTOF_AL", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PrMinusTOF_AL", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hAntiLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hAntiLambdaEta", "", {HistType::kTH1D, {axisEta}});
      histos.add("hAntiLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, axisPt, axisMultiplicity}}});
      histos.add("hLambdas", "", {HistType::kTH1D, {singleCount}});

      histos.add("hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{10, 0, 10}}});
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(1, "Lambda candidates");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(2, "Daughter pt");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(3, "Mass cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(4, "Rapidity cut");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(5, "DCA to PV");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(6, "DCA between daughters");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(7, "V0radius");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(8, "CosPA");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(9, "Proper lifetime");
      histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(10, "Daughter track selection");
    }

    histos.add("hEventCount", "Number of Events;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered event");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseMultCorrCut + 1, "Multiplicity correlation Cut");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseT0AV0ACut + 1, "T0AV0A cut");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseVertexITSTPC + 1, "kIsVertexITSTPC");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseTVXinTRD + 1, "kTVXinTRD");

    histos.add("hTrackCount", "Number of Tracks;; Count", {HistType::kTH1D, {{7, -0.5, 6.5}}});
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kFilteredTracks + 1, "Filtered track");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kUseGlobalTracks + 1, "Global tracks");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kUsePvContributor + 1, "PV contributor");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kItsClustersCut + 1, "ITS clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kHasTpcSignal + 1, "TPC signal");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kTpcClustersCut + 1, "TPC clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kTpcCrossedRowsCut + 1, "TPC crossed rows");

    if (cfgOutputNUAWeights) {
      histos.add<TH3>("NUA/hPhiEtaVtxz_ref", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_k0", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_lambda", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_anlambda", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_phi", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});

      histos.add<TH3>("NUA/hPhiPtCent_ref", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, {20, 0, 100}}});
      histos.add<TH3>("NUA/hPhiPtCent_k0", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, {20, 0, 100}}});
      histos.add<TH3>("NUA/hPhiPtCent_lambda", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, {20, 0, 100}}});
      histos.add<TH3>("NUA/hPhiPtCent_anlambda", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, {20, 0, 100}}});
      histos.add<TH3>("NUA/hPhiPtCent_phi", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, {20, 0, 100}}});

      histos.add<TH3>("NUA/hPhiEtaPt_ref", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_k0", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_lambda", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_anlambda", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_phi", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, axisPt}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    fPhiMassAxis = new TAxis(resoSwitchVals[PHI][kMassBins], resoCutVals[PHI][kMassMin], resoCutVals[PHI][kMassMax]);
    fK0MassAxis = new TAxis(resoSwitchVals[K0][kMassBins], resoCutVals[K0][kMassMin], resoCutVals[K0][kMassMax]);
    fLambdaMassAxis = new TAxis(resoSwitchVals[LAMBDA][kMassBins], resoCutVals[LAMBDA][kMassMin], resoCutVals[LAMBDA][kMassMax]);

    int nPhisPtMassBins = nPtBins * resoSwitchVals[PHI][kMassBins];
    int nK0sPtMassBins = nPtBins * resoSwitchVals[K0][kMassBins];
    int nLambdasPtMassBins = nPtBins * resoSwitchVals[LAMBDA][kMassBins];
    int nPtMassBins;

    //********** Defining the regions  **********
    for (auto i(0); i < regions.GetSize(); ++i) {
      if (regions.GetNames()[i].ends_with("phi")) {
        nPtMassBins = nPhisPtMassBins;
      } else if (regions.GetNames()[i].ends_with("k0")) {
        nPtMassBins = nK0sPtMassBins;
      } else if (regions.GetNames()[i].ends_with("lam") || regions.GetNames()[i].ends_with("antilam")) {
        nPtMassBins = nLambdasPtMassBins;
      } else {
        nPtMassBins = nPtBins;
      }
      fGFW->AddRegion(regions.GetNames()[i], regions.GetEtaMin()[i], regions.GetEtaMax()[i], (regions.GetpTDifs()[i]) ? nPtMassBins + 1 : 1, regions.GetBitmasks()[i]);
    }

    //********** Defining the correlations  ************
    for (auto i = 0; i < configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(configs.GetCorrs()[i], configs.GetHeads()[i], configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();

    // Multiplicity correlation cuts
    if (eventCuts[kUseMultCorrCut]) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }
    if (eventCuts[kUseT0AV0ACut]) {
      fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    // Track density correction
    if (cfgTrackDensityCorrUse) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = cfgTrackDensityP0;
      std::vector<double> f1p1 = cfgTrackDensityP1;
      for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
        funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
        funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc]);
      }
      funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
      funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
      funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
    }
  }

  template <typename P>
  int findComponent(std::vector<std::shared_ptr<P>>& ptr, const std::string& name)
  {
    int nIndex = -1;
    for (int i = 0; i < static_cast<int>(ptr.size()); i++) {
      if (ptr[i]->GetName() == name) {
        nIndex = i;
      }
    }

    return nIndex;
  }

  void fillProfileBoot(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile> profile, const double& cent)
  {
    double dnx, val;
    if (!corrconf.pTDif) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0)
        return;
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        profile->Fill(cent, val, dnx);
      return;
    }
    return;
  }

  void fillProfileBoot3D(const GFW::CorrConfig& corrconf, std::shared_ptr<TProfile3D> profile, const double& cent, TAxis* partaxis)
  {
    double dnx, val;
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      for (int j = 1; j <= partaxis->GetNbins(); j++) {
        dnx = fGFW->Calculate(corrconf, ((i - 1) * partaxis->GetNbins()) + (j - 1), kTRUE).real();
        if (dnx == 0)
          continue;
        val = fGFW->Calculate(corrconf, ((i - 1) * partaxis->GetNbins()) + (j - 1), kFALSE).real() / dnx;
        if (std::fabs(val) < 1)
          profile->Fill(fPtAxis->GetBinCenter(i), partaxis->GetBinCenter(j), cent, val, dnx);
      }
    }
    return;
  }

  // Cosine pointing angle cut
  template <typename TTrack1, typename TTrack2>
  double cosinePointingAngle(const TTrack1& track1, const TTrack2& track2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = track1.pt();
    pt2 = track2.pt();
    pz1 = track1.pz();
    pz2 = track2.pz();
    p1 = track1.p();
    p2 = track2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));

    return angle;
  }

  template <typename TTrack>
  bool isFakeKaon(TTrack const& track)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (std::abs(pglobal - ptpc) > cfgFakeKaonCut) {
      return true;
    }
    return false;
  }

  template <typename TTrack>
  bool isGoodTrack(const TTrack& track)
  {
    histos.fill(HIST("hTrackCount"), kFilteredTracks); // Filtered tracks

    if (cfgUseGlobalTrack && !(track.isGlobalTrack())) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kUseGlobalTracks); // After global track selection

    if (!(track.isPVContributor())) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kUsePvContributor); // After PV contributor selection

    if (!(track.itsNCls() > cfgITScluster)) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kItsClustersCut); // After ITS cluster selection

    if (!(track.hasTPC())) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kHasTpcSignal); // If track has TPC signal

    if (!(track.tpcNClsFound() > cfgTpcCluster)) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kTpcClustersCut); // After TPC cluster selection

    if (!(track.tpcNClsCrossedRows() > cfgTpcCrossRows)) {
      return 0;
    }
    histos.fill(HIST("hTrackCount"), kTpcCrossedRowsCut); // After TPC crossed rows selection

    return 1;
  }

  template <typename TTrack>
  int getNsigmaPIDTpcTof(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr())};
    int pid = -1;
    float nsigma = cfgTpcCut;

    // Choose which nSigma to use
    std::array<float, 3> nSigmaToUse = (track.pt() > cfgTofPtCut && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;
    if (track.pt() > cfgTofPtCut && !track.hasTOF())
      return 0;

    const int numSpecies = 3;
    int pidCount = 0;
    // Select particle with the lowest nsigma
    for (int i = 0; i < numSpecies; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        if (pidCount > 0 && cfgUseStrictPID)
          return 0; // more than one particle with low nsigma

        pidCount++;
        pid = i;
        if (!cfgUseStrictPID)
          nsigma = std::abs(nSigmaToUse[i]);
      }
    }
    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  template <typename TTrack>
  int getNsigmaPIDAssymmetric(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = -1;

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[0] < detectorNsigmaCut[0] && nSigmaToUse[0] > detectorNsigmaCut[0 + 3];
    bool isDetectedKaon = nSigmaToUse[1] < detectorNsigmaCut[1] && nSigmaToUse[1] > detectorNsigmaCut[1 + 3];
    bool isDetectedProton = nSigmaToUse[2] < detectorNsigmaCut[2] && nSigmaToUse[2] > detectorNsigmaCut[2 + 3];

    bool isTofPion = nSigmaTOF[0] < tofNsigmaCut[0] && nSigmaTOF[0] > tofNsigmaCut[0 + 3];
    bool isTofKaon = nSigmaTOF[1] < tofNsigmaCut[1] && nSigmaTOF[1] > tofNsigmaCut[1 + 3];
    bool isTofProton = nSigmaTOF[2] < tofNsigmaCut[2] && nSigmaTOF[2] > tofNsigmaCut[2 + 3];

    if (track.pt() > cfgTofPtCut && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgTofPtCut && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS;
    } else if (isKaon) {
      pid = KAONS;
    } else if (isProton) {
      pid = PROTONS;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      uint64_t timestamp = bc.timestamp();
      mAcceptance.clear();
      mAcceptance.resize(kCount_OutputSpecies);

      mAcceptance[K0] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_k0", timestamp);
      if (mAcceptance[K0])
        LOGF(info, "Loaded acceptance weights from %s_k0 (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[K0]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_k0 (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[K0]);

      mAcceptance[LAMBDA] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_lambda", timestamp);
      if (mAcceptance[LAMBDA])
        LOGF(info, "Loaded acceptance weights from %s_lambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[LAMBDA]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_lambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[LAMBDA]);

      mAcceptance[PHI] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_phi", timestamp);
      if (mAcceptance[PHI])
        LOGF(info, "Loaded acceptance weights from %s_phi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[PHI]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_phi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[PHI]);

      mAcceptance[ANLAMBDA] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_anlambda", timestamp);
      if (mAcceptance[ANLAMBDA])
        LOGF(info, "Loaded acceptance weights from %s_anlambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[ANLAMBDA]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_anlambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[ANLAMBDA]);

      mAcceptance[REF] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ref", timestamp);
      if (mAcceptance[REF])
        LOGF(info, "Loaded acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[REF]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[REF]);
    }

    correctionsLoaded = true;
  }

  template <typename TTrack, typename TCollision>
  double getAcceptance(TTrack track, const TCollision collision, int pid_index_reso)
  { // 0 = k0, 1 = lambda, 2 = phi, 3 = anti-lambda, 4 = ref
    if (pid_index_reso < 0 || pid_index_reso >= kCount_OutputSpecies) {
      return 1;
    }

    double wacc = 1;
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();

    if ((cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaPt && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiEtaPt)) {
      LOGF(fatal, "Only one of the three weight options can be used at a time");
    }
    if (!mAcceptance.empty() && correctionsLoaded) {
      if (!mAcceptance[pid_index_reso]) {
        LOGF(fatal, "Acceptance weights not loaded for pidIndex %d", pid_index_reso);
        return 1;
      }
      if (cfgUseWeightPhiEtaVtxz)
        wacc = mAcceptance[pid_index_reso]->getNUA(track.phi(), track.eta(), vtxz);
      if (cfgUseWeightPhiPtCent)
        wacc = mAcceptance[pid_index_reso]->getNUA(track.phi(), track.pt(), cent);
      if (cfgUseWeightPhiEtaPt)
        wacc = mAcceptance[pid_index_reso]->getNUA(track.phi(), track.eta(), track.pt());
    }
    return wacc;
  }

  template <typename vector, typename TCollision>
  double getAcceptancePhi(vector mom, const TCollision collision, int pid_index_reso)
  { // 0 = k0, 1 = lambda, 2 = phi, 3 = anti-lambda, 4 = ref
    if (pid_index_reso < 0 || pid_index_reso >= kCount_OutputSpecies) {
      return 1;
    }

    double wacc = 1;
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();
    double phi = mom.Phi();
    phi = RecoDecay::constrainAngle(phi, 0.0, 1); // constrain azimuthal angle to [0,2pi]

    if ((cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaPt && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiEtaPt)) {
      LOGF(fatal, "Only one of the three weight options can be used at a time");
    }
    if (!mAcceptance.empty() && correctionsLoaded) {
      if (!mAcceptance[pid_index_reso]) {
        LOGF(fatal, "Acceptance weights not loaded for pidIndex %d", pid_index_reso);
        return 1;
      }
      if (cfgUseWeightPhiEtaVtxz)
        wacc = mAcceptance[pid_index_reso]->getNUA(phi, mom.Eta(), vtxz);
      if (cfgUseWeightPhiPtCent)
        wacc = mAcceptance[pid_index_reso]->getNUA(phi, mom.Pt(), cent);
      if (cfgUseWeightPhiEtaPt)
        wacc = mAcceptance[pid_index_reso]->getNUA(phi, mom.Eta(), mom.Pt());
    }
    return wacc;
  }

  template <typename TTrack, typename TCollision>
  void fillWeights(const TTrack track, const TCollision collision, const int& pid_index_reso)
  {
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();
    double pt = track.pt();
    bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

    if (withinPtRef && pid_index_reso == REF) {
      histos.fill(HIST("NUA/hPhiEtaVtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
      histos.fill(HIST("NUA/hPhiPtCent_ref"), track.phi(), track.pt(), cent);
      histos.fill(HIST("NUA/hPhiEtaPt_ref"), track.phi(), track.eta(), track.pt());
    }

    if (withinPtPOI) {
      switch (pid_index_reso) {
        case K0:
          histos.fill(HIST("NUA/hPhiEtaVtxz_k0"), track.phi(), track.eta(), vtxz); // K0 weights
          histos.fill(HIST("NUA/hPhiPtCent_k0"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_k0"), track.phi(), track.eta(), track.pt());
          break;
        case LAMBDA:
          histos.fill(HIST("NUA/hPhiEtaVtxz_lambda"), track.phi(), track.eta(), vtxz); // Lambda weights
          histos.fill(HIST("NUA/hPhiPtCent_lambda"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_lambda"), track.phi(), track.eta(), track.pt());
          break;
        case ANLAMBDA:
          histos.fill(HIST("NUA/hPhiEtaVtxz_anlambda"), track.phi(), track.eta(), vtxz); // Anti-Lambda weights
          histos.fill(HIST("NUA/hPhiPtCent_anlambda"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_anlambda"), track.phi(), track.eta(), track.pt());
          break;
          // Phi weights are filled in the resurrectPhi function
      }
    }
  }

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid)
  {
    if (!(track.itsNCls() > cfgITScluster))
      return 0;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < cfgTpcCluster)
      return false;
    if (!(track.tpcNClsCrossedRows() > cfgTpcCrossRows))
      return 0;

    if (cfgUseOnlyTPC) {
      if (pid == PIONS && std::abs(track.tpcNSigmaPi()) > cfgTpcCut)
        return false;
      if (pid == KAONS && std::abs(track.tpcNSigmaKa()) > cfgTpcCut)
        return false;
      if (pid == PROTONS && std::abs(track.tpcNSigmaPr()) > cfgTpcCut)
        return false;
    } else {
      int partIndex = cfgUseAsymmetricPID ? getNsigmaPIDAssymmetric(track) : getNsigmaPIDTpcTof(track);
      int pidIndex = partIndex - 1; // 0 = pion, 1 = kaon, 2 = proton
      if (pidIndex != pid)
        return false;
    }

    return true;
  }

  template <typename TTrack, typename vector, char... chars, typename TCollision>
  void resurrectPhi(TTrack trackplus, TTrack trackminus, const TCollision collision, vector plusdaug, vector minusdaug, vector mom, double plusmass, const ConstStr<chars...>& hist)
  {
    for (auto const& [partplus, partminus] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trackplus, trackminus))) {
      histos.fill(HIST("hPhiCount"), 0.5);
      if (!selectionV0Daughter(partplus, KAONS) || !selectionV0Daughter(partminus, KAONS)) // 0 = pion, 1 = kaon, 2 = proton
        continue;
      histos.fill(HIST("hPhiCount"), 1.5);

      if (isFakeKaon(partplus) || isFakeKaon(partminus))
        continue;
      histos.fill(HIST("hPhiCount"), 2.5);

      if (resoSwitchVals[PHI][kUseCosPA] && cosinePointingAngle(partplus, partminus) < resoCutVals[PHI][kCosPA])
        continue;
      histos.fill(HIST("hPhiCount"), 3.5);

      histos.fill(HIST("KaPlusTPC"), partplus.pt(), partplus.tpcNSigmaKa());
      histos.fill(HIST("KaPlusTOF"), partplus.pt(), partplus.tofNSigmaKa());
      histos.fill(HIST("KaMinusTPC"), partminus.pt(), partminus.tpcNSigmaKa());
      histos.fill(HIST("KaMinusTOF"), partminus.pt(), partminus.tofNSigmaKa());

      // Calculation using ROOT vectors
      plusdaug = ROOT::Math::PxPyPzMVector(partplus.px(), partplus.py(), partplus.pz(), plusmass);
      minusdaug = ROOT::Math::PxPyPzMVector(partminus.px(), partminus.py(), partminus.pz(), plusmass);
      mom = plusdaug + minusdaug;

      double pt = mom.Pt();
      double invMass = mom.M();
      double phi = mom.Phi();
      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);

      phi = RecoDecay::constrainAngle(phi, 0.0, 1); // constrain azimuthal angle to [0,2pi]

      if (std::abs(mom.Rapidity()) < resoCutVals[PHI][kRapidity]) {
        histos.fill(HIST("hPhiCount"), 4.5);

        histos.fill(hist, invMass, pt, collision.centFT0C());
        histos.fill(HIST("hPhiPhi"), phi);
        histos.fill(HIST("hPhiEta"), mom.Eta());

        // Fill Phi weights
        if (cfgOutputNUAWeights && withinPtPOI) {
          histos.fill(HIST("NUA/hPhiEtaVtxz_phi"), phi, mom.Eta(), collision.posZ());
          histos.fill(HIST("NUA/hPhiPtCent_phi"), phi, pt, collision.centFT0C());
          histos.fill(HIST("NUA/hPhiEtaPt_phi"), phi, mom.Eta(), pt);
        }
        double weff = 1;
        double waccPOI = getAcceptancePhi(mom, collision, PHI);

        if (withinPtPOI)
          fGFW->Fill(mom.Eta(), ((fPtAxis->FindBin(pt) - 1) * fPhiMassAxis->GetNbins()) + (fPhiMassAxis->FindBin(invMass) - 1), phi, weff * waccPOI, 2);
        if (withinPtPOI && withinPtRef)
          fGFW->Fill(mom.Eta(), ((fPtAxis->FindBin(pt) - 1) * fPhiMassAxis->GetNbins()) + (fPhiMassAxis->FindBin(invMass) - 1), phi, weff * waccPOI, 32);
      }
    } // end of combinations loop
    return;
  }

  template <typename TTrack, char... chars, typename TCollision>
  void likeSignPhi(TTrack track, const TCollision collision, double plusmass, const ConstStr<chars...>& hist)
  {
    ROOT::Math::PxPyPzMVector daug1, daug2, mom;
    for (auto const& [part1, part2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(track, track))) {

      if (!selectionV0Daughter(part1, KAONS) || !selectionV0Daughter(part2, KAONS)) // 0 = pion, 1 = kaon, 2 = proton
        continue;
      if (isFakeKaon(part1) || isFakeKaon(part2))
        continue;

      // Calculation using ROOT vectors
      daug1 = ROOT::Math::PxPyPzMVector(part1.px(), part1.py(), part1.pz(), plusmass);
      daug2 = ROOT::Math::PxPyPzMVector(part2.px(), part2.py(), part2.pz(), plusmass);
      mom = daug1 + daug2;

      double pt = mom.Pt();
      double invMass = mom.M();
      double phi = mom.Phi();

      phi = RecoDecay::constrainAngle(phi, 0.0, 1); // constrain azimuthal angle to [0,2pi]

      if (std::abs(mom.Rapidity()) < resoCutVals[PHI][kRapidity]) {
        histos.fill(hist, invMass, pt, collision.centFT0C());
      }
    } // end of positive combinations loop
    return;
  }

  template <typename TCollision, typename V0>
  bool selectionLambda(TCollision const& collision, V0 const& candidate)
  {
    bool isL = false;  // Is lambda candidate
    bool isAL = false; // Is anti-lambda candidate

    double mlambda = candidate.mLambda();
    double mantilambda = candidate.mAntiLambda();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<AodTracksWithoutBayes>();
    auto negtrack = candidate.template negTrack_as<AodTracksWithoutBayes>();

    histos.fill(HIST("hLambdaCount"), 0.5);
    if (postrack.pt() < resoCutVals[LAMBDA][kPosTrackPt] || negtrack.pt() < resoCutVals[LAMBDA][kNegTrackPt])
      return false;

    histos.fill(HIST("hLambdaCount"), 1.5);
    if (mlambda > resoCutVals[LAMBDA][kMassMin] && mlambda < resoCutVals[LAMBDA][kMassMax])
      isL = true;
    if (mantilambda > resoCutVals[LAMBDA][kMassMin] && mantilambda < resoCutVals[LAMBDA][kMassMax])
      isAL = true;

    if (!isL && !isAL) {
      return false;
    }
    histos.fill(HIST("hLambdaCount"), 2.5);

    // Rapidity correction
    if (candidate.yLambda() > resoCutVals[LAMBDA][kRapidity])
      return false;
    histos.fill(HIST("hLambdaCount"), 3.5);
    // DCA cuts for lambda and antilambda
    if (isL) {
      if (std::abs(candidate.dcapostopv()) < resoCutVals[LAMBDA][kDCAPosToPVMin] || std::abs(candidate.dcanegtopv()) < resoCutVals[LAMBDA][kDCANegToPVMin])
        return false;
    }
    if (isAL) {
      if (std::abs(candidate.dcapostopv()) < resoCutVals[LAMBDA][kDCANegToPVMin] || std::abs(candidate.dcanegtopv()) < resoCutVals[LAMBDA][kDCAPosToPVMin])
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > resoSwitchVals[LAMBDA][kDCABetDaug])
      return false;
    histos.fill(HIST("hLambdaCount"), 5.5);
    // v0 radius cuts
    if (resoSwitchVals[LAMBDA][kUseV0Radius] && (candidate.v0radius() < resoCutVals[LAMBDA][kRadiusMin] || candidate.v0radius() > resoCutVals[LAMBDA][kRadiusMax]))
      return false;
    histos.fill(HIST("hLambdaCount"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < resoCutVals[LAMBDA][kCosPA])
      return false;
    histos.fill(HIST("hLambdaCount"), 7.5);
    // Proper lifetime
    if (resoSwitchVals[LAMBDA][kUseProperLifetime] && candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > resoCutVals[LAMBDA][kLifeTime])
      return false;
    histos.fill(HIST("hLambdaCount"), 8.5);
    if (isL) {
      if (!selectionV0Daughter(postrack, PROTONS) || !selectionV0Daughter(negtrack, PIONS))
        return false;
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, PIONS) || !selectionV0Daughter(negtrack, PROTONS))
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 9.5);
    bool withinPtPOI = (cfgCutPtPOIMin < candidate.pt()) && (candidate.pt() < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < candidate.pt()) && (candidate.pt() < cfgCutPtMax);

    float weff = 1;

    if (isL) {
      if (cfgOutputNUAWeights)
        fillWeights(candidate, collision, LAMBDA);

      double waccPOI = getAcceptance(candidate, collision, LAMBDA);
      if (withinPtPOI)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mlambda) - 1), candidate.phi(), waccPOI * weff, 8);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mlambda) - 1), candidate.phi(), waccPOI * weff, 128);

      histos.fill(HIST("hLambdaMass_sparse"), mlambda, candidate.pt(), collision.centFT0C());
      histos.fill(HIST("hLambdaPhi"), candidate.phi());
      histos.fill(HIST("hLambdaEta"), candidate.eta());
      histos.fill(HIST("PrPlusTPC_L"), postrack.pt(), postrack.tpcNSigmaKa());
      histos.fill(HIST("PrPlusTOF_L"), postrack.pt(), postrack.tofNSigmaKa());
      histos.fill(HIST("PiMinusTPC_L"), negtrack.pt(), negtrack.tpcNSigmaKa());
      histos.fill(HIST("PiMinusTOF_L"), negtrack.pt(), negtrack.tofNSigmaKa());
    }
    if (isAL) {
      if (cfgOutputNUAWeights)
        fillWeights(candidate, collision, ANLAMBDA);

      double waccPOI = getAcceptance(candidate, collision, ANLAMBDA);
      if (withinPtPOI)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mantilambda) - 1), candidate.phi(), waccPOI * weff, 16);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mantilambda) - 1), candidate.phi(), waccPOI * weff, 256);

      histos.fill(HIST("hAntiLambdaMass_sparse"), mantilambda, candidate.pt(), collision.centFT0C());
      histos.fill(HIST("hAntiLambdaPhi"), candidate.phi());
      histos.fill(HIST("hAntiLambdaEta"), candidate.eta());
      histos.fill(HIST("PiPlusTPC_AL"), postrack.pt(), postrack.tpcNSigmaKa());
      histos.fill(HIST("PiPlusTOF_AL"), postrack.pt(), postrack.tofNSigmaKa());
      histos.fill(HIST("PrMinusTPC_AL"), negtrack.pt(), negtrack.tpcNSigmaKa());
      histos.fill(HIST("PrMinusTOF_AL"), negtrack.pt(), negtrack.tofNSigmaKa());
    }
    return true;
  }

  template <typename TCollision, typename V0>
  bool selectionK0(TCollision const& collision, V0 const& candidate)
  {
    double mk0 = candidate.mK0Short();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<AodTracksWithoutBayes>();
    auto negtrack = candidate.template negTrack_as<AodTracksWithoutBayes>();

    histos.fill(HIST("hK0Count"), 0.5);
    if (postrack.pt() < resoCutVals[K0][kPosTrackPt] || negtrack.pt() < resoCutVals[K0][kNegTrackPt])
      return false;
    histos.fill(HIST("hK0Count"), 1.5);
    if (mk0 < resoCutVals[K0][kMassMin] && mk0 > resoCutVals[K0][kMassMax])
      return false;
    histos.fill(HIST("hK0Count"), 2.5);
    // Rapidity correction
    if (candidate.yK0Short() > resoCutVals[K0][kRapidity])
      return false;
    histos.fill(HIST("hK0Count"), 3.5);
    // DCA cuts for K0short
    if (std::abs(candidate.dcapostopv()) < resoCutVals[K0][kDCAPosToPVMin] || std::abs(candidate.dcanegtopv()) < resoCutVals[K0][kDCANegToPVMin])
      return false;
    histos.fill(HIST("hK0Count"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > resoSwitchVals[K0][kDCABetDaug])
      return false;
    histos.fill(HIST("hK0Count"), 5.5);
    // v0 radius cuts
    if (resoSwitchVals[K0][kUseV0Radius] && (candidate.v0radius() < resoCutVals[K0][kRadiusMin] || candidate.v0radius() > resoCutVals[K0][kRadiusMax]))
      return false;
    histos.fill(HIST("hK0Count"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < resoCutVals[K0][kCosPA])
      return false;
    histos.fill(HIST("hK0Count"), 7.5);
    // Proper lifetime
    if (resoSwitchVals[K0][kUseProperLifetime] && candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0Short > resoCutVals[K0][kLifeTime])
      return false;
    histos.fill(HIST("hK0Count"), 8.5);
    if (!selectionV0Daughter(postrack, PIONS) || !selectionV0Daughter(negtrack, PIONS))
      return false;
    histos.fill(HIST("hK0Count"), 9.5);
    bool withinPtPOI = (cfgCutPtPOIMin < candidate.pt()) && (candidate.pt() < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < candidate.pt()) && (candidate.pt() < cfgCutPtMax);

    if (cfgOutputNUAWeights)
      fillWeights(candidate, collision, K0);

    float weff = 1;
    double waccPOI = getAcceptance(candidate, collision, K0);

    if (withinPtPOI)
      fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fK0MassAxis->GetNbins()) + (fK0MassAxis->FindBin(mk0) - 1), candidate.phi(), waccPOI * weff, 4);
    if (withinPtPOI && withinPtRef)
      fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fK0MassAxis->GetNbins()) + (fK0MassAxis->FindBin(mk0) - 1), candidate.phi(), waccPOI * weff, 64);

    histos.fill(HIST("hK0Mass_sparse"), mk0, candidate.pt(), collision.centFT0C());
    histos.fill(HIST("hK0Phi"), candidate.phi());
    histos.fill(HIST("hK0Eta"), candidate.eta());
    histos.fill(HIST("PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaKa());
    histos.fill(HIST("PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaKa());
    histos.fill(HIST("PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaKa());
    histos.fill(HIST("PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaKa());

    return true;
  }

  template <typename TCollision>
  bool selectionEvent(TCollision collision, const int mult, const float cent)
  {
    histos.fill(HIST("hEventCount"), kFilteredEvents);
    if (!collision.sel8()) {
      return 0;
    }
    histos.fill(HIST("hEventCount"), kAfterSel8);

    if (eventCuts[kUseNoTimeFrameBorder] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (eventCuts[kUseNoTimeFrameBorder])
      histos.fill(HIST("hEventCount"), kUseNoTimeFrameBorder);

    if (eventCuts[kUseNoITSROFrameBorder] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (eventCuts[kUseNoITSROFrameBorder])
      histos.fill(HIST("hEventCount"), kUseNoITSROFrameBorder);

    if (eventCuts[kUseNoSameBunchPileup] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (eventCuts[kUseNoSameBunchPileup])
      histos.fill(HIST("hEventCount"), kUseNoSameBunchPileup);

    if (eventCuts[kUseGoodZvtxFT0vsPV] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (eventCuts[kUseGoodZvtxFT0vsPV])
      histos.fill(HIST("hEventCount"), kUseGoodZvtxFT0vsPV);

    if (eventCuts[kUseNoCollInTimeRangeStandard] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (eventCuts[kUseNoCollInTimeRangeStandard])
      histos.fill(HIST("hEventCount"), kUseNoCollInTimeRangeStandard);

    if (eventCuts[kUseGoodITSLayersAll] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return 0;
    }
    if (eventCuts[kUseGoodITSLayersAll])
      histos.fill(HIST("hEventCount"), kUseGoodITSLayersAll);

    if (eventCuts[kUseNoCollInRofStandard] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return 0;
    }
    if (eventCuts[kUseNoCollInRofStandard])
      histos.fill(HIST("hEventCount"), kUseNoCollInRofStandard);

    if (eventCuts[kUseNoHighMultCollInPrevRof] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return 0;
    }
    if (eventCuts[kUseNoHighMultCollInPrevRof])
      histos.fill(HIST("hEventCount"), kUseNoHighMultCollInPrevRof);

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[kUseOccupancy] && (occupancy < cfgCutOccupancyMin || occupancy > cfgCutOccupancyMax)) {
      return 0;
    }
    if (eventCuts[kUseOccupancy])
      histos.fill(HIST("hEventCount"), kUseOccupancy);

    if (eventCuts[kUseMultCorrCut]) {
      if (multNTracksPV < fMultPVCutLow->Eval(cent))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(cent))
        return 0;
      if (mult < fMultCutLow->Eval(cent))
        return 0;
      if (mult > fMultCutHigh->Eval(cent))
        return 0;
    }
    if (eventCuts[kUseMultCorrCut])
      histos.fill(HIST("hEventCount"), kUseMultCorrCut);

    // V0A T0A 5 sigma cut
    if (eventCuts[kUseT0AV0ACut] && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > cfgV0AT0Acut * fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (eventCuts[kUseT0AV0ACut])
      histos.fill(HIST("hEventCount"), kUseT0AV0ACut);

    if (eventCuts[kUseVertexITSTPC] && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return 0;
    if (eventCuts[kUseVertexITSTPC])
      histos.fill(HIST("hEventCount"), kUseVertexITSTPC);

    if (eventCuts[kUseTVXinTRD] && collision.alias_bit(kTVXinTRD)) {
      return 0;
    }
    if (eventCuts[kUseTVXinTRD])
      histos.fill(HIST("hEventCount"), kUseTVXinTRD);

    return 1;
  }

  // using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, o2::aod::evsel::NumTracksInTimeRange>;
  ROOT::Math::PxPyPzMVector phiMom, kaonPlus, kaonMinus;
  double massKaPlus = o2::constants::physics::MassKPlus;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0Short = o2::constants::physics::MassK0Short;

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks, aod::V0Datas const& V0s)
  {
    int nTot = tracks.size();
    if (nTot < 1)
      return;

    float vtxz = collision.posZ();
    const auto cent = collision.centFT0C();

    if (!selectionEvent(collision, nTot, cent))
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), cent);
    fGFW->Clear();

    float weff = 1;

    loadCorrections(bc); // load corrections for the each event

    // Track loop for calculating the Qn angles
    double psi2Est = 0, psi3Est = 0, psi4Est = 0;
    float wEPeff = 1;
    double v2 = 0, v3 = 0, v4 = 0;
    // be cautious, this only works for Pb-Pb
    // esimate the Qn angles and vn for this event
    if (cfgTrackDensityCorrUse) {
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        bool withinPtRef = (cfgCutPtMin < track.pt()) && (track.pt() < cfgCutPtMax); // within RF pT rang
        if (withinPtRef) {
          q2x += std::cos(2 * track.phi());
          q2y += std::sin(2 * track.phi());
          q3x += std::cos(3 * track.phi());
          q3y += std::sin(3 * track.phi());
          q4x += std::cos(4 * track.phi());
          q4y += std::sin(4 * track.phi());
        }
      }
      psi2Est = std::atan2(q2y, q2x) / 2.;
      psi3Est = std::atan2(q3y, q3x) / 3.;
      psi4Est = std::atan2(q4y, q4x) / 4.;
      v2 = funcV2->Eval(cent);
      v3 = funcV3->Eval(cent);
      v4 = funcV4->Eval(cent);
    }

    // Actual track loop
    for (auto const& track : tracks) {
      if (!isGoodTrack(track))
        continue;

      double pt = track.pt();
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);

      weff = 1; // Initializing weff for each track

      if (withinPtRef)
        if (cfgOutputNUAWeights)
          fillWeights(track, collision, REF);

      double waccRef = getAcceptance(track, collision, REF);

      if (cfgTrackDensityCorrUse && withinPtRef) {
        double fphi = v2 * std::cos(2 * (track.phi() - psi2Est)) + v3 * std::cos(3 * (track.phi() - psi3Est)) + v4 * std::cos(4 * (track.phi() - psi4Est));
        fphi = (1 + 2 * fphi);
        int pTBinForEff = hFindPtBin->FindBin(track.pt());
        if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
          wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * tracks.size());
          if (wEPeff > 0.) {
            wEPeff = 1. / wEPeff;
            weff *= wEPeff;
          }
        }
      }

      fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccRef * weff, 1);
    }

    auto posSlicedTracks = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negSlicedTracks = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (resoSwitchVals[PHI][kUseParticle]) {
      resurrectPhi(posSlicedTracks, negSlicedTracks, collision, kaonPlus, kaonMinus, phiMom, massKaPlus, HIST("hPhiMass_sparse"));
    }

    if (cfgUseLsPhi) {
      likeSignPhi(posSlicedTracks, collision, massKaPlus, HIST("hLsPhiMass_sparse"));
      likeSignPhi(negSlicedTracks, collision, massKaPlus, HIST("hLsPhiMass_sparse"));
    }

    // ---------------------- Analyzing the V0s
    for (auto const& v0s : V0s) {
      if (resoSwitchVals[K0][kUseParticle]) {
        if (selectionK0(collision, v0s) == true)
          histos.fill(HIST("hK0s"), 1);
      }
      if (resoSwitchVals[LAMBDA][kUseParticle]) {
        if (selectionLambda(collision, v0s) == true)
          histos.fill(HIST("hLambdas"), 1);
      }
    } // End of v0 loop

    // Filling cumulant profiles
    double r = fRndm->Rndm();
    int bootId = static_cast<int>(r * 10);

    for (auto i = 0; i < static_cast<int>(corrconfigs.size()); ++i) {
      if (resoSwitchVals[PHI][kUseParticle] && corrconfigs.at(i).Head.starts_with("Phi")) {
        int pIndex = findComponent(phiV2, Form("h%spt", corrconfigs.at(i).Head.c_str()));
        fillProfileBoot3D(corrconfigs.at(i), phiV2[pIndex], cent, fPhiMassAxis);

        if (cfgUseBootStrap) {
          fillProfileBoot3D(corrconfigs.at(i), phiBoot[bootId][pIndex], cent, fPhiMassAxis);
        }
      } // end of phi condition

      if (resoSwitchVals[K0][kUseParticle] && corrconfigs.at(i).Head.starts_with("K0")) {
        int pIndex = findComponent(k0V2, Form("h%spt", corrconfigs.at(i).Head.c_str()));
        fillProfileBoot3D(corrconfigs.at(i), k0V2[pIndex], cent, fK0MassAxis);

        if (cfgUseBootStrap) {
          fillProfileBoot3D(corrconfigs.at(i), k0Boot[bootId][pIndex], cent, fK0MassAxis);
        }
      } // end of K0 condition

      if (resoSwitchVals[LAMBDA][kUseParticle] && (corrconfigs.at(i).Head.starts_with("Lam") || corrconfigs.at(i).Head.starts_with("AnLam"))) {
        int pIndex = findComponent(lambdaV2, Form("h%spt", corrconfigs.at(i).Head.c_str()));
        fillProfileBoot3D(corrconfigs.at(i), lambdaV2[pIndex], cent, fLambdaMassAxis);

        if (cfgUseBootStrap) {
          fillProfileBoot3D(corrconfigs.at(i), lambdaBoot[bootId][pIndex], cent, fLambdaMassAxis);
        }
      } // end of lambda condition

      if (configs.GetHeads()[i].starts_with("Ref")) {
        int pIndex = findComponent(refV2, Form("h%s", corrconfigs.at(i).Head.c_str()));
        fillProfileBoot(corrconfigs.at(i), refV2[pIndex], cent);

        if (cfgUseBootStrap) {
          fillProfileBoot(corrconfigs.at(i), refBoot[bootId][pIndex], cent);
        }
      } // end of ref condition
    } // end of loop over correlation configurations
  } // end of processReso
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ResonancesGfwFlow>(cfgc)};
}
