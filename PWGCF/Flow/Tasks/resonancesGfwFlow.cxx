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

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <utility>
#include <array>
#include <string>
#include <memory>

#include "Math/Vector4D.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/GenericFramework/Core/GFWWeightsList.h"

#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"

#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

namespace
{
std::shared_ptr<TProfile> phiC22Boot[10];
std::shared_ptr<TProfile> phiC24Boot[10];
std::shared_ptr<TProfile> k0C22Boot[10];
std::shared_ptr<TProfile> k0C24Boot[10];
std::shared_ptr<TProfile> lambdaC22Boot[10];
std::shared_ptr<TProfile> lambdaC24Boot[10];
std::shared_ptr<TProfile> anLambdaC22Boot[10];
std::shared_ptr<TProfile> anLambdaC24Boot[10];

std::shared_ptr<TProfile3D> phiD22PtBoot[10];
std::shared_ptr<TProfile3D> phiD24PtBoot[10];
std::shared_ptr<TProfile3D> k0D22PtBoot[10];
std::shared_ptr<TProfile3D> k0D24PtBoot[10];
std::shared_ptr<TProfile3D> lambdaD22PtBoot[10];
std::shared_ptr<TProfile3D> lambdaD24PtBoot[10];
std::shared_ptr<TProfile3D> anLambdaD22PtBoot[10];
std::shared_ptr<TProfile3D> anLambdaD24PtBoot[10];
} // namespace

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct ResonancesGfwFlow {
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgTpcNsigmaCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofNsigmaCut, float, 3.0f, "TOF N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 0, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(cfgCutTOFBeta, float, 0.0, "cut TOF beta")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancy, int, 3000, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgUseGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgFakeKaonCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
  O2_DEFINE_CONFIGURABLE(cfgRapidityCut, float, 0.5, "Rapidity cut for the reconstructed particles")
  O2_DEFINE_CONFIGURABLE(cfgUseCosPA, bool, false, "Use Pointing angle for resonances")
  O2_DEFINE_CONFIGURABLE(cfgPhiCosPA, float, 0.04f, "Minimum Pointing angle for Phi")
  O2_DEFINE_CONFIGURABLE(cfgK0CosPA, float, 0.97f, "Minimum Pointing angle for K0")
  O2_DEFINE_CONFIGURABLE(cfgLambdaCosPA, float, 0.995f, "Minimum Pointing angle for Lambda")
  O2_DEFINE_CONFIGURABLE(cfgUseV0Radius, bool, true, "Use V0 radius for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgLambdaRadiusMin, float, 0.5f, "Minimum Lambda radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgLambdaRadiusMax, float, 200.0f, "Maximum Lambda radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgK0RadiusMin, float, 0.5f, "Minimum K0 radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgK0RadiusMax, float, 200.0f, "Maximum K0 radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgUseProperLifetime, bool, false, "Use proper lifetime for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgK0LifeTime, float, 20.0f, "Maximum lifetime for K0 in cm")
  O2_DEFINE_CONFIGURABLE(cfgLambdaLifeTime, float, 30.0f, "Maximum lifetime for Lambda in cm")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgDCALambdaPosToPVMin, float, 0.1f, "minimum DCA to PV for Lambda positive track")
  O2_DEFINE_CONFIGURABLE(cfgDCALambdaNegToPVMin, float, 0.25f, "minimum DCA to PV for Lambda negative track")
  O2_DEFINE_CONFIGURABLE(cfgDCALambdaBetDaug, int, 1, "Maximum DCA between Lambda daughter tracks")
  O2_DEFINE_CONFIGURABLE(cfgDCAK0PosToPVMin, float, 0.06f, "minimum DCA to PV for K0 positive track")
  O2_DEFINE_CONFIGURABLE(cfgDCAK0NegToPVMin, float, 0.06f, "minimum DCA to PV for K0 negative track")
  O2_DEFINE_CONFIGURABLE(cfgDCAK0BetDaug, int, 1, "Maximum DCA between K0 daughter tracks")
  O2_DEFINE_CONFIGURABLE(cfgMassLambdaMin, float, 1.1f, "minimum lambda mass")
  O2_DEFINE_CONFIGURABLE(cfgMassLambdaMax, float, 1.16f, "maximum lambda mass")
  O2_DEFINE_CONFIGURABLE(cfgMassK0Min, float, 0.44f, "minimum K0short mass")
  O2_DEFINE_CONFIGURABLE(cfgMassK0Max, float, 0.56f, "maximum K0short mass")
  O2_DEFINE_CONFIGURABLE(cfgUseMCCLambda, bool, false, "Use mass cross check for lambda")
  O2_DEFINE_CONFIGURABLE(cfgUseMCCK0, bool, false, "Use mass cross check for K0")

  O2_DEFINE_CONFIGURABLE(cfgNPhiMassBins, double, 70, "Invasriant mass bins for phi")
  O2_DEFINE_CONFIGURABLE(cfgNK0MassBins, double, 70, "Invasriant mass bins for K0")
  O2_DEFINE_CONFIGURABLE(cfgNLambdaMassBins, double, 70, "Invasriant mass bins for lambda")

  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")

  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaVtxz, bool, true, "Use Phi, Eta, VertexZ dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiPtCent, bool, false, "Use Phi, Pt, Centrality dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaPt, bool, false, "Use Phi, Eta, Pt dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseBootStrap, bool, true, "Use bootstrap for error estimation")

  // Defining configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};
  ConfigurableAxis axisPhiMass{"axisPhiMass", {cfgNPhiMassBins, 0.99, 1.06}, "axis for invariant mass distibution for Phi"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {cfgNK0MassBins, cfgMassK0Min, cfgMassK0Max}, "axis for invariant mass distibution for K0"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {cfgNLambdaMassBins, cfgMassLambdaMin, cfgMassLambdaMax}, "axis for invariant mass distibution for Lambda"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {10000, 0, 1000}, "axis for TPC signal"};
  ConfigurableAxis axisTOFsignal{"axisTOFsignal", {10000, 0, 1000}, "axis for TOF signal"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  SliceCache cache;
  Partition<AodTracksWithoutBayes> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<AodTracksWithoutBayes> negTracks = aod::track::signed1Pt < 0.0f;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TAxis* fPhiMassAxis;
  TAxis* fK0MassAxis;
  TAxis* fLambdaMassAxis;
  TRandom3* fRndm = new TRandom3(0);

  enum OutputSpecies {
    hRef = 0,
    hK0 = 1,
    hLambda = 2,
    hAnLambda = 3,
    hPhi = 4,
    kCount_OutputSpecies
  };

  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan.value);

    AxisSpec singleCount = {1, 0, 1};

    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});

    histos.add("KaplusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("KaminusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("KaplusTOF", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("KaminusTOF", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("hPhiPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hPhiEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});

    histos.add("PlusTPC_L", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("MinusTPC_L", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("PlusTOF_L", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("MinusTOF_L", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, axisPt, axisMultiplicity}}});
    histos.add("PlusTPC_AL", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("MinusTPC_AL", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("PlusTOF_AL", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("MinusTOF_AL", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("hAntiLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hAntiLambdaEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hAntiLambdaMass_sparse", "", {HistType::kTHnSparseF, {{axisLambdaMass, axisPt, axisMultiplicity}}});
    histos.add("hLambdas", "", {HistType::kTH1D, {singleCount}});

    histos.add("PlusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("MinusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("PlusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("MinusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hK0Mass_sparse", "", {HistType::kTHnSparseF, {{axisK0Mass, axisPt, axisMultiplicity}}});
    histos.add("hK0s", "", {HistType::kTH1D, {singleCount}});

    histos.add("Phic22", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("Phic24", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("Phiv22pt", "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}});
    histos.add("Phiv24pt", "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}});

    histos.add("K0c22", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("K0c24", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("K0v22pt", "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}});
    histos.add("K0v24pt", "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}});

    histos.add("Lambdac22", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("Lambdac24", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("Lambdav22pt", "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
    histos.add("Lambdav24pt", "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});

    histos.add("AnLambdac22", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("AnLambdac24", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("AnLambdav22pt", "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
    histos.add("AnLambdav24pt", "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});

    histos.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{8, 0, 8}}});
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "After sel8");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
    histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(8, "After Occupancy");

    histos.add("hLambdaCount", "Number of Lambda;; Count", {HistType::kTH1D, {{11, 0, 11}}});
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
    histos.get<TH1>(HIST("hLambdaCount"))->GetXaxis()->SetBinLabel(11, "Mass cross check");

    histos.add("hK0Count", "Number of K0;; Count", {HistType::kTH1D, {{11, 0, 11}}});
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
    histos.get<TH1>(HIST("hK0Count"))->GetXaxis()->SetBinLabel(11, "Mass cross check");

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

    if (cfgUseBootStrap) {
      for (int i = 0; i < 10; i++) {
        phiC22Boot[i] = histos.add<TProfile>(Form("BootStrap/Phic22_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        phiC24Boot[i] = histos.add<TProfile>(Form("BootStrap/Phic24_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        k0C22Boot[i] = histos.add<TProfile>(Form("BootStrap/k0c22_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        k0C24Boot[i] = histos.add<TProfile>(Form("BootStrap/k0c24_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        lambdaC22Boot[i] = histos.add<TProfile>(Form("BootStrap/lambdac22_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        lambdaC24Boot[i] = histos.add<TProfile>(Form("BootStrap/lambdac24_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        anLambdaC22Boot[i] = histos.add<TProfile>(Form("BootStrap/anlambdac22_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});
        anLambdaC24Boot[i] = histos.add<TProfile>(Form("BootStrap/anlambdac24_bootstrap_%d", i), "", {HistType::kTProfile, {axisMultiplicity}});

        phiD22PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/Phid22pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}});
        phiD24PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/Phid24pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisPhiMass, axisMultiplicity}});
        k0D22PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/k0d22pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}});
        k0D24PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/k0d24pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisK0Mass, axisMultiplicity}});
        lambdaD22PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/lambdad22pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
        lambdaD24PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/lambdad24pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
        anLambdaD22PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/anlambdad22pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
        anLambdaD24PtBoot[i] = histos.add<TProfile3D>(Form("BootStrap/anlambdad24pt_bootstrap_%d", i), "", {HistType::kTProfile3D, {axisPt, axisLambdaMass, axisMultiplicity}});
      }
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    fPhiMassAxis = new TAxis(cfgNPhiMassBins, 0.99, 1.06);
    fK0MassAxis = new TAxis(cfgNK0MassBins, cfgMassK0Min, cfgMassK0Max);
    fLambdaMassAxis = new TAxis(cfgNLambdaMassBins, cfgMassLambdaMin, cfgMassLambdaMax);

    int nPhisPtMassBins = nPtBins * cfgNPhiMassBins;
    int nK0sPtMassBins = nPtBins * cfgNK0MassBins;
    int nLambdasPtMassBins = nPtBins * cfgNLambdaMassBins;

    //********** Defining the regions  **********
    // reference particles
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);

    // phi
    fGFW->AddRegion("poiNphi", -0.8, -0.4, 1 + nPhisPtMassBins, 2);
    fGFW->AddRegion("poiPphi", -0.8, -0.4, 1 + nPhisPtMassBins, 2);
    fGFW->AddRegion("olNphi", -0.8, -0.4, 1 + nPhisPtMassBins, 32);
    fGFW->AddRegion("olPphi", -0.8, -0.4, 1 + nPhisPtMassBins, 32);

    // kshort
    fGFW->AddRegion("poiNk0", -0.8, -0.4, 1 + nK0sPtMassBins, 4);
    fGFW->AddRegion("poiPk0", -0.8, -0.4, 1 + nK0sPtMassBins, 4);
    fGFW->AddRegion("olNk0", -0.8, -0.4, 1 + nK0sPtMassBins, 64);
    fGFW->AddRegion("olPk0", -0.8, -0.4, 1 + nK0sPtMassBins, 64);

    // lambda
    fGFW->AddRegion("poiNlam", -0.8, -0.4, 1 + nLambdasPtMassBins, 8);
    fGFW->AddRegion("poiPlam", -0.8, -0.4, 1 + nLambdasPtMassBins, 8);
    fGFW->AddRegion("olNlam", -0.8, -0.4, 1 + nLambdasPtMassBins, 128);
    fGFW->AddRegion("olPlam", -0.8, -0.4, 1 + nLambdasPtMassBins, 128);

    // antilambda
    fGFW->AddRegion("poiNantilam", -0.8, -0.4, 1 + nLambdasPtMassBins, 16);
    fGFW->AddRegion("poiPantilam", -0.8, -0.4, 1 + nLambdasPtMassBins, 16);
    fGFW->AddRegion("olNantilam", -0.8, -0.4, 1 + nLambdasPtMassBins, 256);
    fGFW->AddRegion("olPantilam", -0.8, -0.4, 1 + nLambdasPtMassBins, 256);

    //********** Defining the correlations  ************
    // NAMING CONVENTION:
    //    F: Forward --> REF from negative eta + POI from negative eta correlated to REF from positive eta
    //    B: Backward --> REF from negative eta correlated to REF from positive eta + POI from positive eta

    //--------- reference particles
    // Forward and Backward correlations are the same for reference particles
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "PhiF08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "KsF08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "LamF08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "AnLamF08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "PhiF08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "KsF08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "LamF08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "AnLamF08Gap24", kFALSE));

    //--------- pt differential pois
    // Forward
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNphi refN08 | olNphi {2} refP08 {-2}", "PhiF08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk0 refN08 | olNk0 {2} refP08 {-2}", "KsF08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNlam refN08 | olNlam {2} refP08 {-2}", "LamF08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNantilam refN08 | olNantilam {2} refP08 {-2}", "AnLamF08Gap22", kTRUE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNphi refN08 | olNphi {2 2} refP08 {-2 -2}", "PhiF08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk0 refN08 | olNk0 {2 2} refP08 {-2 -2}", "KsF08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNlam refN08 | olNlam {2 2} refP08 {-2 -2}", "LamF08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNantilam refN08 | olNantilam {2 2} refP08 {-2 -2}", "AnLamF08Gap24", kTRUE));

    // Backward
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPphi refP08 | olPphi {2} refN08 {-2}", "PhiB08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPk0 refP08 | olPk0 {2} refN08 {-2}", "KsB08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPlam refP08 | olPlam {2} refN08 {-2}", "LamB08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPantilam refP08 | olPantilam {2} refN08 {-2}", "AnLamB08Gap22", kTRUE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPphi refP08 | olPphi {2 2} refN08 {-2 -2}", "PhiB08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPk0 refP08 | olPk0 {2 2} refN08 {-2 -2}", "KsB08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPlam refP08 | olPlam {2 2} refN08 {-2 -2}", "LamB08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPantilam refP08 | olPantilam {2 2} refN08 {-2 -2}", "AnLamB08Gap24", kTRUE));

    fGFW->CreateRegions();
  }

  template <char... chars>
  void fillResoProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent, TAxis* partaxis)
  {
    double dnx, val;
    if (!corrconf.pTDif) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0)
        return;
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        histos.fill(tarName, cent, val, dnx);
      return;
    }
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      for (int j = 1; j <= partaxis->GetNbins(); j++) {
        dnx = fGFW->Calculate(corrconf, ((i - 1) * partaxis->GetNbins()) + (j - 1), kTRUE).real();
        if (dnx == 0)
          continue;
        val = fGFW->Calculate(corrconf, ((i - 1) * partaxis->GetNbins()) + (j - 1), kFALSE).real() / dnx;
        if (std::fabs(val) < 1)
          histos.fill(tarName, fPtAxis->GetBinCenter(i), partaxis->GetBinCenter(j), cent, val, dnx);
      }
    }
    return;
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
  bool selectionPair(const TTrack1& track1, const TTrack2& track2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = track1.pt();
    pt2 = track2.pt();
    pz1 = track1.pz();
    pz2 = track2.pz();
    p1 = track1.p();
    p2 = track2.p();
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgUseCosPA && angle < cfgPhiCosPA) {
      return false;
    }
    return true;
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
  bool selectionTrack(const TTrack& track)
  {
    if (cfgUseGlobalTrack && !(track.isGlobalTrack() && track.isPVContributor() && track.itsNCls() > cfgITScluster && track.tpcNClsFound() > cfgTpcCluster && track.hasTPC())) {
      return false;
    }
    if (!cfgUseGlobalTrack && !(track.isPVContributor() && track.itsNCls() > cfgITScluster && track.hasTPC())) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr())};
    int pid = -1;
    float nsigma = cfgTpcNsigmaCut;

    // Choose which nSigma to use
    std::array<float, 3> nSigmaToUse = (track.pt() >= cfgTofPtCut && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;
    if (track.pt() >= cfgTofPtCut && !track.hasTOF())
      return -1;

    // Select particle with the lowest nsigma
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        pid = i;
        nsigma = std::abs(nSigmaToUse[i]);
      }
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

      mAcceptance[hRef] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ref", timestamp);
      if (mAcceptance[hRef])
        LOGF(info, "Loaded acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hRef]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hRef]);

      mAcceptance[hK0] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_k0", timestamp);
      if (mAcceptance[hK0])
        LOGF(info, "Loaded acceptance weights from %s_k0 (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hK0]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_k0 (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hK0]);

      mAcceptance[hLambda] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_lambda", timestamp);
      if (mAcceptance[hLambda])
        LOGF(info, "Loaded acceptance weights from %s_lambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hLambda]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_lambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hLambda]);

      mAcceptance[hAnLambda] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_hanlambda", timestamp);
      if (mAcceptance[hAnLambda])
        LOGF(info, "Loaded acceptance weights from %s_anlambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hAnLambda]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_anlambda (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hAnLambda]);

      mAcceptance[hPhi] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_phi", timestamp);
      if (mAcceptance[hPhi])
        LOGF(info, "Loaded acceptance weights from %s_phi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPhi]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_phi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPhi]);
    }

    correctionsLoaded = true;
  }

  template <typename TTrack, typename TCollision>
  double getAcceptance(TTrack track, const TCollision collision, int pid_index_reso)
  { // 0 = ref, 1 = k0, 2 = lambda, 3 = anti-lambda, 4 = phi
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

  template <typename TTrack, typename TCollision>
  void fillWeights(const TTrack track, const TCollision collision, const int& pid_index_reso)
  {
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();
    double pt = track.pt();
    bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

    if (withinPtRef && !pid_index_reso)
      histos.fill(HIST("NUA/hPhiEtaVtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
    histos.fill(HIST("NUA/hPhiPtCent_ref"), track.phi(), track.pt(), cent);
    histos.fill(HIST("NUA/hPhiEtaPt_ref"), track.phi(), track.eta(), track.pt());

    if (withinPtPOI) {
      switch (pid_index_reso) {
        case hK0:
          histos.fill(HIST("NUA/hPhiEtaVtxz_k0"), track.phi(), track.eta(), vtxz); // K0 weights
          histos.fill(HIST("NUA/hPhiPtCent_k0"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_k0"), track.phi(), track.eta(), track.pt());
          break;
        case hLambda:
          histos.fill(HIST("NUA/hPhiEtaVtxz_lambda"), track.phi(), track.eta(), vtxz); // Lambda weights
          histos.fill(HIST("NUA/hPhiPtCent_lambda"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_lambda"), track.phi(), track.eta(), track.pt());
          break;
        case hAnLambda:
          histos.fill(HIST("NUA/hPhiEtaVtxz_anlambda"), track.phi(), track.eta(), vtxz); // Anti-Lambda weights
          histos.fill(HIST("NUA/hPhiPtCent_anlambda"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_anlambda"), track.phi(), track.eta(), track.pt());
          break;
        case hPhi:
          histos.fill(HIST("NUA/hPhiEtaVtxz_phi"), track.phi(), track.eta(), vtxz); // Phi weights
          histos.fill(HIST("NUA/hPhiPtCent_phi"), track.phi(), track.pt(), cent);
          histos.fill(HIST("NUA/hPhiEtaPt_phi"), track.phi(), track.eta(), track.pt());
          break;
      }
    }
  }

  template <typename TTrack, typename vector, char... chars>
  void resurrectPhi(TTrack trackplus, TTrack trackminus, vector plusdaug, vector minusdaug, vector mom, double plusmass, const ConstStr<chars...>& hist, const double cent)
  {
    for (auto const& [partplus, partminus] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trackplus, trackminus))) {
      if (getNsigmaPID(partplus) != 2)
        continue;
      if (getNsigmaPID(partminus) != 2)
        continue;
      if (isFakeKaon(partplus) || isFakeKaon(partminus))
        continue;
      if (!selectionPair(partplus, partminus))
        continue;
      if (!selectionTrack(partplus) || !selectionTrack(partminus))
        continue;

      histos.fill(HIST("KaplusTPC"), partplus.pt(), partplus.tpcNSigmaKa());
      histos.fill(HIST("KaplusTOF"), partplus.pt(), partplus.tofNSigmaKa());
      histos.fill(HIST("KaminusTPC"), partminus.pt(), partminus.tpcNSigmaKa());
      histos.fill(HIST("KaminusTOF"), partminus.pt(), partminus.tofNSigmaKa());

      plusdaug = ROOT::Math::PxPyPzMVector(partplus.px(), partplus.py(), partplus.pz(), plusmass);
      minusdaug = ROOT::Math::PxPyPzMVector(partminus.px(), partminus.py(), partminus.pz(), plusmass);
      mom = plusdaug + minusdaug;

      double pt = mom.Pt();
      double invMass = mom.M();
      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);

      if (std::abs(mom.Rapidity()) < cfgRapidityCut) {
        histos.fill(hist, invMass, pt, cent);
        histos.fill(HIST("hPhiPhi"), mom.Phi());
        histos.fill(HIST("hPhiEta"), mom.Eta());

        if (withinPtPOI)
          fGFW->Fill(mom.Eta(), ((fPtAxis->FindBin(pt) - 1) * fPhiMassAxis->GetNbins()) + (fPhiMassAxis->FindBin(invMass) - 1), mom.Phi(), 1.0, 2);
        if (withinPtPOI && withinPtRef)
          fGFW->Fill(mom.Eta(), ((fPtAxis->FindBin(pt) - 1) * fPhiMassAxis->GetNbins()) + (fPhiMassAxis->FindBin(invMass) - 1), mom.Phi(), 1.0, 32);
      }
    }
    return;
  }

  template <typename TTrack>
  bool selectionV0Daughter(TTrack const& track, int pid) // pid 1: proton, pid 0: pion
  {
    if (track.tpcNClsFound() < cfgTpcCluster)
      return false;
    if (!track.hasTPC())
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPr()) > cfgTpcNsigmaCut)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPi()) > cfgTpcNsigmaCut)
      return false;

    return true;
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
    if (postrack.pt() < 0.15 || negtrack.pt() < 0.15)
      return false;

    histos.fill(HIST("hLambdaCount"), 1.5);
    if (mlambda > cfgMassLambdaMin && mlambda < cfgMassLambdaMax)
      isL = true;
    if (mantilambda > cfgMassLambdaMin && mantilambda < cfgMassLambdaMax)
      isAL = true;

    if (!isL && !isAL) {
      return false;
    }
    histos.fill(HIST("hLambdaCount"), 2.5);

    // Rapidity correction
    if (candidate.yLambda() > 0.5)
      return false;
    histos.fill(HIST("hLambdaCount"), 3.5);
    // DCA cuts for lambda and antilambda
    if (isL) {
      if (std::abs(candidate.dcapostopv()) < cfgDCALambdaPosToPVMin || std::abs(candidate.dcanegtopv()) < cfgDCALambdaNegToPVMin)
        return false;
    }
    if (isAL) {
      if (std::abs(candidate.dcapostopv()) < cfgDCALambdaNegToPVMin || std::abs(candidate.dcanegtopv()) < cfgDCALambdaPosToPVMin)
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > cfgDCALambdaBetDaug)
      return false;
    histos.fill(HIST("hLambdaCount"), 5.5);
    // v0 radius cuts
    if (cfgUseV0Radius && (candidate.v0radius() < cfgLambdaRadiusMin || candidate.v0radius() > cfgLambdaRadiusMax))
      return false;
    histos.fill(HIST("hLambdaCount"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < cfgLambdaCosPA)
      return false;
    histos.fill(HIST("hLambdaCount"), 7.5);
    // Proper lifetime
    if (cfgUseProperLifetime && candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgLambdaLifeTime)
      return false;
    histos.fill(HIST("hLambdaCount"), 8.5);
    if (isL) {
      if (!selectionV0Daughter(postrack, 1) || !selectionV0Daughter(negtrack, 0))
        return false;
    }
    if (isAL) {
      if (!selectionV0Daughter(postrack, 0) || !selectionV0Daughter(negtrack, 1))
        return false;
    }
    histos.fill(HIST("hLambdaCount"), 9.5);
    // Mass cross check
    if (cfgUseMCCLambda && std::abs(massK0Short - 0.497614) < 0.01)
      return false;
    histos.fill(HIST("hLambdaCount"), 10.5);
    bool withinPtPOI = (cfgCutPtPOIMin < candidate.pt()) && (candidate.pt() < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < candidate.pt()) && (candidate.pt() < cfgCutPtMax);

    float weff = 1;

    if (isL) {
      if (cfgOutputNUAWeights)
        fillWeights(candidate, collision, hLambda);

      double waccPOI = getAcceptance(candidate, collision, hLambda);
      if (withinPtPOI)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mlambda) - 1), candidate.phi(), waccPOI * weff, 8);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mlambda) - 1), candidate.phi(), waccPOI * weff, 128);

      histos.fill(HIST("hLambdaMass_sparse"), mlambda, candidate.pt(), collision.centFT0C());
      histos.fill(HIST("hLambdaPhi"), candidate.phi());
      histos.fill(HIST("hLambdaEta"), candidate.eta());
      histos.fill(HIST("PlusTPC_L"), postrack.pt(), postrack.tpcNSigmaKa());
      histos.fill(HIST("PlusTOF_L"), postrack.pt(), postrack.tofNSigmaKa());
      histos.fill(HIST("MinusTPC_L"), negtrack.pt(), negtrack.tpcNSigmaKa());
      histos.fill(HIST("MinusTOF_L"), negtrack.pt(), negtrack.tofNSigmaKa());
    }
    if (isAL) {
      if (cfgOutputNUAWeights)
        fillWeights(candidate, collision, hAnLambda);

      double waccPOI = getAcceptance(candidate, collision, hAnLambda);
      if (withinPtPOI)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mantilambda) - 1), candidate.phi(), waccPOI * weff, 16);
      if (withinPtPOI && withinPtRef)
        fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fLambdaMassAxis->GetNbins()) + (fLambdaMassAxis->FindBin(mantilambda) - 1), candidate.phi(), waccPOI * weff, 256);

      histos.fill(HIST("hAntiLambdaMass_sparse"), mantilambda, candidate.pt(), collision.centFT0C());
      histos.fill(HIST("hAntiLambdaPhi"), candidate.phi());
      histos.fill(HIST("hAntiLambdaEta"), candidate.eta());
      histos.fill(HIST("PlusTPC_AL"), postrack.pt(), postrack.tpcNSigmaKa());
      histos.fill(HIST("PlusTOF_AL"), postrack.pt(), postrack.tofNSigmaKa());
      histos.fill(HIST("MinusTPC_AL"), negtrack.pt(), negtrack.tpcNSigmaKa());
      histos.fill(HIST("MinusTOF_AL"), negtrack.pt(), negtrack.tofNSigmaKa());
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
    if (postrack.pt() < 0.15 || negtrack.pt() < 0.15)
      return false;
    histos.fill(HIST("hK0Count"), 1.5);
    if (mk0 < cfgMassK0Min && mk0 > cfgMassK0Max)
      return false;
    histos.fill(HIST("hK0Count"), 2.5);
    // Rapidity correction
    if (candidate.yK0Short() > 0.5)
      return false;
    histos.fill(HIST("hK0Count"), 3.5);
    // DCA cuts for K0short
    if (std::abs(candidate.dcapostopv()) < cfgDCAK0PosToPVMin || std::abs(candidate.dcanegtopv()) < cfgDCAK0NegToPVMin)
      return false;
    histos.fill(HIST("hK0Count"), 4.5);
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAK0BetDaug)
      return false;
    histos.fill(HIST("hK0Count"), 5.5);
    // v0 radius cuts
    if (cfgUseV0Radius && (candidate.v0radius() < cfgK0RadiusMin || candidate.v0radius() > cfgK0RadiusMax))
      return false;
    histos.fill(HIST("hK0Count"), 6.5);
    // cosine pointing angle cuts
    if (candidate.v0cosPA() < cfgK0CosPA)
      return false;
    histos.fill(HIST("hK0Count"), 7.5);
    // Proper lifetime
    if (cfgUseProperLifetime && candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0Short > cfgK0LifeTime)
      return false;
    histos.fill(HIST("hK0Count"), 8.5);
    if (!selectionV0Daughter(postrack, 0) || !selectionV0Daughter(negtrack, 0))
      return false;
    histos.fill(HIST("hK0Count"), 9.5);
    // Mass cross check
    if (cfgUseMCCK0 && std::abs(massLambda - 1.11568) < 0.005)
      return false;
    if (cfgUseMCCK0 && std::abs(massLambda - 1.11568) < 0.005)
      return false;
    histos.fill(HIST("hK0Count"), 10.5);
    bool withinPtPOI = (cfgCutPtPOIMin < candidate.pt()) && (candidate.pt() < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < candidate.pt()) && (candidate.pt() < cfgCutPtMax);

    if (cfgOutputNUAWeights)
      fillWeights(candidate, collision, hK0);

    float weff = 1;
    double waccPOI = getAcceptance(candidate, collision, hK0);

    if (withinPtPOI)
      fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fK0MassAxis->GetNbins()) + (fK0MassAxis->FindBin(mk0) - 1), candidate.phi(), waccPOI * weff, 4);
    if (withinPtPOI && withinPtRef)
      fGFW->Fill(candidate.eta(), ((fPtAxis->FindBin(candidate.pt()) - 1) * fK0MassAxis->GetNbins()) + (fK0MassAxis->FindBin(mk0) - 1), candidate.phi(), waccPOI * weff, 64);

    histos.fill(HIST("hK0Mass_sparse"), mk0, candidate.pt(), collision.centFT0C());
    histos.fill(HIST("hK0Phi"), candidate.phi());
    histos.fill(HIST("hK0Eta"), candidate.eta());
    histos.fill(HIST("PlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaKa());
    histos.fill(HIST("PlusTOF_K0"), postrack.pt(), postrack.tofNSigmaKa());
    histos.fill(HIST("MinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaKa());
    histos.fill(HIST("MinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaKa());

    return true;
  }

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, o2::aod::evsel::NumTracksInTimeRange>;
  ROOT::Math::PxPyPzMVector phiMom, kaonPlus, kaonMinus;
  double massKaPlus = o2::constants::physics::MassKPlus;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0Short = o2::constants::physics::MassK0Short;

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks, aod::V0Datas const& V0s)
  {
    histos.fill(HIST("hEventCount"), 0.5);
    int nTot = tracks.size();
    if (nTot < 1)
      return;

    if (!collision.sel8())
      return;
    histos.fill(HIST("hEventCount"), 1.5);

    if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return;
    histos.fill(HIST("hEventCount"), 2.5);

    if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return;
    histos.fill(HIST("hEventCount"), 3.5);

    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return;
    histos.fill(HIST("hEventCount"), 4.5);

    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    histos.fill(HIST("hEventCount"), 5.5);

    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return;
    histos.fill(HIST("hEventCount"), 6.5);

    int occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy > cfgCutOccupancy)
      return;
    histos.fill(HIST("hEventCount"), 7.5);

    const auto cent = collision.centFT0C();
    float vtxz = collision.posZ();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), cent);
    fGFW->Clear();

    float weff = 1;

    loadCorrections(bc); // load corrections for the each event

    for (auto const& track : tracks) {
      if (!selectionTrack(track))
        continue;
      double pt = track.pt();
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);

      if (withinPtRef)
        if (cfgOutputNUAWeights)
          fillWeights(track, collision, hRef);

      double waccRef = getAcceptance(track, collision, 0);
      fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccRef * weff, 1);
    }

    auto posSlicedTracks = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negSlicedTracks = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    resurrectPhi(posSlicedTracks, negSlicedTracks, kaonPlus, kaonMinus, phiMom, massKaPlus, HIST("hPhiMass_sparse"), cent);

    for (auto const& v0s : V0s) {
      if (selectionLambda(collision, v0s) == true)
        histos.fill(HIST("hLambdas"), 1);
      if (selectionK0(collision, v0s) == true)
        histos.fill(HIST("hK0s"), 1);
    } // End of v0 loop

    fillResoProfile(corrconfigs.at(0), HIST("Phic22"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(1), HIST("K0c22"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(2), HIST("Lambdac22"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(3), HIST("AnLambdac22"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(4), HIST("Phic24"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(5), HIST("K0c24"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(6), HIST("Lambdac24"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(7), HIST("AnLambdac24"), cent, fLambdaMassAxis);

    fillResoProfile(corrconfigs.at(8), HIST("Phiv22pt"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(9), HIST("K0v22pt"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(10), HIST("Lambdav22pt"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(11), HIST("AnLambdav22pt"), cent, fLambdaMassAxis);

    fillResoProfile(corrconfigs.at(12), HIST("Phiv24pt"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(13), HIST("K0v24pt"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(14), HIST("Lambdav24pt"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(15), HIST("AnLambdav24pt"), cent, fLambdaMassAxis);

    fillResoProfile(corrconfigs.at(16), HIST("Phiv22pt"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(17), HIST("K0v22pt"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(18), HIST("Lambdav22pt"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(19), HIST("AnLambdav22pt"), cent, fLambdaMassAxis);

    fillResoProfile(corrconfigs.at(20), HIST("Phiv24pt"), cent, fPhiMassAxis);
    fillResoProfile(corrconfigs.at(21), HIST("K0v24pt"), cent, fK0MassAxis);
    fillResoProfile(corrconfigs.at(22), HIST("Lambdav24pt"), cent, fLambdaMassAxis);
    fillResoProfile(corrconfigs.at(23), HIST("AnLambdav24pt"), cent, fLambdaMassAxis);

    // bootstraping
    if (cfgUseBootStrap) {
      double r = fRndm->Rndm();
      int bootId = static_cast<int>(r * 10);

      fillProfileBoot(corrconfigs.at(0), phiC22Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(1), k0C22Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(2), lambdaC22Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(3), anLambdaC22Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(4), phiC24Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(5), k0C24Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(6), lambdaC24Boot[bootId], cent);
      fillProfileBoot(corrconfigs.at(7), anLambdaC24Boot[bootId], cent);

      fillProfileBoot3D(corrconfigs.at(8), phiD22PtBoot[bootId], cent, fPhiMassAxis);
      fillProfileBoot3D(corrconfigs.at(9), k0D22PtBoot[bootId], cent, fK0MassAxis);
      fillProfileBoot3D(corrconfigs.at(10), lambdaD22PtBoot[bootId], cent, fLambdaMassAxis);
      fillProfileBoot3D(corrconfigs.at(11), anLambdaD22PtBoot[bootId], cent, fLambdaMassAxis);

      fillProfileBoot3D(corrconfigs.at(12), phiD24PtBoot[bootId], cent, fPhiMassAxis);
      fillProfileBoot3D(corrconfigs.at(13), k0D24PtBoot[bootId], cent, fK0MassAxis);
      fillProfileBoot3D(corrconfigs.at(14), lambdaD24PtBoot[bootId], cent, fLambdaMassAxis);
      fillProfileBoot3D(corrconfigs.at(15), anLambdaD24PtBoot[bootId], cent, fLambdaMassAxis);

      fillProfileBoot3D(corrconfigs.at(16), phiD22PtBoot[bootId], cent, fPhiMassAxis);
      fillProfileBoot3D(corrconfigs.at(17), k0D22PtBoot[bootId], cent, fK0MassAxis);
      fillProfileBoot3D(corrconfigs.at(18), lambdaD22PtBoot[bootId], cent, fLambdaMassAxis);
      fillProfileBoot3D(corrconfigs.at(19), anLambdaD22PtBoot[bootId], cent, fLambdaMassAxis);

      fillProfileBoot3D(corrconfigs.at(20), phiD24PtBoot[bootId], cent, fPhiMassAxis);
      fillProfileBoot3D(corrconfigs.at(21), k0D24PtBoot[bootId], cent, fK0MassAxis);
      fillProfileBoot3D(corrconfigs.at(22), lambdaD24PtBoot[bootId], cent, fLambdaMassAxis);
      fillProfileBoot3D(corrconfigs.at(23), anLambdaD24PtBoot[bootId], cent, fLambdaMassAxis);
    }
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ResonancesGfwFlow>(cfgc)};
}
