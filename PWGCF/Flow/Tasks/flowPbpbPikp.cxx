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

/// \file flowPbpbPikp.cxx
/// \brief PID flow using the generic framework
/// \author Preet Bhanjan Pati <preet.bhanjan.pati@cern.ch>

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWConfig.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/GenericFramework/Core/GFWWeightsList.h"

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
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::analysis::genericframework
{
GFWRegions regions;
GFWCorrConfigs configs;
} // namespace o2::analysis::genericframework

using namespace o2::analysis::genericframework;

struct FlowPbpbPikp {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCrossRows, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputRunByRun, bool, true, "Fill and output NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgTpcCut, float, 2.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")

  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Minimum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Maximum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgUseGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 5, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, false, "Use track density efficiency correction")

  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaVtxz, bool, false, "Use Phi, Eta, VertexZ dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiPtCent, bool, false, "Use Phi, Pt, Centrality dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseWeightPhiEtaPt, bool, true, "Use Phi, Eta, Pt dependent NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgUseStrictPID, bool, true, "Use strict PID cuts for TPC")
  O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  O2_DEFINE_CONFIGURABLE(cfgUseAsymmetricPID, bool, false, "Use asymmetric PID cuts")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, true, "Use Nch for multiplicity selection instead of centrality")

  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTofNsigmaCut{"cfgTofNsigmaCut", std::vector<double>{1.5, 1.5, 1.5, -1.5, -1.5, -1.5}, "TOF n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> cfgItsNsigmaCut{"cfgItsNsigmaCut", std::vector<double>{3, 3, 3, -3, -3, -3}, "ITS n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<double>> cfgTpcNsigmaCut{"cfgTpcNsigmaCut", std::vector<double>{10, 10, 10, -10, -10, -10}, "TOF n-sigma cut for pions_posNsigma, kaons_posNsigma, protons_posNsigma, pions_negNsigma, kaons_negNsigma, protons_negNsigma"};
  Configurable<std::vector<int>> cfgUseEventCuts{"cfgUseEventCuts", std::vector<int>{1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0}, "Switch for various event cuts [Filtered Events, Sel8, kNoTimeFrameBorder, kNoITSROFrameBorder, kNoSameBunchPileup, kIsGoodZvtxFT0vsPV, kNoCollInTimeRangeStandard, kIsGoodITSLayersAll, kNoCollInRofStandard, kNoHighMultCollInPrevRof, Occupancy, Multiplicity correlation, T0AV0A 3 sigma cut, kIsVertexITSTPC, kTVXinTRD]"};

  Configurable<GFWRegions> cfgRegions{"cfgRegions", {{"refN08", "refP08", "full", "poiN", "olN", "poiP", "olP", "poi", "ol", "poiNpi", "olNpi", "poiPpi", "olPpi", "poifullpi", "olfullpi", "poiNka", "olNka", "poiPka", "olPka", "poifullka", "olfullka", "poiNpr", "olNpr", "poiPpr", "olPpr", "poifullpr", "olfullpr"}, {-0.8, 0.4, -0.8, -0.8, -0.8, 0.4, 0.4, -0.8, -0.8, -0.8, -0.8, 0.4, 0.4, -0.8, -0.8, -0.8, -0.8, 0.4, 0.4, -0.8, -0.8, -0.8, -0.8, 0.4, 0.4, -0.8, -0.8}, {-0.4, 0.8, 0.8, -0.4, -0.4, 0.8, 0.8, 0.8, 0.8, -0.4, -0.4, 0.8, 0.8, 0.8, 0.8, -0.4, -0.4, 0.8, 0.8, 0.8, 0.8, -0.4, -0.4, 0.8, 0.8, 0.8, 0.8}, {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 128, 256, 128, 256, 128, 256, 2, 16, 2, 16, 2, 16, 4, 32, 4, 32, 4, 32, 8, 64, 8, 64, 8, 64}}, "Configurations for GFW regions"};
  Configurable<GFWCorrConfigs> cfgCorrConfig{"cfgCorrConfig", {{"full {2 -2}", "full {2 -2}", "full {2 -2}", "full {2 -2}", "refN08 {2} refP08 {-2}", "refN08 {2} refP08 {-2}", "refN08 {2} refP08 {-2}", "refN08 {2} refP08 {-2}", "refP08 {-2} refN08 {2}", "refP08 {-2} refN08 {2}", "refP08 {-2} refN08 {2}", "refP08 {-2} refN08 {2}", "full {2 2 -2 -2}", "full {2 2 -2 -2}", "full {2 2 -2 -2}", "full {2 2 -2 -2}", "poi full | ol {2 -2}", "poifullpi full | olfullpi {2 -2}", "poifullka full | olfullka {2 -2}", "poifullpr full | olfullpr {2 -2}", "poiN refN08 | olN {2} refP08 {-2}", "poiNpi refN08 | olNpi {2} refP08 {-2}", "poiNka refN08 | olNka {2} refP08 {-2}", "poiNpr refN08 | olNpr {2} refP08 {-2}", "poiP refP08 | olP {2} refN08 {-2}", "poiPpi refP08 | olPpi {2} refN08 {-2}", "poiPka refP08 | olPka {2} refN08 {-2}", "poiPpr refP08 | olPpr {2} refN08 {-2}", "poi full | ol {2 2 -2 -2}", "poifullpi full | olfullpi {2 2 -2 -2}", "poifullka full | olfullka {2 2 -2 -2}", "poifullpr full | olfullpr {2 2 -2 -2}", "refN08 {2 2} refP08 {-2 -2}", "refP08 {-2 -2} refN08 {2 2}", "poiNka refN08 | olNka {2 2} refP08 {-2 -2}", "poiPka refP08 | olPka {2 2} refN08 {-2 -2}"}, {"ChFull22", "PiFull22", "KaFull22", "PrFull22", "Ch08FGap22", "Pi08FGap22", "Ka08FGap22", "Pr08FGap22", "Ch08BGap22", "Pi08BGap22", "Ka08BGap22", "Pr08BGap22", "ChFull24", "PiFull24", "KaFull24", "PrFull24", "ChFull22", "PiFull22", "KaFull22", "PrFull22", "Ch08FGap22", "Pi08FGap22", "Ka08FGap22", "Pr08FGap22", "Ch08BGap22", "Pi08BGap22", "Ka08BGap22", "Pr08BGap22", "ChFull24", "PiFull24", "KaFull24", "PrFull24", "Ka08FGap24", "Ka08BGap24", "Ka08FGap24", "Ka08BGap24"}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, "Configurations for each correlation to calculate"};

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {16, -0.8, 0.8}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {10000, 0, 1000}, "axis for TPC signal"};
  ConfigurableAxis axisTOFbeta{"axisTOFbeta", {200, 0, 2}, "axis for TOF beta"};
  ConfigurableAxis axisNch{"axisNch", {200, 2000, 4000}, "N_{ch}"};

  std::vector<double> tofNsigmaCut;
  std::vector<double> itsNsigmaCut;
  std::vector<double> tpcNsigmaCut;
  std::vector<int> eventCuts;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  enum OutputSpecies {
    hRef = 0,
    hCharge,
    hPion,
    hKaon,
    hProton,
    kCount_OutputSpecies
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
    kTpcCrossedRowsCut,
    kNumPions,
    kNumKaons,
    kNumProtons
  };

  int lastRunNumer = -1;
  std::vector<int> runNumbers;
  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  // Local track density correction - Copy from flowTask.cxx
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
    eventCuts = cfgUseEventCuts;

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan.value);

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

    const AxisSpec axisCentForQA{100, 0, 100, "centrality (%)"};

    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    histos.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hPhiWeighted", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hPt", "", {HistType::kTH1D, {axisPt}});
    histos.add("c22_full_ch", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_full_ch_Nch", "", {HistType::kTProfile, {axisNch}});

    histos.add("TpcdEdx", "", {HistType::kTH2D, {axisPt, axisTPCsignal}});
    histos.add("TofBeta", "", {HistType::kTH2D, {axisPt, axisTOFbeta}});

    histos.add("globalTracks_centT0C", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});

    histos.add("TofTpcNsigma_before", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaTPC, axisNsigmaTOF, axisPt}}});
    if (!cfgUseItsPID)
      histos.add("TofTpcNsigma_after", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaTPC, axisNsigmaTOF, axisPt}}});

    histos.add("TofItsNsigma_before", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaITS, axisNsigmaTOF, axisPt}}});
    if (cfgUseItsPID)
      histos.add("TofItsNsigma_after", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaITS, axisNsigmaTOF, axisPt}}});

    histos.add("partCount", "", {HistType::kTHnSparseD, {{axisParticles, axisMultiplicity, axisPt}}});

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

    histos.add("hTrackCount", "Number of Tracks;; Count", {HistType::kTH1D, {{10, -0.5, 9.5}}});
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kFilteredTracks + 1, "Filtered track");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kUseGlobalTracks + 1, "Global tracks");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kUsePvContributor + 1, "PV contributor");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kItsClustersCut + 1, "ITS clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kHasTpcSignal + 1, "TPC signal");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kTpcClustersCut + 1, "TPC clusters");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kTpcCrossedRowsCut + 1, "TPC crossed rows");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kNumPions + 1, "Pions");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kNumKaons + 1, "Kaons");
    histos.get<TH1>(HIST("hTrackCount"))->GetXaxis()->SetBinLabel(kNumProtons + 1, "Protons");

    if (cfgOutputNUAWeights && !cfgOutputRunByRun) {
      histos.add<TH3>("NUA/hPhiEtaVtxz_ref", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_ch", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_pi", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_ka", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_pr", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});

      histos.add<TH3>("NUA/hPhiPtCent_ref", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
      histos.add<TH3>("NUA/hPhiPtCent_ch", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
      histos.add<TH3>("NUA/hPhiPtCent_pi", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
      histos.add<TH3>("NUA/hPhiPtCent_ka", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
      histos.add<TH3>("NUA/hPhiPtCent_pr", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});

      histos.add<TH3>("NUA/hPhiEtaPt_ref", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_ch", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_pi", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_ka", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
      histos.add<TH3>("NUA/hPhiEtaPt_pr", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
    }

    if (!cfgAcceptance.value.empty()) {
      if (cfgUseWeightPhiEtaVtxz) {
        histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_pi_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_ka_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
        histos.add<TH3>("PhiCorrected/hPhiEtaVtxz_pr_corrd", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      }
      if (cfgUseWeightPhiPtCent) {
        histos.add<TH3>("PhiCorrected/hPhiPtCent_pi_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
        histos.add<TH3>("PhiCorrected/hPhiPtCent_ka_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
        histos.add<TH3>("PhiCorrected/hPhiPtCent_pr_corrd", ";#varphi;p_{T};Cent", {HistType::kTH3D, {axisPhi, axisPt, axisMultiplicity}});
      }
      if (cfgUseWeightPhiEtaPt) {
        histos.add<TH3>("PhiCorrected/hPhiEtaPt_pi_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
        histos.add<TH3>("PhiCorrected/hPhiEtaPt_ka_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
        histos.add<TH3>("PhiCorrected/hPhiEtaPt_pr_corrd", ";#varphi;#eta;p_{T}", {HistType::kTH3D, {axisPhi, axisEta, axisPt}});
      }
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    // Defining the regions
    for (auto i(0); i < regions.GetSize(); ++i) {
      fGFW->AddRegion(regions.GetNames()[i], regions.GetEtaMin()[i], regions.GetEtaMax()[i], (regions.GetpTDifs()[i]) ? nPtBins + 1 : 1, regions.GetBitmasks()[i]);
    }

    // Defining the correlators
    for (auto i = 0; i < configs.GetSize(); ++i) {
      corrconfigs.push_back(fGFW->GetCorrelatorConfig(configs.GetCorrs()[i], configs.GetHeads()[i], configs.GetpTDifs()[i]));
    }
    if (corrconfigs.empty())
      LOGF(error, "Configuration contains vectors of different size - check the GFWCorrConfig configurable");
    fGFW->CreateRegions();

    // Defining the flow container
    TObjArray* oba = new TObjArray();
    addConfigObjectsToObjArray(oba, corrconfigs);

    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);

    if (!cfgUseNch) {
      fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    } else {
      fFC->Initialize(oba, axisNch, cfgNbootstrap);
    }

    delete oba;

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

    tofNsigmaCut = cfgTofNsigmaCut;
    itsNsigmaCut = cfgItsNsigmaCut;
    tpcNsigmaCut = cfgTpcNsigmaCut;
  } // End of init()

  enum Particles {
    PIONS,
    KAONS,
    PROTONS
  };

  void addConfigObjectsToObjArray(TObjArray* oba, const std::vector<GFW::CorrConfig>& configs)
  {
    for (auto it = configs.begin(); it != configs.end(); ++it) {
      if (it->pTDif) {
        std::string suffix = "_ptDiff";
        for (auto i = 0; i < fPtAxis->GetNbins(); ++i) {
          std::string index = Form("_pt_%i", i + 1);
          oba->Add(new TNamed(it->Head.c_str() + index, it->Head.c_str() + suffix));
        }
      } else {
        oba->Add(new TNamed(it->Head.c_str(), it->Head.c_str()));
      }
    }
  }

  template <typename TTrack>
  bool selectionTrack(const TTrack& track)
  {
    histos.fill(HIST("hTrackCount"), kFilteredTracks); // Filtered tracks
    if (cfgUseGlobalTrack && !(track.isGlobalTrack())) {
      return 0;
    }
    if (cfgUseGlobalTrack)
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

  template <typename TCollision, typename TTrack>
  void fillQA(const TCollision collision, const TTrack track, int pidIndex, double wacc)
  {
    histos.fill(HIST("partCount"), pidIndex - 1, collision.centFT0C(), track.pt());
    switch (pidIndex) {
      case 1:
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), pidIndex - 1, track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), pidIndex - 1, itsResponse.nSigmaITS<o2::track::PID::Pion>(track), track.tofNSigmaPi(), track.pt());
        histos.fill(HIST("hTrackCount"), kNumPions); // Pion count
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaVtxz)
          histos.fill(HIST("PhiCorrected/hPhiEtaVtxz_pi_corrd"), track.phi(), track.eta(), collision.posZ(), wacc); // pion weights
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiPtCent)
          histos.fill(HIST("PhiCorrected/hPhiPtCent_pi_corrd"), track.phi(), track.pt(), collision.centFT0C(), wacc);
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaPt)
          histos.fill(HIST("PhiCorrected/hPhiEtaPt_pi_corrd"), track.phi(), track.eta(), track.pt(), wacc);
        break;
      case 2:
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), pidIndex - 1, track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), pidIndex - 1, itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), track.tofNSigmaKa(), track.pt());
        histos.fill(HIST("hTrackCount"), kNumKaons); // Kaon count
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaVtxz)
          histos.fill(HIST("PhiCorrected/hPhiEtaVtxz_ka_corrd"), track.phi(), track.eta(), collision.posZ(), wacc); // kaon weights
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiPtCent)
          histos.fill(HIST("PhiCorrected/hPhiPtCent_ka_corrd"), track.phi(), track.pt(), collision.centFT0C(), wacc);
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaPt)
          histos.fill(HIST("PhiCorrected/hPhiEtaPt_ka_corrd"), track.phi(), track.eta(), track.pt(), wacc);
        break;
      case 3:
        if (!cfgUseItsPID)
          histos.fill(HIST("TofTpcNsigma_after"), pidIndex - 1, track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
        if (cfgUseItsPID)
          histos.fill(HIST("TofItsNsigma_after"), pidIndex - 1, itsResponse.nSigmaITS<o2::track::PID::Proton>(track), track.tofNSigmaPr(), track.pt());
        histos.fill(HIST("hTrackCount"), kNumProtons); // Proton count
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaVtxz)
          histos.fill(HIST("PhiCorrected/hPhiEtaVtxz_pr_corrd"), track.phi(), track.eta(), collision.posZ(), wacc); // proton weights
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiPtCent)
          histos.fill(HIST("PhiCorrected/hPhiPtCent_pr_corrd"), track.phi(), track.pt(), collision.centFT0C(), wacc);
        if (!cfgAcceptance.value.empty() && cfgUseWeightPhiEtaPt)
          histos.fill(HIST("PhiCorrected/hPhiEtaPt_pr_corrd"), track.phi(), track.eta(), track.pt(), wacc);
        break;
    } // end of switch
  } // end of fillQA

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

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;            // Choose which nSigma to use: TPC or ITS
    std::vector<double> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

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

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
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
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        histos.fill(tarName, fPtAxis->GetBinCenter(i), val, dnx);
    }
    return;
  }

  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    if (!corrconf.pTDif) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0) {
        return;
      }
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
    }
    return;
  }

  void createRunByRunHistos(int runNumber)
  {
    if (cfgOutputNUAWeights) {
      std::vector<std::shared_ptr<TH3>> tH3s(kCount_OutputSpecies);
      tH3s[hRef] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ref", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[hCharge] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ch", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[hPion] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pi", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[hKaon] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ka", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      tH3s[hProton] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pr", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, axisEta, axisVertex}});
      th3sList.insert(std::make_pair(runNumber, tH3s));
    }
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

      mAcceptance[hCharge] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ch", timestamp);
      if (mAcceptance[hCharge])
        LOGF(info, "Loaded acceptance weights from %s_ch (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hCharge]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ch (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hCharge]);

      mAcceptance[hPion] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_pi", timestamp);
      if (mAcceptance[hPion])
        LOGF(info, "Loaded acceptance weights from %s_pi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPion]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_pi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPion]);

      mAcceptance[hKaon] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ka", timestamp);
      if (mAcceptance[hKaon])
        LOGF(info, "Loaded acceptance weights from %s_ka (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hKaon]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ka (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hKaon]);

      mAcceptance[hProton] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_pr", timestamp);
      if (mAcceptance[hProton])
        LOGF(info, "Loaded acceptance weights from %s_pr (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hProton]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_pr (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hProton]);
    }

    correctionsLoaded = true;
  }

  template <typename TTrack, typename TCollision>
  double getAcceptance(TTrack track, const TCollision collision, int index)
  { // 0 = ref, 1 = ch, 2 = pi, 3 = ka, 4 = pr
    if (index < 0 || index >= kCount_OutputSpecies) {
      return 1;
    }
    double wacc = 1;
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();

    if ((cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaPt && cfgUseWeightPhiPtCent) || (cfgUseWeightPhiEtaVtxz && cfgUseWeightPhiEtaPt)) {
      LOGF(fatal, "Only one of the three weight options can be used at a time");
    }

    if (!mAcceptance.empty() && correctionsLoaded) {
      if (!mAcceptance[index]) {
        LOGF(fatal, "Acceptance weights not loaded for index %d", index);
        return 1;
      }
      if (cfgUseWeightPhiEtaVtxz)
        wacc = mAcceptance[index]->getNUA(track.phi(), track.eta(), vtxz);
      if (cfgUseWeightPhiPtCent)
        wacc = mAcceptance[index]->getNUA(track.phi(), track.pt(), cent);
      if (cfgUseWeightPhiEtaPt)
        wacc = mAcceptance[index]->getNUA(track.phi(), track.eta(), track.pt());
    }
    return wacc;
  }

  template <typename TTrack, typename TCollision>
  void fillWeights(const TTrack track, const TCollision collision, const int& pid_index, const int& run)
  {
    double cent = collision.centFT0C();
    double vtxz = collision.posZ();
    double pt = track.pt();
    bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

    if (cfgOutputRunByRun) {
      if (withinPtRef && !pid_index)
        th3sList[run][hRef]->Fill(track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
      if (withinPtPOI)
        th3sList[run][hCharge + pid_index]->Fill(track.phi(), track.eta(), vtxz); // charged and id'ed particle weights
    } else {
      if (withinPtRef && !pid_index) {
        histos.fill(HIST("NUA/hPhiEtaVtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
        histos.fill(HIST("NUA/hPhiPtCent_ref"), track.phi(), track.pt(), cent);
        histos.fill(HIST("NUA/hPhiEtaPt_ref"), track.phi(), track.eta(), track.pt());
      }

      if (withinPtPOI) {
        switch (pid_index) {
          case 0:
            histos.fill(HIST("NUA/hPhiEtaVtxz_ch"), track.phi(), track.eta(), vtxz); // charged particle weights
            histos.fill(HIST("NUA/hPhiPtCent_ch"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_ch"), track.phi(), track.eta(), track.pt());
            break;
          case 1:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pi"), track.phi(), track.eta(), vtxz); // pion weights
            histos.fill(HIST("NUA/hPhiPtCent_pi"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_pi"), track.phi(), track.eta(), track.pt());
            break;
          case 2:
            histos.fill(HIST("NUA/hPhiEtaVtxz_ka"), track.phi(), track.eta(), vtxz); // kaon weights
            histos.fill(HIST("NUA/hPhiPtCent_ka"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_ka"), track.phi(), track.eta(), track.pt());
            break;
          case 3:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pr"), track.phi(), track.eta(), vtxz); // proton weights
            histos.fill(HIST("NUA/hPhiPtCent_pr"), track.phi(), track.pt(), cent);
            histos.fill(HIST("NUA/hPhiEtaPt_pr"), track.phi(), track.eta(), track.pt());
            break;
        }
      }
    }
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

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks)
  {
    int nTot = tracks.size();
    if (nTot < 1)
      return;

    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    const auto cent = collision.centFT0C();

    if (!selectionEvent(collision, nTot, cent))
      return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    if (cfgOutputRunByRun && runNumber != lastRunNumer) {
      lastRunNumer = runNumber;
      if (std::find(runNumbers.begin(), runNumbers.end(), runNumber) == runNumbers.end()) {
        // if run number is not in the preconfigured list, create new output histograms for this run
        createRunByRunHistos(runNumber);
        runNumbers.push_back(runNumber);
      }
    }

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), cent);
    histos.fill(HIST("globalTracks_centT0C"), cent, nTot);
    fGFW->Clear();

    float weff = 1;
    int pidIndex;
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
        bool withinPtRef = (cfgCutPtMin < track.pt()) && (track.pt() < cfgCutPtMax); // within RF pT range
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
      if (!selectionTrack(track))
        continue;

      double pt = track.pt();
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPt"), pt);

      histos.fill(HIST("TpcdEdx"), pt, track.tpcSignal());
      histos.fill(HIST("TofBeta"), pt, track.beta());

      histos.fill(HIST("TofTpcNsigma_before"), PIONS, track.tpcNSigmaPi(), track.tofNSigmaPi(), pt);
      histos.fill(HIST("TofTpcNsigma_before"), KAONS, track.tpcNSigmaKa(), track.tofNSigmaKa(), pt);
      histos.fill(HIST("TofTpcNsigma_before"), PROTONS, track.tpcNSigmaPr(), track.tofNSigmaPr(), pt);

      histos.fill(HIST("TofItsNsigma_before"), PIONS, itsResponse.nSigmaITS<o2::track::PID::Pion>(track), track.tofNSigmaPi(), pt);
      histos.fill(HIST("TofItsNsigma_before"), KAONS, itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), track.tofNSigmaKa(), pt);
      histos.fill(HIST("TofItsNsigma_before"), PROTONS, itsResponse.nSigmaITS<o2::track::PID::Proton>(track), track.tofNSigmaPr(), pt);

      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

      pidIndex = cfgUseAsymmetricPID ? getNsigmaPIDAssymmetric(track) : getNsigmaPIDTpcTof(track);

      weff = 1; // Initializing weff for each track
      // NUA weights
      if (cfgOutputNUAWeights)
        fillWeights(track, collision, pidIndex, runNumber);

      if (!withinPtPOI && !withinPtRef)
        return;
      double waccRef = getAcceptance(track, collision, 0);
      double waccPOI = withinPtPOI ? getAcceptance(track, collision, pidIndex + 1) : getAcceptance(track, collision, 0);
      if (withinPtRef && withinPtPOI && pidIndex)
        waccRef = waccPOI; // if particle is both (then it's overlap), override ref with POI

      // Track density correction
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
      } // end of track density correction loop

      if (withinPtRef) {
        histos.fill(HIST("hPhiWeighted"), track.phi(), waccRef);
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccRef * weff, 1);
      }
      if (withinPtPOI) {
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 128);
      }
      if (withinPtPOI && withinPtRef) {
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 256);
      }

      if (pidIndex) {
        fillQA(collision, track, pidIndex, waccPOI);
        if (withinPtPOI)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 1 << (pidIndex));
        if (withinPtPOI && withinPtRef)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 1 << (pidIndex + 3));
      }
    } // track loop ends

    // Filling cumulants with ROOT TProfile
    fillProfile(corrconfigs.at(0), HIST("c22_full_ch"), cent);
    fillProfile(corrconfigs.at(0), HIST("c22_full_ch_Nch"), nTot);

    if (!cfgUseNch) {
      for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
        fillFC(corrconfigs.at(l_ind), cent, lRandom);
      }
    } else {
      for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
        fillFC(corrconfigs.at(l_ind), nTot, lRandom);
      }
    }

  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlowPbpbPikp>(cfgc)};
}
