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

/// \file piDeFemtoSystematics.cxx
/// \brief Histogram-only pion-deuteron femtoscopy task for train subwagon cut variations
///
/// The task intentionally produces no derived AOD tables. A train subwagon
/// supplies one complete set of cuts. The default wagon and every systematic
/// subwagon therefore run the same event and pair workflow, while a systematic
/// subwagon changes exactly one configurable cut.

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/GenVector/Boost.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults>;
using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                                  aod::TrackSelection, aod::pidTPCFullDe, aod::pidTOFFullDe,
                                  aod::pidTPCFullPi, aod::pidTOFFullPi,
                                  aod::TOFSignal, aod::TOFEvTime>;

namespace
{
constexpr std::array<double, 6> BetheBlochDeDefault{-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09};
const std::vector<std::string> betheBlochParticleNames{"De"};
const std::vector<std::string> betheBlochParameterNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr std::array<float, 9> tpcRadii{{85.f, 105.f, 125.f, 145.f, 165.f, 185.f, 205.f, 225.f, 245.f}};
constexpr int NumberOfBetheBlochShapeParameters = 5;
constexpr int BetheBlochResolutionIndex = 5;
constexpr double AutomaticMagneticFieldThreshold = -998.;
constexpr float NominalL3Current = 30000.f;
constexpr float InvalidPhiStar = 999.f;
constexpr int ClosePairRadiusModeAtPV = 0;
constexpr int ClosePairRadiusModeAtSpecificRadius = 2;
using PairLorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;

enum SelectionStep {
  kAllTracks = 0,
  kTrackCuts,
  kPidCuts,
  kNSelectionSteps
};

enum PairChannel {
  kLikeSignMatter = 0,
  kLikeSignAntimatter,
  kUnlikeSignMatter,
  kUnlikeSignAntimatter,
  kNPairChannels
};
} // namespace

struct PiDeFemtoSystematics {
  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"event"};
    Configurable<float> vertexZMax{"vertexZMax", 10.f, "Maximum absolute collision z vertex"};
    Configurable<int> numberOfMixedEvents{"numberOfMixedEvents", 5, "Number of mixed events per event"};
  } eventCuts;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pion"};
    Configurable<float> ptMin{"ptMin", 0.14f, "Minimum pion pT"};
    Configurable<float> ptMax{"ptMax", 4.0f, "Maximum pion pT"};
    Configurable<float> etaMax{"etaMax", 0.8f, "Maximum absolute pion eta"};
    Configurable<int> itsInnerBarrelClustersMin{"itsInnerBarrelClustersMin", 3, "Minimum pion ITS inner-barrel clusters"};
    Configurable<int> itsClustersMin{"itsClustersMin", 7, "Minimum pion ITS clusters"};
    Configurable<int> tpcClustersMin{"tpcClustersMin", 80, "Minimum pion found TPC clusters"};
    Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 90, "Minimum pion crossed TPC rows"};
    Configurable<float> dcaXYOffset{"dcaXYOffset", 0.004f, "Pion DCAxy cut offset"};
    Configurable<float> dcaXYPtCoefficient{"dcaXYPtCoefficient", 0.013f, "Pion DCAxy inverse-pT coefficient"};
    Configurable<float> dcaZOffset{"dcaZOffset", 0.004f, "Pion DCAz cut offset"};
    Configurable<float> dcaZPtCoefficient{"dcaZPtCoefficient", 0.013f, "Pion DCAz inverse-pT coefficient"};
    Configurable<float> combinedPidMomentumMin{"combinedPidMomentumMin", 0.5f, "Momentum above which pion TPC+TOF PID is required"};
    Configurable<float> tpcNsigmaMax{"tpcNsigmaMax", 3.f, "Maximum absolute pion TPC n-sigma below the TOF threshold"};
    Configurable<float> combinedNsigmaMax{"combinedNsigmaMax", 3.f, "Maximum pion combined TPC+TOF n-sigma"};
  } pionCuts;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"deuteron"};
    Configurable<float> ptMin{"ptMin", 0.6f, "Minimum deuteron pT"};
    Configurable<float> ptMax{"ptMax", 2.0f, "Maximum deuteron pT"};
    Configurable<float> etaMax{"etaMax", 0.8f, "Maximum absolute deuteron eta"};
    Configurable<int> tpcClustersMin{"tpcClustersMin", 110, "Minimum deuteron found TPC clusters"};
    Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 100, "Minimum deuteron crossed TPC rows"};
    Configurable<int> itsClustersMin{"itsClustersMin", 5, "Minimum deuteron ITS clusters"};
    Configurable<int> itsInnerBarrelClustersMin{"itsInnerBarrelClustersMin", 1, "Minimum deuteron ITS inner-barrel clusters"};
    Configurable<int> sharedTPCClustersMax{"sharedTPCClustersMax", 160, "Maximum deuteron shared TPC clusters"};
    Configurable<float> sharedTPCFractionMax{"sharedTPCFractionMax", 1.f, "Maximum deuteron shared TPC-cluster fraction"};
    Configurable<float> tpcInnerParamMin{"tpcInnerParamMin", 0.f, "Minimum deuteron TPC inner parameter"};
    Configurable<float> tofMomentumMin{"tofMomentumMin", 1.2f, "TPC inner parameter above which combined TPC+TOF PID is required"};
    Configurable<float> nsigmaMax{"nsigmaMax", 2.5f, "Common deuteron TPC, ITS and combined PID threshold"};
    Configurable<bool> requireIndividualNsigma{"requireIndividualNsigma", false, "Also apply individual TPC and TOF cuts in the combined PID branch"};
    Configurable<float> dcaXYOffset{"dcaXYOffset", 0.004f, "Deuteron DCAxy cut offset"};
    Configurable<float> dcaXYPtCoefficient{"dcaXYPtCoefficient", 0.013f, "Deuteron DCAxy inverse-pT coefficient"};
    Configurable<float> dcaZOffset{"dcaZOffset", 0.004f, "Deuteron DCAz cut offset"};
    Configurable<float> dcaZPtCoefficient{"dcaZPtCoefficient", 0.013f, "Deuteron DCAz inverse-pT coefficient"};
  } deuteronCuts;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pair"};
    Configurable<bool> enableClosePairRejection{"enableClosePairRejection", true, "Enable pion-deuteron close-pair rejection"};
    Configurable<float> closePairDeltaEtaMax{"closePairDeltaEtaMax", 0.01f, "CPR ellipse delta-eta radius"};
    Configurable<float> closePairDeltaPhiMax{"closePairDeltaPhiMax", 0.01f, "CPR ellipse delta-phi-star radius"};
    Configurable<int> closePairRadiusMode{"closePairRadiusMode", 1, "CPR mode: 0 PV, 1 average TPC phi-star, 2 one TPC radius"};
    Configurable<float> closePairSpecificRadius{"closePairSpecificRadius", 85.f, "TPC radius in cm for CPR mode 2"};
  } pairCuts;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"analysis"};
    Configurable<bool> useBetheBlochDeuteronNsigma{"useBetheBlochDeuteronNsigma", false, "Compute deuteron TPC n-sigma from the configured Bethe-Bloch parameters"};
    Configurable<float> lowKstarYieldMax{"lowKstarYieldMax", 0.3f, "Upper kstar used for the SE-yield stability check"};
  } analysisSettings;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"zorro"};
    Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Process a Zorro-selected deuteron skim"};
  } zorroSettings;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"ccdb"};
    Configurable<double> magneticField{"magneticField", -999., "Magnetic field; -999 reads it from CCDB"};
    Configurable<std::string> url{"url", "http://alice-ccdb.cern.ch", "CCDB URL"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "GRP object path"};
    Configurable<std::string> grpMagPath{"grpMagPath", "GLO/Config/GRPMagField", "GRP magnetic-field object path"};
    Configurable<int> materialCorrection{"materialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Material correction used by the optional pair DCA fitter"};
  } ccdbSettings;

  struct : ConfigurableGroup {
    // cppcheck-suppress unusedStructMember
    std::string prefix{"pidCalibration"};
    Configurable<LabeledArray<double>> betheBlochParameters{"betheBlochParameters",
                                                            {BetheBlochDeDefault.data(), 1, 6, betheBlochParticleNames, betheBlochParameterNames},
                                                            "Deuteron TPC Bethe-Bloch parameters"};
  } pidCalibration;

  ConfigurableAxis mixingVertexAxis{"mixingVertexAxis", {30, -10., 10.}, "Event-mixing z-vertex bins"};
  ConfigurableAxis mixingCentralityAxis{"mixingCentralityAxis", {40, 0., 100.}, "Event-mixing centrality bins"};
  ConfigurableAxis kstarAxis{"kstarAxis", {300, 0., 3.}, "kstar axis"};
  ConfigurableAxis mtAxis{"mtAxis", {240, 0.8, 3.2}, "pair mT axis"};
  ConfigurableAxis centralityAxis{"centralityAxis", {100, 0., 100.}, "centrality axis"};

  using MixingBinning = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  MixingBinning mixingBinning{{mixingVertexAxis, mixingCentralityAxis}, true};
  SliceCache cache;
  SameKindPair<CollisionsFull, TrackCandidates, MixingBinning> mixedEventPairs{
    mixingBinning, eventCuts.numberOfMixedEvents, -1, &cache};
  Preslice<TrackCandidates> tracksPerCollision = aod::track::collisionId;

  HistogramRegistry registry{"PiDeFemtoSystematics", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  o2::vertexing::DCAFitterN<2> pairFitter;

  std::array<float, 6> betheBlochDe{};
  int currentRunNumber{0};
  float magneticField{0.f};

  void init(InitContext const&)
  {
    const AxisSpec channelAxis{kNPairChannels, -0.5, static_cast<double>(kNPairChannels) - 0.5, "pair channel"};
    registry.add("Pairs/hSE", "Same-event pion-deuteron pairs", HistType::kTHnSparseF,
                 {kstarAxis, mtAxis, centralityAxis, channelAxis});
    registry.add("Pairs/hME", "Mixed-event pion-deuteron pairs", HistType::kTHnSparseF,
                 {kstarAxis, mtAxis, centralityAxis, channelAxis});
    registry.add("Pairs/hSELowKstarYield", "SE low-kstar yield;pair channel;Entries",
                 HistType::kTH1F, {channelAxis});
    registry.add("Pairs/hMELowKstarYield", "ME low-kstar yield;pair channel;Entries",
                 HistType::kTH1F, {channelAxis});
    registry.add("Pairs/hCPRBefore", "CPR before selection;#Delta#eta;#Delta#varphi^{*}",
                 HistType::kTH2F, {{300, -0.15, 0.15}, {400, -0.2, 0.2}});
    registry.add("Pairs/hCPRAfter", "CPR after selection;#Delta#eta;#Delta#varphi^{*}",
                 HistType::kTH2F, {{300, -0.15, 0.15}, {400, -0.2, 0.2}});
    registry.add("Event/hCounter", "Event selection;step;Entries",
                 HistType::kTH1F, {{3, -0.5, 2.5}});
    registry.add("Event/hVertexZ", "Selected collision vertex;z (cm);Entries",
                 HistType::kTH1F, {{200, -20., 20.}});
    registry.add("Event/hCentrality", "Selected collision centrality;FT0C percentile;Entries",
                 HistType::kTH1F, {{100, 0., 100.}});
    registry.add("Pion/hSelection", "Pion selection;step;Entries",
                 HistType::kTH1F, {{kNSelectionSteps, -0.5, static_cast<double>(kNSelectionSteps) - 0.5}});
    registry.add("Pion/hPt", "Selected pions;signed p_{T} (GeV/#it{c});Entries",
                 HistType::kTH1F, {{400, -5., 5.}});
    registry.add("Pion/hEta", "Selected pions;#eta;Entries",
                 HistType::kTH1F, {{200, -1., 1.}});
    registry.add("Pion/hTPCNsigma", "Selected pions;signed p_{T} (GeV/#it{c});n#sigma_{TPC}",
                 HistType::kTH2F, {{400, -5., 5.}, {200, -5., 5.}});
    registry.add("Pion/hTOFNsigma", "Selected pions;signed p_{T} (GeV/#it{c});n#sigma_{TOF}",
                 HistType::kTH2F, {{280, -7., 7.}, {200, -5., 5.}});
    registry.add("Pion/hITSNsigma", "Selected pions;signed p_{T} (GeV/#it{c});n#sigma_{ITS}",
                 HistType::kTH2F, {{280, -7., 7.}, {120, -3., 3.}});
    registry.add("Pion/hTPCTOFNsigma", "Selected pions;signed p_{T} (GeV/#it{c});n#sigma_{TPC+TOF}",
                 HistType::kTH2F, {{280, -7., 7.}, {100, 0., 5.}});
    registry.add("Deuteron/hSelection", "Deuteron selection;step;Entries",
                 HistType::kTH1F, {{kNSelectionSteps, -0.5, static_cast<double>(kNSelectionSteps) - 0.5}});
    registry.add("Deuteron/hPt", "Selected deuterons;signed p_{T} (GeV/#it{c});Entries",
                 HistType::kTH1F, {{400, -5., 5.}});
    registry.add("Deuteron/hEta", "Selected deuterons;#eta;Entries",
                 HistType::kTH1F, {{200, -1., 1.}});
    registry.add("Deuteron/hTPCNsigma", "Selected deuterons;signed p_{T} (GeV/#it{c});n#sigma_{TPC}",
                 HistType::kTH2F, {{400, -5., 5.}, {200, -5., 5.}});
    registry.add("Deuteron/hTOFNsigma", "Selected deuterons;signed p_{T} (GeV/#it{c});n#sigma_{TOF}",
                 HistType::kTH2F, {{280, -7., 7.}, {200, -5., 5.}});
    registry.add("Deuteron/hITSNsigma", "Selected deuterons;signed p_{T} (GeV/#it{c});n#sigma_{ITS}",
                 HistType::kTH2F, {{280, -7., 7.}, {120, -3., 3.}});
    registry.add("Deuteron/hTPCTOFNsigma", "Selected deuterons;signed p_{T} (GeV/#it{c});n#sigma_{TPC+TOF}",
                 HistType::kTH2F, {{280, -7., 7.}, {100, 0., 5.}});
    registry.add("Deuteron/hTPCClusters", "Selected deuterons;TPC clusters;Entries",
                 HistType::kTH1F, {{161, -0.5, 160.5}});
    registry.add("Deuteron/hTPCCrossedRows", "Selected deuterons;TPC crossed rows;Entries",
                 HistType::kTH1F, {{161, -0.5, 160.5}});

    addCutConfigurationHistogram();
    setAxisLabels();

    zorroSummary.setObject(zorro.getZorroSummary());
    ccdb->setURL(ccdbSettings.url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    pairFitter.setPropagateToPCA(true);
    pairFitter.setMaxR(200.);
    pairFitter.setMinParamChange(1.e-3);
    pairFitter.setMinRelChi2Change(0.9);
    pairFitter.setMaxDZIni(1.e9);
    pairFitter.setMaxChi2(1.e9);
    pairFitter.setUseAbsDCA(true);
    pairFitter.setMatCorrType(
      static_cast<o2::base::Propagator::MatCorrType>(ccdbSettings.materialCorrection.value));

    for (int i = 0; i < NumberOfBetheBlochShapeParameters; ++i) {
      betheBlochDe[i] = pidCalibration.betheBlochParameters->get("De", Form("p%i", i));
    }
    betheBlochDe[BetheBlochResolutionIndex] = pidCalibration.betheBlochParameters->get("De", "resolution");
  }

  void setAxisLabels()
  {
    const std::array<const char*, kNPairChannels> channelLabels{
      "LS matter", "LS antimatter", "US matter", "US antimatter"};
    const auto setSparseChannelLabels = [&channelLabels](THnSparse* histogram) {
      for (int channel = 0; channel < kNPairChannels; ++channel) {
        histogram->GetAxis(3)->SetBinLabel(channel + 1, channelLabels[channel]);
      }
    };
    setSparseChannelLabels(registry.get<THnSparse>(HIST("Pairs/hSE")).get());
    setSparseChannelLabels(registry.get<THnSparse>(HIST("Pairs/hME")).get());

    const auto setYieldChannelLabels = [&channelLabels](TH1* histogram) {
      for (int channel = 0; channel < kNPairChannels; ++channel) {
        histogram->GetXaxis()->SetBinLabel(channel + 1, channelLabels[channel]);
      }
    };
    setYieldChannelLabels(registry.get<TH1>(HIST("Pairs/hSELowKstarYield")).get());
    setYieldChannelLabels(registry.get<TH1>(HIST("Pairs/hMELowKstarYield")).get());

    const std::array<const char*, kNSelectionSteps> selectionLabels{"all", "track cuts", "PID"};
    const auto setSelectionLabels = [&selectionLabels](TH1* histogram) {
      for (int step = 0; step < kNSelectionSteps; ++step) {
        histogram->GetXaxis()->SetBinLabel(step + 1, selectionLabels[step]);
      }
    };
    setSelectionLabels(registry.get<TH1>(HIST("Pion/hSelection")).get());
    setSelectionLabels(registry.get<TH1>(HIST("Deuteron/hSelection")).get());
    auto eventCounter = registry.get<TH1>(HIST("Event/hCounter"));
    eventCounter->GetXaxis()->SetBinLabel(1, "all");
    eventCounter->GetXaxis()->SetBinLabel(2, "sel8 and vertex");
    eventCounter->GetXaxis()->SetBinLabel(3, "Zorro");
  }

  void addCutConfigurationHistogram()
  {
    const std::vector<std::pair<std::string, double>> cuts{
      {"pionPtMin", pionCuts.ptMin.value},
      {"pionPtMax", pionCuts.ptMax.value},
      {"pionEtaMax", pionCuts.etaMax.value},
      {"pionITSInnerBarrelClustersMin", pionCuts.itsInnerBarrelClustersMin.value},
      {"pionITSClustersMin", pionCuts.itsClustersMin.value},
      {"pionTPCClustersMin", pionCuts.tpcClustersMin.value},
      {"pionTPCCrossedRowsMin", pionCuts.tpcCrossedRowsMin.value},
      {"pionDCAxyOffset", pionCuts.dcaXYOffset.value},
      {"pionDCAxyPtCoefficient", pionCuts.dcaXYPtCoefficient.value},
      {"pionDCAzOffset", pionCuts.dcaZOffset.value},
      {"pionDCAzPtCoefficient", pionCuts.dcaZPtCoefficient.value},
      {"pionCombinedPidMomentumMin", pionCuts.combinedPidMomentumMin.value},
      {"pionTPCNsigmaMax", pionCuts.tpcNsigmaMax.value},
      {"pionCombinedNsigmaMax", pionCuts.combinedNsigmaMax.value},
      {"deuteronPtMin", deuteronCuts.ptMin.value},
      {"deuteronPtMax", deuteronCuts.ptMax.value},
      {"deuteronEtaMax", deuteronCuts.etaMax.value},
      {"deuteronTPCClustersMin", deuteronCuts.tpcClustersMin.value},
      {"deuteronTPCCrossedRowsMin", deuteronCuts.tpcCrossedRowsMin.value},
      {"deuteronITSClustersMin", deuteronCuts.itsClustersMin.value},
      {"deuteronITSInnerBarrelClustersMin", deuteronCuts.itsInnerBarrelClustersMin.value},
      {"deuteronSharedTPCClustersMax", deuteronCuts.sharedTPCClustersMax.value},
      {"deuteronSharedTPCFractionMax", deuteronCuts.sharedTPCFractionMax.value},
      {"deuteronTPCInnerParamMin", deuteronCuts.tpcInnerParamMin.value},
      {"deuteronTOFMomentumMin", deuteronCuts.tofMomentumMin.value},
      {"deuteronNsigmaMax", deuteronCuts.nsigmaMax.value},
      {"requireIndividualNsigma", deuteronCuts.requireIndividualNsigma.value},
      {"deuteronDCAxyOffset", deuteronCuts.dcaXYOffset.value},
      {"deuteronDCAxyPtCoefficient", deuteronCuts.dcaXYPtCoefficient.value},
      {"deuteronDCAzOffset", deuteronCuts.dcaZOffset.value},
      {"deuteronDCAzPtCoefficient", deuteronCuts.dcaZPtCoefficient.value},
      {"enableCPR", pairCuts.enableClosePairRejection.value},
      {"cprDeltaEtaMax", pairCuts.closePairDeltaEtaMax.value},
      {"cprDeltaPhiMax", pairCuts.closePairDeltaPhiMax.value},
      {"cprRadiusMode", pairCuts.closePairRadiusMode.value},
      {"cprSpecificRadius", pairCuts.closePairSpecificRadius.value},
      {"useBetheBlochDeuteronNsigma", analysisSettings.useBetheBlochDeuteronNsigma.value},
      {"skimmedProcessing", zorroSettings.skimmedProcessing.value},
      {"vertexZMax", eventCuts.vertexZMax.value},
      {"numberOfMixedEvents", eventCuts.numberOfMixedEvents.value}};

    const AxisSpec cutAxis{static_cast<int>(cuts.size()), -0.5,
                           static_cast<double>(cuts.size()) - 0.5, "cut"};
    registry.add("Config/hCutValues", "Configured cut values;cut;value",
                 HistType::kTH1D, {cutAxis});
    auto histogram = registry.get<TH1>(HIST("Config/hCutValues"));
    for (std::size_t index = 0; index < cuts.size(); ++index) {
      histogram->GetXaxis()->SetBinLabel(static_cast<int>(index) + 1, cuts[index].first.c_str());
      histogram->SetBinContent(static_cast<int>(index) + 1, cuts[index].second);
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (currentRunNumber == bc.runNumber()) {
      return;
    }
    if (zorroSettings.skimmedProcessing.value) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDe");
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    auto* grp = ccdb->getForTimeStamp<o2::parameters::GRPObject>(
      ccdbSettings.grpPath.value, bc.timestamp());
    if (grp) {
      o2::base::Propagator::initFieldFromGRP(grp);
      magneticField = ccdbSettings.magneticField.value < AutomaticMagneticFieldThreshold ? grp->getNominalL3Field()
                                                                                         : ccdbSettings.magneticField.value;
    } else {
      auto* grpMag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(
        ccdbSettings.grpMagPath.value, bc.timestamp());
      if (!grpMag) {
        LOG(fatal) << "Could not load magnetic field from " << ccdbSettings.grpPath.value
                   << " or " << ccdbSettings.grpMagPath.value;
      }
      o2::base::Propagator::initFieldFromGRP(grpMag);
      magneticField = ccdbSettings.magneticField.value < AutomaticMagneticFieldThreshold
                        ? std::lround(5.f * grpMag->getL3Current() / NominalL3Current)
                        : ccdbSettings.magneticField.value;
    }
    currentRunNumber = bc.runNumber();
  }

  template <typename Collision>
  bool isCollisionSelected(Collision const& collision, bool fillQA)
  {
    if (fillQA) {
      registry.fill(HIST("Event/hCounter"), 0.);
    }
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    if (!collision.sel8() || std::abs(collision.posZ()) > eventCuts.vertexZMax.value) {
      return false;
    }
    if (zorroSettings.skimmedProcessing.value) {
      if (!zorro.isSelected(bc.globalBC())) {
        return false;
      }
      if (fillQA) {
        registry.fill(HIST("Event/hCounter"), 2.);
      }
    }
    if (fillQA) {
      registry.fill(HIST("Event/hCounter"), 1.);
      registry.fill(HIST("Event/hVertexZ"), collision.posZ());
      registry.fill(HIST("Event/hCentrality"), collision.centFT0C());
    }
    return true;
  }

  template <typename Track>
  bool passPionTrackCuts(Track const& track) const
  {
    const float absPt = std::abs(track.pt());
    if (absPt <= 0.f || absPt < pionCuts.ptMin.value || absPt > pionCuts.ptMax.value ||
        std::abs(track.eta()) > pionCuts.etaMax.value) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < pionCuts.itsInnerBarrelClustersMin.value ||
        track.itsNCls() < pionCuts.itsClustersMin.value ||
        track.tpcNClsFound() < pionCuts.tpcClustersMin.value ||
        track.tpcNClsCrossedRows() < pionCuts.tpcCrossedRowsMin.value) {
      return false;
    }
    const float dcaXYMax = pionCuts.dcaXYOffset.value +
                           pionCuts.dcaXYPtCoefficient.value / absPt;
    const float dcaZMax = pionCuts.dcaZOffset.value +
                          pionCuts.dcaZPtCoefficient.value / absPt;
    return std::abs(track.dcaXY()) <= dcaXYMax &&
           std::abs(track.dcaZ()) <= dcaZMax;
  }

  template <typename Track>
  bool passPionPID(Track const& track) const
  {
    const float tpcNsigma = track.tpcNSigmaPi();
    if (std::abs(track.p()) <= pionCuts.combinedPidMomentumMin.value) {
      return std::abs(tpcNsigma) <= pionCuts.tpcNsigmaMax.value;
    }
    if (!track.hasTOF()) {
      return false;
    }
    const float tofNsigma = track.tofNSigmaPi();
    if (std::hypot(tpcNsigma, tofNsigma) >
        pionCuts.combinedNsigmaMax.value) {
      return false;
    }
    return !deuteronCuts.requireIndividualNsigma.value ||
           (std::abs(tpcNsigma) <= pionCuts.combinedNsigmaMax.value &&
            std::abs(tofNsigma) <= pionCuts.combinedNsigmaMax.value);
  }

  template <typename Track>
  bool passDeuteronTrackCuts(Track const& track) const
  {
    const float absPt = std::abs(track.pt());
    if (absPt <= 0.f || absPt < deuteronCuts.ptMin.value ||
        absPt > deuteronCuts.ptMax.value ||
        std::abs(track.eta()) > deuteronCuts.etaMax.value) {
      return false;
    }
    if (track.tpcNClsFound() < deuteronCuts.tpcClustersMin.value ||
        track.tpcNClsCrossedRows() < deuteronCuts.tpcCrossedRowsMin.value ||
        track.tpcNClsShared() > deuteronCuts.sharedTPCClustersMax.value ||
        track.tpcFractionSharedCls() > deuteronCuts.sharedTPCFractionMax.value ||
        track.itsNCls() < deuteronCuts.itsClustersMin.value ||
        track.itsNClsInnerBarrel() < deuteronCuts.itsInnerBarrelClustersMin.value) {
      return false;
    }
    const float dcaXYMax = deuteronCuts.dcaXYOffset.value +
                           deuteronCuts.dcaXYPtCoefficient.value / absPt;
    const float dcaZMax = deuteronCuts.dcaZOffset.value +
                          deuteronCuts.dcaZPtCoefficient.value / absPt;
    return std::abs(track.dcaXY()) <= dcaXYMax &&
           std::abs(track.dcaZ()) <= dcaZMax;
  }

  template <typename Track>
  float deuteronTPCNsigma(Track const& track) const
  {
    if (!analysisSettings.useBetheBlochDeuteronNsigma.value) {
      return track.tpcNSigmaDe();
    }
    const float expected = o2::common::BetheBlochAleph(
      static_cast<float>(track.tpcInnerParam() /
                         o2::constants::physics::MassDeuteron),
      betheBlochDe[0], betheBlochDe[1], betheBlochDe[2],
      betheBlochDe[3], betheBlochDe[4]);
    return static_cast<float>((track.tpcSignal() - expected) /
                              (expected * betheBlochDe[BetheBlochResolutionIndex]));
  }

  template <typename Track>
  bool passDeuteronPID(Track const& track) const
  {
    const float tpcInnerParam = track.tpcInnerParam();
    if (std::abs(tpcInnerParam) < deuteronCuts.tpcInnerParamMin.value) {
      return false;
    }
    const float tpcNsigma = deuteronTPCNsigma(track);
    if (track.hasTOF() &&
        tpcInnerParam > deuteronCuts.tofMomentumMin.value) {
      const float tofNsigma = track.tofNSigmaDe();
      if (std::hypot(tpcNsigma, tofNsigma) > deuteronCuts.nsigmaMax.value) {
        return false;
      }
      return !deuteronCuts.requireIndividualNsigma.value ||
             (std::abs(tpcNsigma) <= deuteronCuts.nsigmaMax.value &&
              std::abs(tofNsigma) <= deuteronCuts.nsigmaMax.value);
    }
    if (tpcInnerParam <= deuteronCuts.tofMomentumMin.value) {
      if (std::abs(tpcNsigma) > deuteronCuts.nsigmaMax.value) {
        return false;
      }
      o2::aod::ITSResponse itsResponse;
      const float itsNsigma =
        itsResponse.nSigmaITS<o2::track::PID::Deuteron>(
          track.itsClusterSizes(), track.p(), track.eta());
      return std::abs(itsNsigma) <= deuteronCuts.nsigmaMax.value;
    }
    return false;
  }

  template <typename Tracks>
  void fillSingleParticleQA(Tracks const& tracks)
  {
    for (auto const& track : tracks) {
      registry.fill(HIST("Pion/hSelection"), static_cast<float>(kAllTracks));
      if (passPionTrackCuts(track)) {
        registry.fill(HIST("Pion/hSelection"), static_cast<float>(kTrackCuts));
        if (passPionPID(track)) {
          registry.fill(HIST("Pion/hSelection"), static_cast<float>(kPidCuts));
          const float signedPt = track.sign() * track.pt();
          const float tpcNsigma = track.tpcNSigmaPi();
          registry.fill(HIST("Pion/hPt"), signedPt);
          registry.fill(HIST("Pion/hEta"), track.eta());
          registry.fill(HIST("Pion/hTPCNsigma"), signedPt, tpcNsigma);
          o2::aod::ITSResponse itsResponse;
          const float itsNsigma = itsResponse.nSigmaITS<o2::track::PID::Pion>(
            track.itsClusterSizes(), track.p(), track.eta());
          registry.fill(HIST("Pion/hITSNsigma"), signedPt, itsNsigma);
          if (track.hasTOF() &&
              std::abs(track.p()) > pionCuts.combinedPidMomentumMin.value) {
            const float tofNsigma = track.tofNSigmaPi();
            registry.fill(HIST("Pion/hTOFNsigma"), signedPt, tofNsigma);
            registry.fill(HIST("Pion/hTPCTOFNsigma"), signedPt,
                          std::hypot(tpcNsigma, tofNsigma));
          }
        }
      }

      registry.fill(HIST("Deuteron/hSelection"),
                    static_cast<float>(kAllTracks));
      if (passDeuteronTrackCuts(track)) {
        registry.fill(HIST("Deuteron/hSelection"),
                      static_cast<float>(kTrackCuts));
        if (passDeuteronPID(track)) {
          registry.fill(HIST("Deuteron/hSelection"),
                        static_cast<float>(kPidCuts));
          const float signedPt = track.sign() * track.pt();
          const float tpcNsigma = deuteronTPCNsigma(track);
          registry.fill(HIST("Deuteron/hPt"), signedPt);
          registry.fill(HIST("Deuteron/hEta"), track.eta());
          registry.fill(HIST("Deuteron/hTPCNsigma"), signedPt, tpcNsigma);
          if (track.tpcInnerParam() <= deuteronCuts.tofMomentumMin.value) {
            o2::aod::ITSResponse itsResponse;
            const float itsNsigma =
              itsResponse.nSigmaITS<o2::track::PID::Deuteron>(
                track.itsClusterSizes(), track.p(), track.eta());
            registry.fill(HIST("Deuteron/hITSNsigma"), signedPt, itsNsigma);
          } else if (track.hasTOF()) {
            const float tofNsigma = track.tofNSigmaDe();
            registry.fill(HIST("Deuteron/hTOFNsigma"), signedPt, tofNsigma);
            registry.fill(HIST("Deuteron/hTPCTOFNsigma"), signedPt,
                          std::hypot(tpcNsigma, tofNsigma));
          }
          registry.fill(HIST("Deuteron/hTPCClusters"),
                        static_cast<float>(track.tpcNClsFound()));
          registry.fill(HIST("Deuteron/hTPCCrossedRows"),
                        static_cast<float>(track.tpcNClsCrossedRows()));
        }
      }
    }
  }

  template <typename Track>
  float phiAtRadius(Track const& track, float radius) const
  {
    const float absPt = std::abs(track.pt());
    if (absPt <= 0.f) {
      return InvalidPhiStar;
    }
    const float argument =
      0.3f * static_cast<float>(track.sign()) * 0.1f * magneticField *
      radius * 0.01f / (2.f * absPt);
    if (std::abs(argument) >= 1.f) {
      return InvalidPhiStar;
    }
    return track.phi() - std::asin(argument);
  }

  static float wrapDeltaPhi(float deltaPhi)
  {
    return std::atan2(std::sin(deltaPhi), std::cos(deltaPhi));
  }

  template <typename FirstTrack, typename SecondTrack>
  float averagePhiStar(FirstTrack const& first,
                       SecondTrack const& second) const
  {
    float sum = 0.f;
    int entries = 0;
    for (const auto& radius : tpcRadii) {
      const float firstPhi = phiAtRadius(first, radius);
      const float secondPhi = phiAtRadius(second, radius);
      if (firstPhi == InvalidPhiStar || secondPhi == InvalidPhiStar) {
        continue;
      }
      sum += wrapDeltaPhi(firstPhi - secondPhi);
      ++entries;
    }
    return entries > 0 ? sum / static_cast<float>(entries) : InvalidPhiStar;
  }

  template <typename FirstTrack, typename SecondTrack>
  bool isClosePair(FirstTrack const& first, SecondTrack const& second,
                   bool fillQA)
  {
    if (!pairCuts.enableClosePairRejection.value ||
        first.sign() != second.sign()) {
      return false;
    }
    const float deltaEta = first.eta() - second.eta();
    float deltaPhi = averagePhiStar(first, second);
    if (pairCuts.closePairRadiusMode.value == ClosePairRadiusModeAtPV) {
      deltaPhi = wrapDeltaPhi(first.phi() - second.phi());
    } else if (pairCuts.closePairRadiusMode.value == ClosePairRadiusModeAtSpecificRadius) {
      const float firstPhi =
        phiAtRadius(first, pairCuts.closePairSpecificRadius.value);
      const float secondPhi =
        phiAtRadius(second, pairCuts.closePairSpecificRadius.value);
      if (firstPhi == InvalidPhiStar || secondPhi == InvalidPhiStar) {
        return false;
      }
      deltaPhi = wrapDeltaPhi(firstPhi - secondPhi);
    }
    if (deltaPhi == InvalidPhiStar ||
        pairCuts.closePairDeltaEtaMax.value <= 0.f ||
        pairCuts.closePairDeltaPhiMax.value <= 0.f) {
      return false;
    }
    if (fillQA) {
      registry.fill(HIST("Pairs/hCPRBefore"), deltaEta, deltaPhi);
    }
    const bool rejected =
      deltaEta * deltaEta /
          (pairCuts.closePairDeltaEtaMax.value * pairCuts.closePairDeltaEtaMax.value) +
        deltaPhi * deltaPhi /
          (pairCuts.closePairDeltaPhiMax.value * pairCuts.closePairDeltaPhiMax.value) <
      1.f;
    if (fillQA && !rejected) {
      registry.fill(HIST("Pairs/hCPRAfter"), deltaEta, deltaPhi);
    }
    return rejected;
  }

  template <typename FirstTrack, typename SecondTrack>
  bool hasValidPairFit(FirstTrack const& first, SecondTrack const& second)
  {
    try {
      return pairFitter.process(getTrackParCov(first),
                                getTrackParCov(second)) > 0;
    } catch (...) {
      return false;
    }
  }

  float pairKstar(std::array<float, 3> const& pionMomentum,
                  std::array<float, 3> const& deuteronMomentum) const
  {
    const PairLorentzVector pion(
      pionMomentum[0], pionMomentum[1], pionMomentum[2],
      o2::constants::physics::MassPiPlus);
    const PairLorentzVector deuteron(
      deuteronMomentum[0], deuteronMomentum[1], deuteronMomentum[2],
      o2::constants::physics::MassDeuteron);
    const PairLorentzVector pair = pion + deuteron;
    const float beta = pair.Beta();
    const ROOT::Math::Boost boost(
      -beta * std::cos(pair.Phi()) * std::sin(pair.Theta()),
      -beta * std::sin(pair.Phi()) * std::sin(pair.Theta()),
      -beta * std::cos(pair.Theta()));
    const PairLorentzVector relativeMomentum =
      boost(pion) - boost(deuteron);
    return 0.5f * relativeMomentum.P();
  }

  float pairMT(std::array<float, 3> const& pionMomentum,
               std::array<float, 3> const& deuteronMomentum) const
  {
    const PairLorentzVector pion(
      pionMomentum[0], pionMomentum[1], pionMomentum[2],
      o2::constants::physics::MassPiPlus);
    const PairLorentzVector deuteron(
      deuteronMomentum[0], deuteronMomentum[1], deuteronMomentum[2],
      o2::constants::physics::MassDeuteron);
    const float kT = 0.5f * (pion + deuteron).Pt();
    const float averageMass =
      0.5f * (o2::constants::physics::MassPiPlus +
              o2::constants::physics::MassDeuteron);
    return std::hypot(kT, averageMass);
  }

  static int pairChannel(int deuteronSign, int pionSign)
  {
    const bool unlikeSign = deuteronSign * pionSign < 0;
    const bool antimatter = deuteronSign < 0;
    return 2 * static_cast<int>(unlikeSign) +
           static_cast<int>(antimatter);
  }

  template <typename DeuteronTrack, typename PionTrack>
  void fillPair(DeuteronTrack const& deuteron, PionTrack const& pion,
                float centrality, bool mixedEvent)
  {
    const std::array<float, 3> deuteronMomentum{
      deuteron.px(), deuteron.py(), deuteron.pz()};
    const std::array<float, 3> pionMomentum{
      pion.px(), pion.py(), pion.pz()};
    const float kstar = pairKstar(pionMomentum, deuteronMomentum);
    const float mt = pairMT(pionMomentum, deuteronMomentum);
    const auto channel = static_cast<float>(
      pairChannel(deuteron.sign(), pion.sign()));
    if (mixedEvent) {
      registry.fill(HIST("Pairs/hME"), kstar, mt, centrality, channel);
      if (kstar < analysisSettings.lowKstarYieldMax.value) {
        registry.fill(HIST("Pairs/hMELowKstarYield"), channel);
      }
    } else {
      registry.fill(HIST("Pairs/hSE"), kstar, mt, centrality, channel);
      if (kstar < analysisSettings.lowKstarYieldMax.value) {
        registry.fill(HIST("Pairs/hSELowKstarYield"), channel);
      }
    }
  }

  template <typename DeuteronTracks, typename PionTracks>
  void buildPairs(DeuteronTracks const& deuteronTracks,
                  PionTracks const& pionTracks, float centrality,
                  bool mixedEvent)
  {
    for (auto const& deuteron : deuteronTracks) {
      if (!passDeuteronTrackCuts(deuteron) ||
          !passDeuteronPID(deuteron)) {
        continue;
      }
      for (auto const& pion : pionTracks) {
        if (!mixedEvent && deuteron.globalIndex() == pion.globalIndex()) {
          continue;
        }
        if (!passPionTrackCuts(pion) || !passPionPID(pion)) {
          continue;
        }
        if (isClosePair(deuteron, pion, !mixedEvent)) {
          continue;
        }
        if (!mixedEvent && !hasValidPairFit(deuteron, pion)) {
          continue;
        }
        fillPair(deuteron, pion, centrality, mixedEvent);
      }
    }
  }

  void processSameEvent(CollisionsFull const& collisions,
                        TrackCandidates const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    for (auto const& collision : collisions) {
      if (!isCollisionSelected(collision, true)) {
        continue;
      }
      auto collisionTracks =
        tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      collisionTracks.bindExternalIndices(&tracks);
      fillSingleParticleQA(collisionTracks);
      buildPairs(collisionTracks, collisionTracks,
                 collision.centFT0C(), false);
    }
  }
  PROCESS_SWITCH(PiDeFemtoSystematics, processSameEvent,
                 "Fill same-event pion-deuteron pairs", true);

  void processMixedEvent(CollisionsFull const& /*collisions*/,
                         TrackCandidates const& /*tracks*/,
                         aod::BCsWithTimestamps const&)
  {
    for (auto const& [collision1, tracks1, collision2, tracks2] :
         mixedEventPairs) {
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }
      auto bc1 = collision1.template bc_as<aod::BCsWithTimestamps>();
      auto bc2 = collision2.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc1);
      buildPairs(tracks1, tracks2, collision1.centFT0C(), true);
      initCCDB(bc2);
      buildPairs(tracks2, tracks1, collision2.centFT0C(), true);
    }
  }
  PROCESS_SWITCH(PiDeFemtoSystematics, processMixedEvent,
                 "Fill mixed-event pion-deuteron pairs", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PiDeFemtoSystematics>(cfgc)};
}
