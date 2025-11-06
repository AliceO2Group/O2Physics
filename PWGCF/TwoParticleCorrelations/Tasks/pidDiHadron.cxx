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

/// \file pidDiHadron.cxx
/// \brief di-hadron correlation of PID for O-O, Pb-Pb collisions
/// \author Preet Bhanjan Pati (preet.bhanjan.pati@cern.ch), Zhiyong Lu (zhiyong.lu@cern.ch)
/// \since  July/29/2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
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

#include "TF1.h"
#include "TRandom3.h"
#include <TPDGCode.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

template <typename T, typename P>
auto readMatrix(Array2D<T> const& mat, P& array, int rowStart, int rowEnd, int colStart, int colEnd)
{
  for (auto i = rowStart; i < rowEnd; ++i) {
    for (auto j = colStart; j < colEnd; ++j) {
      array[i][j] = mat(i, j);
    }
  }

  return;
}

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
static constexpr int LongArrayInt[3][20] = {{1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {2, 2, 2, -2, -2, -2, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}, {3, 3, 3, -3, -3, -3, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1}};

struct PidDiHadron {
  o2::aod::ITSResponse itsResponse;
  Service<ccdb::BasicCCDBManager> ccdb;

  enum ResoParticles {
    K0 = 0,
    LAMBDA = 1,
    PHI = 2,
    NResoParticles = 3
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
    kRapidity,
    kNParticleCuts
  };
  enum ParticleSwitches {
    kUseCosPA = 0,
    kMassBins,
    kDCABetDaug,
    kUseProperLifetime,
    kUseV0Radius,
    kNParticleSwitches
  };
  enum Particles {
    PIONS = 0,
    KAONS,
    PROTONS
  };
  enum ParticleNsigma {
    kPionUpCut = 0,
    kKaonUpCut,
    kProtonUpCut,
    kPionLowCut,
    kKaonLowCut,
    kProtonLowCut
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
    kUseTVXinTRD,
    kNEventCuts
  };
  enum {
    kCharged = 0,
    kPions,
    kKaons,
    kProtons,
    kK0,
    kLambda,
    kPhi
  };

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "minimum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "maximum accepted track pT")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta cut")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, float, 50.0f, "minimum TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgTpcCrossRows, float, 70.0f, "minimum TPC crossed rows")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgTpcCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, float, 5.0f, "minimum ITS clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMax, int, 2000, "Minimum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancyMin, int, 0, "Maximum occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgFakeKaonCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgCutMerging, float, 0.0, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgSelCollByNch, bool, true, "Select collisions by Nch or centrality")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMin, int, 0, "Minimum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutMultMax, int, 10, "Maximum multiplicity for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMin, float, 60.0f, "Minimum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgCutCentMax, float, 80.0f, "Maximum centrality for collision")
  O2_DEFINE_CONFIGURABLE(cfgMixEventNumMin, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgCentEstimator, int, 0, "0:FT0C; 1:FT0CVariant1; 2:FT0M; 3:FT0A")
  O2_DEFINE_CONFIGURABLE(cfgCentTableUnavailable, bool, false, "if a dataset does not provide centrality information")
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgV0AT0Acut, int, 5, "V0AT0A cut")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, bool, false, "Verbose output")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrder, bool, false, "enable trigger pT < associated pT cut")
  O2_DEFINE_CONFIGURABLE(cfgUsePtOrderInMixEvent, bool, false, "enable trigger pT < associated pT cut in mixed event")
  O2_DEFINE_CONFIGURABLE(cfgUseItsPID, bool, true, "Use ITS PID for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgUseOnlyTPC, bool, true, "Use only TPC PID for daughter selection")
  O2_DEFINE_CONFIGURABLE(cfgPIDParticle, int, 0, "1 = pion, 2 = kaon, 3 = proton, 4 = kshort, 5 = lambda, 6 = phi, 0 for no PID")

  SliceCache cache;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity axis for histograms"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "centrality axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -2.4, 2.4}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.2, 0.5, 1, 1.5, 2, 3, 4, 6, 10}, "pt associated axis for histograms"};
  ConfigurableAxis axisVtxMix{"axisVtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis axisMultMix{"axisMultMix", {VARIABLE_WIDTH, 0, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};

  Configurable<LabeledArray<int>> cfgUseEventCuts{"cfgUseEventCuts", {LongArrayInt[0], 1, 15, {"EvCuts"}, {"Filtered Events", "Sel8", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kNoCollInTimeRangeStandard", "kIsGoodITSLayersAll", "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "Occupancy", "Multcorrelation", "T0AV0ACut", "kIsVertexITSTPC", "kTVXinTRD"}}, "Labeled array (int) for various cuts on resonances"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 3, 6, {"TPC", "TOF", "ITS"}, {"pos_pi", "pos_ka", "pos_pr", "neg_pi", "neg_ka", "neg_pr"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  Configurable<LabeledArray<float>> resonanceCuts{"resonanceCuts", {LongArrayFloat[0], 3, 11, {"K0", "Lambda", "Phi"}, {"cos_PAs", "massMin", "massMax", "PosTrackPt", "NegTrackPt", "DCAPosToPVMin", "DCANegToPVMin", "Lifetime", "RadiusMin", "RadiusMax", "Rapidity"}}, "Labeled array (float) for various cuts on resonances"};
  Configurable<LabeledArray<int>> resonanceSwitches{"resonanceSwitches", {LongArrayInt[0], 3, 5, {"K0", "Lambda", "Phi"}, {"UseCosPA", "NMassBins", "DCABetDaug", "UseProperLifetime", "UseV0Radius"}}, "Labeled array (int) for various cuts on resonances"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {1, 0, 1}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {1, 0, 1}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {1, 0, 1}, "pt axis for efficiency histograms"};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using FilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>>;
  using FilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using V0TrackCandidate = aod::V0Datas;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::V0Datas> perCollisionV0 = aod::v0::collisionId;

  // Corrections
  TH3D* mEfficiency = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  OutputObj<CorrelationContainer> sameReso{"sameEventReso"};
  OutputObj<CorrelationContainer> mixedReso{"mixedEventReso"};
  HistogramRegistry histos{"histos"};

  // define global variables
  TRandom3* gRandom = new TRandom3();
  enum CentEstimators {
    kCentFT0C = 0,
    kCentFT0CVariant1,
    kCentFT0M,
    kCentFV0A,
    // Count the total number of enum
    kCount_CentEstimators
  };
  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };
  enum EventCutType {
    kEvCut1 = 0,
    kNEvCutTypes = 1
  };

  std::array<std::array<int, 15>, 1> eventCuts;
  std::array<std::array<float, 11>, 3> resoCutVals;
  std::array<std::array<int, 6>, 3> resoSwitchVals;
  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;

  // persistent caches
  std::vector<float> efficiencyAssociatedCache;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  void init(InitContext&)
  {
    // projectMatrix(nSigmas->getData(), tpcNsigmaCut, tofNsigmaCut, itsNsigmaCut);
    readMatrix(resonanceCuts->getData(), resoCutVals, K0, NResoParticles, kCosPA, kNParticleCuts);
    readMatrix(resonanceSwitches->getData(), resoSwitchVals, K0, NResoParticles, kUseCosPA, kNParticleSwitches);
    readMatrix(cfgUseEventCuts->getData(), eventCuts, kEvCut1, kNEvCutTypes, kFilteredEvents, kNEventCuts);

    tpcNsigmaCut[kPionUpCut] = nSigmas->getData()[kTPC][kPionUpCut];
    tpcNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTPC][kKaonUpCut];
    tpcNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTPC][kProtonUpCut];
    tpcNsigmaCut[kPionLowCut] = nSigmas->getData()[kTPC][kPionLowCut];
    tpcNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTPC][kKaonLowCut];
    tpcNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTPC][kProtonLowCut];

    tofNsigmaCut[kPionUpCut] = nSigmas->getData()[kTOF][kPionUpCut];
    tofNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTOF][kKaonUpCut];
    tofNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTOF][kProtonUpCut];
    tofNsigmaCut[kPionLowCut] = nSigmas->getData()[kTOF][kPionLowCut];
    tofNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTOF][kKaonLowCut];
    tofNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTOF][kProtonLowCut];

    itsNsigmaCut[kPionUpCut] = nSigmas->getData()[kITS][kPionUpCut];
    itsNsigmaCut[kKaonUpCut] = nSigmas->getData()[kITS][kKaonUpCut];
    itsNsigmaCut[kProtonUpCut] = nSigmas->getData()[kITS][kProtonUpCut];
    itsNsigmaCut[kPionLowCut] = nSigmas->getData()[kITS][kPionLowCut];
    itsNsigmaCut[kKaonLowCut] = nSigmas->getData()[kITS][kKaonLowCut];
    itsNsigmaCut[kProtonLowCut] = nSigmas->getData()[kITS][kProtonLowCut];

    AxisSpec axisK0Mass = {resoSwitchVals[K0][kMassBins], resoCutVals[K0][kMassMin], resoCutVals[K0][kMassMax]};
    AxisSpec axisLambdaMass = {resoSwitchVals[LAMBDA][kMassBins], resoCutVals[LAMBDA][kMassMin], resoCutVals[LAMBDA][kMassMax]};
    AxisSpec axisPhiMass = {resoSwitchVals[PHI][kMassBins], resoCutVals[PHI][kMassMin], resoCutVals[PHI][kMassMax]};

    if (cfgCentTableUnavailable && !cfgSelCollByNch) {
      LOGF(fatal, "Centrality table is unavailable, cannot select collisions by centrality");
    }
    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    LOGF(info, "Starting init");

    // Creating mass axis depending on particle - 4 = kshort, 5 = lambda, 6 = phi
    AxisSpec massAxisReso = {10, 0, 1, "mass"};
    if (cfgPIDParticle == kK0)
      massAxisReso = {resoSwitchVals[K0][kMassBins], resoCutVals[K0][kMassMin], resoCutVals[K0][kMassMax], "M_{#pi^{+}#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDParticle == kLambda)
      massAxisReso = {resoSwitchVals[LAMBDA][kMassBins], resoCutVals[LAMBDA][kMassMin], resoCutVals[LAMBDA][kMassMax], "M_{p#pi^{-}} (GeV/c^{2})"};
    if (cfgPIDParticle == kPhi)
      massAxisReso = {resoSwitchVals[PHI][kMassBins], resoCutVals[PHI][kMassMin], resoCutVals[PHI][kMassMax], "M_{K^{+}K^{-}} (GeV/c^{2})"};

    // Event Counter
    if ((doprocessSame || doprocessSameReso) && cfgUseAdditionalEventCut) {
      histos.add("hEventCount", "Number of Events;; Count", {HistType::kTH1D, {{15, -0.5, 14.5}}});
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kFilteredEvents + 1, "Filtered event");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kAfterSel8 + 1, "After sel8");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoTimeFrameBorder + 1, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoITSROFrameBorder + 1, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoSameBunchPileup + 1, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseGoodITSLayersAll + 1, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoCollInRofStandard + 1, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseNoHighMultCollInPrevRof + 1, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseOccupancy + 1, "Occupancy Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseMultCorrCut + 1, "Multiplicity correlation Cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseT0AV0ACut + 1, "T0AV0A cut");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseVertexITSTPC + 1, "kIsVertexITSTPC");
      histos.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(kUseTVXinTRD + 1, "kTVXinTRD");
    }

    if (cfgPIDParticle == kK0) { // For K0
      histos.add("PiPlusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiPlusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_K0", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hK0Phi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hK0Eta", "", {HistType::kTH1D, {axisEta}});

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
    if (cfgPIDParticle == kLambda) { // For Lambda
      histos.add("PrPlusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PiMinusTPC_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTPC}}});
      histos.add("PrPlusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("PiMinusTOF_L", "", {HistType::kTH2D, {{axisPt, axisNsigmaTOF}}});
      histos.add("hLambdaPhi", "", {HistType::kTH1D, {axisPhi}});
      histos.add("hLambdaEta", "", {HistType::kTH1D, {axisEta}});

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
    if (cfgPIDParticle == kPhi) { // For Phi
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

    // Multiplicity correlation cuts
    if (eventCuts[kEvCut1][kUseMultCorrCut]) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }
    if (eventCuts[kEvCut1][kUseT0AV0ACut]) {
      fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    std::string hCentTitle = "Centrality distribution, Estimator " + std::to_string(cfgCentEstimator);
    // Make histograms to check the distributions after cuts
    if (doprocessSame || doprocessSameReso) {
      histos.add("deltaEta_deltaPhi_same", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      histos.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      histos.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      histos.add("Centrality", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}});
      histos.add("Centrality_used", hCentTitle.c_str(), {HistType::kTH1D, {axisCentrality}}); // histogram to see how many events are in the same and mixed event
      histos.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      histos.add("zVtx_used", "zVtx_used", {HistType::kTH1D, {axisVertex}});

      if (cfgPIDParticle == kCharged || cfgPIDParticle == kPions || cfgPIDParticle == kKaons || cfgPIDParticle == kProtons) {
        histos.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
        histos.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
        histos.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
        histos.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
        histos.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
        histos.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }

      if (cfgPIDParticle == kK0 || cfgPIDParticle == kLambda || cfgPIDParticle == kPhi) {
        histos.add("Trig_histReso", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger, massAxisReso}}});
      }
    }
    if (doprocessMixed || doprocessMixedReso) {
      histos.add("deltaEta_deltaPhi_mixed", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
    }

    histos.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");
    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};

    std::vector<AxisSpec> corrAxisReso = {{axisSample, "Sample"},
                                          {axisVertex, "z-vtx (cm)"},
                                          {axisPtTrigger, "p_{T} (GeV/c)"},
                                          {massAxisReso},
                                          {axisDeltaPhi, "#Delta#varphi (rad)"},
                                          {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };
    std::vector<AxisSpec> userAxis;

    if (cfgPIDParticle == kCharged || cfgPIDParticle == kPions || cfgPIDParticle == kKaons || cfgPIDParticle == kProtons) {
      same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));
    }

    if (cfgPIDParticle == kK0 || cfgPIDParticle == kLambda || cfgPIDParticle == kPhi) {
      sameReso.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxisReso, effAxis, userAxis));
      mixedReso.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxisReso, effAxis, userAxis));
    }
    LOGF(info, "End of init");
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision>
  float getCentrality(TCollision const& collision)
  {
    float cent;
    switch (cfgCentEstimator) {
      case kCentFT0C:
        cent = collision.centFT0C();
        break;
      case kCentFT0CVariant1:
        cent = collision.centFT0CVariant1();
        break;
      case kCentFT0M:
        cent = collision.centFT0M();
        break;
      case kCentFV0A:
        cent = collision.centFV0A();
        break;
      default:
        cent = collision.centFT0C();
    }
    return cent;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgTpcCluster) && (track.tpcNClsCrossedRows() >= cfgTpcCrossRows) && (track.itsNCls() >= cfgITScluster));
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
      int partIndex = getNsigmaPID(track);
      int pidIndex = partIndex - 1; // 0 = pion, 1 = kaon, 2 = proton
      if (pidIndex != pid)
        return false;
    }

    return true;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = 0; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[PIONS] < detectorNsigmaCut[kPionUpCut] && nSigmaToUse[PIONS] > detectorNsigmaCut[kPionLowCut];
    bool isDetectedKaon = nSigmaToUse[KAONS] < detectorNsigmaCut[kKaonUpCut] && nSigmaToUse[KAONS] > detectorNsigmaCut[kKaonLowCut];
    bool isDetectedProton = nSigmaToUse[PROTONS] < detectorNsigmaCut[kProtonUpCut] && nSigmaToUse[PROTONS] > detectorNsigmaCut[kProtonLowCut];

    bool isTofPion = nSigmaTOF[PIONS] < tofNsigmaCut[kPionUpCut] && nSigmaTOF[PIONS] > tofNsigmaCut[kPionLowCut];
    bool isTofKaon = nSigmaTOF[KAONS] < tofNsigmaCut[kKaonUpCut] && nSigmaTOF[KAONS] > tofNsigmaCut[kKaonLowCut];
    bool isTofProton = nSigmaTOF[PROTONS] < tofNsigmaCut[kProtonUpCut] && nSigmaTOF[PROTONS] > tofNsigmaCut[kProtonLowCut];

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
      pid = PIONS + 1;
    } else if (isKaon) {
      pid = KAONS + 1;
    } else if (isProton) {
      pid = PROTONS + 1;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = reinterpret_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    float weff1 = 1;
    float vtxz = collision.posZ();
    for (auto const& track1 : tracks) {
      if (!trackSelected(track1))
        continue;
      if (cfgPIDParticle && getNsigmaPID(track1) != cfgPIDParticle)
        continue; // if PID is selected, check if the track has the right PID
      if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), vtxz))
        continue;
      histos.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      histos.fill(HIST("Eta"), track1.eta());
      histos.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      histos.fill(HIST("pT"), track1.pt());
      histos.fill(HIST("pTCorrected"), track1.pt(), weff1);
    }
  }

  float getDPhiStar(float charge1, float charge2, float phi1, float phi2, float pt1, float pt2, float radius, int magField)
  {
    int fbSign = (magField > 0) ? 1 : -1;

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelations(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // Cache efficiency for particles (too many FindBin lookups)
    if (mEfficiency) {
      efficiencyAssociatedCache.clear();
      efficiencyAssociatedCache.reserve(tracks2.size());
      for (const auto& track2 : tracks2) {
        float weff = 1.;
        getEfficiencyCorrection(weff, track2.eta(), track2.pt(), posZ);
        efficiencyAssociatedCache.push_back(weff);
      }
    }

    if (system == SameEvent) {
      if (!cfgCentTableUnavailable)
        histos.fill(HIST("Centrality_used"), cent);
      histos.fill(HIST("Nch_used"), tracks1.size());
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;
      if (cfgPIDParticle && getNsigmaPID(track1) != cfgPIDParticle)
        continue; // if PID is selected, check if the track has the right PID
      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;
      if (system == SameEvent) {
        histos.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;
        if (mEfficiency) {
          associatedWeight = efficiencyAssociatedCache[track2.filteredIndex()];
        }

        if (!cfgUsePtOrder && track1.globalIndex() == track2.globalIndex())
          continue; // For pt-differential correlations, skip if the trigger and associate are the same track
        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1.sign(), track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {

          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {

          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TV0Tracks, typename TTracksAssoc>
  void fillCorrelationsReso(TV0Tracks tracks1, TTracksAssoc tracks2, float posZ, float posY, float posX, int system, int magneticField, float cent, float eventWeight) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {
    // Cache efficiency for particles (too many FindBin lookups)
    if (mEfficiency) {
      efficiencyAssociatedCache.clear();
      efficiencyAssociatedCache.reserve(tracks2.size());
      for (const auto& track2 : tracks2) {
        float weff = 1.;
        getEfficiencyCorrection(weff, track2.eta(), track2.pt(), posZ);
        efficiencyAssociatedCache.push_back(weff);
      }
    }

    if (system == SameEvent) {
      if (!cfgCentTableUnavailable)
        histos.fill(HIST("Centrality_used"), cent);
      histos.fill(HIST("Nch_used"), tracks2.size()); // Taking Nch from tracks2 since tracks1 are V0s
    }

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;
    float associatedWeight = 1.0f;
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      double resoMass = -1;

      // 4 = kshort, 5 = lambda, 6 = phi
      if (cfgPIDParticle == kK0) {
        if (!selectionK0(track1, posZ, posY, posX))
          continue; // Reject if called for K0 but V0 is not K0

        resoMass = track1.mK0Short();
      }

      if (cfgPIDParticle == kLambda) {
        if (!selectionLambda(track1, posZ, posY, posX))
          continue; // Reject if called for Lambda but V0 is not lambda

        resoMass = track1.mLambda();
      }

      if (system == SameEvent) {
        histos.fill(HIST("Trig_histReso"), fSampleIndex, posZ, track1.pt(), resoMass, eventWeight * triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (!trackSelected(track2))
          continue;
        if (track2.pt() < cfgCutPtMin || track2.pt() > cfgCutPtMax) // Select associated particles in the pt range 0.2 - 3.0 GeV/c
          continue;
        if (mEfficiency) {
          associatedWeight = efficiencyAssociatedCache[track2.filteredIndex()];
        }

        if (cfgUsePtOrder && system == SameEvent && track1.pt() <= track2.pt())
          continue; // Without pt-differential correlations, skip if the trigger pt is less than the associate pt
        if (cfgUsePtOrder && system == MixedEvent && cfgUsePtOrderInMixEvent && track1.pt() <= track2.pt())
          continue; // For pt-differential correlations in mixed events, skip if the trigger pt is less than the associate pt

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (std::abs(deltaEta) < cfgCutMerging) {

          double dPhiStarHigh = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgCutMerging;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(0, track2.sign(), track1.phi(), track2.phi(), track1.pt(), track2.pt(), rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {
          sameReso->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        } else if (system == MixedEvent) {
          mixedReso->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), resoMass, deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * triggerWeight * associatedWeight);
        }
      }
    }
  }

  template <typename TCollision>
  bool selectionEvent(TCollision collision, const int mult, const float cent, const bool fillCounter)
  {
    if (fillCounter)
      histos.fill(HIST("hEventCount"), kFilteredEvents);
    if (!collision.sel8()) {
      return 0;
    }
    if (fillCounter)
      histos.fill(HIST("hEventCount"), kAfterSel8);

    if (eventCuts[kEvCut1][kUseNoTimeFrameBorder] && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoTimeFrameBorder])
      histos.fill(HIST("hEventCount"), kUseNoTimeFrameBorder);

    if (eventCuts[kEvCut1][kUseNoITSROFrameBorder] && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoITSROFrameBorder])
      histos.fill(HIST("hEventCount"), kUseNoITSROFrameBorder);

    if (eventCuts[kEvCut1][kUseNoSameBunchPileup] && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoSameBunchPileup])
      histos.fill(HIST("hEventCount"), kUseNoSameBunchPileup);

    if (eventCuts[kEvCut1][kUseGoodZvtxFT0vsPV] && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseGoodZvtxFT0vsPV])
      histos.fill(HIST("hEventCount"), kUseGoodZvtxFT0vsPV);

    if (eventCuts[kEvCut1][kUseNoCollInTimeRangeStandard] && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoCollInTimeRangeStandard])
      histos.fill(HIST("hEventCount"), kUseNoCollInTimeRangeStandard);

    if (eventCuts[kEvCut1][kUseGoodITSLayersAll] && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseGoodITSLayersAll])
      histos.fill(HIST("hEventCount"), kUseGoodITSLayersAll);

    if (eventCuts[kEvCut1][kUseNoCollInRofStandard] && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoCollInRofStandard])
      histos.fill(HIST("hEventCount"), kUseNoCollInRofStandard);

    if (eventCuts[kEvCut1][kUseNoHighMultCollInPrevRof] && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseNoHighMultCollInPrevRof])
      histos.fill(HIST("hEventCount"), kUseNoHighMultCollInPrevRof);

    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();

    if (eventCuts[kEvCut1][kUseOccupancy] && (occupancy < cfgCutOccupancyMin || occupancy > cfgCutOccupancyMax)) {
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseOccupancy])
      histos.fill(HIST("hEventCount"), kUseOccupancy);

    if (eventCuts[kEvCut1][kUseMultCorrCut]) {
      if (multNTracksPV < fMultPVCutLow->Eval(cent))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(cent))
        return 0;
      if (mult < fMultCutLow->Eval(cent))
        return 0;
      if (mult > fMultCutHigh->Eval(cent))
        return 0;
    }

    if (fillCounter && eventCuts[kEvCut1][kUseMultCorrCut])
      histos.fill(HIST("hEventCount"), kUseMultCorrCut);

    // V0A T0A 5 sigma cut
    if (eventCuts[kEvCut1][kUseT0AV0ACut] && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > cfgV0AT0Acut * fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && eventCuts[kEvCut1][kUseT0AV0ACut])
      histos.fill(HIST("hEventCount"), kUseT0AV0ACut);

    if (eventCuts[kEvCut1][kUseVertexITSTPC] && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return 0;
    if (fillCounter && eventCuts[kEvCut1][kUseVertexITSTPC])
      histos.fill(HIST("hEventCount"), kUseVertexITSTPC);

    if (eventCuts[kEvCut1][kUseTVXinTRD] && collision.alias_bit(kTVXinTRD)) {
      return 0;
    }
    if (fillCounter && eventCuts[kEvCut1][kUseTVXinTRD])
      histos.fill(HIST("hEventCount"), kUseTVXinTRD);

    return 1;
  }

  void processSame(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    if (!cfgCentTableUnavailable)
      cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true))
      return;

    if (!cfgCentTableUnavailable)
      histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    loadEfficiency(bc.timestamp());
    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);

    same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()), cent, 1.0f);
  }
  PROCESS_SWITCH(PidDiHadron, processSame, "Process same event", true);

  // the process for filling the mixed events
  void processMixed(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<FilteredCollisions, FilteredTracks, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadEfficiency(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()), cent1, eventWeight);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixed, "Process mixed events", true);

  template <typename V0>
  bool selectionK0(V0 const& candidate, float posZ, float posY, float posX)
  {
    double mk0 = candidate.mK0Short();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

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
    if (resoSwitchVals[K0][kUseProperLifetime] && candidate.distovertotmom(posX, posY, posZ) * massK0Short > resoCutVals[K0][kLifeTime])
      return false;
    histos.fill(HIST("hK0Count"), 8.5);
    if (!selectionV0Daughter(postrack, PIONS) || !selectionV0Daughter(negtrack, PIONS))
      return false;
    histos.fill(HIST("hK0Count"), 9.5);

    histos.fill(HIST("hK0Phi"), candidate.phi());
    histos.fill(HIST("hK0Eta"), candidate.eta());
    histos.fill(HIST("PiPlusTPC_K0"), postrack.pt(), postrack.tpcNSigmaPi());
    histos.fill(HIST("PiPlusTOF_K0"), postrack.pt(), postrack.tofNSigmaPi());
    histos.fill(HIST("PiMinusTPC_K0"), negtrack.pt(), negtrack.tpcNSigmaPi());
    histos.fill(HIST("PiMinusTOF_K0"), negtrack.pt(), negtrack.tofNSigmaPi());

    return true;
  }

  template <typename V0>
  bool selectionLambda(V0 const& candidate, float posZ, float posY, float posX)
  {
    bool isL = false;  // Is lambda candidate
    bool isAL = false; // Is anti-lambda candidate

    double mlambda = candidate.mLambda();
    double mantilambda = candidate.mAntiLambda();

    // separate the positive and negative V0 daughters
    auto postrack = candidate.template posTrack_as<FilteredTracks>();
    auto negtrack = candidate.template negTrack_as<FilteredTracks>();

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
    if (resoSwitchVals[LAMBDA][kUseProperLifetime] && candidate.distovertotmom(posX, posY, posZ) * massLambda > resoCutVals[LAMBDA][kLifeTime])
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

    if (isAL) { // Reject the track if it is antilambda
      return false;
    }

    if (isL) {
      histos.fill(HIST("hLambdaPhi"), candidate.phi());
      histos.fill(HIST("hLambdaEta"), candidate.eta());
      histos.fill(HIST("PrPlusTPC_L"), postrack.pt(), postrack.tpcNSigmaPr());
      histos.fill(HIST("PrPlusTOF_L"), postrack.pt(), postrack.tofNSigmaPr());
      histos.fill(HIST("PiMinusTPC_L"), negtrack.pt(), negtrack.tpcNSigmaPi());
      histos.fill(HIST("PiMinusTOF_L"), negtrack.pt(), negtrack.tofNSigmaPi());
    }

    return true;
  }

  double massKaPlus = o2::constants::physics::MassKPlus;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0Short = o2::constants::physics::MassK0Short;

  void processSameReso(FilteredCollisions::iterator const& collision, FilteredTracks const& tracks, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float cent = -1.;
    if (!cfgCentTableUnavailable)
      cent = getCentrality(collision);
    if (cfgUseAdditionalEventCut && !selectionEvent(collision, tracks.size(), cent, true))
      return;

    if (!cfgCentTableUnavailable)
      histos.fill(HIST("Centrality"), cent);
    histos.fill(HIST("Nch"), tracks.size());
    histos.fill(HIST("zVtx"), collision.posZ());

    if (cfgSelCollByNch && (tracks.size() < cfgCutMultMin || tracks.size() >= cfgCutMultMax)) {
      return;
    }
    if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent < cfgCutCentMin || cent >= cfgCutCentMax)) {
      return;
    }

    loadEfficiency(bc.timestamp());
    histos.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    sameReso->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
    fillCorrelationsReso<CorrelationContainer::kCFStepReconstructed>(V0s, tracks, collision.posZ(), collision.posY(), collision.posX(), SameEvent, getMagneticField(bc.timestamp()), cent, 1.0f);
  }
  PROCESS_SWITCH(PidDiHadron, processSameReso, "Process same event for resonances", true);

  void processMixedReso(FilteredCollisions const& collisions, FilteredTracks const& tracks, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s)
  {
    auto getTracksSize = [&tracks, this](FilteredCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxMix, axisMultMix}, true};

    auto tracksTuple = std::make_tuple(V0s, tracks);
    Pair<FilteredCollisions, aod::V0Datas, FilteredTracks, MixedBinning> pairs{binningOnVtxAndMult, cfgMixEventNumMin, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, v0s1, collision2, tracks2] = *it;
      if (!collision1.sel8() || !collision2.sel8())
        continue;

      auto tracks1 = tracks.sliceByCached(o2::aod::track::collisionId, collision1.globalIndex(), this->cache);
      if (cfgSelCollByNch && (tracks1.size() < cfgCutMultMin || tracks1.size() >= cfgCutMultMax))
        continue;

      if (cfgSelCollByNch && (tracks2.size() < cfgCutMultMin || tracks2.size() >= cfgCutMultMax))
        continue;

      float cent1 = -1;
      float cent2 = -1;
      if (!cfgCentTableUnavailable) {
        cent1 = getCentrality(collision1);
        cent2 = getCentrality(collision2);
      }
      if (cfgUseAdditionalEventCut && !selectionEvent(collision1, tracks1.size(), cent1, false))
        continue;
      if (cfgUseAdditionalEventCut && !selectionEvent(collision2, tracks2.size(), cent2, false))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent1 < cfgCutCentMin || cent1 >= cfgCutCentMax))
        continue;

      if (!cfgSelCollByNch && !cfgCentTableUnavailable && (cent2 < cfgCutCentMin || cent2 >= cfgCutCentMax))
        continue;

      histos.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadEfficiency(bc.timestamp());
      float eventWeight = 1.0f;
      if (cfgUseEventWeights) {
        eventWeight = 1.0f / it.currentWindowNeighbours();
      }

      fillCorrelationsReso<CorrelationContainer::kCFStepReconstructed>(v0s1, tracks2, collision1.posZ(), collision1.posY(), collision1.posX(), MixedEvent, getMagneticField(bc.timestamp()), cent1, eventWeight);
    }
  }
  PROCESS_SWITCH(PidDiHadron, processMixedReso, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidDiHadron>(cfgc),
  };
}
