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

/// \file dndetaHi.cxx
/// \brief Task for dN/deta analysis in heavy ion collisions
/// \author hyungjun lee <leehy@cern.ch>, joonsuk bae <jbae@cern.ch>

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGMM/Mult/DataModel/Index.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod;

using BCsRun3Cent = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse, aod::BCCentFT0Ms>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults>;
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using DaughterTrack = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, DaughterTrack>;
using FiTracks = soa::Filtered<ExTracks>;
using ExDataTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiDataTracks = soa::Filtered<ExDataTracks>;
using Particles = soa::Join<aod::McParticles, aod::ParticlesToTracks>;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiLTracks = soa::Filtered<LabeledTracksEx>;
using DaughterTracks = soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using BoolVector = std::vector<bool>;
using MyCollisionsMcLabels = soa::Join<MyCollisions, aod::McCollisionLabels>;

enum EventCategory {
  EvtBegin = 0,
  EvtData,
  EvtInel,
  EvtInelg0,
  EvtDd,
  EvtSd,
  EvtNd,
  // Append only: preserve existing category indices for nominal results.
  // Combined with TrigSel8g0 these give MB>0 process-resolved corrections
  // with the same generated INEL>0 target as the nominal correction.
  EvtInelg0Dd,
  EvtInelg0Sd,
  EvtInelg0Nd,
  EvtEnd
};

enum TriggerClass {
  TrigBegin = 0,
  TrigSel8,
  TrigSel8g0,
  TrigEnd
};

enum Species {
  SpecBegin = 0,
  SpecK0short,
  SpecLambda,
  SpecAntilambda,
  SpecEnd
};

enum TrackSign {
  SignBegin = 0,
  SignPositive,
  SignNegative,
  SignEnd
};

enum AnalysisStep {
  StepBegin = 0,
  StepAll,
  StepBasiccut,
  StepMasscut,
  StepEnd
};

enum ParticleType {
  PartBegin = 0,
  PartData,
  PartMotherStrange,
  PartBkg,
  PartNotPrimary,
  PartPion,
  PartKaon,
  PartProton,
  PartOther,
  PartEnd
};

enum PtVariation {
  PtVarBegin = 0,
  PtVarNone,
  PtVarUp,
  PtVarDown,
  PtVarEnd
};

namespace
{

constexpr float CutZ = 10.0f;
constexpr float CutEta1 = 1.0f;
constexpr float CutEta08 = 0.8f;
constexpr float CutEta05 = 0.5f;
constexpr float CutEta09 = 0.9f;
constexpr int CutCharge = 3;
constexpr float CutPt = 0.1f;

constexpr int ProcIdND = 101;
constexpr int ProcIdSD1 = 103;
constexpr int ProcIdSD2 = 104;
constexpr int ProcIdDD1 = 105;
constexpr int ProcIdDD2 = 106;

// Linear pT-extrapolation weights applied when pt < CutPt: weight = slope * pt + intercept
constexpr float PtWeightSlopeUp = -10.0f;
constexpr float PtWeightInterceptUp = 2.0f;
constexpr float PtWeightSlopeDown = 5.0f;
constexpr float PtWeightInterceptDown = 0.5f;
// Sentinel value used as centrality when no centrality information is available
constexpr float NoCentrality = -1.0f;

// Source used to assign the detector-level centrality of an MC collision.
// The BC fallback is deliberately restricted to generated collisions with no
// reconstructed collision.  All other nominal entries retain the official
// collision-level (vertex-z equalized) FT0M percentile.
enum CentralitySource {
  CentSourceRecoFT0M = 1,
  CentSourceZeroRecoBCFT0M,
  CentSourceZeroRecoNoFT0,
  CentSourceZeroRecoInvalidBCFT0M,
  CentSourceZeroRecoNoBC,
  CentSourceInvalidRecoFT0M
};

enum ZeroRecoCentralityClass {
  ZeroRecoCent0To1 = 1,
  ZeroRecoCent1To5,
  ZeroRecoCent5To10,
  ZeroRecoCent10To15,
  ZeroRecoCent15To20,
  ZeroRecoCent20To30,
  ZeroRecoCent30To40,
  ZeroRecoCent40To50,
  ZeroRecoCent50To70,
  ZeroRecoCent70To100,
  ZeroRecoCentNoFT0,
  ZeroRecoCentInvalid,
  ZeroRecoCentNoBC
};

ZeroRecoCentralityClass getZeroRecoCentralityClass(float cent)
{
  if (cent < 1.f) {
    return ZeroRecoCent0To1;
  }
  if (cent < 5.f) {
    return ZeroRecoCent1To5;
  }
  if (cent < 10.f) {
    return ZeroRecoCent5To10;
  }
  if (cent < 15.f) {
    return ZeroRecoCent10To15;
  }
  if (cent < 20.f) {
    return ZeroRecoCent15To20;
  }
  if (cent < 30.f) {
    return ZeroRecoCent20To30;
  }
  if (cent < 40.f) {
    return ZeroRecoCent30To40;
  }
  if (cent < 50.f) {
    return ZeroRecoCent40To50;
  }
  if (cent < 70.f) {
    return ZeroRecoCent50To70;
  }
  return ZeroRecoCent70To100;
}

// Bin indices of the "Selection" counter histogram — values must match SetBinLabel calls in init()
enum SelectionStep {
  SelBegin = 0,
  SelAll,
  SelSel8,
  SelSel8z10,
  SelGeneratedInelGt0,
  SelGoodBCs,
  SelBCsWithCollisions,
  SelBCsWithPileup,
  SelInelSel8Mcz10,
  SelInel,
  SelInelg0,
  SelInelGenz10,
  SelInelg0Genz10,
  SelInelSel8,
  SelInelg0Sel8g0,
  SelInelSel8Recz10,
  SelInelg0Sel8g0Recz10,
  SelK0sz10Eta05
};

static constexpr TrackSelectionFlags::flagtype TrackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype TrackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype TrackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
} // namespace

struct DndetaHi {
  SliceCache cache;
  Preslice<Particles> perMCCol = aod::mcparticle::mcCollisionId;
  Preslice<FiDataTracks> perDataCol = aod::track::collisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 2.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isPbPb{"isPbPb", false, "Is Pb-Pb"};
  Configurable<bool> bookFactorizedCorrections{"bookFactorizedCorrections", true, "Additional event/signal-loss inputs in MCCounting; leaves nominal corrections unchanged"};
  Configurable<bool> bookMultiplicityCalibration{"bookMultiplicityCalibration", false, "Independent raw generated-FT0 multiplicity calibration and correction inputs; no change to nominal classes"};
  ConfigurableAxis calibrationRecoCentrality{"calibrationRecoCentrality", {1000, 0., 100.}, "Data FT0M percentile axis for offline class boundaries"};

  Configurable<float> dcaV0Dau{"dcaV0Dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcaNegToPV{"dcaNegToPV", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPV{"dcaPosToPV", 0.06, "DCA Pos To PV"};
  Configurable<double> v0CosPA{"v0CosPA", 0.97, "V0 CosPA"};
  Configurable<float> v0Radius{"v0Radius", 0.5, "Radius"};
  Configurable<float> etaDau{"etaDau", 4, "Eta Daughters"};
  Configurable<float> v0Rapidity{"v0Rapidity", 0.5, "V0 rapidity"};
  ConfigurableAxis centrality{"centrality", {VARIABLE_WIDTH, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100}, ""};
  Configurable<bool> zvtxcut{"zvtxcut", false, "z vtx cut < 10cm"};
  // Only used by processMCUnfolding. Neither the nominal centrality nor any
  // nominal histogram depends on these options or on McCentFT0Ms.
  ConfigurableAxis unfoldingCentrality{"unfoldingCentrality", {VARIABLE_WIDTH, 0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100}, "Independent generator/reconstructed FT0M response binning"};
  Configurable<float> unfoldingZMax{"unfoldingZMax", 10.f, "Truth AND measured |z| limit for the separate unfolding sample"};

  HistogramRegistry registry{
    "registry",
    {{"Selection", ";status;events", {HistType::kTH1F, {{17, 0.5, 17.5}}}}}};

  AxisSpec zAxis = {{-30, -20, -15, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30}, "Z (cm)", "zaxis"};
  AxisSpec deltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
  AxisSpec dcaAxis = {601, -3.01, 3.01, "", "DCA axis"};
  AxisSpec etaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
  AxisSpec v0EtaAxis = {20, -1.0, 1.0, "#etav0", "eta axis"};
  AxisSpec phiAxis = {629, 0, 2 * o2::constants::math::PI, "Rad", "phi axis"};
  AxisSpec ptVarAxis = {static_cast<int>(PtVarEnd) - 1, static_cast<int>(PtVarBegin) + 0.5f, static_cast<int>(PtVarEnd) - 0.5f, "", "ptvar axis"};
  AxisSpec evtClassAxis = {static_cast<int>(EvtEnd) - 1, static_cast<int>(EvtBegin) + 0.5f, static_cast<int>(EvtEnd) - 0.5f, "", "event class"};
  AxisSpec trigClassAxis = {static_cast<int>(TrigEnd) - 1, static_cast<int>(TrigBegin) + 0.5f, static_cast<int>(TrigEnd) - 0.5f, "", "trigger class"};
  AxisSpec particleTypeAxis = {static_cast<int>(PartEnd) - 1, static_cast<int>(PartBegin) + 0.5f, static_cast<int>(PartEnd) - 0.5f, "", "Particle type"};
  AxisSpec speciesAxis = {static_cast<int>(SpecEnd) - 1, static_cast<int>(SpecBegin) + 0.5f, static_cast<int>(SpecEnd) - 0.5f, "", "species class"};
  AxisSpec massAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};
  AxisSpec signAxis = {static_cast<int>(SignEnd) - 1, static_cast<int>(SignBegin) + 0.5f, static_cast<int>(SignEnd) - 0.5f, "", "sign"};
  AxisSpec stepAxis = {static_cast<int>(StepEnd) - 1, static_cast<int>(StepBegin) + 0.5f, static_cast<int>(StepEnd) - 0.5f, "", "step"};
  AxisSpec testAxis = {101, -0.5, 100.5, "", "test"};
  AxisSpec multAxis = {1001, -0.5, 1000.5, "", "Ntrks"};
  AxisSpec statusCodeAxis = {3, -1.5, 2.5, "", "StatusCode"};
  AxisSpec processCodeAxis = {45, -1.5, 44.5, "", "StatusCode"};
  AxisSpec phiBin = {{0, o2::constants::math::PI / 2, o2::constants::math::PI, o2::constants::math::PI * 3. / 2, 2 * o2::constants::math::PI}, "#phi", "phi bin"};

  expressions::Filter trackSelectionProper =
    ((aod::track::trackCutFlag & TrackSelectionITS) == TrackSelectionITS) &&
    ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
           ncheckbit(aod::track::trackCutFlag, TrackSelectionTPC),
           true) &&
    ncheckbit(aod::track::trackCutFlag, TrackSelectionDCA);

  expressions::Filter preFilterV0 =
    nabs(aod::v0data::dcapostopv) > dcaPosToPV&& nabs(aod::v0data::dcanegtopv) > dcaNegToPV&& aod::v0data::dcaV0daughters < dcaV0Dau;

  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Partition<FiDataTracks> tSampleData = nabs(aod::track::eta) < estimatorEta;

  void init(InitContext&)
  {
    if (doprocessMCCounting && doprocessMCCountingPVCalibrated) {
      LOGP(fatal, "Choose either detector-level MCCounting or PV-calibrated MCCounting, never both");
    }
    if (doprocessMCCountingPVCalibrated) {
      if (isPbPb || doprocessGen || doprocessMCUnfolding) {
        LOGP(fatal, "PV-matched classes support pp corrections only; disable processGen and quantile unfolding");
      }
      registry.add("pvMatched/hGenClassINELgt0", "PV-matched class, generated INEL>0;class;events", HistType::kTH1D, {{centrality, "PV-matched FT0M class (%)"}});
    }
    if (bookMultiplicityCalibration) {
      if (isPbPb) {
        LOGP(fatal, "bookMultiplicityCalibration currently supports pp FT0M only");
      }
      const AxisSpec nForwardAxis{5001, -0.5, 5000.5, "N_{FT0M,gen} (charged physical primaries)"};
      const AxisSpec nPVAxis{1001, -0.5, 1000.5, "Reconstructed PV contributors"};
      const AxisSpec calibCentAxis{calibrationRecoCentrality, "Reconstructed FT0M percentile (%)"};
      const AxisSpec recoCountAxis{11, -0.5, 10.5, "Number of reconstructed collisions"};
      // Raw integer N_forward, NOT a preassigned percentile. Retain all event
      // classes (including non-INEL>0) so event contamination can be audited.
      if (doprocessCounting) {
        registry.add("multCalib/hDataEvents", "Selected data;trigger;z_rec;cent_rec;N_PV", HistType::kTHnSparseD, {trigClassAxis, zAxis, calibCentAxis, nPVAxis});
        registry.add("multCalib/hDataTracks", "Selected data tracks;trigger;z_rec;eta_rec;cent_rec", HistType::kTHnSparseD, {trigClassAxis, zAxis, etaAxis, calibCentAxis});
      }
      if (doprocessMCCounting || doprocessMCCountingPVCalibrated) {
        registry.add("multCalib/hGenEvents", "All truth including zero reco;eventClass;z_gen;N_forward;N_reco", HistType::kTHnSparseD, {evtClassAxis, zAxis, nForwardAxis, recoCountAxis});
        registry.add("multCalib/hGenParticles", "All charged primaries;eventClass;z_gen;eta_gen;N_forward", HistType::kTHnSparseD, {evtClassAxis, zAxis, etaAxis, nForwardAxis});
        registry.add("multCalib/hSelectedEvents", "Selected best reco;eventClass;trigger;z_gen;z_rec;N_forward;N_PV", HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, zAxis, nForwardAxis, nPVAxis});
        registry.add("multCalib/hSelectedGenParticles", "All charged primaries in selected event;eventClass;trigger;z_gen;z_rec;eta_gen;N_forward", HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, zAxis, etaAxis, nForwardAxis});
        registry.add("multCalib/hSelectedRecTracks", "All filtered tracks in selected best reco;eventClass;trigger;z_gen;z_rec;eta_rec;N_forward", HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, zAxis, etaAxis, nForwardAxis});
        registry.add("multCalib/hClassCorrelation", "Selected best reco;trigger;z_gen;z_rec;N_forward;cent_rec", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, nForwardAxis, calibCentAxis});
      }
    }
    if (doprocessMCUnfoldingCalib) {
      const AxisSpec trueMultAxis{5001, -0.5, 5000.5, "Generated charged-primary FT0A+FT0C multiplicity"};
      registry.add("unfoldingCalib/hFT0MInelGt0", "Generated INEL>0, all z", HistType::kTH1D, {trueMultAxis});
      registry.add("unfoldingCalib/hFT0MInelGt0Z10", "Generated INEL>0, |z_gen|<10 cm", HistType::kTH1D, {trueMultAxis});
    }
    if (doprocessMCCounting && doprocessGen) {
      LOGP(fatal, "processMCCounting and processGen both fill the nominal generator histograms; enable only processMCCounting for MC corrections");
    }
    if (doprocessMCUnfolding) {
      if (isPbPb || !std::isfinite(unfoldingZMax.value) || unfoldingZMax <= 0.f || unfoldingZMax > 30.f ||
          !std::isfinite(estimatorEta.value) || estimatorEta <= 0.f || estimatorEta > 4.f) {
        LOGP(fatal, "processMCUnfolding requires pp, 0 < unfoldingZMax <= 30 cm and 0 < estimatorEta <= 4");
      }
      const AxisSpec responseCentAxis = {unfoldingCentrality, "FT0M percentile (%)"};
      const AxisSpec extraReasonAxis = {9, 0.5, 9.5, "1:no collision label;2:outside truth event sample;3:split non-best;4:invalid reco percentile;5:no particle label;6:other MC collision;7:not charged primary;8:outside truth eta;9:duplicate track"};
      // Full truth denominators include zero-reco and failed-selection events.
      // Primary particles use unit weights; this is an independent nominal
      // unfolding input, not a replacement for the legacy pT-variation maps.
      registry.add("unfolding/hGenEvents", "Generated INEL>0;z_gen;cent_gen", HistType::kTHnSparseD, {zAxis, responseCentAxis});
      registry.add("unfolding/hGenParticles", "Generated charged primaries;z_gen;eta_gen;cent_gen;species", HistType::kTHnSparseD, {zAxis, etaAxis, responseCentAxis, particleTypeAxis});
      registry.add("unfolding/hRecEvents", "All selected MC reconstructed collisions;trigger;z_rec;cent_rec", HistType::kTHnSparseD, {trigClassAxis, zAxis, responseCentAxis});
      registry.add("unfolding/hRecParticles", "All selected MC tracks;trigger;z_rec;eta_rec;cent_rec", HistType::kTHnSparseD, {trigClassAxis, zAxis, etaAxis, responseCentAxis});
      registry.add("unfolding/hEventResponse", "One best collision per generated event;trigger;z_gen;z_rec;cent_gen;cent_rec", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, responseCentAxis, responseCentAxis});
      registry.add("unfolding/hParticleResponse", "Unique primary track matches;trigger;z_gen;z_rec;eta_gen;eta_rec;cent_gen;cent_rec;species", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, etaAxis, etaAxis, responseCentAxis, responseCentAxis, particleTypeAxis});
      registry.add("unfolding/hMissEvents", "Truth events without a selected match;trigger;z_gen;cent_gen", HistType::kTHnSparseD, {trigClassAxis, zAxis, responseCentAxis});
      registry.add("unfolding/hMissParticles", "Truth particles without a selected match;trigger;z_gen;eta_gen;cent_gen;species", HistType::kTHnSparseD, {trigClassAxis, zAxis, etaAxis, responseCentAxis, particleTypeAxis});
      registry.add("unfolding/hExtraEvents", "Measured collisions outside the one-to-one signal response;trigger;z_rec;cent_rec;reason", HistType::kTHnSparseD, {trigClassAxis, zAxis, responseCentAxis, extraReasonAxis});
      registry.add("unfolding/hExtraParticles", "Measured tracks outside the one-to-one signal response;trigger;z_rec;eta_rec;cent_rec;reason", HistType::kTHnSparseD, {trigClassAxis, zAxis, etaAxis, responseCentAxis, extraReasonAxis});
    }
    AxisSpec centAxis = {centrality, "", "centrality"};
    AxisSpec centAxisPbPb = {centrality, "", "centrality"};
    if ((doprocessMCCounting || doprocessMCCountingPVCalibrated) && bookFactorizedCorrections) {
      // Independent inputs for detector correction followed by f_event/f_signal.
      // All truth includes zero-reco. Selected truth counts ALL charged primaries
      // in the selected best collision's MC event, not only track-matched ones.
      // Keep both vertex coordinates to permit explicit, consistent acceptance.
      const AxisSpec lossCentAxis = {centrality, "Nominal detector percentile (%)"};
      registry.add("factorized/hGenEvents", "INEL>0 before selection;z_gen;cent", HistType::kTHnSparseD, {zAxis, lossCentAxis});
      registry.add("factorized/hGenParticles", "Charged primaries before selection;z_gen;eta_gen;cent", HistType::kTHnSparseD, {zAxis, etaAxis, lossCentAxis});
      registry.add("factorized/hSelectedEvents", "Unique INEL>0 events with selected best reco;trigger;z_gen;z_rec;cent", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, lossCentAxis});
      registry.add("factorized/hSelectedGenParticles", "ALL charged primaries in selected events;trigger;z_gen;z_rec;eta_gen;cent", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, etaAxis, lossCentAxis});
      registry.add("factorized/hSelectedRecParticles", "ALL filtered tracks in selected events;trigger;z_gen;z_rec;eta_rec;cent", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, etaAxis, lossCentAxis});
    }
    if (doprocessMCCounting || doprocessMCCountingPVCalibrated) {
      registry.add("selectionLoss/hAllEvents", "All INEL>0;z_gen;cent", HistType::kTHnSparseD, {zAxis, centAxis});
      registry.add("selectionLoss/hSelectedEvents", "Selected INEL>0;trigger;z_gen;z_rec;cent", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, centAxis});
      registry.add("selectionLoss/hAllParticles", "All INEL>0 primaries;z_gen;eta;cent;species;ptvar;phi", HistType::kTHnSparseD, {zAxis, etaAxis, centAxis, particleTypeAxis, ptVarAxis, phiBin});
      registry.add("selectionLoss/hSelectedParticles", "All primaries in selected INEL>0 events;trigger;z_gen;z_rec;eta;cent;species;ptvar;phi", HistType::kTHnSparseD, {trigClassAxis, zAxis, zAxis, etaAxis, centAxis, particleTypeAxis, ptVarAxis, phiBin});
    }
    registry.add({"hetaresponse", ";etaresponse", {HistType::kTH2D, {{80, -4, 4}, {80, -4, 4}}}});
    registry.add({"hft0multiplicity", ";multiplicity", {HistType::kTH1D, {{10000, 0, 100000}}}});
    registry.add({"hcentrality", isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH1F, {{10000, 0, 100}}}});
    registry.add({"hcentralityvscentraldndeta", isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}}});
    registry.add({"hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, phiBin}}});
    registry.add({"hreczvtx", "evntclass; triggerclass;  zvtex", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hphieta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, phiAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hrecdndetamissing", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hgendndeta", "evntclass;  zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, ptVarAxis, phiBin}}});
    registry.add({"hgenzvtx", "evntclass; zvtex", {HistType::kTHnSparseD, {evtClassAxis, zAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hv0mass", "etaaxis; invmass", {HistType::kTHnSparseD, {isPbPb ? centAxisPbPb : centAxis, speciesAxis, v0EtaAxis, massAxis}}});
    registry.add({"hv0k0s", "invmass", {HistType::kTH1D, {{100, 0.4, 0.6}}}});

    registry.add({"recetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0Sel8genz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0genz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});

    registry.add({"reczINELg0Sel8", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"genzINELg0Sel8", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"genzINELg0", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"hcentmult", ";status;events", {HistType::kTH1D, {{100, 0, 100}}}});

    // MC-only bookkeeping for the zero-reconstructed-vertex FT0M fallback.
    // These histograms are filled for generated INEL>0 collisions with
    // |z_gen| < 10 cm, i.e. the event sample relevant to the correction.
    registry.add({"hMCRecoCollisionsPerGenINELg0", ";N_{reco collisions};events", {HistType::kTH1D, {{6, -0.5, 5.5}}}});
    registry.add({"hMCCentralitySourceINELg0", ";centrality source;events", {HistType::kTH1D, {{6, 0.5, 6.5}}}});
    registry.add({"hZeroRecoBCFT0MPercentileINELg0", ";BC FT0M percentile;events", {HistType::kTH1D, {{103, -1.5, 101.5}}}});
    registry.add({"hZeroRecoBCFT0MClassINELg0", ";BC FT0M class;events", {HistType::kTH1D, {{13, 0.5, 13.5}}}});
    registry.add({"hBCFT0MvsRecoFT0MINELg0", ";collision FT0M percentile;BC FT0M percentile", {HistType::kTH2D, {{100, 0, 100}, {102, 0, 102}}}});

    constexpr int NLogBins = 50;
    constexpr double PtLogLow = 0.01;
    constexpr double PtLogHigh = 10.0;
    std::vector<double> logBins(NLogBins + 1, 0.0);
    double logBinWidth = (std::log(PtLogHigh) - std::log(PtLogLow)) / NLogBins;
    for (int i = 0; i <= NLogBins; i++) {
      logBins[i] = PtLogLow * std::exp(i * logBinWidth);
    }
    AxisSpec ptBins = {logBins, "pT (GeV/c)", "pt bin"};

    if (doprocessMCCounting || doprocessMCCountingPVCalibrated) {
      std::vector<double> trackingPtEdges;
      for (int i = 0; i <= 100; ++i) {
        trackingPtEdges.push_back(0.01 * std::pow(200.0 / 0.01, i / 100.0));
      }
      const AxisSpec trackingPtAxis{trackingPtEdges, "p_{T}^{gen} (GeV/c)"};
      registry.add("trackingPt/hAllGenerated", "Charged physical primaries;eventClass;FT0M;pt_gen", HistType::kTHnSparseD, {evtClassAxis, isPbPb ? centAxisPbPb : centAxis, trackingPtAxis});
      registry.add("trackingPt/hSelectedGenerated", "Denominator: same selected best collision;eventClass;trigger;FT0M;pt_gen", HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, isPbPb ? centAxisPbPb : centAxis, trackingPtAxis});
      registry.add("trackingPt/hMatchedPrimary", "Numerator: unique matched physical primaries;eventClass;trigger;FT0M;pt_gen", HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, isPbPb ? centAxisPbPb : centAxis, trackingPtAxis});
    }

    registry.add({"hrecdndpt", " pt", {HistType::kTH1D, {ptBins}}});
    registry.add({"hdndptefficiency", " pt", {HistType::kTH1D, {ptBins}}});
    registry.add({"hgendndpt", " pt", {HistType::kTH1D, {ptBins}}});
    registry.add({"hgendndpt2", " pt", {HistType::kTH1D, {ptBins}}});
    registry.add({"hgendndpt05", " pt", {HistType::kTH1D, {ptBins}}});

    auto hstat = registry.get<TH1>(HIST("Selection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(SelAll, "All");
    x->SetBinLabel(SelSel8, "Sel8");
    x->SetBinLabel(SelSel8z10, "Sel8z10");
    x->SetBinLabel(SelGeneratedInelGt0, "Generated INEL>0");
    x->SetBinLabel(SelGoodBCs, "Good BCs");
    x->SetBinLabel(SelBCsWithCollisions, "BCs with collisions");
    x->SetBinLabel(SelBCsWithPileup, "BCs with pile-up/splitting");
    x->SetBinLabel(SelInelSel8Mcz10, "INEL&Sel8&mcz10");
    x->SetBinLabel(SelInel, "INEL");
    x->SetBinLabel(SelInelg0, "INELg0");
    x->SetBinLabel(SelInelGenz10, "INELgenz10");
    x->SetBinLabel(SelInelg0Genz10, "INELg0genz10");
    x->SetBinLabel(SelInelSel8, "INELSel8");
    x->SetBinLabel(SelInelg0Sel8g0, "INELg0Sel8g0");
    x->SetBinLabel(SelInelSel8Recz10, "INELSel8recz10");
    x->SetBinLabel(SelInelg0Sel8g0Recz10, "INELg0Sel8g0recz10");
    x->SetBinLabel(SelK0sz10Eta05, "K0Sz10eta05");

    auto* centSourceAxis = registry.get<TH1>(HIST("hMCCentralitySourceINELg0"))->GetXaxis();
    centSourceAxis->SetBinLabel(CentSourceRecoFT0M, "reco FT0M");
    centSourceAxis->SetBinLabel(CentSourceZeroRecoBCFT0M, "zero reco: BC FT0M");
    centSourceAxis->SetBinLabel(CentSourceZeroRecoNoFT0, "zero reco: no FT0");
    centSourceAxis->SetBinLabel(CentSourceZeroRecoInvalidBCFT0M, "zero reco: invalid BC FT0M");
    centSourceAxis->SetBinLabel(CentSourceZeroRecoNoBC, "zero reco: no BC");
    centSourceAxis->SetBinLabel(CentSourceInvalidRecoFT0M, "invalid reco FT0M");

    auto* zeroRecoClassAxis = registry.get<TH1>(HIST("hZeroRecoBCFT0MClassINELg0"))->GetXaxis();
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent0To1, "0-1%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent1To5, "1-5%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent5To10, "5-10%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent10To15, "10-15%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent15To20, "15-20%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent20To30, "20-30%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent30To40, "30-40%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent40To50, "40-50%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent50To70, "50-70%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCent70To100, "70-100%");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCentNoFT0, "no FT0");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCentInvalid, "invalid percentile");
    zeroRecoClassAxis->SetBinLabel(ZeroRecoCentNoBC, "no associated BC");
  }

  void processEventStat(
    const FullBCs& bcs,
    const soa::Join<aod::Collisions, aod::EvSels>& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (const auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(o2::aod::evsel::kIsBBT0A) &&
                        bc.selection_bit(o2::aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("Selection"), SelGoodBCs);
        cols.clear();
        for (const auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("Selection"), SelBCsWithCollisions);
          if (cols.size() > 1) {
            registry.fill(HIST("Selection"), SelBCsWithPileup);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaHi, processEventStat, "Collect event sample stats", false);

  template <typename C>
  void runCounting(
    const C& collisions,
    const FiDataTracks& /*tracks*/)
  {
    for (const auto& collision : collisions) {
      BoolVector btrigc(TrigEnd, false);
      registry.fill(HIST("Selection"), SelAll);
      auto z = collision.posZ();
      auto pertracks = tSampleData->sliceBy(perDataCol, collision.globalIndex());
      auto nTrk = 0;

      if (collision.sel8()) {
        btrigc[TrigSel8] = true;
        registry.fill(HIST("Selection"), SelSel8);
        if (std::abs(z) < CutZ) {
          registry.fill(HIST("Selection"), SelSel8z10);
        }
      }
      if (btrigc[TrigSel8] && std::abs(z) < CutZ) {
        registry.fill(HIST("hft0multiplicity"), collision.multFT0C());
      }

      for (const auto& track : pertracks) {
        [[maybe_unused]] int dummy = track.globalIndex();
        if (std::abs(track.eta()) < CutEta1) {
          nTrk++;
        }
      }
      if (nTrk > 0) {
        if (btrigc[TrigSel8]) {
          btrigc[TrigSel8g0] = true;
        }
      }
      if (btrigc[TrigSel8g0]) {
        registry.fill(HIST("reczINELg0Sel8"), z);
      }

      auto cent = -1.f;
      if (isPbPb) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          cent = collision.centFT0C();
        }
      } else {
        if constexpr (C::template contains<aod::CentFT0Ms>()) {
          cent = collision.centFT0M();
        }
      }

      if (bookMultiplicityCalibration) {
        const auto calibCent = std::isfinite(cent) && cent >= 0.f && cent <= 100.f ? std::min(cent, std::nextafter(100.f, 0.f)) : -1.f;
        for (auto trigger = 1u; trigger < static_cast<unsigned>(TrigEnd); ++trigger) {
          if (!btrigc[trigger]) {
            continue;
          }
          registry.fill(HIST("multCalib/hDataEvents"), trigger, z, calibCent, collision.multNTracksPV());
          for (const auto& track : pertracks) {
            registry.fill(HIST("multCalib/hDataTracks"), trigger, z, track.eta(), calibCent);
          }
        }
      }

      if (isPbPb) {
        if (std::abs(z) < CutZ && btrigc[TrigSel8]) {
          registry.fill(HIST("hcentrality"), cent);
        }
      } else {
        if (std::abs(z) < CutZ && btrigc[TrigSel8g0]) {
          registry.fill(HIST("hcentrality"), cent);
          registry.fill(HIST("hcentralityvscentraldndeta"), cent, nTrk);
        }
      }
      for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); itrigc++) {
        if (btrigc[itrigc]) {
          registry.fill(HIST("hreczvtx"), EvtData, itrigc, z, cent);
        }
      }

      for (const auto& track : pertracks) {
        if (btrigc[TrigSel8] && std::abs(track.eta()) < CutEta08 && std::abs(z) < CutZ) {
          registry.fill(HIST("hrecdndpt"), track.pt());
        }
        if (btrigc[TrigSel8]) {
          registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);
        }
        for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); itrigc++) {
          if (btrigc[itrigc]) {
            registry.fill(HIST("hphieta"), EvtData, itrigc, track.phi(), track.eta(), cent);
            registry.fill(HIST("hrecdndeta"), EvtData, itrigc, z, track.eta(), cent, PartData, track.phi());
          }
        }
      }
    }
  }

  PresliceUnsorted<soa::Join<MyCollisions, aod::McCollisionLabels>> perMcCol = o2::aod::mccollisionlabel::mcCollisionId;
  Preslice<aod::McParticles> perMCColparticles = aod::mcparticle::mcCollisionId;
  Preslice<FiLTracks> perColFiLTracks = aod::track::collisionId;
  // McCollsExtra provides one deterministic reconstructed collision for split
  // MC collisions (the collision with the largest number of contributors).
  // Select exactly one correction mode through the process switches below:
  // processMCCounting: official reconstructed FT0M + zero-reco BC fallback.
  // processMCCountingPVCalibrated: McCentFT0Ms for both truth and reco.
  // Both modes always fill ND/SD/DD INEL>0 categories, selectionLoss and
  // trackingPt inputs; no separate booking switches are required.
  using MCex = soa::Join<aod::McCollisions, aod::HepMCXSections, aod::McCollsExtra>;

  template <typename MCCollisions>
  void runMCCounting(
    const MCCollisions& mcCollisions,
    const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
    const BCsRun3Cent& /*bcs*/,
    const aod::FT0s& /*ft0s*/,
    const Particles& mcParticles,
    const FiLTracks& tracks)
  {
    for (const auto& mcCollision : mcCollisions) {
      BoolVector bevtc(EvtEnd, false);
      bevtc[EvtInel] = true;
      auto procId = mcCollision.processId();
      if (procId == ProcIdND) {
        bevtc[EvtNd] = true;
      } else if (procId == ProcIdSD1 || procId == ProcIdSD2) {
        bevtc[EvtSd] = true;
      } else if (procId == ProcIdDD1 || procId == ProcIdDD2) {
        bevtc[EvtDd] = true;
      }
      registry.fill(HIST("Selection"), SelAll);

      auto mcz = mcCollision.posZ();
      auto genz = mcz;
      auto nTrkGen = 0;
      int nForwardGen = 0;
      auto particles = mcParticles.sliceBy(perMCCol, mcCollision.globalIndex());

      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= CutCharge) {
            if (bookMultiplicityCalibration && ((particle.eta() > 3.5f && particle.eta() < 4.9f) || (particle.eta() > -3.3f && particle.eta() < -2.1f))) {
              ++nForwardGen;
            }
            if (std::abs(particle.eta()) < CutEta1) {
              nTrkGen++;
            }
          }
        }
      }
      if (nTrkGen > 0) {
        registry.fill(HIST("Selection"), SelGeneratedInelGt0);
        bevtc[EvtInelg0] = true;
        // Mandatory process-resolved INEL>0 correction categories in BOTH modes.
        bevtc[EvtInelg0Dd] = bevtc[EvtDd];
        bevtc[EvtInelg0Sd] = bevtc[EvtSd];
        bevtc[EvtInelg0Nd] = bevtc[EvtNd];
      }
      if (bevtc[EvtInel]) {
        registry.fill(HIST("Selection"), SelInel);
      }
      if (bevtc[EvtInel] && std::abs(mcz) < CutZ) {
        registry.fill(HIST("Selection"), SelInelGenz10);
      }
      if (bevtc[EvtInelg0]) {
        registry.fill(HIST("Selection"), SelInelg0);
      }
      if (bevtc[EvtInelg0] && std::abs(mcz) < CutZ) {
        registry.fill(HIST("Selection"), SelInelg0Genz10);
      }

      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= CutCharge) {
            if (bevtc[EvtInel] && std::abs(particle.eta()) < CutEta08 && std::abs(mcz) < CutZ) {
              registry.fill(HIST("hgendndpt"), particle.pt());
              if (particle.pt() < CutPt) {
                registry.fill(HIST("hgendndpt2"), particle.pt(), PtWeightSlopeUp * particle.pt() + PtWeightInterceptUp);
                registry.fill(HIST("hgendndpt05"), particle.pt(), PtWeightSlopeDown * particle.pt() + PtWeightInterceptDown);
              } else {
                registry.fill(HIST("hgendndpt2"), particle.pt());
                registry.fill(HIST("hgendndpt05"), particle.pt());
              }
            }
          }
        }
      }

      auto collisionsample = collisions.sliceBy(perMcCol, mcCollision.globalIndex());
      const auto nRecoCollisions = collisionsample.size();
      const bool hasRecoCollision = nRecoCollisions > 0;
      const bool isCorrectionSample = bevtc[EvtInelg0] && std::abs(mcz) < CutZ;

      if (isCorrectionSample) {
        registry.fill(HIST("hMCRecoCollisionsPerGenINELg0"), nRecoCollisions);
      }

      // Detector-level definition:
      //   * >=1 reconstructed collision: official collision FT0M percentile;
      //   * 0 reconstructed collisions: use valid FT0M information from the
      //     associated BC, so these generated events and particles remain in
      //     a detector-estimator class even though reconstruction failed.
      //     If no usable BC FT0M exists, retain underflow and record the reason.
      // BC FT0M and collision FT0M have different vertex-equalization inputs;
      // the fallback applies only to zero-reco events and is monitored below.
      // For split vertices only the deterministic best reconstructed collision
      // defines the nominal class and contributes to reconstructed histograms.
      auto centForCorrection = NoCentrality;
      auto centSource = CentSourceInvalidRecoFT0M;
      auto zeroRecoClass = ZeroRecoCentInvalid;
      auto rawZeroRecoBCCent = NoCentrality;

      if (hasRecoCollision && mcCollision.bestCollisionIndex() >= 0) {
        for (const auto& collision : collisionsample) {
          if (collision.globalIndex() == mcCollision.bestCollisionIndex()) {
            if (isPbPb) {
              if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>()) {
                centForCorrection = collision.centFT0C();
                centSource = CentSourceRecoFT0M;
              }
            } else {
              if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>()) {
                const auto centBest = collision.centFT0M();
                if (std::isfinite(centBest) && centBest >= 0.f && centBest <= 100.f) {
                  centForCorrection = std::min(centBest, std::nextafter(100.f, 0.f));
                  centSource = CentSourceRecoFT0M;
                } else {
                  LOGP(warn, "Invalid reconstructed FT0M percentile {} for best collision {}; keeping this event in the centrality underflow", centBest, collision.globalIndex());
                }
              }
            }
            break;
          }
        }
      } else if (!hasRecoCollision && !isPbPb) {
        // As in the Run-2 event-level estimator approach, a failed vertex
        // reconstruction must not discard otherwise available detector
        // percentile information from the generated correction denominator.
        // Use the BC estimator only when it has a matched FT0 and a valid
        // percentile; do not invent a class for missing detector information.
        centSource = CentSourceZeroRecoNoBC;
        zeroRecoClass = ZeroRecoCentNoBC;
        if (mcCollision.has_bc()) {
          const auto mcBC = mcCollision.template bc_as<BCsRun3Cent>();
          rawZeroRecoBCCent = mcBC.centFT0M();
          if (!mcBC.has_foundFT0()) {
            centSource = CentSourceZeroRecoNoFT0;
            zeroRecoClass = ZeroRecoCentNoFT0;
          } else if (std::isfinite(rawZeroRecoBCCent) && rawZeroRecoBCCent >= 0.f && rawZeroRecoBCCent <= 100.f) {
            centForCorrection = std::min(rawZeroRecoBCCent, std::nextafter(100.f, 0.f));
            centSource = CentSourceZeroRecoBCFT0M;
            zeroRecoClass = getZeroRecoCentralityClass(rawZeroRecoBCCent);
          } else {
            centSource = CentSourceZeroRecoInvalidBCFT0M;
            zeroRecoClass = ZeroRecoCentInvalid;
          }
        }
      }

      if (isCorrectionSample) {
        registry.fill(HIST("hMCCentralitySourceINELg0"), centSource);
        if (!hasRecoCollision) {
          if (std::isfinite(rawZeroRecoBCCent)) {
            registry.fill(HIST("hZeroRecoBCFT0MPercentileINELg0"), rawZeroRecoBCCent);
          }
          registry.fill(HIST("hZeroRecoBCFT0MClassINELg0"), zeroRecoClass);
        } else if (!isPbPb && centSource == CentSourceRecoFT0M && mcCollision.has_bc()) {
          const auto mcBC = mcCollision.template bc_as<BCsRun3Cent>();
          const auto bcCent = mcBC.centFT0M();
          if (mcBC.has_foundFT0() && std::isfinite(bcCent) && bcCent >= 0.f && bcCent <= 100.f) {
            registry.fill(HIST("hBCFT0MvsRecoFT0MINELg0"), centForCorrection, std::min(bcCent, std::nextafter(100.f, 0.f)));
          }
        }
      }

      // Preserve detector QA and raw calibration correlation. In the optional
      // PV-matched branch, ALL correction inputs use one MC-table class,
      // including events with no reconstructed collision. Nominal path unchanged.
      const auto detectorCentrality = centForCorrection;
      if constexpr (MCCollisions::template contains<aod::McCentFT0Ms>()) {
        const auto calibrated = mcCollision.centFT0M();
        if (!std::isfinite(calibrated) || calibrated < 0.f || calibrated >= 100.f) {
          LOGP(fatal, "Invalid PV-matched MC centrality {} for MC collision {}", calibrated, mcCollision.globalIndex());
        }
        centForCorrection = calibrated;
        if (bevtc[EvtInelg0]) {
          registry.fill(HIST("pvMatched/hGenClassINELgt0"), calibrated);
        }
      }

      for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
        if (bevtc[ievtc]) {
          registry.fill(HIST("hgenzvtx"), ievtc, genz, centForCorrection);
        }
      }

      int pid = 0;
      struct LossParticle {
        double eta, phi, pt;
        int species;
      };
      std::vector<LossParticle> lossParticles;
      std::vector<double> particleetas;
      for (const auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (std::abs(particle.pdgCode()) == kK0Short && std::abs(particle.eta()) < CutEta05 && std::abs(genz) < CutZ) {
          registry.fill(HIST("Selection"), SelK0sz10Eta05);
        }
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        switch (std::abs(particle.pdgCode())) {
          case kPiPlus:
            pid = PartPion;
            break;
          case kKPlus:
            pid = PartKaon;
            break;
          case kProton:
            pid = PartProton;
            break;
          default:
            pid = PartOther;
            break;
        }

        if (p != nullptr) {
          if (std::abs(p->Charge()) >= CutCharge) {
            for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
              if (bevtc[ievtc]) {
                if (std::abs(genz) < CutZ && std::abs(particle.eta()) < CutEta08) {
                  registry.fill(HIST("trackingPt/hAllGenerated"), ievtc, centForCorrection, particle.pt());
                }
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centForCorrection, pid, PtVarNone, particle.phi());
                if (particle.pt() < CutPt) {
                  const auto weightUp = PtWeightSlopeUp * particle.pt() + PtWeightInterceptUp;
                  const auto weightDown = PtWeightSlopeDown * particle.pt() + PtWeightInterceptDown;
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centForCorrection, pid, PtVarUp, particle.phi(), weightUp);
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centForCorrection, pid, PtVarDown, particle.phi(), weightDown);
                } else {
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centForCorrection, pid, PtVarUp, particle.phi());
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centForCorrection, pid, PtVarDown, particle.phi());
                }
              }
            }
            if (pid >= PartPion && pid <= PartOther) {
              particleetas.push_back(particle.eta());
              if (bevtc[EvtInelg0]) {
                lossParticles.push_back({particle.eta(), particle.phi(), particle.pt(), pid});
                registry.fill(HIST("selectionLoss/hAllParticles"), genz, particle.eta(), centForCorrection, pid, PtVarNone, particle.phi());
                const double up = particle.pt() < CutPt ? PtWeightSlopeUp * particle.pt() + PtWeightInterceptUp : 1.;
                const double down = particle.pt() < CutPt ? PtWeightSlopeDown * particle.pt() + PtWeightInterceptDown : 1.;
                registry.fill(HIST("selectionLoss/hAllParticles"), genz, particle.eta(), centForCorrection, pid, PtVarUp, particle.phi(), up);
                registry.fill(HIST("selectionLoss/hAllParticles"), genz, particle.eta(), centForCorrection, pid, PtVarDown, particle.phi(), down);
              }
            }
          }
        }
      }

      // Full generated INEL>0 denominator. This must not depend on whether the
      // event has zero, one, or several reconstructed collisions.
      if (bookMultiplicityCalibration) {
        for (auto eventClass = 1u; eventClass < static_cast<unsigned>(EvtEnd); ++eventClass) {
          if (!bevtc[eventClass]) {
            continue;
          }
          registry.fill(HIST("multCalib/hGenEvents"), eventClass, genz, nForwardGen, nRecoCollisions);
          for (const auto& eta : particleetas) {
            registry.fill(HIST("multCalib/hGenParticles"), eventClass, genz, eta, nForwardGen);
          }
        }
      }
      if (bevtc[EvtInelg0]) {
        for (const auto& eta : particleetas) {
          registry.fill(HIST("genetaINELg0genz10"), eta, mcz);
        }
        registry.fill(HIST("genzINELg0"), genz);
        registry.fill(HIST("selectionLoss/hAllEvents"), genz, centForCorrection);
        if (bookFactorizedCorrections) {
          registry.fill(HIST("factorized/hGenEvents"), genz, centForCorrection);
          for (const auto& eta : particleetas) {
            registry.fill(HIST("factorized/hGenParticles"), genz, eta, centForCorrection);
          }
        }
      }

      for (const auto& collision : collisionsample) {
        // A generated collision contributes at most one reconstructed
        // collision to the correction sample.
        const bool isBestCollision = collision.globalIndex() == mcCollision.bestCollisionIndex();
        if (!isBestCollision) {
          continue;
        }
        // The loop contains only the best reconstructed collision, so use the
        // same sanitized official collision centrality as the MC denominator.
        const auto centRec = centForCorrection;

        BoolVector btrigc(TrigEnd, false);
        auto z = collision.posZ();
        if (collision.sel8()) {
          btrigc[TrigSel8] = true;
          registry.fill(HIST("Selection"), SelSel8);
          if (std::abs(z) < CutZ) {
            registry.fill(HIST("Selection"), SelSel8z10);
          }
        }
        if (bevtc[EvtInel] && btrigc[TrigSel8] && std::abs(z) < CutZ) {
          registry.fill(HIST("hft0multiplicity"), collision.multFT0C());
          registry.fill(HIST("Selection"), SelInelSel8Mcz10);
        }
        if (isBestCollision && bevtc[EvtInelg0] && btrigc[TrigSel8]) {
          for (const auto& eta : particleetas) {
            registry.fill(HIST("genetaINELg0Sel8recz10"), eta, z);
            registry.fill(HIST("genetaINELg0Sel8genz10"), eta, mcz);
          }
          registry.fill(HIST("reczINELg0Sel8"), z);
          registry.fill(HIST("genzINELg0Sel8"), genz);
        }

        auto nTrkRec = 0;
        auto trackspart = tracks.sliceBy(perColFiLTracks, collision.globalIndex());
        for (const auto& track : trackspart) {
          if (std::abs(track.eta()) < CutEta1) {
            nTrkRec++;
          }
        }

        if (nTrkRec > 0) {
          if (btrigc[TrigSel8]) {
            btrigc[TrigSel8g0] = true;
          }
        }
        if (bookMultiplicityCalibration) {
          for (auto trigger = 1u; trigger < static_cast<unsigned>(TrigEnd); ++trigger) {
            if (!btrigc[trigger]) {
              continue;
            }
            registry.fill(HIST("multCalib/hClassCorrelation"), trigger, genz, z, nForwardGen, detectorCentrality);
            for (auto eventClass = 1u; eventClass < static_cast<unsigned>(EvtEnd); ++eventClass) {
              if (!bevtc[eventClass]) {
                continue;
              }
              registry.fill(HIST("multCalib/hSelectedEvents"), eventClass, trigger, genz, z, nForwardGen, collision.multNTracksPV());
              for (const auto& eta : particleetas) {
                registry.fill(HIST("multCalib/hSelectedGenParticles"), eventClass, trigger, genz, z, eta, nForwardGen);
              }
              for (const auto& track : trackspart) {
                registry.fill(HIST("multCalib/hSelectedRecTracks"), eventClass, trigger, genz, z, track.eta(), nForwardGen);
              }
            }
          }
        }

        if (bevtc[EvtInelg0]) {
          for (unsigned trigger = 1; trigger < static_cast<unsigned>(TrigEnd); ++trigger) {
            if (!btrigc[trigger])
              continue;
            registry.fill(HIST("selectionLoss/hSelectedEvents"), trigger, genz, z, centRec);
            for (const auto& p : lossParticles) {
              registry.fill(HIST("selectionLoss/hSelectedParticles"), trigger, genz, z, p.eta, centRec, p.species, PtVarNone, p.phi);
              const double up = p.pt < CutPt ? PtWeightSlopeUp * p.pt + PtWeightInterceptUp : 1.;
              const double down = p.pt < CutPt ? PtWeightSlopeDown * p.pt + PtWeightInterceptDown : 1.;
              registry.fill(HIST("selectionLoss/hSelectedParticles"), trigger, genz, z, p.eta, centRec, p.species, PtVarUp, p.phi, up);
              registry.fill(HIST("selectionLoss/hSelectedParticles"), trigger, genz, z, p.eta, centRec, p.species, PtVarDown, p.phi, down);
            }
          }
        }

        // Independent diagnostic: identical selected event and truth fiducial
        // denominator for every reconstructed numerator particle. Deduplicate
        // only these new counts; do not change nominal track categorization.
        if (std::abs(genz) < CutZ && std::abs(z) < CutZ) {
          std::unordered_set<int64_t> matchedPrimaryIds;
          for (const auto& track : trackspart) {
            if (!track.has_mcParticle() || std::abs(track.eta()) >= CutEta08) {
              continue;
            }
            const auto particle = track.template mcParticle_as<Particles>();
            if (particle.mcCollisionId() != mcCollision.globalIndex() || !particle.isPhysicalPrimary() || std::abs(particle.eta()) >= CutEta08) {
              continue;
            }
            const auto kp = pdg->GetParticle(particle.pdgCode());
            if (kp != nullptr && std::abs(kp->Charge()) >= CutCharge) {
              matchedPrimaryIds.insert(particle.globalIndex());
            }
          }
          for (const auto& particle : particles) {
            if (!particle.isPhysicalPrimary() || std::abs(particle.eta()) >= CutEta08) {
              continue;
            }
            const auto kp = pdg->GetParticle(particle.pdgCode());
            if (kp == nullptr || std::abs(kp->Charge()) < CutCharge) {
              continue;
            }
            const bool matched = matchedPrimaryIds.count(particle.globalIndex()) != 0;
            for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ++ievtc) {
              if (!bevtc[ievtc])
                continue;
              for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); ++itrigc) {
                if (!btrigc[itrigc])
                  continue;
                registry.fill(HIST("trackingPt/hSelectedGenerated"), ievtc, itrigc, centRec, particle.pt());
                if (matched) {
                  registry.fill(HIST("trackingPt/hMatchedPrimary"), ievtc, itrigc, centRec, particle.pt());
                }
              }
            }
          }
        }
        if (bookFactorizedCorrections && bevtc[EvtInelg0]) {
          // Best-collision guard above ensures each generated event is counted
          // at most once per trigger. No additional vertex/centrality cuts here.
          for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); ++itrigc) {
            if (!btrigc[itrigc]) {
              continue;
            }
            registry.fill(HIST("factorized/hSelectedEvents"), itrigc, genz, z, centRec);
            for (const auto& eta : particleetas) {
              registry.fill(HIST("factorized/hSelectedGenParticles"), itrigc, genz, z, eta, centRec);
            }
            // Include backgrounds, duplicates and unlabelled tracks, as in raw
            // data. Their net effect belongs in the detector correction ratio.
            for (const auto& track : trackspart) {
              registry.fill(HIST("factorized/hSelectedRecParticles"), itrigc, genz, z, track.eta(), centRec);
            }
          }
        }
        if (isPbPb) {
          if (std::abs(z) < CutZ && btrigc[TrigSel8]) {
            registry.fill(HIST("hcentrality"), centRec);
          }
        } else {
          if (std::abs(z) < CutZ && btrigc[TrigSel8g0]) {
            registry.fill(HIST("hcentrality"), centRec);
          }
        }

        if (bevtc[EvtInel] && btrigc[TrigSel8]) {
          registry.fill(HIST("Selection"), SelInelSel8);
        }
        if (bevtc[EvtInel] && btrigc[TrigSel8] && std::abs(z) < CutZ) {
          registry.fill(HIST("Selection"), SelInelSel8Recz10);
        }
        if (bevtc[EvtInelg0] && btrigc[TrigSel8g0]) {
          registry.fill(HIST("Selection"), SelInelg0Sel8g0);
        }
        if (bevtc[EvtInelg0] && btrigc[TrigSel8g0] && std::abs(z) < CutZ) {
          registry.fill(HIST("Selection"), SelInelg0Sel8g0Recz10);
        }

        for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
          for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); itrigc++) {
            if (bevtc[ievtc] && btrigc[itrigc]) {
              registry.fill(HIST("hreczvtx"), ievtc, itrigc, z, centRec);
            }
          }
        }

        std::vector<int> mclabels;
        for (const auto& track : trackspart) {
          if (track.has_mcParticle()) {
            int trackPid = PartBkg;
            auto particle = track.template mcParticle_as<Particles>();
            if (particle.isPhysicalPrimary()) {
              switch (std::abs(particle.pdgCode())) {
                case kPiPlus:
                  trackPid = PartPion;
                  break;
                case kKPlus:
                  trackPid = PartKaon;
                  break;
                case kProton:
                  trackPid = PartProton;
                  break;
                default:
                  trackPid = PartOther;
                  break;
              }
            } else {
              trackPid = PartNotPrimary;
            }
            for (auto motherIds = particle.mothersIds().front(); motherIds <= particle.mothersIds().back(); motherIds++) {
              auto mother = mcParticles.rawIteratorAt(motherIds);
              auto pdgMother = mother.pdgCode();
              if (pdgMother == kK0Short || std::abs(pdgMother) == kLambda0) {
                trackPid = PartMotherStrange;
              }
            }
            if (std::find(mclabels.begin(), mclabels.end(), track.mcParticleId()) != mclabels.end()) {
              trackPid = PartBkg;
            }
            mclabels.push_back(track.mcParticleId());
            registry.fill(HIST("hetaresponse"), particle.eta(), track.eta(), centRec);
            if (bevtc[EvtInel] && btrigc[TrigSel8] &&
                std::abs(track.eta()) < CutEta08 && std::abs(z) < CutZ &&
                trackPid != PartBkg && trackPid != PartNotPrimary) {
              registry.fill(HIST("hdndptefficiency"), particle.pt());
            }
            if (btrigc[TrigSel8] && bevtc[EvtInelg0]) {
              registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);
            }
            for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
              for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), ievtc, itrigc, z, particle.eta(), centRec, trackPid, particle.phi());
                  registry.fill(HIST("hphieta"), ievtc, itrigc, track.phi(), track.eta(), centRec);
                }
              }
            }
          } else {
            for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
              for (auto itrigc = 1u; itrigc < static_cast<unsigned>(TrigEnd); itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), ievtc, itrigc, z, track.eta(), centRec, PartBkg, track.phi());
                }
              }
            }
          }
        }
      }
    }
  }

  void processMCCounting(const MCex& mcCollisions,
                         const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
                         const BCsRun3Cent& bcs, const aod::FT0s& ft0s,
                         const Particles& particles, const FiLTracks& tracks)
  {
    runMCCounting(mcCollisions, collisions, bcs, ft0s, particles, tracks);
  }
  PROCESS_SWITCH(DndetaHi, processMCCounting, "MC Count tracks (detector FT0M with BC fallback)", false);

  using MCWithPVCalibration = soa::Join<MCex, aod::McCentFT0Ms>;
  void processMCCountingPVCalibrated(const MCWithPVCalibration& mcCollisions,
                                     const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
                                     const BCsRun3Cent& bcs, const aod::FT0s& ft0s,
                                     const Particles& particles, const FiLTracks& tracks)
  {
    runMCCounting(mcCollisions, collisions, bcs, ft0s, particles, tracks);
  }
  PROCESS_SWITCH(DndetaHi, processMCCountingPVCalibrated, "MC corrections using PV-matched McCentFT0Ms for both generated and reconstructed entries", false);

  // Optional second-pass branch. Enabling this alongside processMCCounting
  // records BOTH methods' inputs without filling any nominal histogram twice.
  // Only this signature requests the generator-percentile table/calibration.
  using MCWithUnfoldingCent = soa::Join<MCex, aod::McCentFT0Ms>;

  // First pass: independent of reconstructed vertices and percentile tables.
  void processMCUnfoldingCalib(aod::McCollisions const& collisions, aod::McParticles const& particles)
  {
    std::unordered_map<int64_t, std::array<int, 2>> counts;
    for (const auto& collision : collisions) {
      counts.emplace(collision.globalIndex(), std::array<int, 2>{0, 0});
    }
    for (const auto& particle : particles) {
      auto event = counts.find(particle.mcCollisionId());
      if (event == counts.end() || !particle.isPhysicalPrimary()) {
        continue;
      }
      const auto* definition = pdg->GetParticle(particle.pdgCode());
      if (!definition || std::abs(definition->Charge()) < CutCharge) {
        continue;
      }
      const auto eta = particle.eta();
      if (std::abs(eta) < CutEta1) {
        ++event->second[0];
      }
      if ((eta > 3.5f && eta < 4.9f) || (eta > -3.3f && eta < -2.1f)) {
        ++event->second[1];
      }
    }
    for (const auto& collision : collisions) {
      const auto& count = counts.at(collision.globalIndex());
      if (count[0] > 0) {
        registry.fill(HIST("unfoldingCalib/hFT0MInelGt0"), count[1]);
        if (std::abs(collision.posZ()) < CutZ) {
          registry.fill(HIST("unfoldingCalib/hFT0MInelGt0Z10"), count[1]);
        }
      }
    }
  }
  PROCESS_SWITCH(DndetaHi, processMCUnfoldingCalib, "First pass: generated INEL>0 FT0 multiplicity without percentile calibration", false);

  void processMCUnfolding(
    const MCWithUnfoldingCent& mcCollisions,
    const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
    const Particles& mcParticles,
    const FiLTracks& tracks)
  {
    struct TruthEvent {
      float z = 0.f;
      float cent = NoCentrality;
      int64_t bestCollision = -1;
      bool inelg0 = false;
      bool fiducial = false;
      std::array<bool, TrigEnd> matched{};
    };
    std::unordered_map<int64_t, TruthEvent> truth;
    truth.reserve(mcCollisions.size());
    std::array<std::unordered_set<int64_t>, TrigEnd> matchedParticles;

    auto isChargedPrimary = [this](const auto& particle) {
      if (!particle.isPhysicalPrimary()) {
        return false;
      }
      const auto* definition = pdg->GetParticle(particle.pdgCode());
      return definition != nullptr && std::abs(definition->Charge()) >= CutCharge;
    };
    auto species = [](const auto& particle) {
      switch (std::abs(particle.pdgCode())) {
        case kPiPlus:
          return PartPion;
        case kKPlus:
          return PartKaon;
        case kProton:
          return PartProton;
        default:
          return PartOther;
      }
    };
    auto validPercentile = [](float percentile) {
      return std::isfinite(percentile) && percentile >= 0.f && percentile <= 100.f;
    };
    auto insidePercentile = [](float percentile) {
      return std::min(percentile, std::nextafter(100.f, 0.f));
    };

    // Use global IDs, not row numbers: input slices may have nonzero offsets.
    for (const auto& collision : mcCollisions) {
      truth.emplace(collision.globalIndex(), TruthEvent{collision.posZ(), collision.centFT0M(), collision.bestCollisionIndex()});
    }
    for (const auto& particle : mcParticles) {
      auto event = truth.find(particle.mcCollisionId());
      if (event != truth.end() && isChargedPrimary(particle) && std::abs(particle.eta()) < CutEta1) {
        event->second.inelg0 = true;
      }
    }
    for (auto& [id, event] : truth) {
      event.fiducial = event.inelg0 && std::abs(event.z) < unfoldingZMax;
      if (!event.fiducial) {
        continue;
      }
      if (!validPercentile(event.cent)) {
        LOGP(fatal, "Invalid generated INEL>0 percentile {} for MC collision {}: processMCUnfolding requires the FIRST-PASS FT0M INEL>0 calibration, not reconstructed/BC fallback percentiles", event.cent, id);
      }
      event.cent = insidePercentile(event.cent);
      registry.fill(HIST("unfolding/hGenEvents"), event.z, event.cent);
    }

    // Traverse ALL measured collisions, including non-best split collisions
    // and collisions with no MC label. Data cannot apply the MC best label.
    for (const auto& collision : collisions) {
      const auto z = collision.posZ();
      if (!collision.sel8() || !(std::abs(z) < unfoldingZMax)) {
        continue;
      }
      const auto perTracks = tracks.sliceBy(perColFiLTracks, collision.globalIndex());
      bool hasTrackEta1 = false;
      for (const auto& track : perTracks) {
        if (std::abs(track.eta()) < estimatorEta && std::abs(track.eta()) < CutEta1) {
          hasTrackEta1 = true;
          break;
        }
      }
      const auto rawCent = collision.centFT0M();
      const bool validRecoCent = validPercentile(rawCent);
      const auto centRec = validRecoCent ? insidePercentile(rawCent) : NoCentrality;
      auto event = truth.find(collision.mcCollisionId());
      int eventExtraReason = 0;
      if (!collision.has_mcCollision() || event == truth.end()) {
        eventExtraReason = 1;
      } else if (!event->second.fiducial) {
        eventExtraReason = 2;
      } else if (collision.globalIndex() != event->second.bestCollision) {
        eventExtraReason = 3;
      } else if (!validRecoCent) {
        eventExtraReason = 4;
      }
      for (int trigger = TrigSel8; trigger < TrigEnd; ++trigger) {
        if (trigger == TrigSel8g0 && !hasTrackEta1) {
          continue;
        }
        registry.fill(HIST("unfolding/hRecEvents"), trigger, z, centRec);
        if (eventExtraReason != 0) {
          registry.fill(HIST("unfolding/hExtraEvents"), trigger, z, centRec, eventExtraReason);
        } else {
          event->second.matched[trigger] = true;
          registry.fill(HIST("unfolding/hEventResponse"), trigger, event->second.z, z, event->second.cent, centRec);
        }
        for (const auto& track : perTracks) {
          // Identical measured eta acceptance to runCounting/tSampleData.
          if (!(std::abs(track.eta()) < estimatorEta)) {
            continue;
          }
          registry.fill(HIST("unfolding/hRecParticles"), trigger, z, track.eta(), centRec);
          int extraReason = eventExtraReason;
          if (extraReason == 0 && !track.has_mcParticle()) {
            extraReason = 5;
          }
          if (extraReason == 0) {
            const auto particle = track.template mcParticle_as<Particles>();
            if (particle.mcCollisionId() != collision.mcCollisionId()) {
              extraReason = 6;
            } else if (!isChargedPrimary(particle)) {
              extraReason = 7;
            } else if (!(std::abs(particle.eta()) < estimatorEta)) {
              extraReason = 8;
            } else if (!matchedParticles[trigger].insert(particle.globalIndex()).second) {
              extraReason = 9;
            } else {
              registry.fill(HIST("unfolding/hParticleResponse"), trigger, event->second.z, z, particle.eta(), track.eta(), event->second.cent, centRec, species(particle));
            }
          }
          if (extraReason != 0) {
            registry.fill(HIST("unfolding/hExtraParticles"), trigger, z, track.eta(), centRec, extraReason);
          }
        }
      }
    }

    // Explicit miss histograms make the identities testable for EACH trigger:
    // truth = response projected on truth + misses;
    // measured MC = response projected on reco + extras.
    for (const auto& [id, event] : truth) {
      if (!event.fiducial) {
        continue;
      }
      for (int trigger = TrigSel8; trigger < TrigEnd; ++trigger) {
        if (!event.matched[trigger]) {
          registry.fill(HIST("unfolding/hMissEvents"), trigger, event.z, event.cent);
        }
      }
    }
    for (const auto& particle : mcParticles) {
      const auto event = truth.find(particle.mcCollisionId());
      if (event == truth.end() || !event->second.fiducial || !isChargedPrimary(particle) || !(std::abs(particle.eta()) < estimatorEta)) {
        continue;
      }
      const auto pid = species(particle);
      registry.fill(HIST("unfolding/hGenParticles"), event->second.z, particle.eta(), event->second.cent, pid);
      for (int trigger = TrigSel8; trigger < TrigEnd; ++trigger) {
        if (matchedParticles[trigger].find(particle.globalIndex()) == matchedParticles[trigger].end()) {
          registry.fill(HIST("unfolding/hMissParticles"), trigger, event->second.z, particle.eta(), event->second.cent, pid);
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaHi, processMCUnfolding, "Also store independent generated-INEL>0/reconstructed-FT0M unfolding inputs (requires McCentFT0Ms)", false);

  // void processTrackEfficiencyGeneral(
  //   MyCollisionsMcLabels::iterator const& collision,
  //   const aod::McCollisions& /*mcCollisions*/,
  //   const Particles& particles,
  //   const FiLTracks& /*tracks*/)
  // {
  //   if (!collision.sel8()) {
  //     return;
  //   }
  //   if (!collision.has_mcCollision()) {
  //     return;
  //   }

  //   auto mcCollision = collision.mcCollision();
  //   [[maybe_unused]] auto particlesPerCol = particles.sliceBy(perMCCol, mcCollision.globalIndex());
  // }
  // PROCESS_SWITCH(DndetaHi, processTrackEfficiencyGeneral, "MC Count tracks", false);

  void processCounting(
    const MyCollisionsCent& collisions,
    const FiDataTracks& tracks)
  {
    runCounting(collisions, tracks);
  }
  PROCESS_SWITCH(DndetaHi, processCounting, "Count tracks with Centrality", false);

  void processGen(
    const aod::McCollisions::iterator& mcCollision,
    const Particles& mcParticles)
  {
    auto genz = mcCollision.posZ();
    BoolVector bevtc(EvtEnd, false);
    bevtc[EvtInel] = true;
    for (const auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= CutCharge) {
          if (std::abs(particle.eta()) < CutEta1) {
            bevtc[EvtInelg0] = true;
          }
        }
      }
    }
    for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
      if (bevtc[ievtc]) {
        registry.fill(HIST("hgenzvtx"), ievtc, genz, NoCentrality);
      }
    }
    int pid = 0;
    for (const auto& particle : mcParticles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (std::abs(particle.pdgCode()) == kK0Short && std::abs(particle.eta()) < CutEta05 && std::abs(genz) < CutZ) {
        registry.fill(HIST("Selection"), SelK0sz10Eta05);
      }
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      switch (std::abs(particle.pdgCode())) {
        case kPiPlus:
          pid = PartPion;
          break;
        case kKPlus:
          pid = PartKaon;
          break;
        case kProton:
          pid = PartProton;
          break;
        default:
          pid = PartOther;
          break;
      }

      if (p != nullptr) {
        if (std::abs(p->Charge()) >= CutCharge) {
          for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
            if (bevtc[ievtc]) {
              registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), NoCentrality, pid, PtVarNone, particle.phi());
              if (particle.pt() < CutPt) {
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), NoCentrality, pid, PtVarUp, particle.phi(), PtWeightInterceptUp);
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), NoCentrality, pid, PtVarDown, particle.phi(), PtWeightInterceptDown);
              } else {
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), NoCentrality, pid, PtVarUp, particle.phi());
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), NoCentrality, pid, PtVarDown, particle.phi());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(DndetaHi, processGen, "Process generator-level info", false);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;

  void processV0Counting(
    const MyCollisionsCent& collisions,
    const aod::V0Datas& fullV0s,
    const FiTracks& /*tracks*/,
    const DaughterTracks& /*dauTrks*/)
  {
    for (const auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      auto z = collision.posZ();
      auto cent = -1.f;
      if (isPbPb) {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>()) {
          cent = collision.centFT0C();
        }
      } else {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>()) {
          cent = collision.centFT0M();
        }
      }

      auto v0sPerColl = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      for (const auto& v0 : v0sPerColl) {
        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();
        if (std::abs(z) < CutZ) {
          if (v0.v0radius() > v0Radius) {
            continue;
          }
          if (v0.dcapostopv() > dcaPosToPV) {
            continue;
          }
          if (v0.dcanegtopv() > dcaNegToPV) {
            continue;
          }
          if (v0.v0cosPA() < v0CosPA) {
            continue;
          }
          if (std::fabs(pTrack.eta()) > CutEta09) {
            continue;
          }
          if (std::fabs(nTrack.eta()) > CutEta09) {
            continue;
          }
          if (std::fabs(v0.eta()) < CutEta05) {
            registry.fill(HIST("hv0k0s"), v0.mK0Short());
          }
          registry.fill(HIST("hv0mass"), cent, SpecK0short, v0.eta(), v0.mK0Short());
          registry.fill(HIST("hv0mass"), cent, SpecLambda, v0.eta(), v0.mLambda());
          registry.fill(HIST("hv0mass"), cent, SpecAntilambda, v0.eta(), v0.mAntiLambda());
        }
      }
    }
  }
  PROCESS_SWITCH(DndetaHi, processV0Counting, "MC V0Count tracks without Centrality", false);

  // void processMCV0Counting(
  //   const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
  //   const aod::McCollisions& /*mcCollisions*/,
  //   const Particles& /*mcParticles*/,
  //   const soa::Filtered<aod::V0Datas>& /*fullV0s*/,
  //   const soa::Filtered<LabeledTracksEx>& /*tracks*/,
  //   const DaughterTracks& /*dauTrks*/)
  // {
  //   for (const auto& collision : collisions) {
  //     if (useEvSel && !collision.sel8()) {
  //       continue;
  //     }
  //     if (!collision.has_mcCollision()) {
  //       continue;
  //     }
  //   }
  // }
  // PROCESS_SWITCH(DndetaHi, processMCV0Counting, "MC V0Count tracks without Centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaHi>(cfgc)};
}
