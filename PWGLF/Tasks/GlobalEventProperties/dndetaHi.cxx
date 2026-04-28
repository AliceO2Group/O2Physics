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
#include "PWGMM/Mult/DataModel/Index.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/InitContext.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TH1.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults>;
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using DaughterTrack = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, DaughterTrack>;
using FiTracks = soa::Filtered<ExTracks>;
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
  Preslice<FiTracks> perCol = aod::track::collisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 2.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isPbPb{"isPbPb", false, "Is Pb-Pb"};

  Configurable<float> dcaV0Dau{"dcaV0Dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcaNegToPV{"dcaNegToPV", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPosToPV{"dcaPosToPV", 0.06, "DCA Pos To PV"};
  Configurable<double> v0CosPA{"v0CosPA", 0.97, "V0 CosPA"};
  Configurable<float> v0Radius{"v0Radius", 0.5, "Radius"};
  Configurable<float> etaDau{"etaDau", 4, "Eta Daughters"};
  Configurable<float> v0Rapidity{"v0Rapidity", 0.5, "V0 rapidity"};
  ConfigurableAxis centrality{"centrality", {VARIABLE_WIDTH, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100}, ""};
  Configurable<bool> zvtxcut{"zvtxcut", false, "z vtx cut < 10cm"};

  HistogramRegistry registry{
    "registry",
    {{"Selection", ";status;events", {HistType::kTH1F, {{17, 0.5, 17.5}}}}}};

  AxisSpec zAxis = {{-30, -20, -15, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30}, "Z (cm)", "zaxis"};
  AxisSpec deltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
  AxisSpec dcaAxis = {601, -3.01, 3.01, "", "DCA axis"};
  AxisSpec etaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
  AxisSpec v0EtaAxis = {20, -1.0, 1.0, "#etav0", "eta axis"};
  AxisSpec phiAxis = {629, 0, 2 * o2::constants::math::PI, "Rad", "phi axis"};
  AxisSpec ptVarAxis = {PtVarEnd - 1, PtVarBegin + 0.5, PtVarEnd - 0.5, "", "ptvar axis"};
  AxisSpec evtClassAxis = {EvtEnd - 1, EvtBegin + 0.5, EvtEnd - 0.5, "", "event class"};
  AxisSpec trigClassAxis = {TrigEnd - 1, TrigBegin + 0.5, TrigEnd - 0.5, "", "trigger class"};
  AxisSpec particleTypeAxis = {PartEnd - 1, PartBegin + 0.5, PartEnd - 0.5, "", "Particle type"};
  AxisSpec speciesAxis = {SpecEnd - 1, SpecBegin + 0.5, SpecEnd - 0.5, "", "species class"};
  AxisSpec massAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};
  AxisSpec signAxis = {SignEnd - 1, SignBegin + 0.5, SignEnd - 0.5, "", "sign"};
  AxisSpec stepAxis = {StepEnd - 1, StepBegin + 0.5, StepEnd - 0.5, "", "step"};
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
    nabs(aod::v0data::dcapostopv) > dcaPosToPV &&
    nabs(aod::v0data::dcanegtopv) > dcaNegToPV &&
    aod::v0data::dcaV0daughters < dcaV0Dau;

  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Partition<FiTracks> tSample3 = nabs(aod::track::eta) < estimatorEta;

  void init(InitContext&)
  {
    AxisSpec centAxis = {centrality, "", "centrality"};
    AxisSpec centAxisPbPb = {centrality, "", "centrality"};
    registry.add({"hetaresponse", ";etaresponse", {HistType::kTH2D, {{80, -4, 4}, {80, -4, 4}}}});
    registry.add({"hft0multiplicity", ";multiplicity", {HistType::kTH1D, {{10000, 0, 100000}}}});
    registry.add({"hcentrality", isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH1F, {{10000, 0, 100}}}});
    registry.add({"hcentralityvscentraldndeta",
                  isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M",
                  {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}}});
    registry.add({"hrecdndeta", "evntclass; triggerclass; zvtex, eta",
                  {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, phiBin}}});
    registry.add({"hreczvtx", "evntclass; triggerclass;  zvtex",
                  {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hphieta", "; #varphi; #eta; tracks",
                  {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, phiAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hrecdndetamissing", "evntclass; triggerclass; zvtex, eta",
                  {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hgendndeta", "evntclass;  zvtex, eta",
                  {HistType::kTHnSparseD, {evtClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, ptVarAxis, phiBin}}});
    registry.add({"hgenzvtx", "evntclass; zvtex",
                  {HistType::kTHnSparseD, {evtClassAxis, zAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hv0mass", "etaaxis; invmass",
                  {HistType::kTHnSparseD, {isPbPb ? centAxisPbPb : centAxis, speciesAxis, v0EtaAxis, massAxis}}});
    registry.add({"hv0k0s", "invmass", {HistType::kTH1D, {{100, 0.4, 0.6}}}});

    registry.add({"recetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0Sel8genz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});
    registry.add({"genetaINELg0genz10", ";etaresponse", {HistType::kTH2D, {etaAxis, zAxis}}});

    registry.add({"reczINELg0Sel8", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"genzINELg0Sel8", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"genzINELg0", ";z", {HistType::kTH1D, {zAxis}}});
    registry.add({"hcentmult", ";status;events", {HistType::kTH1D, {{100, 0, 100}}}});

    constexpr int NLogBins = 50;
    constexpr double PtLogLow = 0.01;
    constexpr double PtLogHigh = 10.0;
    std::vector<double> logBins(NLogBins + 1, 0.0);
    double logBinWidth = (std::log(PtLogHigh) - std::log(PtLogLow)) / NLogBins;
    for (int i = 0; i <= NLogBins; i++) {
      logBins[i] = PtLogLow * std::exp(i * logBinWidth);
    }
    AxisSpec ptBins = {logBins, "pT (GeV/c)", "pt bin"};

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
    const FiTracks& /*tracks*/)
  {
    for (const auto& collision : collisions) {
      BoolVector btrigc(TrigEnd, false);
      registry.fill(HIST("Selection"), SelAll);
      auto z = collision.posZ();
      auto pertracks = tSample3->sliceBy(perCol, collision.globalIndex());
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
  using MCex = soa::Join<aod::McCollisions, aod::HepMCXSections>;

  void processMCCounting(
    const MCex& mcCollisions,
    const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
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
      auto particles = mcParticles.sliceBy(perMCCol, mcCollision.globalIndex());

      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= CutCharge) {
            if (std::abs(particle.eta()) < CutEta1) {
              nTrkGen++;
            }
          }
        }
      }
      if (nTrkGen > 0) {
        registry.fill(HIST("Selection"), SelGeneratedInelGt0);
        bevtc[EvtInelg0] = true;
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
      auto centMC = -1.f;
      if (collisionsample.size() != 1) {
        centMC = -1.0;
      } else {
        for (const auto& collision : collisionsample) {
          if (isPbPb) {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>()) {
              centMC = collision.centFT0C();
            }
          } else {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>()) {
              centMC = collision.centFT0M();
            }
          }
        }
      }

      for (auto ievtc = 1u; ievtc < static_cast<unsigned>(EvtEnd); ievtc++) {
        if (bevtc[ievtc]) {
          registry.fill(HIST("hgenzvtx"), ievtc, genz, centMC);
        }
      }

      int pid = 0;
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
                registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centMC, pid, PtVarNone, particle.phi());
                if (particle.pt() < CutPt) {
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centMC, pid, PtVarUp, particle.phi(), PtWeightSlopeUp * particle.pt() + PtWeightInterceptUp);
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centMC, pid, PtVarDown, particle.phi(), PtWeightSlopeDown * particle.pt() + PtWeightInterceptDown);
                } else {
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centMC, pid, PtVarUp, particle.phi());
                  registry.fill(HIST("hgendndeta"), ievtc, genz, particle.eta(), centMC, pid, PtVarDown, particle.phi());
                }
              }
            }
            if (pid >= PartPion && pid <= PartOther) {
              particleetas.push_back(particle.eta());
            }
          }
        }
      }

      for (const auto& collision : collisionsample) {
        auto centRec = -1.f;
        if (isPbPb) {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>()) {
            centRec = collision.centFT0C();
          }
        } else {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>()) {
            centRec = collision.centFT0M();
          }
        }

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
        }
        if (collisionsample.size() == 1 && bevtc[EvtInelg0] && btrigc[TrigSel8]) {
          for (const auto& eta : particleetas) {
            registry.fill(HIST("genetaINELg0Sel8recz10"), eta, z);
            registry.fill(HIST("genetaINELg0Sel8genz10"), eta, mcz);
          }
          registry.fill(HIST("reczINELg0Sel8"), z);
          registry.fill(HIST("genzINELg0Sel8"), genz);
        }
        if (collisionsample.size() == 1 && bevtc[EvtInelg0]) {
          for (const auto& eta : particleetas) {
            registry.fill(HIST("genetaINELg0genz10"), eta, mcz);
          }
          registry.fill(HIST("genzINELg0"), genz);
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

  PROCESS_SWITCH(DndetaHi, processMCCounting, "MC Count tracks", false);

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
    const FiTracks& tracks)
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
