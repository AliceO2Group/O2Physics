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
#include <TMath.h>
#include <TPDGCode.h>

#include <RtypesCore.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <vector>

#include <math.h>

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
// using FiTracks = ExTracks;
// using Particles = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
using Particles = soa::Join<aod::McParticles, aod::ParticlesToTracks>;
// using Particle = Particles;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiLTracks = soa::Filtered<LabeledTracksEx>;
using DaughterTracks = soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using Bool1d = std::vector<bool>;

// Alias for const-ref-in-process linter issue
using MyCollisionsMcLabels = soa::Join<MyCollisions, aod::McCollisionLabels>;

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kINELg0,
  kDD,
  kSD,
  kND,
  kECend
};
enum {
  kTrigbegin = 0,
  kSel8 = 1,
  kSel8g0,
  kTrigend
};
enum {
  kSpeciesbegin = 0,
  kK0short = 1,
  kLambda,
  kAntilambda,
  kSpeciesend
};
enum {
  kSignbegin = 0,
  kPositive = 1,
  kNegative,
  kSignend
};
enum {
  kStepbegin = 0,
  kAll = 1,
  kBasiccut,
  kMasscut,
  kStepend
};

enum {
  kParTypebegin = 0,
  kParDATA = 1,
  kMotherStrange,
  kBkg,
  kNotPrimary,
  kPion,
  kKaon,
  kProtonMy,
  kOPar,
  kParTypeend
};
enum {
  kPtVarbegin = 0,
  kNoPtVar = 1,
  kPtUp,
  kPtDw,
  kPtVarend
};

namespace
{
template <typename T>
static constexpr bool hasCent()
{
  if constexpr (!soa::is_soa_join_v<T>) {
    return false;
  } else if (T::template contains<aod::HepMCHeavyIons>()) {
    return true;
  } else {
    return false;
  }
}

// UpperCamelCase for all constexpr constants to pass linter
constexpr float CutZ = 10.0f;
constexpr float CutEta1 = 1.0f;
constexpr float CutEta08 = 0.8f;
constexpr float CutEta05 = 0.5f;
constexpr float CutEta09 = 0.9f;
constexpr int CutCharge = 3;
constexpr float CutPt = 0.1f;

// Process IDs magic numbers
constexpr int ProcIdND = 101;
constexpr int ProcIdSD1 = 103;
constexpr int ProcIdSD2 = 104;
constexpr int ProcIdDD1 = 105;
constexpr int ProcIdDD2 = 106;

// Weights and Fill values magic numbers
constexpr float W1 = -10.0f;
constexpr float W2 = 2.0f;
constexpr float W3 = 5.0f;
constexpr float W4 = 0.5f;
constexpr float Fill1 = 1.0f;
constexpr float Fill2 = 2.0f;
constexpr float Fill3 = 3.0f;
constexpr float Fill4 = 4.0f;
constexpr float Fill5 = 5.0f;
constexpr float Fill6 = 6.0f;
constexpr float Fill7 = 7.0f;
constexpr float Fill9 = 9.0f;
constexpr float Fill10 = 10.0f;
constexpr float Fill11 = 11.0f;
constexpr float Fill12 = 12.0f;
constexpr float Fill13 = 13.0f;
constexpr float Fill14 = 14.0f;
constexpr float Fill15 = 15.0f;
constexpr float Fill16 = 16.0f;
constexpr float Fill17 = 17.0f;
constexpr float FillM1 = -1.0f;
} // namespace

AxisSpec zAxis = {{-30, -20, -15, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30}, "Z (cm)", "zaxis"};
AxisSpec deltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
AxisSpec dcaAxis = {601, -3.01, 3.01, "", "DCA axis"};
AxisSpec etaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
AxisSpec v0EtaAxis = {20, -1.0, 1.0, "#etav0", "eta axis"};
AxisSpec phiAxis = {629, 0, 2 * o2::constants::math::PI, "Rad", "phi axis"};
AxisSpec ptVarAxis = {kPtVarend - 1, +kPtVarbegin + 0.5, +kPtVarend - 0.5, "", "ptvar axis"};
AxisSpec evtClassAxis = {kECend - 1, +kECbegin + 0.5, +kECend - 0.5, "", "event class"};
AxisSpec trigClassAxis = {kTrigend - 1, +kTrigbegin + 0.5, +kTrigend - 0.5, "", "trigger class"};
AxisSpec particleTypeAxis = {kParTypeend - 1, +kParTypebegin + 0.5, +kParTypeend - 0.5, "", "Particle type"};
AxisSpec speciesAxis = {kSpeciesend - 1, +kSpeciesbegin + 0.5, +kSpeciesend - 0.5, "", "species class"};
AxisSpec massAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};
AxisSpec signAxis = {kSignend - 1, +kSignbegin + 0.5, +kSignend - 0.5, "", "sign"};
AxisSpec stepAxis = {kStepend - 1, +kStepbegin + 0.5, +kStepend - 0.5, "", "step"};
AxisSpec testAxis = {101, -0.5, 100.5, "", "test"};
AxisSpec multAxis = {1001, -0.5, 1000.5, "", "Ntrks"};
AxisSpec statusCodeAxis = {3, -1.5, 2.5, "", "StatusCode"};
AxisSpec processCodeAxis = {45, -1.5, 44.5, "", "StatusCode"};

AxisSpec phibin = {{0, o2::constants::math::PI / 2, o2::constants::math::PI, o2::constants::math::PI * 3. / 2, 2 * o2::constants::math::PI}, "#phi", "phi bin"};

static constexpr TrackSelectionFlags::flagtype TrackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype TrackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype TrackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

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

  // Configurable<bool> mftanalysis{"mftanalysis", false, "mft analysis switch"};
  Configurable<bool> zvtxcut{"zvtxcut", false, "z vtx cut < 10cm"};

  HistogramRegistry registry{
    "registry",
    {{"Selection", ";status;events", {HistType::kTH1F, {{17, 0.5, 17.5}}}}

    }};

  std::vector<int> usedTracksIds;
  void init(InitContext&)
  {

    AxisSpec centAxis = {centrality, "", "centrality"};
    AxisSpec centAxisPbPb = {centrality, "", "centrality"};
    registry.add({"hetaresponse", ";etaresponse", {HistType::kTH2D, {{80, -4, 4}, {80, -4, 4}}}});
    registry.add({"hft0multiplicity", ";multiplicity", {HistType::kTH1D, {{10000, 0, 100000}}}});
    registry.add({"hcentrality", isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH1F, {{10000, 0, 100}}}});
    registry.add({"hcentralityvscentraldndeta", isPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH2F, {
                                                                                                                              {100, 0, 100},
                                                                                                                              {100, 0, 100},
                                                                                                                            }}});
    registry.add({"hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, phibin}}});
    registry.add({"hreczvtx", "evntclass; triggerclass;  zvtex", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hphieta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, phiAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hrecdndetamissing", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, trigClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis}}});
    registry.add({"hgendndeta", "evntclass;  zvtex, eta", {HistType::kTHnSparseD, {evtClassAxis, zAxis, etaAxis, isPbPb ? centAxisPbPb : centAxis, particleTypeAxis, ptVarAxis, phibin}}});
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
    const int nbins = 50;
    std::vector<double> logbins(nbins + 1, 0);
    double low = 0.01;
    double high = 10;
    double logbw = (std::log(high) - std::log(low)) / nbins;
    for (int ij = 0; ij <= nbins; ij++) {
      logbins[ij] = low * std::exp(ij * logbw);
    }
    AxisSpec ptbins2 = {logbins, "pT (GeV/c)", "pt bin"};

    registry.add({"hrecdndpt", " pt", {HistType::kTH1D, {ptbins2}}});
    registry.add({"hdndptefficiency", " pt", {HistType::kTH1D, {ptbins2}}});
    registry.add({"hgendndpt", " pt", {HistType::kTH1D, {{ptbins2}}}});
    registry.add({"hgendndpt2", " pt", {HistType::kTH1D, {{ptbins2}}}});
    registry.add({"hgendndpt05", " pt", {HistType::kTH1D, {{ptbins2}}}});

    auto hstat = registry.get<TH1>(HIST("Selection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Sel8");
    x->SetBinLabel(3, "Sel8z10");
    x->SetBinLabel(4, "Generated INEL>0");
    x->SetBinLabel(5, "Good BCs");
    x->SetBinLabel(6, "BCs with collisions");
    x->SetBinLabel(7, "BCs with pile-up/splitting");
    x->SetBinLabel(8, "INEL&Sel8&mcz10");
    x->SetBinLabel(9, "INEL");
    x->SetBinLabel(10, "INELg0");
    x->SetBinLabel(11, "INELgenz10");
    x->SetBinLabel(12, "INELg0genz10");
    x->SetBinLabel(13, "INELSel8");
    x->SetBinLabel(14, "INELg0Sel8g0");
    x->SetBinLabel(15, "INELSel8recz10");
    x->SetBinLabel(16, "INELg0Sel8g0recz10");
    x->SetBinLabel(17, "K0Sz10eta05");
  }

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & TrackSelectionITS) == TrackSelectionITS) && ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC), ncheckbit(aod::track::trackCutFlag, TrackSelectionTPC), true) && ncheckbit(aod::track::trackCutFlag, TrackSelectionDCA);
  expressions::Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcaPosToPV&& nabs(aod::v0data::dcanegtopv) > dcaNegToPV&& aod::v0data::dcaV0daughters < dcaV0Dau;

  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Partition<FiTracks> tSample3 = nabs(aod::track::eta) < estimatorEta;
  void processEventStat(
    const FullBCs& bcs,
    const soa::Join<aod::Collisions, aod::EvSels>& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (const auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(o2::aod::evsel::kIsBBT0A) &&
                        bc.selection_bit(o2::aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("Selection"), Fill5);
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
          registry.fill(HIST("Selection"), Fill6);
          if (cols.size() > 1) {
            registry.fill(HIST("Selection"), Fill7);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaHi, processEventStat, "Collect event sample stats", false);

  std::vector<double> tracketas;
  template <typename C>
  void runCounting(
    const C& collisions,
    const FiTracks& /*tracks*/)
  {

    for (const auto& collision : collisions) {
      Bool1d btrigc(kTrigend, false);
      registry.fill(HIST("Selection"), Fill1);
      auto z = collision.posZ();
      auto pertracks = tSample3->sliceBy(perCol, collision.globalIndex());
      auto nTrk = 0;

      // if (collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      // if (collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      if (collision.sel8()) {
        btrigc[kSel8] = true;
        registry.fill(HIST("Selection"), Fill2);
        if (std::abs(z) < CutZ) {
          registry.fill(HIST("Selection"), Fill3);
        }
      }
      if (btrigc[kSel8] && std::abs(z) < CutZ)
        registry.fill(HIST("hft0multiplicity"), collision.multFT0C());

      for (const auto& track : pertracks) {
        [[maybe_unused]] int dummy = track.globalIndex();
        if (std::abs(track.eta()) < CutEta1)
          nTrk++; // charged track check
      }
      if (nTrk > 0) {
        // registry.fill(HIST("Selection"), 3.);
        if (btrigc[kSel8])
          btrigc[kSel8g0] = true;
      }
      if (btrigc[kSel8g0])
        registry.fill(HIST("reczINELg0Sel8"), z);

      auto cent = -1.f;
      if (isPbPb) {
        if constexpr (C::template contains<aod::CentFT0Cs>())
          cent = collision.centFT0C();

      } else {
        if constexpr (C::template contains<aod::CentFT0Ms>())
          cent = collision.centFT0M();
      }

      if (isPbPb) {
        if (std::abs(z) < CutZ && btrigc[kSel8])
          registry.fill(HIST("hcentrality"), cent);
      } else {
        if (std::abs(z) < CutZ && btrigc[kSel8g0])
          registry.fill(HIST("hcentrality"), cent);
        if (std::abs(z) < CutZ && btrigc[kSel8g0])
          registry.fill(HIST("hcentralityvscentraldndeta"), cent, nTrk);
      }
      for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
        if (btrigc[itrigc])
          registry.fill(HIST("hreczvtx"), double(kDATA), double(itrigc), z, cent);
      }

      for (const auto& track : pertracks) {
        if (btrigc[kSel8] && std::abs(track.eta()) < CutEta08 && std::abs(z) < CutZ)
          registry.fill(HIST("hrecdndpt"), track.pt());
        if (btrigc[kSel8])
          registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);

        for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
          if (btrigc[itrigc]) {
            registry.fill(HIST("hphieta"), double(kDATA), double(itrigc), track.phi(), track.eta(), cent);
            registry.fill(HIST("hrecdndeta"), double(kDATA), double(itrigc), z, track.eta(), cent, double(kParDATA), track.phi());
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
    const MCex& mcCollisions, const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions, const Particles& mcParticles,
    const FiLTracks& tracks)
  {
    for (const auto& mcCollision : mcCollisions) {
      Bool1d bevtc(kECend, false);
      bevtc[kINEL] = true;
      auto procId = mcCollision.processId();
      if (procId == ProcIdND) {
        bevtc[kND] = true;
      } else if (procId == ProcIdSD1 || procId == ProcIdSD2) {
        bevtc[kSD] = true;
      } else if (procId == ProcIdDD1 || procId == ProcIdDD2) {
        bevtc[kDD] = true;
      }
      registry.fill(HIST("Selection"), Fill1);

      auto mcz = mcCollision.posZ();
      auto genz = mcz;

      auto nTrkGen = 0;
      auto particles = mcParticles.sliceBy(perMCCol, mcCollision.globalIndex());
      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
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
        registry.fill(HIST("Selection"), Fill4);
        bevtc[kINELg0] = true;
      }
      if (bevtc[kINEL])
        registry.fill(HIST("Selection"), Fill9);
      if (bevtc[kINEL] && std::abs(mcz) < CutZ)
        registry.fill(HIST("Selection"), Fill11);
      if (bevtc[kINELg0])
        registry.fill(HIST("Selection"), Fill10);
      if (bevtc[kINELg0] && std::abs(mcz) < CutZ)
        registry.fill(HIST("Selection"), Fill12);
      for (const auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= CutCharge) {
            if (bevtc[kINEL] && std::abs(particle.eta()) < CutEta08 && std::abs(mcz) < CutZ) {
              registry.fill(HIST("hgendndpt"), particle.pt());
              if (particle.pt() < CutPt) {
                registry.fill(HIST("hgendndpt2"), particle.pt(), W1 * particle.pt() + W2);
                registry.fill(HIST("hgendndpt05"), particle.pt(), W3 * particle.pt() + W4);
              } else {
                registry.fill(HIST("hgendndpt2"), particle.pt());
                registry.fill(HIST("hgendndpt05"), particle.pt());
              }
            }
          }
        }
      }

      auto collisionsample = collisions.sliceBy(perMcCol, mcCollision.globalIndex());
      auto cent = -1.f;

      if (collisionsample.size() != 1) {
        cent = -1.0;
      } else {
        for (const auto& collision : collisionsample) {
          if (isPbPb) {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
              cent = collision.centFT0C();
          } else {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
              cent = collision.centFT0M();

            // auto nTrkRec = 0;
            // auto trackspart = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            // for (const auto& track : trackspart) {
            //   if (std::abs(track.eta()) < 1) {
            //     nTrkRec++;
            //   }
            // }
          }
        }
      }
      for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
        if (bevtc[ievtc])
          registry.fill(HIST("hgenzvtx"), double(ievtc), genz, cent);
      }
      int pid = 0;
      std::vector<double> particleetas;
      for (const auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (std::abs(particle.pdgCode()) == kK0Short && std::abs(particle.eta()) < CutEta05 && std::abs(genz) < CutZ)
          registry.fill(HIST("Selection"), Fill17);
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        switch (std::abs(particle.pdgCode())) {
          case 211:
            pid = kPion;
            break;
          case 321:
            pid = kKaon;
            break;
          case 2212:
            pid = kProtonMy;
            break;
          default:
            pid = kOPar;
            break;
        }

        if (p != nullptr) {
          if (std::abs(p->Charge()) >= CutCharge) {
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              if (bevtc[ievtc]) {
                registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), cent, double(pid), kNoPtVar, particle.phi());
                if (particle.pt() < CutPt) {
                  registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), cent, double(pid), kPtUp, particle.phi(), W1 * particle.pt() + W2);
                  registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), cent, double(pid), kPtDw, particle.phi(), W3 * particle.pt() + W4);

                } else {
                  registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), cent, double(pid), kPtUp, particle.phi());
                  registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), cent, double(pid), kPtDw, particle.phi());
                }
              }
            }
            if (pid >= kPion && pid <= kOPar)
              particleetas.push_back(particle.eta());
          }
        }
      }

      for (const auto& collision : collisionsample) {
        auto cent = -1.f;
        if (isPbPb) {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
            cent = collision.centFT0C();
        } else {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
            cent = collision.centFT0M();
          // auto nTrkRec = 0;
          // auto trackspart = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
          // for (const auto& track : trackspart) {
          //   if (std::abs(track.eta()) < 1) {
          //     nTrkRec++;
          //   }
          // }
        }

        Bool1d btrigc(kTrigend, false);
        auto z = collision.posZ();
        // if (collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        // if (collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        if (collision.sel8()) {
          btrigc[kSel8] = true;
          registry.fill(HIST("Selection"), Fill2);
          if (std::abs(z) < CutZ) {
            registry.fill(HIST("Selection"), Fill3);
          }
        }
        if (bevtc[kINEL] && btrigc[kSel8] && std::abs(z) < CutZ)
          registry.fill(HIST("hft0multiplicity"), collision.multFT0C());
        if (collisionsample.size() == 1 && bevtc[kINELg0] && btrigc[kSel8]) {
          for (const auto& eta : particleetas) {
            registry.fill(HIST("genetaINELg0Sel8recz10"), eta, z);
            registry.fill(HIST("genetaINELg0Sel8genz10"), eta, mcz);
          }
          registry.fill(HIST("reczINELg0Sel8"), z);
          registry.fill(HIST("genzINELg0Sel8"), genz);
        }
        if (collisionsample.size() == 1 && bevtc[kINELg0]) {
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
          if (btrigc[kSel8])
            btrigc[kSel8g0] = true;
        }
        if (isPbPb) {
          if (std::abs(z) < CutZ && btrigc[kSel8])
            registry.fill(HIST("hcentrality"), cent);
        } else {
          if (std::abs(z) < CutZ && btrigc[kSel8g0])
            registry.fill(HIST("hcentrality"), cent);
        }

        if (bevtc[kINEL] && btrigc[kSel8])
          registry.fill(HIST("Selection"), Fill13);
        if (bevtc[kINEL] && btrigc[kSel8] && std::abs(z) < CutZ)
          registry.fill(HIST("Selection"), Fill15);
        if (bevtc[kINELg0] && btrigc[kSel8g0])
          registry.fill(HIST("Selection"), Fill14);
        if (bevtc[kINELg0] && btrigc[kSel8g0] && std::abs(z) < CutZ)
          registry.fill(HIST("Selection"), Fill16);

        for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
          for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
            if (bevtc[ievtc] && btrigc[itrigc]) {
              registry.fill(HIST("hreczvtx"), double(ievtc), double(itrigc), z, cent);
            }
          }
        }
        std::vector<int> mclabels;
        for (const auto& track : trackspart) {
          if (track.has_mcParticle()) {
            int pid = kBkg;
            auto particle = track.template mcParticle_as<Particles>();
            if (particle.isPhysicalPrimary()) {
              switch (std::abs(particle.pdgCode())) {
                case 211:
                  pid = kPion;
                  break;
                case 321:
                  pid = kKaon;
                  break;
                case 2212:
                  pid = kProtonMy;
                  break;
                default:
                  pid = kOPar;
                  break;
              }
            } else {
              pid = kNotPrimary;
            }
            for (auto motherIds = particle.mothersIds().front(); motherIds <= particle.mothersIds().back(); motherIds++) {
              auto mother = mcParticles.rawIteratorAt(motherIds);
              auto pdg_mother = mother.pdgCode();
              if (pdg_mother == kK0Short || std::abs(pdg_mother) == kLambda0) {
                pid = kMotherStrange;
              }
            }
            if (find(mclabels.begin(), mclabels.end(), track.mcParticleId()) != mclabels.end())
              pid = kBkg;
            mclabels.push_back(track.mcParticleId());
            registry.fill(HIST("hetaresponse"), particle.eta(), track.eta(), cent);
            if (bevtc[kINEL] && btrigc[kSel8] && std::abs(track.eta()) < CutEta08 && std::abs(z) < CutZ && pid != kBkg && pid != kNotPrimary)
              registry.fill(HIST("hdndptefficiency"), particle.pt());
            if (btrigc[kSel8] && bevtc[kINELg0])
              registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), double(ievtc), double(itrigc), z, particle.eta(), cent, double(pid), particle.phi());
                  registry.fill(HIST("hphieta"), double(ievtc), double(itrigc), track.phi(), track.eta(), cent);
                }
              }
            }

          } else {
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), double(ievtc), double(itrigc), z, track.eta(), cent, double(kBkg), track.phi());
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(DndetaHi, processMCCounting, "MC Count tracks", false);

  void processTrackEfficiencyGeneral(
    MyCollisionsMcLabels::iterator const& collision,
    const aod::McCollisions& /*mcCollisions*/, const Particles& particles,
    const FiLTracks& /*tracks*/)
  {

    if (!collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    auto mcCollision = collision.mcCollision();
    auto particlesPerCol = particles.sliceBy(perMCCol, mcCollision.globalIndex());
  }

  PROCESS_SWITCH(DndetaHi, processTrackEfficiencyGeneral, "MC Count tracks", false);

  void processCounting(
    const MyCollisionsCent& collisions,
    const FiTracks& tracks)
  {
    runCounting(collisions, tracks);
  }
  PROCESS_SWITCH(DndetaHi, processCounting, "Count tracks with Centrality", false);

  void processGen(
    const aod::McCollisions::iterator& mcCollision, const Particles& mcParticles)
  {

    auto genz = mcCollision.posZ();
    Bool1d bevtc(kECend, false);
    bevtc[kINEL] = true;
    for (const auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= CutCharge) {
          if (std::abs(particle.eta()) < CutEta1)
            bevtc[kINELg0] = true;
        }
      }
    }
    for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
      if (bevtc[ievtc])
        registry.fill(HIST("hgenzvtx"), double(ievtc), genz, FillM1);
    }
    int pid = 0;
    for (const auto& particle : mcParticles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (std::abs(particle.pdgCode()) == kK0Short && std::abs(particle.eta()) < CutEta05 && std::abs(genz) < CutZ)
        registry.fill(HIST("Selection"), Fill17);
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      switch (std::abs(particle.pdgCode())) {
        case 211:
          pid = kPion;
          break;
        case 321:
          pid = kKaon;
          break;
        case 2212:
          pid = kProtonMy;
          break;
        default:
          pid = kOPar;
          break;
      }

      if (p != nullptr) {
        if (std::abs(p->Charge()) >= CutCharge) {
          for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
            if (bevtc[ievtc]) {
              registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), FillM1, double(pid), kNoPtVar, particle.phi());
              if (particle.pt() < CutPt) {
                registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), FillM1, double(pid), kPtUp, particle.phi(), W2);
                registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), FillM1, double(pid), kPtDw, particle.phi(), W4);

              } else {
                registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), FillM1, double(pid), kPtUp, particle.phi());
                registry.fill(HIST("hgendndeta"), double(ievtc), genz, particle.eta(), FillM1, double(pid), kPtDw, particle.phi());
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
    const DaughterTracks& /*Dautrks*/)
  {
    for (const auto& collision : collisions) {
      if (!collision.sel8())
        continue;
      auto z = collision.posZ();

      auto cent = -1.f;
      if (isPbPb) {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
          cent = collision.centFT0C();
      } else {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
          cent = collision.centFT0M();
      }

      auto v0s_per_coll = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      for (const auto& v0 : v0s_per_coll) {

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();
        if (std::abs(z) < CutZ) {
          if (v0.v0radius() > v0Radius)
            continue;
          if (v0.dcapostopv() > dcaPosToPV)
            continue;
          if (v0.dcanegtopv() > dcaNegToPV)
            continue;
          if (v0.v0cosPA() < v0CosPA)
            continue;
          if (std::fabs(pTrack.eta()) > CutEta09)
            continue;
          if (std::fabs(nTrack.eta()) > CutEta09)
            continue;

          if (std::fabs(v0.eta()) < CutEta05)
            registry.fill(HIST("hv0k0s"), v0.mK0Short());
          registry.fill(HIST("hv0mass"), cent, double(kK0short), v0.eta(), double(v0.mK0Short()));
          registry.fill(HIST("hv0mass"), cent, double(kLambda), v0.eta(), double(v0.mLambda()));
          registry.fill(HIST("hv0mass"), cent, double(kAntilambda), v0.eta(), double(v0.mAntiLambda()));
        }
      }
    }
  }
  PROCESS_SWITCH(DndetaHi, processV0Counting, "MC V0Count tracks without Centrality", false);

  void processMCV0Counting(
    const soa::Join<MyCollisionsCent, aod::McCollisionLabels>& collisions,
    const aod::McCollisions& /*mcCollisions*/,
    const Particles& /*mcParticles*/,
    const soa::Filtered<aod::V0Datas>& fullV0s,
    const soa::Filtered<LabeledTracksEx>& /*tracks*/,
    const DaughterTracks& /*dauTrks*/)
  {
    for (const auto& collision : collisions) {
      auto cent = -1.f;

      if (isPbPb) {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
          cent = collision.centFT0C();
      } else {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
          cent = collision.centFT0M();
      }

      if (useEvSel && !collision.sel8()) // event selection cut
        continue;
      if (!collision.has_mcCollision()) // check mc particle
        continue;

      for (const auto& v0 : fullV0s) {

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();

        if (0 && v0.v0radius() > v0Radius &&
            v0.v0cosPA() > v0CosPA &&
            std::abs(pTrack.eta()) < etaDau &&
            std::abs(nTrack.eta()) < etaDau) {

          registry.fill(HIST("hv0mass"), cent, double(kK0short), v0.eta(), double(v0.mK0Short()));
          registry.fill(HIST("hv0mass"), cent, double(kLambda), v0.eta(), double(v0.mLambda()));
          registry.fill(HIST("hv0mass"), cent, double(kAntilambda), v0.eta(), double(v0.mAntiLambda()));
        }
      }
    }
  }
  PROCESS_SWITCH(DndetaHi, processMCV0Counting, "MC V0Count tracks without Centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DndetaHi>(cfgc)};
}
