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

#include "Index.h"
#include "bestCollisionTable.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

using namespace o2::aod;
using namespace o2::analysis;

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
typedef std::vector<Bool_t> Bool_1d;

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kINELg0,
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
} // namespace

AxisSpec ZAxis = {{-30, -20, -15, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30}, "Z (cm)", "zaxis"};
AxisSpec DeltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
AxisSpec DCAAxis = {601, -3.01, 3.01, "", "DCA axis"};
AxisSpec EtaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
AxisSpec V0EtaAxis = {20, -1.0, 1.0, "#etav0", "eta axis"};
AxisSpec PhiAxis = {629, 0, 2 * M_PI, "Rad", "phi axis"};
AxisSpec PtVarAxis = {kPtVarend - 1, +kPtVarbegin + 0.5, +kPtVarend - 0.5, "", "ptvar axis"};
AxisSpec EvtClassAxis = {kECend - 1, +kECbegin + 0.5, +kECend - 0.5, "", "event class"};
AxisSpec TrigClassAxis = {kTrigend - 1, +kTrigbegin + 0.5, +kTrigend - 0.5, "", "trigger class"};
AxisSpec ParticleTypeAxis = {kParTypeend - 1, +kParTypebegin + 0.5, +kParTypeend - 0.5, "", "Particle type"};
AxisSpec SpeciesAxis = {kSpeciesend - 1, +kSpeciesbegin + 0.5, +kSpeciesend - 0.5, "", "species class"};
AxisSpec MassAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};
AxisSpec SignAxis = {kSignend - 1, +kSignbegin + 0.5, +kSignend - 0.5, "", "sign"};
AxisSpec StepAxis = {kStepend - 1, +kStepbegin + 0.5, +kStepend - 0.5, "", "step"};
AxisSpec testAxis = {101, -0.5, 100.5, "", "test"};
AxisSpec multAxis = {1001, -0.5, 1000.5, "", "Ntrks"};
AxisSpec StatusCodeAxis = {3, -1.5, 2.5, "", "StatusCode"};
AxisSpec ProcessCodeAxis = {45, -1.5, 44.5, "", "StatusCode"};

auto pi = TMath::Pi();
AxisSpec phibin = {{0, pi / 2, pi, pi * 3. / 2, 2 * pi}, "#phi", "phi bin"};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

struct MultiplicityCounter {
  SliceCache cache;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 2.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> IsPbPb{"IsPbPb", false, "Is Pb-Pb"};

  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> v0radius{"v0radius", 0.5, "Radius"};
  Configurable<float> etadau{"etadau", 4, "Eta Daughters"};
  Configurable<float> rapidity{"v0rapidity", 0.5, "V0 rapidity"};
  ConfigurableAxis centBinning{"centrality", {VARIABLE_WIDTH, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100}, ""};

  // Configurable<bool> mftanalysis{"mftanalysis", false, "mft analysis switch"};
  Configurable<bool> zvtxcut{"zvtxcut", false, "z vtx cut < 10cm"};

  HistogramRegistry registry{
    "registry",
    {{"Selection", ";status;events", {HistType::kTH1F, {{17, 0.5, 17.5}}}}

    }};

  std::vector<int> usedTracksIds;
  void init(InitContext&)
  {

    AxisSpec CentAxis = {centBinning, "", "centrality"};
    AxisSpec CentAxisPbPb = {centBinning, "", "centrality"};
    registry.add({"hetaresponse", ";etaresponse", {HistType::kTH2D, {{80, -4, 4}, {80, -4, 4}}}});
    registry.add({"hft0multiplicity", ";multiplicity", {HistType::kTH1D, {{10000, 0, 100000}}}});
    registry.add({"hcentrality", IsPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH1F, {{10000, 0, 100}}}});
    registry.add({"hcentralityvscentraldndeta", IsPbPb ? " ; centrality_FT0C (%) " : "; centrality_FT0M", {HistType::kTH2F, {
                                                                                                                              {100, 0, 100},
                                                                                                                              {100, 0, 100},
                                                                                                                            }}});
    registry.add({"hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis, IsPbPb ? CentAxisPbPb : CentAxis, ParticleTypeAxis, phibin}}});
    registry.add({"hreczvtx", "evntclass; triggerclass;  zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, IsPbPb ? CentAxisPbPb : CentAxis}}});
    registry.add({"hphieta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, PhiAxis, EtaAxis, IsPbPb ? CentAxisPbPb : CentAxis}}});
    registry.add({"hrecdndetamissing", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis, IsPbPb ? CentAxisPbPb : CentAxis}}});
    registry.add({"hgendndeta", "evntclass;  zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis, EtaAxis, IsPbPb ? CentAxisPbPb : CentAxis, ParticleTypeAxis, PtVarAxis, phibin}}});
    registry.add({"hgenzvtx", "evntclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis, IsPbPb ? CentAxisPbPb : CentAxis}}});
    registry.add({"hv0mass", "etaaxis; invmass", {HistType::kTHnSparseD, {IsPbPb ? CentAxisPbPb : CentAxis, SpeciesAxis, V0EtaAxis, MassAxis}}});
    registry.add({"hv0k0s", "invmass", {HistType::kTH1D, {{100, 0.4, 0.6}}}});

    registry.add({"recetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {EtaAxis, ZAxis}}});
    registry.add({"genetaINELg0Sel8recz10", ";etaresponse", {HistType::kTH2D, {EtaAxis, ZAxis}}});
    registry.add({"genetaINELg0Sel8genz10", ";etaresponse", {HistType::kTH2D, {EtaAxis, ZAxis}}});
    registry.add({"genetaINELg0genz10", ";etaresponse", {HistType::kTH2D, {EtaAxis, ZAxis}}});

    registry.add({"reczINELg0Sel8", ";z", {HistType::kTH1D, {ZAxis}}});
    registry.add({"genzINELg0Sel8", ";z", {HistType::kTH1D, {ZAxis}}});
    registry.add({"genzINELg0", ";z", {HistType::kTH1D, {ZAxis}}});
    registry.add({"hcentmult", ";status;events", {HistType::kTH1D, {{100, 0, 100}}}});
    const int nbins = 50;
    std::vector<Double_t> logbins(nbins + 1, 0);
    Double_t low = 0.01;
    Double_t high = 10;
    Double_t logbw = (log(high) - log(low)) / nbins;
    for (int ij = 0; ij <= nbins; ij++) {
      logbins[ij] = low * exp(ij * logbw);
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

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) && ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC), ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true) && ncheckbit(aod::track::trackCutFlag, trackSelectionDCA);
  expressions::Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Partition<FiTracks> tSample3 = nabs(aod::track::eta) < estimatorEta;
  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(o2::aod::evsel::kIsBBT0A) &&
                        bc.selection_bit(o2::aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("Selection"), 5.);
        cols.clear();
        for (auto& collision : collisions) {
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
          registry.fill(HIST("Selection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("Selection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStat, "Collect event sample stats", false);

  std::vector<Double_t> tracketas;
  template <typename C>
  void runCounting(
    C const& collisions,
    FiTracks const& /*tracks*/)
  {

    for (auto& collision : collisions) {
      Bool_1d btrigc(kTrigend, false);
      registry.fill(HIST("Selection"), 1.);
      auto z = collision.posZ();
      auto pertracks = tSample3->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto Ntrk = 0;

      // if (collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      // if (collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      if (collision.sel8()) {
        btrigc[kSel8] = true;
        registry.fill(HIST("Selection"), 2.);
        if (std::abs(z) < 10) {
          registry.fill(HIST("Selection"), 3.);
        }
      }
      if (btrigc[kSel8] && std::abs(z) < 10)
        registry.fill(HIST("hft0multiplicity"), collision.multFT0C());

      for (auto& track : pertracks) {
        [[maybe_unused]] int dummy = track.globalIndex();
        if (std::abs(track.eta()) < 1)
          Ntrk++; // charged track check
      }
      if (Ntrk > 0) {
        // registry.fill(HIST("Selection"), 3.);
        if (btrigc[kSel8])
          btrigc[kSel8g0] = true;
      }
      if (btrigc[kSel8g0])
        registry.fill(HIST("reczINELg0Sel8"), z);

      auto cent = -1.f;
      if (IsPbPb) {
        if constexpr (C::template contains<aod::CentFT0Cs>())
          cent = collision.centFT0C();

      } else {
        if constexpr (C::template contains<aod::CentFT0Ms>())
          cent = collision.centFT0M();
      }

      if (IsPbPb) {
        if (std::abs(z) < 10 && btrigc[kSel8])
          registry.fill(HIST("hcentrality"), cent);
      } else {
        if (std::abs(z) < 10 && btrigc[kSel8g0])
          registry.fill(HIST("hcentrality"), cent);
        if (std::abs(z) < 10 && btrigc[kSel8g0])
          registry.fill(HIST("hcentralityvscentraldndeta"), cent, Ntrk);
      }
      for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
        if (btrigc[itrigc])
          registry.fill(HIST("hreczvtx"), Double_t(kDATA), Double_t(itrigc), z, cent);
      }

      for (auto& track : pertracks) {
        if (btrigc[kSel8] && std::abs(track.eta()) < 0.8 && std::abs(z) < 10)
          registry.fill(HIST("hrecdndpt"), track.pt());
        if (btrigc[kSel8])
          registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);

        for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
          if (btrigc[itrigc]) {
            registry.fill(HIST("hphieta"), Double_t(kDATA), Double_t(itrigc), track.phi(), track.eta(), cent);
            registry.fill(HIST("hrecdndeta"), Double_t(kDATA), Double_t(itrigc), z, track.eta(), cent, Double_t(kParDATA), track.phi());
          }
        }
      }
    }
  }

  PresliceUnsorted<soa::Join<MyCollisions, aod::McCollisionLabels>> perMcCol = o2::aod::mccollisionlabel::mcCollisionId;
  Preslice<aod::McParticles> perMCColparticles = aod::mcparticle::mcCollisionId;
  void processMCCounting(
    aod::McCollisions const& mcCollisions, soa::Join<MyCollisionsCent, aod::McCollisionLabels> const& collisions, Particles const& mcParticles,
    FiLTracks const& tracks)
  {
    for (auto& mcCollision : mcCollisions) {
      Bool_1d bevtc(kECend, false);
      bevtc[kINEL] = true;
      registry.fill(HIST("Selection"), 1.);

      auto mcz = mcCollision.posZ();
      auto genz = mcz;

      auto Ntrk_gen = 0;
      auto particles = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      for (auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= 3) {
            if (std::abs(particle.eta()) < 1) {
              Ntrk_gen++;
            }
          }
        }
      }
      if (Ntrk_gen > 0) {
        registry.fill(HIST("Selection"), 4.);
        bevtc[kINELg0] = true;
      }
      if (bevtc[kINEL])
        registry.fill(HIST("Selection"), 9);
      if (bevtc[kINEL] && std::abs(mcz) < 10)
        registry.fill(HIST("Selection"), 11);
      if (bevtc[kINELg0])
        registry.fill(HIST("Selection"), 10);
      if (bevtc[kINELg0] && std::abs(mcz) < 10)
        registry.fill(HIST("Selection"), 12);
      for (auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= 3) {
            if (bevtc[kINEL] && std::abs(particle.eta()) < 0.8 && std::abs(mcz) < 10) {
              registry.fill(HIST("hgendndpt"), particle.pt());
              if (particle.pt() < 0.1) {
                registry.fill(HIST("hgendndpt2"), particle.pt(), -10.0 * particle.pt() + 2);
                registry.fill(HIST("hgendndpt05"), particle.pt(), 5.0 * particle.pt() + 0.5);
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
        for (auto& collision : collisionsample) {
          if (IsPbPb) {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
              cent = collision.centFT0C();
          } else {
            if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
              cent = collision.centFT0M();

            // auto Ntrk_rec = 0;
            // auto trackspart = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
            // for (auto& track : trackspart) {
            //   if (std::abs(track.eta()) < 1) {
            //     Ntrk_rec++;
            //   }
            // }
          }
        }
      }
      for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
        if (bevtc[ievtc])
          registry.fill(HIST("hgenzvtx"), Double_t(ievtc), genz, cent);
      }
      Int_t pid = 0;
      std::vector<Double_t> particleetas;
      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (std::abs(particle.pdgCode()) == 310 && std::abs(particle.eta()) < 0.5 && std::abs(genz) < 10)
          registry.fill(HIST("Selection"), 17.);
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
          if (std::abs(p->Charge()) >= 3) {
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              if (bevtc[ievtc]) {
                registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), cent, Double_t(pid), kNoPtVar, particle.phi());
                if (particle.pt() < 0.1) {
                  registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), cent, Double_t(pid), kPtUp, particle.phi(), -10.0 * particle.pt() + 2);
                  registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), cent, Double_t(pid), kPtDw, particle.phi(), 5.0 * particle.pt() + 0.5);

                } else {
                  registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), cent, Double_t(pid), kPtUp, particle.phi());
                  registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), cent, Double_t(pid), kPtDw, particle.phi());
                }
              }
            }
            if (pid >= kPion && pid <= kOPar)
              particleetas.push_back(particle.eta());
          }
        }
      }

      for (auto& collision : collisionsample) {
        auto cent = -1.f;
        if (IsPbPb) {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
            cent = collision.centFT0C();
        } else {
          if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
            cent = collision.centFT0M();
          // auto Ntrk_rec = 0;
          // auto trackspart = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
          // for (auto& track : trackspart) {
          //   if (std::abs(track.eta()) < 1) {
          //     Ntrk_rec++;
          //   }
          // }
        }

        Bool_1d btrigc(kTrigend, false);
        auto z = collision.posZ();
        // if (collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        // if (collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        if (collision.sel8()) {
          btrigc[kSel8] = true;
          registry.fill(HIST("Selection"), 2.);
          if (std::abs(z) < 10) {
            registry.fill(HIST("Selection"), 3.);
          }
        }
        if (bevtc[kINEL] && btrigc[kSel8] && std::abs(z) < 10)
          registry.fill(HIST("hft0multiplicity"), collision.multFT0C());
        if (collisionsample.size() == 1 && bevtc[kINELg0] && btrigc[kSel8]) {
          for (auto eta : particleetas) {
            registry.fill(HIST("genetaINELg0Sel8recz10"), eta, z);
            registry.fill(HIST("genetaINELg0Sel8genz10"), eta, mcz);
          }
          registry.fill(HIST("reczINELg0Sel8"), z);
          registry.fill(HIST("genzINELg0Sel8"), genz);
        }
        if (collisionsample.size() == 1 && bevtc[kINELg0]) {
          for (auto eta : particleetas) {
            registry.fill(HIST("genetaINELg0genz10"), eta, mcz);
          }
          registry.fill(HIST("genzINELg0"), genz);
        }

        auto Ntrk_rec = 0;
        auto trackspart = tracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        for (auto& track : trackspart) {
          if (std::abs(track.eta()) < 1) {
            Ntrk_rec++;
          }
        }

        if (Ntrk_rec > 0) {
          if (btrigc[kSel8])
            btrigc[kSel8g0] = true;
        }
        if (IsPbPb) {
          if (std::abs(z) < 10 && btrigc[kSel8])
            registry.fill(HIST("hcentrality"), cent);
        } else {
          if (std::abs(z) < 10 && btrigc[kSel8g0])
            registry.fill(HIST("hcentrality"), cent);
        }

        if (bevtc[kINEL] && btrigc[kSel8])
          registry.fill(HIST("Selection"), 13);
        if (bevtc[kINEL] && btrigc[kSel8] && std::abs(z) < 10)
          registry.fill(HIST("Selection"), 15);
        if (bevtc[kINELg0] && btrigc[kSel8g0])
          registry.fill(HIST("Selection"), 14);
        if (bevtc[kINELg0] && btrigc[kSel8g0] && std::abs(z) < 10)
          registry.fill(HIST("Selection"), 16);

        for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
          for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
            if (bevtc[ievtc] && btrigc[itrigc]) {
              registry.fill(HIST("hreczvtx"), Double_t(ievtc), Double_t(itrigc), z, cent);
            }
          }
        }
        std::vector<Int_t> mclabels;
        for (auto& track : trackspart) {
          if (track.has_mcParticle()) {
            Int_t pid = kBkg;
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
            for (auto MotherIDs = particle.mothersIds().front(); MotherIDs <= particle.mothersIds().back(); MotherIDs++) {
              auto mother = mcParticles.rawIteratorAt(MotherIDs);
              auto pdg_mother = mother.pdgCode();
              if (pdg_mother == 310 || std::abs(pdg_mother) == 3122) {
                pid = kMotherStrange;
              }
            }
            if (find(mclabels.begin(), mclabels.end(), track.mcParticleId()) != mclabels.end())
              pid = kBkg;
            mclabels.push_back(track.mcParticleId());
            registry.fill(HIST("hetaresponse"), particle.eta(), track.eta(), cent);
            if (bevtc[kINEL] && btrigc[kSel8] && std::abs(track.eta()) < 0.8 && std::abs(z) < 10 && pid != kBkg && pid != kNotPrimary)
              registry.fill(HIST("hdndptefficiency"), particle.pt());
            if (btrigc[kSel8] && bevtc[kINELg0])
              registry.fill(HIST("recetaINELg0Sel8recz10"), track.eta(), z);
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), Double_t(ievtc), Double_t(itrigc), z, particle.eta(), cent, Double_t(pid), particle.phi());
                  registry.fill(HIST("hphieta"), Double_t(ievtc), Double_t(itrigc), track.phi(), track.eta(), cent);
                }
              }
            }

          } else {
            for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
              for (auto itrigc = 1u; itrigc < kTrigend; itrigc++) {
                if (bevtc[ievtc] && btrigc[itrigc]) {
                  registry.fill(HIST("hrecdndeta"), Double_t(ievtc), Double_t(itrigc), z, track.eta(), cent, Double_t(kBkg), track.phi());
                }
              }
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processMCCounting, "MC Count tracks", false);

  void processTrackEfficiencyGeneral(
    typename soa::Join<MyCollisions, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& /*mcCollisions*/, Particles const& particles,
    FiLTracks const& /*tracks*/)
  {

    if (!collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }

    auto mcCollision = collision.mcCollision();
    auto particlesPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyGeneral, "MC Count tracks", false);

  void processCounting(
    MyCollisionsCent const& collisions,
    FiTracks const& tracks)
  {
    runCounting(collisions, tracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks with Centrality", false);

  void processGen(
    aod::McCollisions::iterator const& mcCollision, Particles const& mcParticles)
  {

    auto genz = mcCollision.posZ();
    Bool_1d bevtc(kECend, false);
    bevtc[kINEL] = true;
    for (auto& particle : mcParticles) {
      if (!particle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          if (std::abs(particle.eta()) < 1)
            bevtc[kINELg0] = true;
        }
      }
    }
    for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
      if (bevtc[ievtc])
        registry.fill(HIST("hgenzvtx"), Double_t(ievtc), genz, -1.0);
    }
    Int_t pid = 0;
    for (auto& particle : mcParticles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (std::abs(particle.pdgCode()) == 310 && std::abs(particle.eta()) < 0.5 && std::abs(genz) < 10)
        registry.fill(HIST("Selection"), 17.);
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
        if (std::abs(p->Charge()) >= 3) {
          for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
            if (bevtc[ievtc]) {
              registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), -1.0, Double_t(pid), kNoPtVar, particle.phi());
              if (particle.pt() < 0.1) {
                registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), -1.0, Double_t(pid), kPtUp, particle.phi(), 2.0);
                registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), -1.0, Double_t(pid), kPtDw, particle.phi(), 0.5);

              } else {
                registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), -1.0, Double_t(pid), kPtUp, particle.phi());
                registry.fill(HIST("hgendndeta"), Double_t(ievtc), genz, particle.eta(), -1.0, Double_t(pid), kPtDw, particle.phi());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);
  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  void processV0Counting(
    MyCollisionsCent const& collisions,
    aod::V0Datas const& fullV0s,
    FiTracks const& /*tracks*/,
    DaughterTracks const& /*Dautrks*/)
  {
    for (auto& collision : collisions) {
      if (!collision.sel8())
        continue;
      auto z = collision.posZ();

      auto cent = -1.f;
      if (IsPbPb) {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Cs>())
          cent = collision.centFT0C();
      } else {
        if constexpr (MyCollisionsCent::template contains<aod::CentFT0Ms>())
          cent = collision.centFT0M();
      }

      auto v0s_per_coll = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      for (auto& v0 : v0s_per_coll) {

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();
        if (std::abs(z) < 10) {
          if (v0.v0radius() > v0radius)
            continue;
          if (v0.dcapostopv() > dcapostopv)
            continue;
          if (v0.dcanegtopv() > dcanegtopv)
            continue;
          if (v0.v0cosPA() < v0cospa)
            continue;
          if (fabs(pTrack.eta()) > 0.9)
            continue;
          if (fabs(nTrack.eta()) > 0.9)
            continue;

          if (fabs(v0.eta()) < 0.5)
            registry.fill(HIST("hv0k0s"), v0.mK0Short());
          registry.fill(HIST("hv0mass"), cent, Double_t(kK0short), v0.eta(), Double_t(v0.mK0Short()));
          registry.fill(HIST("hv0mass"), cent, Double_t(kLambda), v0.eta(), Double_t(v0.mLambda()));
          registry.fill(HIST("hv0mass"), cent, Double_t(kAntilambda), v0.eta(), Double_t(v0.mAntiLambda()));
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processV0Counting, "MC V0Count tracks without Centrality", false);

  void processMCV0Counting(
    soa::Join<MyCollisionsCent, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& /*MCCollisions*/,
    Particles const& /*mcParticles*/,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    soa::Filtered<LabeledTracksEx> const& /*tracks*/,
    DaughterTracks const& /*Dautrks*/)
  {
    for (auto& collision : collisions) {
      auto cent = -1.f;

      if (IsPbPb) {
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

      for (auto& v0 : fullV0s) {

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();

        if (0 && v0.v0radius() > v0radius &&
            v0.v0cosPA() > v0cospa &&
            abs(pTrack.eta()) < etadau &&
            abs(nTrack.eta()) < etadau) {

          registry.fill(HIST("hv0mass"), cent, Double_t(kK0short), v0.eta(), Double_t(v0.mK0Short()));
          registry.fill(HIST("hv0mass"), cent, Double_t(kLambda), v0.eta(), Double_t(v0.mLambda()));
          registry.fill(HIST("hv0mass"), cent, Double_t(kAntilambda), v0.eta(), Double_t(v0.mAntiLambda()));
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCV0Counting, "MC V0Count tracks without Centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
