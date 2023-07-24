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

#include <Math/Vector4D.h>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "bestCollisionTable.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Index.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_cand_bplus;
using namespace o2::analysis::hf_cuts_bplus_to_d0_pi;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using DaughterTrack = soa::Join<aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, DaughterTrack>;
using FiTracks = soa::Filtered<ExTracks>;
using Particles = soa::Filtered<aod::McParticles>;
using Particle = Particles::iterator;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using DaughterTracks = soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

enum {
  kECbegin = 0,
  kDATA = 1,
  kDATAbg0,
  kINEL,
  kINELbg0,
  kECend
};
enum {
  kTrigbegin = 0,
  kSel8 = 1,
  kBackground,
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

AxisSpec ZAxis = {60, -30, 30, "Z (cm)", "zaxis"};
AxisSpec DeltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
AxisSpec DCAAxis = {601, -3.01, 3.01, "", "DCA axis"};
AxisSpec EtaAxis = {80, -4.0, 4.0, "#eta", "eta axis"};
AxisSpec PhiAxis = {629, 0, 2 * M_PI, "Rad", "phi axis"};
AxisSpec PtAxis = {2401, -0.005, 24.005, "#it{p}_{T} (GeV/c)", "P_{T} axis"};
AxisSpec EvtClassAxis = {kECend - 1, kECbegin + 0.5, kECend - 0.5, "", "event class"};
AxisSpec TrigClassAxis = {kTrigend - 1, kTrigbegin + 0.5, kTrigend - 0.5, "", "trigger class"};
std::vector<double> centBinning = {0, 10., 20., 30., 40., 50., 60., 70., 80., 100};
AxisSpec CentAxis = {centBinning, "", "centrality"};
AxisSpec SpeciesAxis = {kSpeciesend - 1, kSpeciesbegin + 0.5, kSpeciesend - 0.5, "", "species class"};
AxisSpec MassAxis = {600, 0.3f, 1.3f, "Mass (GeV/c^{2})", "Inv. Mass (GeV/c^{2})"};
AxisSpec SignAxis = {kSignend - 1, kSignbegin + 0.5, kSignend - 0.5, "", "sign"};
AxisSpec StepAxis = {kStepend - 1, kStepbegin + 0.5, kStepend - 0.5, "", "step"};
AxisSpec testAxis = {101, -0.5, 100.5, "", "test"};
AxisSpec multAxis = {1001, -0.5, 1000.5, "", "Ntrks"};
AxisSpec StatusCodeAxis = {3, -1.5, 2.5, "", "StatusCode"};
AxisSpec ProcessCodeAxis = {45, -1.5, 44.5, "", "StatusCode"};

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

  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> v0radius{"v0radius", 0.5, "Radius"};
  Configurable<float> etadau{"etadau", 4, "Eta Daughters"};
  Configurable<float> rapidity{"v0rapidity", 0.5, "V0 rapidity"};

  Configurable<bool> mftanalysis{"mftanalysis", false, "mft analysis switch"};
  Configurable<bool> zvtxcut{"zvtxcut", false, "z vtx cut < 10cm"};

  HistogramRegistry registry{
    "registry",
    {{"Events/Selection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}}}};

  std::vector<int> usedTracksIds;
  void init(InitContext&)
  {
    if (doprocessCountingWithCent) {
      registry.add({"Tracks/ProcessCounting/Centrality/Centrality", " ; centrality_FT0C (%) ", {HistType::kTH1F, {CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hrecpt", " eventclass;  pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hreczvtx", "evntclass; triggerclass;  zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis, CentAxis}}});

      registry.add({"Tracks/ProcessCounting/Centrality/Ntrk_rec", " ; Ntrk_rec", {HistType::kTHnSparseD, {EvtClassAxis, multAxis, CentAxis}}});
    }
    if (doprocessCountingWithoutCent) {
      registry.add({"Tracks/ProcessCounting/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessCounting/hrecpt", " eventclass; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis}}});
      registry.add({"Tracks/ProcessCounting/hreczvtx", "evntclass; triggerclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis}}});
      registry.add({"Tracks/ProcessCounting/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessCounting/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessCounting/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});

      registry.add({"Tracks/ProcessCounting/Ntrk_rec", " ; Ntrk_rec", {HistType::kTHnSparseD, {EvtClassAxis, multAxis}}});
    }

    if (doprocessV0CountingWithCent) {
      registry.add({"Tracks/ProcessV0Counting/Centrality/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis, CentAxis}}});
      registry.add({"Tracks/ProcessV0Counting/Centrality/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessV0Counting/Centrality/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis, CentAxis}}});
    }
    if (doprocessV0CountingWithoutCent) {
      registry.add({"Tracks/ProcessV0Counting/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis}}});
      registry.add({"Tracks/ProcessV0Counting/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessV0Counting/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis}}});
    }

    if (doprocessMCCountingWithoutCent) {
      registry.add({"Tracks/ProcessMCCounting/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hreczvtx", "evntclass; triggerclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hrecpt", " eventclass; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hgenpt", " eventclass; centrality; pt;  ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis}}});
      registry.add({"Tracks/ProcessMCCounting/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCCounting/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessMCCounting/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});

      registry.add({"Tracks/ProcessMCCounting/Multiplicity", " ; Ntrk_rec; Ntrk_gen", {HistType::kTHnSparseD, {EvtClassAxis, multAxis, multAxis}}});
      registry.add({"Tracks/ProcessMCCounting/Zvtx", " ; Ntrk_rec; Ntrk_gen", {HistType::kTHnSparseD, {ZAxis, ZAxis}}});
      registry.add({"Tracks/ProcessMCCounting/Ntrk_missing", " ; Ntrk_missing", {HistType::kTHnSparseD, {EvtClassAxis, multAxis}}});
      registry.add({"Tracks/ProcessMCCounting/Ntrk_fake", " ; Ntrk_fake", {HistType::kTHnSparseD, {EvtClassAxis, multAxis}}});
      registry.add({"Tracks/ProcessMCCounting/Ntrk_rec", " ; Ntrk_rec", {HistType::kTHnSparseD, {EvtClassAxis, multAxis}}});
      registry.add({"Tracks/ProcessMCCounting/Ntrk_gen", " ; Ntrk_gen", {HistType::kTHnSparseD, {EvtClassAxis, multAxis}}});

      registry.add({"Tracks/ProcessMCCounting/hStatusCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, StatusCodeAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hMCStatusCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, StatusCodeAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hProcessCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, ProcessCodeAxis}}});
    }

    if (doprocessMCV0CountingWithoutCent) {
      registry.add({"Tracks/ProcessMCV0Counting/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis}}});
      registry.add({"Tracks/ProcessMCV0Counting/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCV0Counting/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis}}});
      registry.add({"Tracks/ProcessMCV0Counting/hMotherV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis}}});
    }
    if (doprocessGen) {
      registry.add({"Tracks/ProcessGen/hgendndeta", "evntclass;  zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessGen/hgenzvtx", "evntclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis}}});
    }

    auto hstat = registry.get<TH1>(HIST("Events/Selection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Sel8");
    x->SetBinLabel(3, "Selected INEL>0");
    x->SetBinLabel(4, "Generated INEL>0");
    x->SetBinLabel(5, "Good BCs");
    x->SetBinLabel(6, "BCs with collisions");
    x->SetBinLabel(7, "BCs with pile-up/splitting");
  }

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                                    (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                                    true) &&
                                             ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA);
  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  expressions::Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Partition<aod::MFTTracks> tSample2 = (aod::fwdtrack::eta < -2.f) && (aod::fwdtrack::eta > -4.f);
  Partition<FiTracks> tSample3 = nabs(aod::track::eta) < estimatorEta;
  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;

  Preslice<FiTracks> perCol = aod::track::collisionId;
  Preslice<soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>> mcmfttracks_slice = o2::aod::fwdtrack::collisionId;
  Preslice<soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>> mctracks_slice = aod::track::collisionId;
  Preslice<aod::McParticles> mcparticle_slice = o2::aod::mcparticle::mcCollisionId;
  Preslice<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>> tracks_slice = aod::track::collisionId;
  Preslice<aod::MFTTracks> mfttracks_slice = o2::aod::fwdtrack::collisionId;

  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(kIsBBT0A) &
                        bc.selection_bit(kIsBBT0C)) != 0) {
        registry.fill(HIST("Events/Selection"), 5.);
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
          registry.fill(HIST("Events/Selection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("Events/Selection"), 7.);
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
    FiTracks const& tracks,
    aod::MFTTracks const& mfttracks)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("Events/Selection"), 1.);
      auto z = collision.posZ();
      auto permfttracks = tSample2->sliceByCached(aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto pertracks = tSample3->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto Ntrk = 0;

      if (useEvSel && !collision.sel8())
        continue; // event selection cut
      registry.fill(HIST("Events/Selection"), 2.);
      float cent = 0;
      for (auto& track : pertracks) {
        [[maybe_unused]] int dummy = track.globalIndex();
        if (std::abs(track.eta()) < 0.5)
          Ntrk++; // charged track check
      }
      if (Ntrk > 0)
        registry.fill(HIST("Events/Selection"), 3.);

      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        // with centrality
        cent = collision.centFT0C();
        registry.fill(HIST("Tracks/ProcessCounting/Centrality/Centrality"), cent);
        if (Ntrk > 0)
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/hreczvtx"), Double_t(kDATAbg0), Double_t(kSel8), z, cent);
        registry.fill(HIST("Tracks/ProcessCounting/Centrality/hreczvtx"), Double_t(kDATA), Double_t(kSel8), z, cent);

        for (auto& track : pertracks) {
          if (Ntrk > 0) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/PhiEta"), Double_t(kDATAbg0), track.phi(), track.eta(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAXY"), Double_t(kDATAbg0), track.dcaXY(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAZ"), Double_t(kDATAbg0), track.dcaZ(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecpt"), Double_t(kDATAbg0), -1, track.pt(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecdndeta"), Double_t(kDATAbg0), Double_t(kSel8), z, track.eta(), cent);
          }
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/PhiEta"), Double_t(kDATA), track.phi(), track.eta(), cent);
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAXY"), Double_t(kDATA), track.dcaXY(), cent);
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAZ"), Double_t(kDATA), track.dcaZ(), cent);
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecpt"), Double_t(kDATA), -1, track.pt(), cent);
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecdndeta"), Double_t(kDATA), Double_t(kSel8), z, track.eta(), cent);
        }
        if (mftanalysis) {
          for (auto& track : permfttracks) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecdndeta"), Double_t(kDATA), Double_t(kSel8), z, track.eta(), cent);
          }
        }

        if (Ntrk > 0) {
          if (zvtxcut && std::abs(z) < 10) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/Ntrk_rec"), Double_t(kDATAbg0), Ntrk, cent);
          } else if (!zvtxcut) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/Ntrk_rec"), Double_t(kDATAbg0), Ntrk, cent);
          }
        }

        registry.fill(HIST("Tracks/ProcessCounting/Centrality/Ntrk_rec"), Double_t(kDATA), Ntrk, cent);
      } else {
        // without centrality
        if (Ntrk > 0)
          registry.fill(HIST("Tracks/ProcessCounting/hreczvtx"), Double_t(kDATAbg0), Double_t(kSel8), z);
        registry.fill(HIST("Tracks/ProcessCounting/hreczvtx"), Double_t(kDATA), Double_t(kSel8), z);

        for (auto& track : pertracks) {
          if (Ntrk > 0) {
            registry.fill(HIST("Tracks/ProcessCounting/PhiEta"), Double_t(kDATAbg0), track.phi(), track.eta());
            registry.fill(HIST("Tracks/ProcessCounting/DCAXY"), Double_t(kDATAbg0), track.dcaXY());
            registry.fill(HIST("Tracks/ProcessCounting/DCAZ"), Double_t(kDATAbg0), track.dcaZ());
            registry.fill(HIST("Tracks/ProcessCounting/hrecpt"), Double_t(kDATAbg0), -1, track.pt());
            registry.fill(HIST("Tracks/ProcessCounting/hrecdndeta"), Double_t(kDATAbg0), Double_t(kSel8), z, track.eta());
          }
          registry.fill(HIST("Tracks/ProcessCounting/PhiEta"), Double_t(kDATA), track.phi(), track.eta());
          registry.fill(HIST("Tracks/ProcessCounting/DCAXY"), Double_t(kDATA), track.dcaXY());
          registry.fill(HIST("Tracks/ProcessCounting/DCAZ"), Double_t(kDATA), track.dcaZ());
          registry.fill(HIST("Tracks/ProcessCounting/hrecpt"), Double_t(kDATA), -1, track.pt());
          registry.fill(HIST("Tracks/ProcessCounting/hrecdndeta"), Double_t(kDATA), Double_t(kSel8), z, track.eta());
        }
        if (mftanalysis) {
          for (auto& track : permfttracks) {
            registry.fill(HIST("Tracks/ProcessCounting/hrecdndeta"), Double_t(kDATA), Double_t(kSel8), z, track.eta());
          }
        }
        if (Ntrk > 0) {
          if (zvtxcut && std::abs(z) < 10) {
            registry.fill(HIST("Tracks/ProcessCounting/Ntrk_rec"), Double_t(kDATAbg0), Ntrk);
          } else if (!zvtxcut) {
            registry.fill(HIST("Tracks/ProcessCounting/Ntrk_rec"), Double_t(kDATAbg0), Ntrk);
          }
        }
        registry.fill(HIST("Tracks/ProcessCounting/Ntrk_rec"), Double_t(kDATA), Ntrk);
      }
    }
  }

  template <typename C, typename MC>
  void runMCCounting(typename soa::Join<C, aod::McCollisionLabels> const& collisions,
                     MC const&,
                     Particles const& mcParticles,
                     soa::Filtered<LabeledTracksEx> const&,
                     DaughterTracks const&,
                     soa::SmallGroups<aod::ReassignedTracksCore> const& atracks,
                     soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks)
  {
    [[maybe_unused]] constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>() || hasCent<MC>();
    for (auto& collision : collisions) {
      registry.fill(HIST("Events/Selection"), 1.);
      auto z = collision.posZ();
      auto mcCollision = collision.mcCollision();
      auto mcz = mcCollision.posZ();

      auto Ntrk_rec = 0;
      auto Ntrk_gen = 0;

      if (useEvSel && !collision.sel8()) // event selection cut
        continue;
      registry.fill(HIST("Events/Selection"), 2.);
      if (!collision.has_mcCollision()) // check mc particle
        continue;

      registry.fill(HIST("Tracks/ProcessMCCounting/hreczvtx"), Double_t(kINEL), Double_t(kSel8), z);

      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      tracks.bindExternalIndices(&mcParticles);
      auto permcmfttracks = mfttracks.sliceBy(mcmfttracks_slice, collision.globalIndex());

      usedTracksIds.clear();

      for (auto& track : atracks) {
        auto ktrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ktrack.globalIndex());
        if (ktrack.has_mcParticle()) {
          if (std::abs(ktrack.eta()) < 0.5)
            Ntrk_rec++;
        }
      }

      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end())
          continue;
        if (track.has_mcParticle()) {
          if (std::abs(track.eta()) < 0.5)
            Ntrk_rec++;
        }
      }

      if (Ntrk_rec > 0)
        registry.fill(HIST("Events/Selection"), 3.);
      usedTracksIds.clear();

      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          if (Ntrk_rec > 0) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINELbg0), Double_t(kSel8), z, ttrack.mcParticle_as<Particles>().eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINELbg0), ttrack.mcParticle_as<Particles>().pt(), ttrack.pt());
            registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINELbg0), ttrack.phi(), ttrack.eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINELbg0), ttrack.dcaXY());
            registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINELbg0), ttrack.dcaZ());
          }
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kSel8), z, ttrack.mcParticle_as<Particles>().eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINEL), ttrack.mcParticle_as<Particles>().pt(), ttrack.pt());
          registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINEL), ttrack.phi(), ttrack.eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINEL), ttrack.dcaXY());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINEL), ttrack.dcaZ());

          registry.fill(HIST("Tracks/ProcessMCCounting/hStatusCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getGenStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hMCStatusCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getHepMCStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hProcessCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getProcess());

        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end())
          continue;
        if (track.has_mcParticle()) {
          if (Ntrk_rec > 0) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINELbg0), Double_t(kSel8), z, track.template mcParticle_as<Particles>().eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINELbg0), track.template mcParticle_as<Particles>().pt(), track.pt());
            registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINELbg0), track.phi(), track.eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINELbg0), track.dcaXY());
            registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINELbg0), track.dcaZ());
          }
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kSel8), z, track.template mcParticle_as<Particles>().eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINEL), track.template mcParticle_as<Particles>().pt(), track.pt());
          registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINEL), track.phi(), track.eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINEL), track.dcaXY());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINEL), track.dcaZ());

          registry.fill(HIST("Tracks/ProcessMCCounting/hStatusCode"), Double_t(kINEL), Double_t(kAll), track.template mcParticle_as<Particles>().getGenStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hMCStatusCode"), Double_t(kINEL), Double_t(kAll), track.template mcParticle_as<Particles>().getHepMCStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hProcessCode"), Double_t(kINEL), Double_t(kAll), track.template mcParticle_as<Particles>().getProcess());
        } else {
          // when secondary
        }
      }
      if (mftanalysis) {
        for (auto& track : permcmfttracks) {
          if (track.has_mcParticle()) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kSel8), z, track.template mcParticle_as<Particles>().eta());
          }
        }
      }
      for (auto& particle : particles) {
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= 3) {
            if (std::abs(particle.eta()) < 0.5)
              Ntrk_gen++;
          }
        }
      }

      if (Ntrk_gen > 0)
        registry.fill(HIST("Events/Selection"), 4.);

      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          if (std::abs(p->Charge()) >= 3) {
            if (Ntrk_gen > 0)
              registry.fill(HIST("Tracks/ProcessMCCounting/hgenpt"), Double_t(kINELbg0), particle.pt());
            registry.fill(HIST("Tracks/ProcessMCCounting/hgenpt"), Double_t(kINEL), particle.pt());
          }
        }
      }
      if (Ntrk_rec > 0 && Ntrk_gen > 0) {
        if (std::abs(z) < 10 && std::abs(mcz) < 10)
          registry.fill(HIST("Tracks/ProcessMCCounting/Multiplicity"), Double_t(kINELbg0), Ntrk_rec, Ntrk_gen); // multiplicity matrix
        if (std::abs(z) > 10 && std::abs(mcz) < 10)
          registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_missing"), Double_t(kINELbg0), Ntrk_gen); // missing trk
        if (std::abs(z) < 10 && std::abs(mcz) > 10)
          registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_fake"), Double_t(kINELbg0), Ntrk_rec); // fake trk
        if (std::abs(z) < 10)
          registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_rec"), Double_t(kINELbg0), Ntrk_rec); // MC rec trk
        if (std::abs(mcz) < 10)
          registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_gen"), Double_t(kINELbg0), Ntrk_gen); // MC gen trk
      }
      registry.fill(HIST("Tracks/ProcessMCCounting/Multiplicity"), Double_t(kINEL), Ntrk_rec, Ntrk_gen); // multiplicity matrix
      registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_missing"), Double_t(kINEL), Ntrk_gen);           // missing trk
      registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_fake"), Double_t(kINEL), Ntrk_rec);              // fake trk
      registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_rec"), Double_t(kINEL), Ntrk_rec);               // MC rec trk
      registry.fill(HIST("Tracks/ProcessMCCounting/Ntrk_gen"), Double_t(kINEL), Ntrk_gen);               // MC gen trk
      registry.fill(HIST("Tracks/ProcessMCCounting/Zvtx"), z, mcz);                                      // zvtx matrix
    }
  }
  template <typename C>
  void runV0Counting(
    C const& collisions,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    FiTracks const& tracks)
  {
    for (auto& collision : collisions) {
      if (useEvSel && !collision.sel8())
        continue;
      float cent = 0;
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        cent = collision.centFT0C();
        for (auto& v0 : fullV0s) {
          registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kAll), cent);
          registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kAll), cent);
          registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kAll), cent);

          auto pTrack = v0.template posTrack_as<FiTracks>();
          auto nTrack = v0.template negTrack_as<FiTracks>();

          if (v0.v0radius() > v0radius &&
              v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
              abs(pTrack.eta()) < etadau &&
              abs(nTrack.eta()) < etadau) {
            registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kBasiccut), cent);
            registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kBasiccut), cent);
            registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kBasiccut), cent);

            if (abs(v0.yK0Short()) < rapidity) {
              registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short(), cent);
              if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509) {
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kK0short), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kK0short), nTrack.eta(), cent);
              }
            }

            if (abs(v0.yLambda()) < rapidity) {
              registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kLambda), v0.mLambda(), cent);
              if (1.11 < v0.mLambda() && v0.mLambda() < 1.12) {
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kLambda), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kLambda), nTrack.eta(), cent);
              }
              registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), v0.mAntiLambda(), cent);
              if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12) {
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessV0Counting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta(), cent);
              }
            }
          }
        }
      } else {
        for (auto& v0 : fullV0s) {
          registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kAll));
          registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kAll));
          registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kAll));

          auto pTrack = v0.template posTrack_as<FiTracks>();
          auto nTrack = v0.template negTrack_as<FiTracks>();

          if (v0.v0radius() > v0radius &&
              v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
              abs(pTrack.eta()) < etadau &&
              abs(nTrack.eta()) < etadau) {
            registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kBasiccut));
            registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kBasiccut));
            registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kBasiccut));

            if (abs(v0.yK0Short()) < rapidity) {
              registry.fill(HIST("Tracks/ProcessV0Counting/hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short());
              if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509) {
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kK0short), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kK0short), nTrack.eta());
              }
            }

            if (abs(v0.yLambda()) < rapidity) {
              registry.fill(HIST("Tracks/ProcessV0Counting/hV0Mass"), Double_t(kDATA), Double_t(kLambda), v0.mLambda());
              if (1.11 < v0.mLambda() && v0.mLambda() < 1.12) {
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kLambda), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kLambda), nTrack.eta());
              }
              registry.fill(HIST("Tracks/ProcessV0Counting/hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), v0.mAntiLambda());
              if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12) {
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessV0Counting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta());
              }
            }
          }
        }
      }
    }
  }

  template <typename C, typename MC>
  void runMCV0Counting(
    C const& collisions,
    MC const&,
    Particles const& mcParticles,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    soa::Filtered<LabeledTracksEx> const&,
    DaughterTracks const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    [[maybe_unused]] constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>() || hasCent<MC>();
    for (auto& collision : collisions) {
      auto mcCollision = collision.mcCollision();

      if (useEvSel && !collision.sel8()) // event selection cut
        continue;
      if (!collision.has_mcCollision()) // check mc particle
        continue;

      for (auto& v0 : fullV0s) {
        registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kAll));
        registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kAll));
        registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kAll));

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();

        if (v0.v0radius() > v0radius &&
            v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
            abs(pTrack.eta()) < etadau &&
            abs(nTrack.eta()) < etadau) {
          registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kBasiccut));
          registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kBasiccut));
          registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kBasiccut));

          if (abs(v0.yK0Short()) < rapidity) {
            registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Mass"), Double_t(kINEL), Double_t(kK0short), v0.mK0Short());
            if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kMasscut));
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kK0short), pTrack.eta());
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kK0short), nTrack.eta());
            }
          }

          if (abs(v0.yLambda()) < rapidity) {
            registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Mass"), Double_t(kINEL), Double_t(kLambda), v0.mLambda());
            if (1.11 < v0.mLambda() && v0.mLambda() < 1.12) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kMasscut));
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kLambda), pTrack.eta());
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kLambda), nTrack.eta());
            }
            registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Mass"), Double_t(kINEL), Double_t(kAntilambda), v0.mAntiLambda());
            if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kMasscut));
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta());
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta());
            }
          }
        }
      }

      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          for (auto MotherIDs = ttrack.mcParticle_as<Particles>().mothersIds().front(); MotherIDs <= ttrack.mcParticle_as<Particles>().mothersIds().back(); MotherIDs++) {
            auto mother = mcParticles.rawIteratorAt(MotherIDs);
            auto pdg_mother = mother.pdgCode();
            if (pdg_mother == 310) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kK0short));
            }
            if (pdg_mother == 3122) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kLambda));
            }
            if (pdg_mother == -3122) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kAntilambda));
            }
          }
        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end())
          continue;

        if (track.has_mcParticle()) {
          for (auto MotherIDs = track.template mcParticle_as<Particles>().mothersIds().front(); MotherIDs <= track.template mcParticle_as<Particles>().mothersIds().back(); MotherIDs++) {
            auto mother = mcParticles.rawIteratorAt(MotherIDs);
            auto pdg_mother = mother.pdgCode();
            if (pdg_mother == 310) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kK0short));
            }
            if (pdg_mother == 3122) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kLambda));
            }
            if (pdg_mother == -3122) {
              registry.fill(HIST("Tracks/ProcessMCV0Counting/hMotherV0Count"), Double_t(kINEL), Double_t(kAntilambda));
            }
          }
        } else {
          // when secondary
        }
      }
    }
  }

  void processCountingWithCent(
    MyCollisionsCent const& collisions,
    FiTracks const& tracks,
    aod::MFTTracks const& mfttracks)
  {
    // runCounting(collisions, tracks, fullV0s, mfttracks);
    runCounting(collisions, tracks, mfttracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithCent, "Count tracks with Centrality", false);

  void processCountingWithoutCent(
    MyCollisions const& collisions,
    FiTracks const& tracks,
    aod::MFTTracks const& mfttracks)
  {
    // runCounting(collisions, tracks, fullV0s, mfttracks);
    runCounting(collisions, tracks, mfttracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithoutCent, "Count tracks with No Centrality", false);

  void processV0CountingWithCent(
    MyCollisionsCent const& collisions,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    FiTracks const& tracks)
  {
    runV0Counting(collisions, fullV0s, tracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processV0CountingWithCent, "V0Count tracks with Centrality", false);

  void processV0CountingWithoutCent(
    MyCollisions const& collisions,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    FiTracks const& tracks)
  {
    runV0Counting(collisions, fullV0s, tracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processV0CountingWithoutCent, "V0Count tracks without Centrality", false);

  void processMCV0CountingWithoutCent(
    soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& MCCollisions,
    Particles const& mcParticles,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    soa::Filtered<LabeledTracksEx> const& tracks,
    DaughterTracks const& Dautrks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    runMCV0Counting(collisions, MCCollisions, mcParticles, fullV0s, tracks, Dautrks, atracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCV0CountingWithoutCent, "MC V0Count tracks without Centrality", false);

  void processMCCountingWithoutCent(
    soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const& MCCollisions,
    Particles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const& Fitrks,
    DaughterTracks const& Dautrks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks,
    soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& mfttracks)
  {
    // runMCCounting(collisions, MCCollisions, fullV0s, mcParticles, Fitrks, Dautrks, atracks, mfttracks);
    runMCCounting(collisions, MCCollisions, mcParticles, Fitrks, Dautrks, atracks, mfttracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCCountingWithoutCent, "MC Count tracks", false);

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
    Particles const& mcParticles, FiTracks const& tracks)
  {
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto genz = mcCollision.posZ();
    registry.fill(HIST("Tracks/ProcessGen/hgenzvtx"), Double_t(kINEL), genz);
    for (auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          registry.fill(HIST("Tracks/ProcessGen/hgendndeta"), Double_t(kINEL), genz, particle.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
