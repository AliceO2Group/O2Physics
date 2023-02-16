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
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiTracks = soa::Join<soa::Filtered<ExTracks>, DaughterTrack>;
using Particles = soa::Filtered<aod::McParticles>;
using Particle = Particles::iterator;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using DaughterTracks = soa::Join<aod::Tracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kECend
};
enum {
  kTrigbegin = 0,
  kMBAND = 1,
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
AxisSpec testAxis = {4, 0.5, 4.5, "", "test"};
AxisSpec particleIDAxis = {10000, 0.5, 10000.5, "", "particleID"};
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
  Service<O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isMC{"isMC", false, "check if MC"};

  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  AxisSpec MultAxis = {multBinning, "N"};

  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> v0radius{"v0radius", 0.5, "Radius"};
  Configurable<float> etadau{"etadau", 4, "Eta Daughters"};
  Configurable<float> rapidity{"v0rapidity", 0.5, "V0 rapidity"};

  HistogramRegistry registry{
    "registry",
    {{"Events/Selection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}}}};

  std::vector<int> usedTracksIds;
  void init(InitContext&)
  {
    if (doprocessCountingWithCent) {
      registry.add({"Tracks/ProcessCounting/Centrality/Centrality", " ; centrality_FT0C (%) ", {HistType::kTH1F, {CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hrecpt", " eventclass;  pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hreczvtx", "evntclass; triggerclass;  zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/ProcessCounting/Centrality/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis, CentAxis}}});
    }
    if (doprocessCountingWithoutCent) {
      registry.add({"Tracks/ProcessCounting/Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis}}});
      registry.add({"Tracks/ProcessCounting/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessCounting/hrecpt", " eventclass; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis}}});
      registry.add({"Tracks/ProcessCounting/hreczvtx", "evntclass; triggerclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis}}});
      registry.add({"Tracks/ProcessCounting/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessCounting/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessCounting/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessCounting/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis}}});
      registry.add({"Tracks/ProcessCounting/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessCounting/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis}}});
    }
    if (doprocessMCCounting) {
      registry.add({"Tracks/ProcessMCCounting/Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hrecdndeta", "evntclass; triggerclass; zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hreczvtx", "evntclass; triggerclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, ZAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hrecpt", " eventclass; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis, PtAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hgenpt", " eventclass; centrality; pt;  ", {HistType::kTHnSparseD, {EvtClassAxis, PtAxis}}});
      registry.add({"Tracks/ProcessMCCounting/PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCCounting/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessMCCounting/DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, StepAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hV0DauEta", "", {HistType::kTHnSparseD, {EvtClassAxis, SignAxis, SpeciesAxis, EtaAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis, MassAxis}}});

      registry.add({"Tracks/ProcessMCCounting/hStatusCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, StatusCodeAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hMCStatusCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, StatusCodeAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hProcessCode", "", {HistType::kTHnSparseD, {EvtClassAxis, StepAxis, ProcessCodeAxis}}});
      registry.add({"Tracks/ProcessMCCounting/hMotherV0Count", "", {HistType::kTHnSparseD, {EvtClassAxis, SpeciesAxis}}});
    }
    if (doprocessGen) {
      registry.add({"Tracks/ProcesGen/hgendndeta", "evntclass;  zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis, EtaAxis}}});
      registry.add({"Tracks/ProcesGen/hgenzvtx", "evntclass; zvtex", {HistType::kTHnSparseD, {EvtClassAxis, ZAxis}}});
    }
  }
  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection()[kIsBBT0A] &
                        bc.selection()[kIsBBT0C]) != 0) {
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

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                                    (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                                    true) &&
                                             ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA);
  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));

  std::vector<Double_t> tracketas;
  template <typename CollisionTypes>
  void runCounting(typename CollisionTypes::iterator const& collision, BCsRun3 const& bcs, aod::FT0s const& ft0s, aod::FV0As const& fv0as, FiTracks const& tracks, soa::Filtered<aod::V0Datas> const& fullV0s, soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {

    const auto& foundBC = collision.template foundBC_as<BCsRun3>();
    float multT0A = 0;
    float multT0C = 0;
    float multV0A = 0;

    registry.fill(HIST("Events/Selection"), 1.);
    auto z = collision.posZ();
    if (!useEvSel || collision.sel8()) { // event selection cut
      if (std::abs(z) < 10) {            // z-vtx cut
        if (foundBC.has_ft0()) {
          for (auto amplitude : foundBC.ft0().amplitudeA()) {
            multT0A += amplitude;
          }
          for (auto amplitude : foundBC.ft0().amplitudeC()) {
            multT0C += amplitude;
          }
        } else {
          multT0A = multT0C = -999;
        }

        if (foundBC.has_fv0a()) {
          for (auto amplitude : foundBC.fv0a().amplitude()) {
            multV0A += amplitude;
          }
        } else {
          multV0A = -999;
        }
        float cent = 0;
        if constexpr (CollisionTypes::template contains<aod::CentFT0Cs>()) { // with centrality
          cent = collision.centFT0C();
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/Centrality"), cent);

          registry.fill(HIST("Events/Selection"), 2.);
          registry.fill(HIST("Tracks/ProcessCounting/Centrality/Multiplicity"), multV0A, multT0A, multT0C, cent);

          registry.fill(HIST("Tracks/ProcessCounting/Centrality/hreczvtx"), Double_t(kDATA), Double_t(kMBAND), z, cent);

          usedTracksIds.clear();
          tracketas.clear();

          for (auto& track : atracks) {
            auto otrack = track.template track_as<FiTracks>();
            tracketas.push_back(otrack.eta());
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/PhiEta"), Double_t(kDATA), otrack.phi(), otrack.eta(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAXY"), Double_t(kDATA), otrack.dcaXY(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAZ"), Double_t(kDATA), otrack.dcaZ(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecpt"), Double_t(kDATA), -1, otrack.pt(), cent);
          }

          for (auto& track : tracks) {
            if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
              continue;
            }
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/PhiEta"), Double_t(kDATA), track.phi(), track.eta(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAXY"), Double_t(kDATA), track.dcaXY(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/DCAZ"), Double_t(kDATA), track.dcaZ(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecpt"), Double_t(kDATA), -1, track.pt(), cent);
            tracketas.push_back(track.eta());
          }

          for (auto& eta : tracketas) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), z, eta, cent);
          }

          for (auto& v0 : fullV0s) {
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kAll), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kAll), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kAll), cent);

            auto pTrack = v0.template posTrack_as<FiTracks>();
            auto nTrack = v0.template negTrack_as<FiTracks>();

            if (v0.v0radius() > v0radius &&
                v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
                abs(pTrack.eta()) < etadau &&
                abs(nTrack.eta()) < etadau) {
              registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kBasiccut), cent);
              registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kBasiccut), cent);
              registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kBasiccut), cent);

              if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509 && abs(v0.yK0Short()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kK0short), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kK0short), nTrack.eta(), cent);
              }

              if (1.11 < v0.mLambda() && v0.mLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kLambda), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kLambda), nTrack.eta(), cent);
              }

              if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kMasscut), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta(), cent);
                registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta(), cent);
              }
            }

            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kLambda), v0.mLambda(), cent);
            registry.fill(HIST("Tracks/ProcessCounting/Centrality/hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), v0.mAntiLambda(), cent);
          }
        } else {
          // without centrality
          registry.fill(HIST("Events/Selection"), 2.);
          registry.fill(HIST("Tracks/ProcessCounting/Multiplicity"), multV0A, multT0A, multT0C);

          registry.fill(HIST("Tracks/ProcessCounting/hreczvtx"), Double_t(kDATA), Double_t(kMBAND), z);

          usedTracksIds.clear();
          tracketas.clear();

          for (auto& track : atracks) {
            auto otrack = track.template track_as<FiTracks>();
            tracketas.push_back(otrack.eta());
            registry.fill(HIST("Tracks/ProcessCounting/PhiEta"), Double_t(kDATA), otrack.phi(), otrack.eta());
            registry.fill(HIST("Tracks/ProcessCounting/DCAXY"), Double_t(kDATA), otrack.dcaXY());
            registry.fill(HIST("Tracks/ProcessCounting/DCAZ"), Double_t(kDATA), otrack.dcaZ());
            registry.fill(HIST("Tracks/ProcessCounting/hrecpt"), Double_t(kDATA), -1, otrack.pt());
          }

          for (auto& track : tracks) {
            if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
              continue;
            }
            registry.fill(HIST("Tracks/ProcessCounting/PhiEta"), Double_t(kDATA), track.phi(), track.eta());
            registry.fill(HIST("Tracks/ProcessCounting/DCAXY"), Double_t(kDATA), track.dcaXY());
            registry.fill(HIST("Tracks/ProcessCounting/DCAZ"), Double_t(kDATA), track.dcaZ());
            registry.fill(HIST("Tracks/ProcessCounting/hrecpt"), Double_t(kDATA), -1, track.pt());
            tracketas.push_back(track.eta());
          }

          for (auto& eta : tracketas) {
            registry.fill(HIST("Tracks/ProcessCounting/hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), z, eta);
          }

          for (auto& v0 : fullV0s) {
            registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kAll));
            registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kAll));
            registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kAll));

            auto pTrack = v0.template posTrack_as<FiTracks>();
            auto nTrack = v0.template negTrack_as<FiTracks>();

            if (v0.v0radius() > v0radius &&
                v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
                abs(pTrack.eta()) < etadau &&
                abs(nTrack.eta()) < etadau) {
              registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kBasiccut));
              registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kBasiccut));
              registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kBasiccut));

              if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509 && abs(v0.yK0Short()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kK0short), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kK0short), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kK0short), nTrack.eta());
              }

              if (1.11 < v0.mLambda() && v0.mLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kLambda), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kLambda), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kLambda), nTrack.eta());
              }

              if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
                registry.fill(HIST("Tracks/ProcessCounting/hV0Count"), Double_t(kDATA), Double_t(kAntilambda), Double_t(kMasscut));
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta());
                registry.fill(HIST("Tracks/ProcessCounting/hV0DauEta"), Double_t(kDATA), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta());
              }
            }

            registry.fill(HIST("Tracks/ProcessCounting/hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short());
            registry.fill(HIST("Tracks/ProcessCounting/hV0Mass"), Double_t(kDATA), Double_t(kLambda), v0.mLambda());
            registry.fill(HIST("Tracks/ProcessCounting/hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), v0.mAntiLambda());
          }
        }
      }
    }
  }

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  void processCountingWithCent(
    MyCollisionsCent::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    runCounting<MyCollisionsCent>(collision, bcs, ft0s, fv0as, tracks, fullV0s, atracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithCent, "Count tracks with Centrality", false);

  void processCountingWithoutCent(
    MyCollisions::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    runCounting<MyCollisions>(collision, bcs, ft0s, fv0as, tracks, fullV0s, atracks);
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingWithoutCent, "Count tracks with No Centrality", false);

  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<aod::Tracks> tSample = nabs(aod::track::eta) < estimatorEta;
  Preslice<FiTracks> perCol = aod::track::collisionId;
  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;

  void processMCCounting(
    soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&,
    soa::Filtered<aod::V0Datas> const& fullV0s,
    Particles const& mcParticles,
    // aod::McParticles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const&,
    DaughterTracks const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    for (auto& collision : collisions) {
      auto z = collision.posZ();
      if (useEvSel && !collision.sel8()) { // event selection cut
        continue;
      }
      if (!collision.has_mcCollision()) { // check mc particle
        continue;
      }
      if (std::abs(z) > 10) { // z-vtx cut
        continue;
      }
      for (auto& v0 : fullV0s) {
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kAll));
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kAll));
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kAll));

        auto pTrack = v0.template posTrack_as<DaughterTracks>();
        auto nTrack = v0.template negTrack_as<DaughterTracks>();

        if (v0.v0radius() > v0radius &&
            v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
            abs(pTrack.eta()) < etadau &&
            abs(nTrack.eta()) < etadau) {
          registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kBasiccut));
          registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kBasiccut));
          registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kBasiccut));

          if (0.482 < v0.mK0Short() && v0.mK0Short() < 0.509 && abs(v0.yK0Short()) < rapidity) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kK0short), Double_t(kMasscut));
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kK0short), pTrack.eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kK0short), nTrack.eta());
          }
          if (1.11 < v0.mLambda() && v0.mLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kLambda), Double_t(kMasscut));
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kLambda), pTrack.eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kLambda), nTrack.eta());
          }
          if (1.11 < v0.mAntiLambda() && v0.mAntiLambda() < 1.12 && abs(v0.yLambda()) < rapidity) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0Count"), Double_t(kINEL), Double_t(kAntilambda), Double_t(kMasscut));
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kPositive), Double_t(kAntilambda), pTrack.eta());
            registry.fill(HIST("Tracks/ProcessMCCounting/hV0DauEta"), Double_t(kINEL), Double_t(kNegative), Double_t(kAntilambda), nTrack.eta());
          }
        }
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Mass"), Double_t(kINEL), Double_t(kK0short), v0.mK0Short());
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Mass"), Double_t(kINEL), Double_t(kLambda), v0.mLambda());
        registry.fill(HIST("Tracks/ProcessMCCounting/hV0Mass"), Double_t(kINEL), Double_t(kAntilambda), v0.mAntiLambda());
      }

      registry.fill(HIST("Tracks/ProcessMCCounting/hreczvtx"), Double_t(kINEL), Double_t(kMBAND), z);
      auto mcCollision = collision.mcCollision();
      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), z, ttrack.mcParticle_as<Particles>().eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINEL), ttrack.mcParticle_as<Particles>().pt(), ttrack.pt());
          registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINEL), ttrack.phi(), ttrack.eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINEL), ttrack.dcaXY());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINEL), ttrack.dcaZ());

          registry.fill(HIST("Tracks/ProcessMCCounting/hStatusCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getGenStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hMCStatusCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getHepMCStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hProcessCode"), Double_t(kINEL), Double_t(kAll), ttrack.mcParticle_as<Particles>().getProcess());

          if (ttrack.mcParticle_as<Particles>().getProcess() != 0) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kBackground), z, ttrack.mcParticle_as<Particles>().eta());
          }

          for (auto MotherIDs = ttrack.mcParticle_as<Particles>().mothersIds().front(); MotherIDs <= ttrack.mcParticle_as<Particles>().mothersIds().back(); MotherIDs++) {
            auto mother = mcParticles.rawIteratorAt(MotherIDs);
            auto pdg_mother = mother.pdgCode();
            if (pdg_mother == 310) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kK0short));
            }
            if (pdg_mother == 3122) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kLambda));
            }
            if (pdg_mother == -3122) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kAntilambda));
            }
          }
        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (track.has_mcParticle()) {
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), z, track.mcParticle_as<Particles>().eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/hrecpt"), Double_t(kINEL), track.mcParticle_as<Particles>().pt(), track.pt());
          registry.fill(HIST("Tracks/ProcessMCCounting/PhiEta"), Double_t(kINEL), track.phi(), track.eta());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAXY"), Double_t(kINEL), track.dcaXY());
          registry.fill(HIST("Tracks/ProcessMCCounting/DCAZ"), Double_t(kINEL), track.dcaZ());

          registry.fill(HIST("Tracks/ProcessMCCounting/hStatusCode"), Double_t(kINEL), Double_t(kAll), track.mcParticle_as<Particles>().getGenStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hMCStatusCode"), Double_t(kINEL), Double_t(kAll), track.mcParticle_as<Particles>().getHepMCStatusCode());
          registry.fill(HIST("Tracks/ProcessMCCounting/hProcessCode"), Double_t(kINEL), Double_t(kAll), track.mcParticle_as<Particles>().getProcess());

          if (track.mcParticle_as<Particles>().getProcess() != 0) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hrecdndeta"), Double_t(kINEL), Double_t(kBackground), z, track.mcParticle_as<Particles>().eta());
          }

          for (auto MotherIDs = track.mcParticle_as<Particles>().mothersIds().front(); MotherIDs <= track.mcParticle_as<Particles>().mothersIds().back(); MotherIDs++) {
            auto mother = mcParticles.rawIteratorAt(MotherIDs);
            auto pdg_mother = mother.pdgCode();
            if (pdg_mother == 310) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kK0short));
            }
            if (pdg_mother == 3122) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kLambda));
            }
            if (pdg_mother == -3122) {
              registry.fill(HIST("Tracks/ProcessMCCounting/hMotherV0Count"), Double_t(kINEL), Double_t(kAntilambda));
            }
          }
        } else {
          // when secondary
        }
      }
      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          if (std::abs(p->Charge()) >= 3) {
            registry.fill(HIST("Tracks/ProcessMCCounting/hgenpt"), Double_t(kINEL), particle.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCCounting, "MC Count tracks", false);

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
    Particles const& mcParticles, FiTracks const& tracks)
  {
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto genz = mcCollision.posZ();
    registry.fill(HIST("Tracks/ProcesGen/hgenzvtx"), Double_t(kINEL), genz);
    for (auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          registry.fill(HIST("Tracks/ProcesGen/hgendndeta"), Double_t(kINEL), genz, particle.eta());
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
