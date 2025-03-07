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
///
/// \file   flowEfficiencyCasc.cxx
/// \author Fuchun Cui(fcui@cern.ch)
/// \since  Feb/21/2025
/// \brief  This task is to calculate V0s and cascades local density efficiency

#include <vector>
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/cascqaanalysis.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowEfficiencyCasc {
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
  // topological cut for V0
  O2_DEFINE_CONFIGURABLE(cfgv0_radius, float, 5.0f, "minimum decay radius")
  O2_DEFINE_CONFIGURABLE(cfgv0_v0cospa, float, 0.995f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgv0_dcadautopv, float, 0.1f, "minimum daughter DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgv0_dcav0dau, float, 0.5f, "maximum DCA among V0 daughters")
  O2_DEFINE_CONFIGURABLE(cfgv0_mk0swindow, float, 0.1f, "Invariant mass window of K0s")
  O2_DEFINE_CONFIGURABLE(cfgv0_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
  O2_DEFINE_CONFIGURABLE(cfgv0_ArmPodocut, float, 0.2f, "Armenteros Podolski cut for K0")
  // topological cut for cascade
  O2_DEFINE_CONFIGURABLE(cfgcasc_radius, float, 0.5f, "minimum decay radius")
  O2_DEFINE_CONFIGURABLE(cfgcasc_casccospa, float, 0.999f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgcasc_v0cospa, float, 0.998f, "minimum cosine of pointing angle")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0topv, float, 0.01f, "minimum daughter DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcabachtopv, float, 0.01f, "minimum bachelor DCA to PV")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcacascdau, float, 0.3f, "maximum DCA among cascade daughters")
  O2_DEFINE_CONFIGURABLE(cfgcasc_dcav0dau, float, 1.0f, "maximum DCA among V0 daughters")
  O2_DEFINE_CONFIGURABLE(cfgcasc_mlambdawindow, float, 0.04f, "Invariant mass window of lambda")
  // track quality and type selections
  O2_DEFINE_CONFIGURABLE(cfgtpcclusters, int, 0, "minimum number of TPC clusters requirement")
  O2_DEFINE_CONFIGURABLE(cfgitsclusters, int, 0, "minimum number of ITS clusters requirement")
  O2_DEFINE_CONFIGURABLE(cfgtpcclufindable, int, 0, "minimum number of findable TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgtpccrossoverfindable, int, 0, "minimum number of Ratio crossed rows over findable clusters")
  O2_DEFINE_CONFIGURABLE(cfgcheckDauTPC, bool, true, "check daughter tracks TPC or not")
  O2_DEFINE_CONFIGURABLE(cfgcheckDauTOF, bool, false, "check daughter tracks TOF or not")
  O2_DEFINE_CONFIGURABLE(cfgCasc_rapidity, float, 0.5, "rapidity")

  O2_DEFINE_CONFIGURABLE(cfgNSigmatpctof, std::vector<float>, (std::vector<float>{3, 3, 3}), "tpc and tof NSigma for Pion Kaon Proton")

  ConfigurableAxis cfgaxisPt{"cfgaxisPt", {VARIABLE_WIDTH, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtXi{"cfgaxisPtXi", {VARIABLE_WIDTH, 0, 0.1, 0.5, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtOmega{"cfgaxisPtOmega", {VARIABLE_WIDTH, 0, 0.1, 0.5, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisPtV0{"cfgaxisPtV0", {VARIABLE_WIDTH, 0, 0.1, 0.5, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9}, "pt (GeV)"};
  ConfigurableAxis cfgaxisMultiplicity{"cfgaxisMultiplicity", {1000, 0, 5000}, "Nch"};

  using MyCollisions = soa::Join<aod::StraCollisions, aod::StraEvSels>;
  using MyMcCollisions = soa::Join<aod::StraMCCollisions, aod::StraMCCollMults>;
  using CascMCCandidates = soa::Join<aod::CascCollRefs, aod::CascCores, aod::CascExtras, aod::CascBBs, aod::CascCoreMCLabels>;
  using V0MCCandidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras, aod::V0CoreMCLabels>;
  using DaughterTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

  // Define the output
  HistogramRegistry registry{"registry"};

  std::vector<float> cfgNSigma = cfgNSigmatpctof;

  void init(InitContext const&)
  {
    const AxisSpec axisCounter{1, 0, +1, ""};
    // create histograms
    registry.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    registry.add("mcEventCounter", "Monte Carlo Truth EventCounter", kTH1F, {axisCounter});

    registry.add("h2DGenK0s", "", {HistType::kTH2D, {cfgaxisPtV0, cfgaxisMultiplicity}});
    registry.add("h2DGenLambda", "", {HistType::kTH2D, {cfgaxisPtV0, cfgaxisMultiplicity}});
    registry.add("h2DGenXi", "", {HistType::kTH2D, {cfgaxisPtXi, cfgaxisMultiplicity}});
    registry.add("h2DGenOmega", "", {HistType::kTH2D, {cfgaxisPtOmega, cfgaxisMultiplicity}});
    registry.add("h2DRecK0s", "", {HistType::kTH2D, {cfgaxisPtV0, cfgaxisMultiplicity}});
    registry.add("h2DRecLambda", "", {HistType::kTH2D, {cfgaxisPtV0, cfgaxisMultiplicity}});
    registry.add("h2DRecXi", "", {HistType::kTH2D, {cfgaxisPtXi, cfgaxisMultiplicity}});
    registry.add("h2DRecOmega", "", {HistType::kTH2D, {cfgaxisPtOmega, cfgaxisMultiplicity}});
  }

  void processRec(MyCollisions::iterator const& collision, V0MCCandidates const& V0s, CascMCCandidates const& Cascades, DaughterTracks const&, soa::Join<aod::CascMCCores, aod::CascMCCollRefs> const&, soa::Join<aod::V0MCCores, aod::V0MCCollRefs> const&)
  {
    registry.fill(HIST("eventCounter"), 0.5);
    if (!collision.sel8())
      return;
    int rectracknum = collision.multNTracksGlobal();
    for (const auto& casc : Cascades) {
      if (!casc.has_cascMCCore())
        continue;
      auto cascMC = casc.cascMCCore_as<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>();
      auto negdau = casc.negTrackExtra_as<DaughterTracks>();
      auto posdau = casc.posTrackExtra_as<DaughterTracks>();
      auto bachelor = casc.bachTrackExtra_as<DaughterTracks>();
      // track quality check
      if (bachelor.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (posdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (negdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (bachelor.itsNCls() < cfgitsclusters)
        continue;
      if (posdau.itsNCls() < cfgitsclusters)
        continue;
      if (negdau.itsNCls() < cfgitsclusters)
        continue;
      // topological cut
      if (casc.cascradius() < cfgcasc_radius)
        continue;
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_casccospa)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cfgcasc_v0cospa)
        continue;
      if (std::fabs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < cfgcasc_dcav0topv)
        continue;
      if (std::fabs(casc.dcabachtopv()) < cfgcasc_dcabachtopv)
        continue;
      if (casc.dcacascdaughters() > cfgcasc_dcacascdau)
        continue;
      if (casc.dcaV0daughters() > cfgcasc_dcav0dau)
        continue;
      if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cfgcasc_mlambdawindow)
        continue;
      // Omega and antiOmega
      int pdgCode{cascMC.pdgCode()};
      if (std::abs(pdgCode) == kOmegaMinus && std::abs(cascMC.pdgCodeV0()) == kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == kKPlus) {
        if (casc.sign() < 0 && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71) && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
            (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
          registry.fill(HIST("h2DRecOmega"), casc.pt(), rectracknum);
        } else if (casc.sign() < 0 && (casc.mOmega() > 1.63) && (casc.mOmega() < 1.71) && std::fabs(casc.yOmega()) < cfgCasc_rapidity &&
                   (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaKa()) < cfgNSigma[2] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
          registry.fill(HIST("h2DRecOmega"), casc.pt(), rectracknum);
        }
      }
      // Xi and antiXi
      if (std::abs(pdgCode) == kXiMinus && std::abs(cascMC.pdgCodeV0()) == kLambda0 && std::abs(cascMC.pdgCodeBachelor()) == kPiPlus) {
        if (casc.sign() < 0 && (casc.mXi() > 1.30) && (casc.mXi() < 1.37) && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
            (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(negdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
          registry.fill(HIST("h2DRecXi"), casc.pt(), rectracknum);
        } else if (casc.sign() < 0 && (casc.mXi() > 1.30) && (casc.mXi() < 1.37) && std::fabs(casc.yXi()) < cfgCasc_rapidity &&
                   (!cfgcheckDauTPC || (std::fabs(bachelor.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(posdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
          registry.fill(HIST("h2DRecXi"), casc.pt(), rectracknum);
        }
      }
    }

    for (const auto& v0 : V0s) {
      if (!v0.has_v0MCCore())
        continue;
      auto v0MC = v0.v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();
      auto v0negdau = v0.negTrackExtra_as<DaughterTracks>();
      auto v0posdau = v0.posTrackExtra_as<DaughterTracks>();

      // track quality check
      if (v0posdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0negdau.tpcNClsFound() < cfgtpcclusters)
        continue;
      if (v0posdau.tpcNClsFindable() < cfgtpcclufindable)
        continue;
      if (v0negdau.tpcNClsFindable() < cfgtpcclufindable)
        continue;
      if (v0posdau.tpcCrossedRowsOverFindableCls() < cfgtpccrossoverfindable)
        continue;
      if (v0posdau.itsNCls() < cfgitsclusters)
        continue;
      if (v0negdau.itsNCls() < cfgitsclusters)
        continue;
      // topological cut
      if (v0.v0radius() < cfgv0_radius)
        continue;
      if (v0.v0cosPA() < cfgv0_v0cospa)
        continue;
      if (v0.dcaV0daughters() > cfgv0_dcav0dau)
        continue;
      if (std::fabs(v0.dcapostopv()) < cfgv0_dcadautopv)
        continue;
      if (std::fabs(v0.dcanegtopv()) < cfgv0_dcadautopv)
        continue;

      int pdgCode{v0MC.pdgCode()};
      // K0short
      if (std::abs(pdgCode) == kK0Short && v0MC.pdgCodePositive() == kPiPlus && v0MC.pdgCodeNegative() == kPiMinus) {
        if (v0.qtarm() / std::fabs(v0.alpha()) > cfgv0_ArmPodocut && std::fabs(v0.y()) < 0.5 && std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < cfgv0_mk0swindow &&
            (!cfgcheckDauTPC || (std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
          registry.fill(HIST("h2DRecK0s"), v0.pt(), rectracknum);
        }
      }
      // Lambda and antiLambda
      if (std::fabs(v0.y()) < 0.5 && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < cfgv0_mlambdawindow &&
          (!cfgcheckDauTPC || (std::fabs(v0posdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0negdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
        if (std::abs(pdgCode) == kLambda0 && v0MC.pdgCodePositive() == kProton && v0MC.pdgCodeNegative() == kPiMinus)
          registry.fill(HIST("h2DRecLambda"), v0.pt(), rectracknum);
      } else if (std::fabs(v0.y()) < 0.5 && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda) < cfgv0_mlambdawindow &&
                 (!cfgcheckDauTPC || (std::fabs(v0negdau.tpcNSigmaPr()) < cfgNSigma[1] && std::fabs(v0posdau.tpcNSigmaPi()) < cfgNSigma[0]))) {
        if (std::abs(pdgCode) == kLambda0 && v0MC.pdgCodePositive() == kPiPlus && v0MC.pdgCodeNegative() == kProtonBar)
          registry.fill(HIST("h2DRecLambda"), v0.pt(), rectracknum);
      }
    }
  }
  PROCESS_SWITCH(FlowEfficiencyCasc, processRec, "process reconstructed information", true);

  void processGen(MyMcCollisions::iterator const&, soa::SmallGroups<soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCollLabels>> const& coll, const soa::SmallGroups<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>& cascMCs, const soa::SmallGroups<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>& v0MCs)
  {
    registry.fill(HIST("mcEventCounter"), 0.5);
    int rectracknum = 0;
    for (const auto& col : coll) {
      rectracknum = col.multNTracksGlobal();
    }
    for (auto const& cascmc : cascMCs) {
      if (std::abs(cascmc.pdgCode()) == kXiMinus) {
        if (std::fabs(cascmc.yMC()) < cfgCasc_rapidity)
          registry.fill(HIST("h2DGenXi"), cascmc.ptMC(), rectracknum);
      }
      if (std::abs(cascmc.pdgCode()) == kOmegaMinus) {
        if (std::fabs(cascmc.yMC()) < cfgCasc_rapidity)
          registry.fill(HIST("h2DGenOmega"), cascmc.ptMC(), rectracknum);
      }
    }
    for (auto const& v0mc : v0MCs) {
      if (std::abs(v0mc.pdgCode()) == kK0Short) {
        if (std::fabs(v0mc.yMC()) < cfgCasc_rapidity)
          registry.fill(HIST("h2DGenK0s"), v0mc.ptMC(), rectracknum);
      }
      if (std::abs(v0mc.pdgCode()) == kLambda0) {
        if (std::fabs(v0mc.yMC()) < cfgCasc_rapidity)
          registry.fill(HIST("h2DGenLambda"), v0mc.ptMC(), rectracknum);
      }
    }
  }
  PROCESS_SWITCH(FlowEfficiencyCasc, processGen, "process gen information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEfficiencyCasc>(cfgc)};
}
