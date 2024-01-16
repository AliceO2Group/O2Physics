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
//
// Example V0 analysis task
// ========================
//
// This code loops over a V0Data table and produces some
// standard analysis output. It is meant to be run over
// derived data.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

struct derivedlambdakzeroanalysis {

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Selection criteria: acceptance
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
  Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

  // Standard 5 topological criteria
  Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
  Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
  Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
  Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};

  // PID (TPC)
  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  Configurable<bool> doQA{"doQA", true, "do topological variable QA histograms"};
  Configurable<float> massWindowQA{"massWindowQA", 0.005f, "mass window for QA plots"};

  static constexpr float defaultLifetimeCuts[1][2] = {{30., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0.0f, 100.0f}, ""};

  // topological variable QA axes
  ConfigurableAxis axisDCAtoPV{"axisDCAtoPV", {100, -0.5f, 0.5f}, "DCA (cm)"};
  ConfigurableAxis axisDCAdau{"axisDCAdau", {200, 0.0f, 2.0f}, "DCA (cm)"};
  ConfigurableAxis axisPointingAngle{"axisPointingAngle", {200, 0.0f, 2.0f}, "pointing angle (rad)"};
  ConfigurableAxis axisV0Radius{"axisV0Radius", {200, 0.0f, 60.0f}, "V0 2D radius (cm)"};

  enum species { spK0Short = 0,
                 spLambda,
                 spAntiLambda };

  void init(InitContext const&)
  {
    // Event Counters
    histos.add("hEventSelection", "hEventSelection", kTH1F, {{2, -0.5f, +2.5f}});
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");

    histos.add("hEventCentrality", "hEventCentrality", kTH1F, {{100, 0.0f, +100.0f}});

    // histograms versus mass
    histos.add("h3dMassK0Short", "h3dMassK0Short", kTH3F, {axisCentrality, axisPt, axisK0Mass});
    histos.add("h3dMassLambda", "h3dMassLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});
    histos.add("h3dMassAntiLambda", "h3dMassAntiLambda", kTH3F, {axisCentrality, axisPt, axisLambdaMass});

    // QA histograms if requested
    if (doQA) {
      // initialize for K0short...
      histos.add("K0Short/h3dPosDCAxy", "h3dPosDCAxy", kTH3F, {axisCentrality, axisPt, axisDCAtoPV});
      histos.add("K0Short/h3dNegDCAxy", "h3dNegDCAxy", kTH3F, {axisCentrality, axisPt, axisDCAtoPV});
      histos.add("K0Short/h3dDCADaughters", "h3dDCADaughters", kTH3F, {axisCentrality, axisPt, axisDCAdau});
      histos.add("K0Short/h3dPointingAngle", "h3dPointingAngle", kTH3F, {axisCentrality, axisPt, axisPointingAngle});
      histos.add("K0Short/h3dV0Radius", "h3dV0Radius", kTH3F, {axisCentrality, axisPt, axisV0Radius});

      // ...clone for Lambda and AntiLambda
      histos.addClone("K0Short/", "Lambda/");
      histos.addClone("K0Short/", "AntiLambda/");
    }
  }

  template <typename TV0>
  bool compatibleTPC(TV0 v0, int sp)
  {
    float pidPos = TMath::Abs(v0.template posTrackExtra_as<dauTracks>().tpcNSigmaPi());
    float pidNeg = TMath::Abs(v0.template negTrackExtra_as<dauTracks>().tpcNSigmaPi());

    if (sp == spLambda)
      pidPos = TMath::Abs(v0.template posTrackExtra_as<dauTracks>().tpcNSigmaPr());
    if (sp == spAntiLambda)
      pidNeg = TMath::Abs(v0.template negTrackExtra_as<dauTracks>().tpcNSigmaPr());

    if (pidPos < TpcPidNsigmaCut && pidNeg < TpcPidNsigmaCut)
      return true;

    // if not, then not
    return false;
  }

  template <typename TV0>
  void fillQAHistograms(TV0 v0, int sp, float centrality)
  {
    float massRef = o2::constants::physics::MassKaonNeutral;

    if (sp != spK0Short)
      massRef = o2::constants::physics::MassLambda;

    if (std::abs(v0.mK0Short() - massRef) < massWindowQA && sp == spK0Short) {
      histos.fill(HIST("K0Short/h3dPosDCAxy"), centrality, v0.pt(), v0.dcapostopv());
      histos.fill(HIST("K0Short/h3dNegDCAxy"), centrality, v0.pt(), v0.dcanegtopv());
      histos.fill(HIST("K0Short/h3dDCADaughters"), centrality, v0.pt(), v0.dcaV0daughters());
      histos.fill(HIST("K0Short/h3dPointingAngle"), centrality, v0.pt(), TMath::ACos(v0.v0cosPA()));
      histos.fill(HIST("K0Short/h3dV0Radius"), centrality, v0.pt(), v0.v0radius());
    }
    if (std::abs(v0.mLambda() - massRef) < massWindowQA && sp == spLambda) {
      histos.fill(HIST("Lambda/h3dPosDCAxy"), centrality, v0.pt(), v0.dcapostopv());
      histos.fill(HIST("Lambda/h3dNegDCAxy"), centrality, v0.pt(), v0.dcanegtopv());
      histos.fill(HIST("Lambda/h3dDCADaughters"), centrality, v0.pt(), v0.dcaV0daughters());
      histos.fill(HIST("Lambda/h3dPointingAngle"), centrality, v0.pt(), TMath::ACos(v0.v0cosPA()));
      histos.fill(HIST("Lambda/h3dV0Radius"), centrality, v0.pt(), v0.v0radius());
    }
    if (std::abs(v0.mAntiLambda() - massRef) < massWindowQA && sp == spAntiLambda) {
      histos.fill(HIST("AntiLambda/h3dPosDCAxy"), centrality, v0.pt(), v0.dcapostopv());
      histos.fill(HIST("AntiLambda/h3dNegDCAxy"), centrality, v0.pt(), v0.dcanegtopv());
      histos.fill(HIST("AntiLambda/h3dDCADaughters"), centrality, v0.pt(), v0.dcaV0daughters());
      histos.fill(HIST("AntiLambda/h3dPointingAngle"), centrality, v0.pt(), TMath::ACos(v0.v0cosPA()));
      histos.fill(HIST("AntiLambda/h3dV0Radius"), centrality, v0.pt(), v0.v0radius());
    }
  }

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters<dcav0dau && aod::v0data::v0cosPA> v0cospa;

  void process(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels>::iterator const& collision, soa::Filtered<soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>> const& fullV0s, dauTracks const&)
  {
    histos.fill(HIST("hEventSelection"), 0.5 /* all collisions */);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 1.5 /* sel8 collisions */);

    if (std::abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hEventSelection"), 2.5 /* vertex-Z selected */);
    histos.fill(HIST("hEventCentrality"), collision.centFT0C());

    for (auto& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > daughterEtaCut || std::abs(v0.positiveeta()) > daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC

      if (v0.v0radius() > v0radius) {
        // ___________________________________
        // Analysis 1: Lambda
        if (TMath::Abs(v0.yLambda()) < rapidityCut) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 < lifetimecut->get("lifetimecutLambda")) {
            // Lambda daughter PID
            if (compatibleTPC(v0, spLambda)) {
              histos.fill(HIST("h3dMassLambda"), collision.centFT0C(), v0.pt(), v0.mLambda());
              if (doQA)
                fillQAHistograms(v0, spLambda, collision.centFT0C());
            }
            // AntiLambda daughter PID
            if (compatibleTPC(v0, spAntiLambda)) {
              histos.fill(HIST("h3dMassAntiLambda"), collision.centFT0C(), v0.pt(), v0.mAntiLambda());
              if (doQA)
                fillQAHistograms(v0, spAntiLambda, collision.centFT0C());
            }
          }
        }
        // ___________________________________
        // Analysis 2: K0Short
        if (TMath::Abs(v0.yK0Short()) < rapidityCut) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short < lifetimecut->get("lifetimecutK0S")) {
            if (compatibleTPC(v0, spK0Short)) {
              histos.fill(HIST("h3dMassK0Short"), collision.centFT0C(), v0.pt(), v0.mK0Short());
              if (doQA)
                fillQAHistograms(v0, spK0Short, collision.centFT0C());
            }
          }
        }
      } // end radius check
    }   // env V0 loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<derivedlambdakzeroanalysis>(cfgc)};
}
