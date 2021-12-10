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
// Example cascade analysis task
// =============================
//
// This code loops over a CascData table and produces some
// standard analysis output. It requires either
// the cascadefinder or the cascadeproducer tasks
// to have been executed in the workflow (before).
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
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

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

struct cascadeQA {
  //Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},

      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f}}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
    },
  };

  void process(aod::Collision const& collision, aod::CascDataExt const& Cascades)
  {
    for (auto& casc : Cascades) {
      if (casc.sign() < 0) { //FIXME: could be done better...
        registry.fill(HIST("hMassXiMinus"), casc.mXi());
        registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
      } else {
        registry.fill(HIST("hMassXiPlus"), casc.mXi());
        registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
      }
      //The basic eleven!
      registry.fill(HIST("hV0Radius"), casc.v0radius());
      registry.fill(HIST("hCascRadius"), casc.cascradius());
      registry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv());
      registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv());
      registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters());
      registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters());
      registry.fill(HIST("hLambdaMass"), casc.mLambda());
    }
  }
};

struct cascadeanalysis {
  HistogramRegistry registry{
    "registry",
    {
      {"h3dMassXiMinus", "h3dMassXiMinus", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.222f, 1.422f}}}},
      {"h3dMassXiPlus", "h3dMassXiPlus", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.222f, 1.422f}}}},
      {"h3dMassOmegaMinus", "h3dMassOmegaMinus", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.572f, 1.772f}}}},
      {"h3dMassOmegaPlus", "h3dMassOmegaPlus", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.572f, 1.772f}}}},
    },
  };

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.999, "V0 CosPA"};       //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> casccospa{"casccospa", 0.999, "Casc CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcacascdau{"dcacascdau", .3, "DCA Casc Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", .1, "DCA Bach To PV"};
  Configurable<float> dcav0topv{"dcav0topv", .1, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 2.0, "v0radius"};
  Configurable<float> cascradius{"cascradius", 1.0, "cascradius"};
  Configurable<float> v0masswindow{"v0masswindow", 0.008, "v0masswindow"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  Filter preFilterV0 =
    nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > dcabachtopv&& aod::cascdata::dcaV0daughters < dcav0dau&& aod::cascdata::dcacascdaughters < dcacascdau;

  // void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades)
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades)
  {
    if (eventSelection && !collision.sel8()) {
      return;
    }
    for (auto& casc : Cascades) {
      //FIXME: dynamic columns cannot be filtered on?
      if (casc.v0radius() > v0radius &&
          casc.cascradius() > cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa &&
          casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) {
        if (casc.sign() < 0) { //FIXME: could be done better...
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h3dMassXiMinus"), 0., casc.pt(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h3dMassOmegaMinus"), 0., casc.pt(), casc.mOmega());
          }
        } else {
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h3dMassXiPlus"), 0., casc.pt(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h3dMassOmegaPlus"), 0., casc.pt(), casc.mOmega());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(cascadeanalysis, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades)
  {
    if (eventSelection && !collision.alias()[kINT7]) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    for (auto& casc : Cascades) {
      //FIXME: dynamic columns cannot be filtered on?
      if (casc.v0radius() > v0radius &&
          casc.cascradius() > cascradius &&
          casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa &&
          casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) {
        if (casc.sign() < 0) { //FIXME: could be done better...
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h3dMassXiMinus"), collision.centV0M(), casc.pt(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h3dMassOmegaMinus"), collision.centV0M(), casc.pt(), casc.mOmega());
          }
        } else {
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h3dMassXiPlus"), collision.centV0M(), casc.pt(), casc.mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h3dMassOmegaPlus"), collision.centV0M(), casc.pt(), casc.mOmega());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(cascadeanalysis, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeanalysis>(cfgc),
    adaptAnalysisTask<cascadeQA>(cfgc)};
}
