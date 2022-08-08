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
// the cascade finder or the cascade builder tasks
// to have been executed in the workflow (before).
//
// Note, however, that only the cascade builder does
// time-frame-aware cascade building at this time. This
// is deliberate, as time-frame information is available
// during reconstruction and collision association is better
// done at that stage.
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

//use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

struct cascadeQa {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2}²)"}}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "Inv. Mass (GeV/c^{2})"}}}},
    },
  };

  void process(aod::Collision const& collision, aod::CascDataExt const& Cascades)
  {
    for (auto& casc : Cascades) {
      if (casc.sign() < 0) { // FIXME: could be done better...
        registry.fill(HIST("hMassXiMinus"), casc.mXi());
        registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
      } else {
        registry.fill(HIST("hMassXiPlus"), casc.mXi());
        registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
      }
      // The basic eleven!
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

struct cascadeAnalysis {
  HistogramRegistry registry{
    "registry",
    {},
  };

  void init(InitContext const&)
  {
    AxisSpec ptAxis = {200, 0.0f, 10.0f, "it{p}_{T} (GeV/c)"};
    AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
    AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});
    registry.add("h2dMassXiMinus", "h2dMassXiMinus", {HistType::kTH2F, {ptAxis, massAxisXi}});
    registry.add("h2dMassXiPlus", "h2dMassXiPlus", {HistType::kTH2F, {ptAxis, massAxisXi}});
    registry.add("h2dMassOmegaMinus", "h2dMassOmegaMinus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
    registry.add("h2dMassOmegaPlus", "h2dMassOmegaPlus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
  }

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"};       // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<double> casccospa{"casccospa", 0.95, "Casc CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA Casc Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .05, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .05, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", .05, "DCA Bach To PV"};
  Configurable<float> dcav0topv{"dcav0topv", .05, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 0.9, "v0radius"};
  Configurable<float> cascradius{"cascradius", 0.5, "cascradius"};
  Configurable<float> v0masswindow{"v0masswindow", 0.008, "v0masswindow"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  Configurable<int> tpcClusters{"tpcClusters", 70, "minimum number of TPC clusters requirement"};
  Configurable<int> itsClusters{"itsClusters", 4, "minimum number of ITS clusters requirement for ITSSA tracks"};
  Configurable<bool> allowITSSAbachelor{"allowITSSAbachelor", true, "allow for bachelor <- cascade track to be via ITS tracking only"};
  Configurable<bool> allowITSSAproton{"allowITSSAproton", true, "allow for proton <- lambda track to be via ITS tracking only"};
  Configurable<bool> allowITSSApion{"allowITSSApion", false, "allow for pion <- lambda track to be via ITS tracking only "};

  Filter preFilter =
    nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > dcabachtopv&& aod::cascdata::dcaV0daughters < dcav0dau&& aod::cascdata::dcacascdaughters < dcacascdau;

  Partition<LabeledCascades> associatedCascades = (aod::mccasclabel::mcParticleId > -1);

  template <class TCascTracksTo>
  void fillCascadeOutput(soa::Filtered<aod::CascDataExt> const& cascades, float const& pvx, float const& pvy, float const& pvz)
  //function to process cascades and generate corresponding invariant mass distributions
  {
    for (auto& casc : cascades) {
      registry.fill(HIST("hCandidateCounter"), 0.5); //all candidates
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        continue; //skip those cascades for which V0 doesn't exist
      }
      registry.fill(HIST("hCandidateCounter"), 1.5); //v0data exists
      auto v0data = v0.v0Data();                     // de-reference index to correct v0data in case it exists
      auto bachTrackCast = casc.bachelor_as<TCascTracksTo>();
      auto posTrackCast = v0data.posTrack_as<TCascTracksTo>();
      auto negTrackCast = v0data.negTrack_as<TCascTracksTo>();

      //track-level selections
      Bool_t lEnoughTPCNClsBac = kTRUE;
      Bool_t lEnoughTPCNClsPos = kTRUE;
      Bool_t lEnoughTPCNClsNeg = kTRUE;
      Bool_t lEnoughITSNClsBac = kTRUE;
      Bool_t lEnoughITSNClsPos = kTRUE;
      Bool_t lEnoughITSNClsNeg = kTRUE;

      if (bachTrackCast.tpcNClsFound() < tpcClusters)
        lEnoughTPCNClsBac = kFALSE;
      if (posTrackCast.tpcNClsFound() < tpcClusters)
        lEnoughTPCNClsPos = kFALSE;
      if (negTrackCast.tpcNClsFound() < tpcClusters)
        lEnoughTPCNClsNeg = kFALSE;
      if (bachTrackCast.itsNCls() < itsClusters)
        lEnoughITSNClsBac = kFALSE;
      if (posTrackCast.itsNCls() < itsClusters)
        lEnoughITSNClsPos = kFALSE;
      if (negTrackCast.itsNCls() < itsClusters)
        lEnoughITSNClsNeg = kFALSE;

      //Logic: either you have enough TPC clusters, OR you enabled ITSSA and have enough ITS clusters as requested
      //N.B.: This will require dedicated studies!

      Bool_t lGoodCandidate = kFALSE;
      if (casc.sign() < 0) {
        if (
          (lEnoughTPCNClsBac || (allowITSSAbachelor && lEnoughITSNClsBac)) && //bachelor conditional
          (lEnoughTPCNClsPos || (allowITSSAproton && lEnoughITSNClsPos)) &&   //bachelor conditional
          (lEnoughTPCNClsNeg || (allowITSSApion && lEnoughITSNClsNeg))        //bachelor conditional
        ) {
          lGoodCandidate = kTRUE;
        }
      }
      if (casc.sign() > 0) {
        if (
          (lEnoughTPCNClsBac || (allowITSSAbachelor && lEnoughITSNClsBac)) && //bachelor conditional
          (lEnoughTPCNClsPos || (allowITSSApion && lEnoughITSNClsPos)) &&     //bachelor conditional
          (lEnoughTPCNClsNeg || (allowITSSAproton && lEnoughITSNClsNeg))      //bachelor conditional
        ) {
          lGoodCandidate = kTRUE;
        }
      }
      if (!lGoodCandidate)
        continue;
      registry.fill(HIST("hCandidateCounter"), 2.5); //okay track quality

      if (casc.v0radius() > v0radius &&
          casc.cascradius() > cascradius &&
          casc.v0cosPA(pvx, pvy, pvz) > v0cospa &&
          casc.casccosPA(pvx, pvy, pvz) > casccospa &&
          casc.dcav0topv(pvx, pvy, pvz) > dcav0topv &&
          TMath::Abs(casc.mLambda() - 1.115683) < v0masswindow) {
        registry.fill(HIST("hCandidateCounter"), 3.5); //pass cascade selections
        if (casc.sign() < 0) {                         // FIXME: could be done better...
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h2dMassXiMinus"), casc.template pt(), casc.template mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h2dMassOmegaMinus"), casc.template pt(), casc.template mOmega());
          }
        } else {
          if (TMath::Abs(casc.yXi()) < 0.5) {
            registry.fill(HIST("h2dMassXiPlus"), casc.template pt(), casc.template mXi());
          }
          if (TMath::Abs(casc.yOmega()) < 0.5) {
            registry.fill(HIST("h2dMassOmegaPlus"), casc.template pt(), casc.template mOmega());
          }
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIU const&)
  //process function subscribing to Run 3-like analysis objects
  {
    //Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    //fill cascade information with tracksIU typecast (Run 3)
    fillCascadeOutput<FullTracksExtIU>(Cascades, collision.posX(), collision.posY(), collision.posZ());
  }
  PROCESS_SWITCH(cascadeAnalysis, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExt const&)
  //process function subscribing to Run 3-like analysis objects
  {
    //Run 2 event selection criteria
    if (eventSelection && !collision.alias()[kINT7]) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    //fill cascade information with tracks typecast (Run 2)
    fillCascadeOutput<FullTracksExt>(Cascades, collision.posX(), collision.posY(), collision.posZ());
  }
  PROCESS_SWITCH(cascadeAnalysis, processRun2, "Process Run 2 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeAnalysis>(cfgc),
    adaptAnalysisTask<cascadeQa>(cfgc)};
}
