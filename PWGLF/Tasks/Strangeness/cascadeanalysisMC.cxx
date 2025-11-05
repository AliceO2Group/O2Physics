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
// this task makes use of labels provided by the cascade
// builder to process cascades quickly.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

struct cascadeGenerated {
  // Basic generated histograms
  HistogramRegistry registry{
    "registry",
    {},
  };

  Configurable<float> maxPt{"maxPt", 20.0, "max generated pT"};
  Configurable<int> nPtBins{"nPtBins", 200, "number of pT bins"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "max (absolute) rapidity of generated cascade"};
  Configurable<int> nRapidityBins{"nRapidityBins", 200, "number of pT bins"};
  Configurable<bool> requirePhysicalPrimary{"requirePhysicalPrimary", false, "require the generated cascade to be a physical primary"};

  void init(InitContext const&)
  {
    AxisSpec ptAxis = {nPtBins, 0.0f, maxPt, "#it{p}_{T} (GeV/c)"};
    AxisSpec rapidityAxis = {nRapidityBins, -rapidityCut, rapidityCut, "y"};

    registry.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});
    registry.add("hPtXiMinus", "hPtXiMinus", {HistType::kTH1F, {ptAxis}});
    registry.add("hPtXiPlus", "hPtXiPlus", {HistType::kTH1F, {ptAxis}});
    registry.add("hPtOmegaMinus", "hPtOmegaMinus", {HistType::kTH1F, {ptAxis}});
    registry.add("hPtOmegaPlus", "hPtOmegaPlus", {HistType::kTH1F, {ptAxis}});
    registry.add("h2DXiMinus", "h2DXiMinus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("h2DXiPlus", "h2DXiPlus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("h2DOmegaMinus", "h2DOmegaMinus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("h2DOmegaPlus", "h2DOmegaPlus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
  }

  void process(aod::McCollision const& /*collision*/, aod::McParticles const& mcparts)
  {
    // Count monte carlo events
    // WARNING: MC collision <-> real collision association has to be understood
    registry.fill(HIST("hEventCounter"), 0.5f);

    // Count all generated MC particles
    // WARNING: event-level losses have to be understood too
    for (auto& particle : mcparts) {
      if (TMath::Abs(particle.y()) > rapidityCut)
        continue;
      if (requirePhysicalPrimary && !particle.isPhysicalPrimary())
        continue;
      if (particle.pdgCode() == 3312) {
        registry.fill(HIST("hPtXiMinus"), particle.pt());
        registry.fill(HIST("h2DXiMinus"), particle.pt(), particle.y());
      }
      if (particle.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlus"), particle.pt());
        registry.fill(HIST("h2DXiPlus"), particle.pt(), particle.y());
      }
      if (particle.pdgCode() == 3334) {
        registry.fill(HIST("hPtOmegaMinus"), particle.pt());
        registry.fill(HIST("h2DOmegaMinus"), particle.pt(), particle.y());
      }
      if (particle.pdgCode() == -3334) {
        registry.fill(HIST("hPtOmegaPlus"), particle.pt());
        registry.fill(HIST("h2DOmegaPlus"), particle.pt(), particle.y());
      }
    }
  }
};

struct cascadeQaMC {
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

  Configurable<bool> assocMC{"assocMC", true, "fill histograms only for MC associated candidates"};

  void process(aod::Collision const& collision, LabeledCascades const& Cascades, aod::McParticles const&)
  {
    for (auto& casc : Cascades) {
      Int_t lPDG = 0;
      if (assocMC) {
        if (!casc.has_mcParticle())
          continue;
        auto cascmc = casc.mcParticle();
        if (TMath::Abs(cascmc.pdgCode()) == 3312 || TMath::Abs(cascmc.pdgCode()) == 3334)
          lPDG = cascmc.pdgCode();
      }

      if (casc.sign() < 0) { // FIXME: could be done better...
        if (!assocMC || lPDG == 3312)
          registry.fill(HIST("hMassXiMinus"), casc.mXi());
        if (!assocMC || lPDG == 3334)
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
      } else {
        if (!assocMC || lPDG == -3312)
          registry.fill(HIST("hMassXiPlus"), casc.mXi());
        if (!assocMC || lPDG == -3334)
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

struct cascadeAnalysisMC {
  HistogramRegistry registry{
    "registry",
    {},
  };

  Configurable<float> rapidityCut{"rapidityCut", 0.5, "max (absolute) rapidity of cascade candidate"};
  Configurable<int> nRapidityBins{"nRapidityBins", 200, "number of pT bins"};

  void init(InitContext const&)
  {
    AxisSpec ptAxis = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"};
    AxisSpec rapidityAxis = {nRapidityBins, -rapidityCut, rapidityCut, "y"};
    AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
    AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});
    registry.add("h2dMassXiMinus", "h2dMassXiMinus", {HistType::kTH2F, {ptAxis, massAxisXi}});
    registry.add("h2dMassXiPlus", "h2dMassXiPlus", {HistType::kTH2F, {ptAxis, massAxisXi}});
    registry.add("h2dMassOmegaMinus", "h2dMassOmegaMinus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
    registry.add("h2dMassOmegaPlus", "h2dMassOmegaPlus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
    registry.add("hPtYMassXiMinus", "hPtYMassXiMinus", {HistType::kTH3F, {ptAxis, rapidityAxis, massAxisXi}});
    registry.add("hPtYMassXiPlus", "hPtYMassXiPlus", {HistType::kTH3F, {ptAxis, rapidityAxis, massAxisXi}});
    registry.add("hPtYMassOmegaMinus", "hPtYMassOmegaMinus", {HistType::kTH3F, {ptAxis, rapidityAxis, massAxisOmega}});
    registry.add("hPtYMassOmegaPlus", "hPtYMassOmegaPlus", {HistType::kTH3F, {ptAxis, rapidityAxis, massAxisOmega}});
  }

  // Selection criteria
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Selection criteria - compatible with core wagon autodetect
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.95, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};

  Configurable<int> tpcClusters{"tpcClusters", 70, "minimum number of TPC clusters requirement"};
  Configurable<int> itsClusters{"itsClusters", 4, "minimum number of ITS clusters requirement for ITSSA tracks"};
  Configurable<bool> allowITSSAbachelor{"allowITSSAbachelor", true, "allow for bachelor <- cascade track to be via ITS tracking only"};
  Configurable<bool> allowITSSAproton{"allowITSSAproton", true, "allow for proton <- lambda track to be via ITS tracking only"};
  Configurable<bool> allowITSSApion{"allowITSSApion", false, "allow for pion <- lambda track to be via ITS tracking only "};
  Configurable<bool> assocMC{"assocMC", true, "fill histograms only for MC associated candidates"};

  // Track identification configurables
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 4, "TPC NSigma bachelor (>10 is no cut)"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 4, "TPC NSigma proton <- lambda (>10 is no cut)"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 4, "TPC NSigma pion <- lambda (>10 is no cut)"};

  // Switch for centrality
  Configurable<bool> doCentralityStudy{"doCentralityStudy", false, "do centrality percentile selection (yes/no)"};

  Filter preFilter =
    nabs(aod::cascdata::dcapostopv) > v0setting_dcapostopv&& nabs(aod::cascdata::dcanegtopv) > v0setting_dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < v0setting_dcav0dau&& aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau;

  template <class TCascTracksTo, typename TCascade>
  int checkCascadeTPCPID(TCascade& lCascade)
  // function to check PID of a certain cascade candidate for a hypothesis
  {
    bool lConsistentWithLambda = true;
    bool lConsistentWithXi = true;
    bool lConsistentWithOm = true;

    auto bachTrack = lCascade.template bachelor_as<TCascTracksTo>();
    auto posTrack = lCascade.template posTrack_as<TCascTracksTo>();
    auto negTrack = lCascade.template negTrack_as<TCascTracksTo>();

    // Bachelor: depends on type
    if (TMath::Abs(bachTrack.tpcNSigmaPi()) > tpcNsigmaBachelor && tpcNsigmaBachelor < 9.99)
      lConsistentWithXi = false;
    if (TMath::Abs(bachTrack.tpcNSigmaKa()) > tpcNsigmaBachelor && tpcNsigmaBachelor < 9.99)
      lConsistentWithOm = false;

    // Proton check: depends on cascade sign
    if (lCascade.sign() < 0 && TMath::Abs(posTrack.tpcNSigmaPr()) > tpcNsigmaProton && tpcNsigmaProton < 9.99)
      lConsistentWithLambda = false;
    if (lCascade.sign() > 0 && TMath::Abs(negTrack.tpcNSigmaPr()) > tpcNsigmaProton && tpcNsigmaProton < 9.99)
      lConsistentWithLambda = false;

    // Pion check: depends on cascade sign
    if (lCascade.sign() < 0 && TMath::Abs(negTrack.tpcNSigmaPi()) > tpcNsigmaPion && tpcNsigmaPion < 9.99)
      lConsistentWithLambda = false;
    if (lCascade.sign() > 0 && TMath::Abs(posTrack.tpcNSigmaPi()) > tpcNsigmaPion && tpcNsigmaPion < 9.99)
      lConsistentWithLambda = false;

    // bit-packing (first bit -> consistent with Xi, second bit -> consistent with Omega)
    return lConsistentWithLambda * lConsistentWithXi + 2 * lConsistentWithLambda * lConsistentWithOm;
  }

  template <class TCascTracksTo, typename TCascade>
  void processCascadeCandidate(TCascade const& casc, float const& pvx, float const& pvy, float const& pvz, float lPercentile = 999.0f, int lPIDvalue = 3)
  // function to process cascades and generate corresponding invariant mass distributions
  {
    registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates

    // check mc association if requested
    Int_t lPDG = 0;
    if (assocMC) {
      if (!casc.has_mcParticle())
        return;
      auto cascmc = casc.mcParticle();
      if (TMath::Abs(cascmc.pdgCode()) == 3312 || TMath::Abs(cascmc.pdgCode()) == 3334)
        lPDG = cascmc.pdgCode();
    }

    registry.fill(HIST("hCandidateCounter"), 1.5); // v0data exists - deprecated
    auto bachTrackCast = casc.template bachelor_as<TCascTracksTo>();
    auto posTrackCast = casc.template posTrack_as<TCascTracksTo>();
    auto negTrackCast = casc.template negTrack_as<TCascTracksTo>();

    // track-level selections
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

    // Logic: either you have enough TPC clusters, OR you enabled ITSSA and have enough ITS clusters as requested
    // N.B.: This will require dedicated studies!

    Bool_t lGoodCandidate = kFALSE;
    if (casc.sign() < 0) {
      if (
        (lEnoughTPCNClsBac || (allowITSSAbachelor && lEnoughITSNClsBac)) && // bachelor conditional
        (lEnoughTPCNClsPos || (allowITSSAproton && lEnoughITSNClsPos)) &&   // bachelor conditional
        (lEnoughTPCNClsNeg || (allowITSSApion && lEnoughITSNClsNeg))        // bachelor conditional
      ) {
        lGoodCandidate = kTRUE;
      }
    }
    if (casc.sign() > 0) {
      if (
        (lEnoughTPCNClsBac || (allowITSSAbachelor && lEnoughITSNClsBac)) && // bachelor conditional
        (lEnoughTPCNClsPos || (allowITSSApion && lEnoughITSNClsPos)) &&     // bachelor conditional
        (lEnoughTPCNClsNeg || (allowITSSAproton && lEnoughITSNClsNeg))      // bachelor conditional
      ) {
        lGoodCandidate = kTRUE;
      }
    }
    if (!lGoodCandidate)
      return;
    registry.fill(HIST("hCandidateCounter"), 2.5); // okay track quality

    // assign TPC PID compatibility booleans
    bool lCompatiblePID_Xi = (lPIDvalue >> 0 & 1);
    bool lCompatiblePID_Om = (lPIDvalue >> 1 & 1);

    if (casc.v0radius() > v0setting_radius &&
        casc.cascradius() > cascadesetting_cascradius &&
        casc.v0cosPA(pvx, pvy, pvz) > v0setting_cospa &&
        casc.casccosPA(pvx, pvy, pvz) > cascadesetting_cospa &&
        casc.dcav0topv(pvx, pvy, pvz) > cascadesetting_mindcav0topv &&
        TMath::Abs(casc.mLambda() - 1.115683) < cascadesetting_v0masswindow) {
      registry.fill(HIST("hCandidateCounter"), 3.5); // pass cascade selections
      if (casc.sign() < 0) {                         // FIXME: could be done better...
        if (TMath::Abs(casc.yXi()) < rapidityCut && lCompatiblePID_Xi && ((!assocMC) || (lPDG == 3312))) {
          registry.fill(HIST("hPtYMassXiMinus"), casc.pt(), casc.yXi(), casc.mXi());
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassXiMinus"), casc.pt(), casc.mXi());
          } else {
            registry.fill(HIST("h3dMassXiMinus"), lPercentile, casc.pt(), casc.mXi());
          }
        }
        if (TMath::Abs(casc.yOmega()) < rapidityCut && lCompatiblePID_Om && ((!assocMC) || (lPDG == 3334))) {
          registry.fill(HIST("hPtYMassOmegaMinus"), casc.pt(), casc.yOmega(), casc.mOmega());
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassOmegaMinus"), casc.pt(), casc.mOmega());
          } else {
            registry.fill(HIST("h3dMassOmegaMinus"), lPercentile, casc.pt(), casc.mOmega());
          }
        }
      } else {
        if (TMath::Abs(casc.yXi()) < rapidityCut && lCompatiblePID_Xi && ((!assocMC) || (lPDG == -3312))) {
          registry.fill(HIST("hPtYMassXiPlus"), casc.pt(), casc.yXi(), casc.mXi());
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassXiPlus"), casc.pt(), casc.mXi());
          } else {
            registry.fill(HIST("h3dMassXiPlus"), lPercentile, casc.pt(), casc.mXi());
          }
        }
        if (TMath::Abs(casc.yOmega()) < rapidityCut && lCompatiblePID_Om && ((!assocMC) || (lPDG == -3334))) {
          registry.fill(HIST("hPtYMassOmegaPlus"), casc.pt(), casc.yOmega(), casc.mOmega());
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassOmegaPlus"), casc.pt(), casc.mOmega());
          } else {
            registry.fill(HIST("h3dMassOmegaPlus"), lPercentile, casc.pt(), casc.mOmega());
          }
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtIU const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      processCascadeCandidate<FullTracksExtIU>(casc, collision.posX(), collision.posY(), collision.posZ());
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExt const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 2 event selection criteria
    if (eventSelection && !collision.alias_bit(kINT7)) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    // fill cascade information with tracks typecast (Run 2)
    for (auto& casc : Cascades) {
      processCascadeCandidate<FullTracksExt>(casc, collision.posX(), collision.posY(), collision.posZ());
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun2, "Process Run 2 data", false);

  void processRun3VsMultiplicity(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtIU const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      processCascadeCandidate<FullTracksExtIU>(casc, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M());
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun3VsMultiplicity, "Process Run 3 data vs multiplicity", false);

  void processRun2VsMultiplicity(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExt const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 2 event selection criteria
    if (eventSelection && !collision.alias_bit(kINT7)) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    // fill cascade information with tracks typecast (Run 2)
    for (auto& casc : Cascades) {
      processCascadeCandidate<FullTracksExt>(casc, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M());
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun2VsMultiplicity, "Process Run 2 data vs multiplicity", false);

  void processRun3WithPID(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtIUWithPID const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      int lPIDvalue = checkCascadeTPCPID<FullTracksExtIUWithPID>(casc);
      processCascadeCandidate<FullTracksExtIUWithPID>(casc, collision.posX(), collision.posY(), collision.posZ(), -999, lPIDvalue);
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun3WithPID, "Process Run 3 data  with PID", false);

  void processRun2WithPID(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtWithPID const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 2 event selection criteria
    if (eventSelection && !collision.alias_bit(kINT7)) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    // fill cascade information with tracks typecast (Run 2)
    for (auto& casc : Cascades) {
      int lPIDvalue = checkCascadeTPCPID<FullTracksExtWithPID>(casc);
      processCascadeCandidate<FullTracksExtWithPID>(casc, collision.posX(), collision.posY(), collision.posZ(), -999, lPIDvalue);
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun2WithPID, "Process Run 2 data  with PID", false);

  void processRun3VsMultiplicityWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtIUWithPID const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      int lPIDvalue = checkCascadeTPCPID<FullTracksExtIUWithPID>(casc);
      processCascadeCandidate<FullTracksExtIUWithPID>(casc, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), lPIDvalue);
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun3VsMultiplicityWithPID, "Process Run 3 data vs multiplicity with PID", false);

  void processRun2VsMultiplicityWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<LabeledCascades> const& Cascades, FullTracksExtWithPID const&, aod::McParticles const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 2 event selection criteria
    if (eventSelection && !collision.alias_bit(kINT7)) {
      return;
    }
    if (eventSelection && !collision.sel7()) {
      return;
    }
    // fill cascade information with tracks typecast (Run 2)
    for (auto& casc : Cascades) {
      int lPIDvalue = checkCascadeTPCPID<FullTracksExtWithPID>(casc);
      processCascadeCandidate<FullTracksExtWithPID>(casc, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), lPIDvalue);
    }
  }
  PROCESS_SWITCH(cascadeAnalysisMC, processRun2VsMultiplicityWithPID, "Process Run 2 data vs multiplicity with PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeGenerated>(cfgc),
    adaptAnalysisTask<cascadeAnalysisMC>(cfgc),
    adaptAnalysisTask<cascadeQaMC>(cfgc)};
}
