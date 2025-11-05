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
  // This "basic QA" plot a few distributions of topological variables for inspection.
  // it is not meant to be a complete study (yet) - it could be enhanced significantly.
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
    AxisSpec centAxis = {100, 0.0f, 100.0f, "mult percentile"};
    AxisSpec ptAxis = {200, 0.0f, 10.0f, "it{p}_{T} (GeV/c)"};
    AxisSpec massAxisXi = {200, 1.222f, 1.422f, "Inv. Mass (GeV/c^{2})"};
    AxisSpec massAxisOmega = {200, 1.572f, 1.772f, "Inv. Mass (GeV/c^{2})"};

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{10, 0.0f, 10.0f}}});

    // have registrey with 2d histograms (no centrality selection)
    if (!doCentralityStudy) {
      registry.add("h2dMassXiMinus", "h2dMassXiMinus", {HistType::kTH2F, {ptAxis, massAxisXi}});
      registry.add("h2dMassXiPlus", "h2dMassXiPlus", {HistType::kTH2F, {ptAxis, massAxisXi}});
      registry.add("h2dMassOmegaMinus", "h2dMassOmegaMinus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
      registry.add("h2dMassOmegaPlus", "h2dMassOmegaPlus", {HistType::kTH2F, {ptAxis, massAxisOmega}});
    } else {
      registry.add("h3dMassXiMinus", "h3dMassXiMinus", {HistType::kTH3F, {centAxis, ptAxis, massAxisXi}});
      registry.add("h3dMassXiPlus", "h3dMassXiPlus", {HistType::kTH3F, {centAxis, ptAxis, massAxisXi}});
      registry.add("h3dMassOmegaMinus", "h3dMassOmegaMinus", {HistType::kTH3F, {centAxis, ptAxis, massAxisOmega}});
      registry.add("h3dMassOmegaPlus", "h3dMassOmegaPlus", {HistType::kTH3F, {centAxis, ptAxis, massAxisOmega}});
    }
  }

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

  // Track quality and type selections
  Configurable<int> tpcClusters{"tpcClusters", 70, "minimum number of TPC clusters requirement"};
  Configurable<int> itsClusters{"itsClusters", 4, "minimum number of ITS clusters requirement for ITSSA tracks"};
  Configurable<bool> allowITSSAbachelor{"allowITSSAbachelor", true, "allow for bachelor <- cascade track to be via ITS tracking only"};
  Configurable<bool> allowITSSAproton{"allowITSSAproton", true, "allow for proton <- lambda track to be via ITS tracking only"};
  Configurable<bool> allowITSSApion{"allowITSSApion", false, "allow for pion <- lambda track to be via ITS tracking only "};

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
        if (TMath::Abs(casc.yXi()) < 0.5 && lCompatiblePID_Xi) {
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassXiMinus"), casc.pt(), casc.mXi());
          } else {
            registry.fill(HIST("h3dMassXiMinus"), lPercentile, casc.pt(), casc.mXi());
          }
        }
        if (TMath::Abs(casc.yOmega()) < 0.5 && lCompatiblePID_Om) {
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassOmegaMinus"), casc.pt(), casc.mOmega());
          } else {
            registry.fill(HIST("h3dMassOmegaMinus"), lPercentile, casc.pt(), casc.mOmega());
          }
        }
      } else {
        if (TMath::Abs(casc.yXi()) < 0.5 && lCompatiblePID_Xi) {
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassXiPlus"), casc.pt(), casc.mXi());
          } else {
            registry.fill(HIST("h3dMassXiPlus"), lPercentile, casc.pt(), casc.mXi());
          }
        }
        if (TMath::Abs(casc.yOmega()) < 0.5 && lCompatiblePID_Om) {
          if (!doCentralityStudy) {
            registry.fill(HIST("h2dMassOmegaPlus"), casc.pt(), casc.mOmega());
          } else {
            registry.fill(HIST("h3dMassOmegaPlus"), lPercentile, casc.pt(), casc.mOmega());
          }
        }
      }
    }
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtIU const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun3, "Process Run 3 data", true);

  void processRun2(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExt const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun2, "Process Run 2 data", false);

  void processRun3VsMultiplicity(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtIU const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      processCascadeCandidate<FullTracksExtIU>(casc, collision.posX(), collision.posY(), collision.posZ(), collision.centFT0M());
    }
  }
  PROCESS_SWITCH(cascadeAnalysis, processRun3VsMultiplicity, "Process Run 3 data vs multiplicity", false);

  void processRun2VsMultiplicity(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExt const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun2VsMultiplicity, "Process Run 2 data vs multiplicity", false);

  void processRun3WithPID(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtIUWithPID const&)
  // process function subscribing to Run 3-like analysis objects
  {
    // Run 3 event selection criteria
    if (eventSelection && !collision.sel8()) {
      return;
    }
    // fill cascade information with tracksIU typecast (Run 3)
    for (auto& casc : Cascades) {
      int lPIDvalue = checkCascadeTPCPID<FullTracksExtWithPID>(casc);
      processCascadeCandidate<FullTracksExtIUWithPID>(casc, collision.posX(), collision.posY(), collision.posZ(), -999, lPIDvalue);
    }
  }
  PROCESS_SWITCH(cascadeAnalysis, processRun3WithPID, "Process Run 3 data  with PID", false);

  void processRun2WithPID(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtWithPID const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun2WithPID, "Process Run 2 data  with PID", false);

  void processRun3VsMultiplicityWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtIUWithPID const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun3VsMultiplicityWithPID, "Process Run 3 data vs multiplicity with PID", false);

  void processRun2VsMultiplicityWithPID(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, FullTracksExtWithPID const&)
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
  PROCESS_SWITCH(cascadeAnalysis, processRun2VsMultiplicityWithPID, "Process Run 2 data vs multiplicity with PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeAnalysis>(cfgc),
    adaptAnalysisTask<cascadeQa>(cfgc)};
}
