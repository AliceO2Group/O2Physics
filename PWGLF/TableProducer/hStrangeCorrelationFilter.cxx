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
/// \brief This task pre-filters tracks, V0s and cascades to do h-strangeness
///        correlations with an analysis task.
/// \author Kai Cui
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFHStrangeCorrelationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct hstrangecorrelationfilter {
  // Trigger particle selections in phase space
  Configurable<float> triggerEtaMin{"triggerEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> triggerEtaMax{"triggerEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> triggerPtCutMin{"triggerPtCutMin", 1, "triggerptmin"};
  Configurable<float> triggerPtCutMax{"triggerPtCutMax", 3, "triggerptmax"};

  // Track quality
  Configurable<int> minTPCNCrossedRows{"minTPCNCrossedRows", 50, "Minimum TPC crossed rows"};
  Configurable<bool> triggerRequireITS{"triggerRequireITS", true, "require ITS signal in trigger tracks"};

  // Associated particle selections in phase space
  Configurable<float> assocEtaMin{"assocEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> assocEtaMax{"assocEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> assocPtCutMin{"assocPtCutMin", 1, "assocptmin"};
  Configurable<float> assocPtCutMax{"assocPtCutMax", 3, "assocptmax"};

  // V0 selections
  Configurable<double> v0Cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcaV0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcaNegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0RadiusMin{"v0radiusmin", 0.5, "v0radius"};
  Configurable<float> v0RadiusMax{"v0radiusmax", 200, "v0radius"};

  // cascade selections
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.1, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};

  // Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daus only!
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcaPostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > dcaNegtopv&& aod::v0data::dcaV0daughters < dcaV0dau;
  Filter preFilterCascade =
    nabs(aod::cascdata::dcapostopv) > dcaPostopv&& nabs(aod::cascdata::dcanegtopv) > dcaNegtopv&& nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv&& aod::cascdata::dcaV0daughters < dcaV0dau&& aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau;

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {{"correlationHadronHadron", "correlationHadronHadron", {HistType::kTH1F, {{40, -0.5 * M_PI, 1.5 * M_PI, "#Phi"}}}},
     {"correlationHadronV0", "correlationHadronV0", {HistType::kTH1F, {{40, -0.5 * M_PI, 1.5 * M_PI, "#Phi"}}}},
     {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{100, -15., 15.}}}},
     {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{250, 0, 250}}}},
     {"hV0Eta", "hV0Eta", {HistType::kTH1F, {{200, -1, 1, "#Eta"}}}},
     {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -1, 1, "#Eta"}}}},
     {"hTrackSign", "hTrackSign", {HistType::kTH1F, {{5, -2, 2}}}},
     {"hV0dauDCA", "hV0dauDCA", {HistType::kTH1F, {{200, -1, 1}}}},
     {"hID", "hID", {HistType::kTH1F, {{20000, 0, 20000}}}},
     {"hV0CPA", "hV0CPA", {HistType::kTH1F, {{100, 0, 1}}}},
     {"hPosDCAtoPV", "hPosDCAtoPV", {HistType::kTH1F, {{400, 0.05, 0.45}}}},
     {"hNegDCAtoPV", "hNegDCAtoPV", {HistType::kTH1F, {{400, 0.05, 0.45}}}},
     {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.450f, 0.550f}}}},
     {"hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}},
     {"hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}},
     {"hMassXiMinus", "hMassXiMinus", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}},
     {"hMassXiPlus", "hMassXiPlus", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}},
     {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}},
     {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH1F, {{200, 1.016f, 1.216f}}}}}};
  using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;

  Produces<aod::TriggerTracks> triggerTrack;
  Produces<aod::AssocV0s> assocV0;
  Produces<aod::AssocCascades> assocCascades;

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults>::iterator const& collision, DauTracks const& tracks, soa::Filtered<aod::V0Datas> const& V0s, soa::Filtered<aod::CascDatas> const& Cascades, aod::V0sLinked const&)
  {
    // Perform basic event selection
    if (!collision.sel8()) {
      return;
    }
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    // No need to correlate stuff that's in far collisions
    if (TMath::Abs(collision.posZ()) > 10.0) {
      return;
    }

    /// _________________________________________________
    /// Step 1: Populate table with trigger tracks
    for (auto const& track : tracks) {
      if (track.eta() > triggerEtaMax || track.eta() < triggerEtaMin) {
        continue;
      }
      // if (track.sign()= 1 ) {continue;}
      if (track.pt() > triggerPtCutMax || track.pt() < triggerPtCutMin) {
        continue;
      }
      if (track.tpcNClsCrossedRows() < minTPCNCrossedRows) {
        continue; // crossed rows
      }
      if (!track.hasITS() && triggerRequireITS) {
        continue; // skip, doesn't have ITS signal (skips lots of TPC-only!)
      }
      triggerTrack(
        track.collisionId(),
        track.globalIndex());

      registry.fill(HIST("hTrackEta"), track.eta());
      registry.fill(HIST("hTrackSign"), track.sign());
      registry.fill(HIST("hID"), track.collisionId());
    }

    /// _________________________________________________
    /// Step 2: Populate table with associated V0s
    for (auto const& v0 : V0s) {
      if (v0.v0radius() < v0RadiusMin || v0.v0radius() > v0RadiusMax || v0.eta() > assocEtaMax || v0.eta() < assocEtaMin || v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0Cospa) {
        continue;
      }
      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0Eta"), v0.eta());
      registry.fill(HIST("hV0dauDCA"), v0.dcaV0daughters());
      registry.fill(HIST("hPosDCAtoPV"), v0.dcapostopv());
      registry.fill(HIST("hNegDCAtoPV"), v0.dcanegtopv());
      registry.fill(HIST("hV0CPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

      // check dE/dx compatibility
      bool compatibleK0Short = false;
      bool compatibleLambda = false;
      bool compatibleAntiLambda = false;

      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();

      if (negdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posdau.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      if (TMath::Abs(posdau.tpcNSigmaPi()) < 5 && TMath::Abs(negdau.tpcNSigmaPi()) < 5) {
        registry.fill(HIST("hMassK0Short"), v0.mK0Short());
        compatibleK0Short = true;
      }
      if (TMath::Abs(posdau.tpcNSigmaPr()) < 5 && TMath::Abs(negdau.tpcNSigmaPi()) < 5) {
        registry.fill(HIST("hMassLambda"), v0.mLambda());
        compatibleLambda = true;
      }
      if (TMath::Abs(posdau.tpcNSigmaPi()) < 5 && TMath::Abs(negdau.tpcNSigmaPr()) < 5) {
        registry.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
        compatibleAntiLambda = true;
      }
      assocV0(v0.collisionId(), v0.globalIndex(), compatibleK0Short, compatibleLambda, compatibleAntiLambda);
    }

    /// _________________________________________________
    /// Step 3: Populate table with associated Cascades
    for (auto const& casc : Cascades) {
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        return; // skip those cascades for which V0 doesn't exist
      }
      auto v0data = v0.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachTrackCast = casc.bachelor_as<DauTracks>();
      auto posTrackCast = v0data.posTrack_as<DauTracks>();
      auto negTrackCast = v0data.negTrack_as<DauTracks>();

      // minimum TPC crossed rows
      if (bachTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (posTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;
      if (negTrackCast.tpcNClsCrossedRows() < minTPCNCrossedRows)
        continue;

      // check dE/dx compatibility
      bool compatibleXiMinus = false;
      bool compatibleXiPlus = false;
      bool compatibleOmegaMinus = false;
      bool compatibleOmegaPlus = false;

      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < 5) {
        registry.fill(HIST("hMassXiMinus"), casc.mXi());
        compatibleXiMinus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaPi()) < 5) {
        registry.fill(HIST("hMassXiPlus"), casc.mXi());
        compatibleXiPlus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < 5) {
        registry.fill(HIST("hMassOmegaMinus"), casc.mOmega());
        compatibleOmegaMinus = true;
      }
      if (TMath::Abs(posTrackCast.tpcNSigmaPi()) < 5 && TMath::Abs(negTrackCast.tpcNSigmaPr()) < 5 && TMath::Abs(bachTrackCast.tpcNSigmaKa()) < 5) {
        registry.fill(HIST("hMassOmegaPlus"), casc.mOmega());
        compatibleOmegaPlus = true;
      }

      assocCascades(casc.collisionId(), casc.globalIndex(), compatibleXiMinus, compatibleXiPlus, compatibleOmegaMinus, compatibleOmegaPlus);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hstrangecorrelationfilter>(cfgc)};
}
