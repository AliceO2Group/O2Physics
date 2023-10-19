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
// 2-particle correlations for cascades
// =============================
//
// Author: Rik Spijkers (rik.spijkers@cern.ch)
//

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

// Add a column to the cascdataext table: IsSelected.
// 0 = not selected, 1 = Xi, 2 = Omega, 3 = both
namespace o2::aod
{
namespace cascadeflags
{
DECLARE_SOA_COLUMN(IsSelected, isSelected, int); //~!
} // namespace cascadeflags
DECLARE_SOA_TABLE(CascadeFlags, "AOD", "CASCADEFLAGS", //!
                  cascadeflags::IsSelected);
using CascDataExtSelected = soa::Join<CascDataExt, CascadeFlags>;
} // namespace o2::aod

struct cascadeSelector {
  Produces<aod::CascadeFlags> cascflags;

  // Configurables
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 3, "TPC NSigma bachelor (>10 is no cut)"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 3, "TPC NSigma proton <- lambda (>10 is no cut)"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 3, "TPC NSigma pion <- lambda (>10 is no cut)"};

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::CascDataExt const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIUWithPID const&)
  {
    for (auto& casc : Cascades) {
      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        cascflags(0);
        continue; // reject if no v0data
      }
      auto v0data = v0.v0Data();

      // TODO: inv mass selection, other cuts.

      // Let's try to do some PID
      // these are the tracks:
      auto bachTrack = casc.bachelor_as<FullTracksExtIUWithPID>();
      auto posTrack = v0data.posTrack_as<FullTracksExtIUWithPID>();
      auto negTrack = v0data.negTrack_as<FullTracksExtIUWithPID>();

      // PID
      // Lambda check
      if (casc.sign() < 0) {
        // Proton check:
        if (TMath::Abs(posTrack.tpcNSigmaPr()) > tpcNsigmaProton && tpcNsigmaProton < 9.99) {
          cascflags(0);
          continue;
        }
        // Pion check:
        if (TMath::Abs(negTrack.tpcNSigmaPi()) > tpcNsigmaPion) {
          cascflags(0);
          continue;
        }
      } else {
        // Proton check:
        if (TMath::Abs(negTrack.tpcNSigmaPr()) > tpcNsigmaProton) {
          cascflags(0);
          continue;
        }
        // Pion check:
        if (TMath::Abs(posTrack.tpcNSigmaPi()) > tpcNsigmaPion) {
          cascflags(0);
          continue;
        }
      }
      // Bachelor check
      if (TMath::Abs(bachTrack.tpcNSigmaPi()) < tpcNsigmaBachelor) {
        if (TMath::Abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor) {
          // consistent with both!
          cascflags(3);
          continue;
        }
        cascflags(1);
        continue;
      } else if (TMath::Abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor) {
        cascflags(2);
        continue;
      }
      // if we reach here, the bachelor was neither pion nor kaon
      cascflags(0);
    } // cascade loop
  }   // process
};    // struct

struct cascadeCorrelations {
  AxisSpec invMassAxis = {3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"};
  AxisSpec deltaPhiAxis = {100, -PI / 2, 1.5 * PI, "#Delta#varphi"};
  AxisSpec deltaEtaAxis = {40, -2, 2, "#Delta#eta"};
  AxisSpec ptAxis = {200, 0, 15, "#it{p}_{T}"};
  AxisSpec selectionFlagAxis = {4, -0.05f, 3.5f, "Selection flag of casc candidate"};
  AxisSpec vertexAxis = {1000, -10.0f, 10.0f, "cm"};

  HistogramRegistry registry{
    "registry",
    {
      // inv mass
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH2F, {invMassAxis, ptAxis}}},

      // basic selection variables
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{1000, 0.95f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hSelectionFlag", "hSelectionFlag", {HistType::kTH1I, {selectionFlagAxis}}},
      {"hAutoCorrelation", "hAutoCorrelation", {HistType::kTH1I, {{4, -0.05f, 3.5f, "Types of autocorrelation"}}}},
      {"hPhi", "hPhi", {HistType::kTH1F, {{100, 0, 2 * PI, "#varphi"}}}},
      {"hEta", "hEta", {HistType::kTH1F, {{100, -2, 2, "#eta"}}}},

      // correlation histos
      {"hDeltaPhiSS", "hDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"hDeltaPhiOS", "hDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},
      // THnSparses containing all relevant dimensions, to be extended with e.g. multiplicity
      {"hSparseOS", "hSparseOS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaEtaAxis, ptAxis, ptAxis, selectionFlagAxis, selectionFlagAxis, vertexAxis}}},
      {"hSparseSS", "hSparseSS", {HistType::kTHnSparseF, {deltaPhiAxis, deltaEtaAxis, ptAxis, ptAxis, selectionFlagAxis, selectionFlagAxis, vertexAxis}}},
    },
  };

  Filter Selector = aod::cascadeflags::isSelected > 0;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::CascDataExtSelected> const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, FullTracksExtIU const&)
  {
    // Some QA on the cascades
    for (auto& casc : Cascades) {

      auto v0 = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0.has_v0Data())) {
        continue; // reject if no v0data
      }

      if (casc.isSelected() != 2) { // not exclusively an Omega --> consistent with Xi or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt());
        }
      }
      if (casc.isSelected() >= 2) { // consistent with Omega or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt());
        }
      }
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

      registry.fill(HIST("hSelectionFlag"), casc.isSelected());
      registry.fill(HIST("hPhi"), casc.phi());
      registry.fill(HIST("hEta"), casc.eta());
    } // casc loop

    for (auto& [c0, c1] : combinations(Cascades, Cascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT. As we can't swap the cascade tables themselves, we swap the addresses and later dereference them
      auto* triggerAddress = &c0;
      auto* assocAddress = &c1;
      if (assocAddress->pt() > triggerAddress->pt()) {
        std::swap(triggerAddress, assocAddress);
      }
      auto trigger = *triggerAddress;
      auto assoc = *assocAddress;

      auto lambdaTrigg = trigger.v0_as<o2::aod::V0sLinked>();
      auto lambdaAssoc = assoc.v0_as<o2::aod::V0sLinked>();
      if (!(lambdaTrigg.has_v0Data()) || !(lambdaAssoc.has_v0Data())) {
        continue; // reject if no v0data in either of the lambda's
      }
      auto v0dataTrigg = lambdaTrigg.v0Data();
      auto v0dataAssoc = lambdaAssoc.v0Data();

      // calculate angular correlations
      double deta = trigger.eta() - assoc.eta();
      double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -0.5 * PI);

      // Fill the correct histograms based on same-sign or opposite-sign
      if (trigger.sign() * assoc.sign() < 0) { // opposite-sign
        registry.fill(HIST("hDeltaPhiOS"), dphi);
        registry.fill(HIST("hSparseOS"), dphi, deta, trigger.pt(), assoc.pt(), trigger.isSelected(), assoc.isSelected(), collision.posZ());
      } else { // same-sign
        // make sure to check for autocorrelations - only possible in same-sign correlations
        if (v0dataTrigg.v0Id() == v0dataAssoc.v0Id()) {
          // LOGF(info, "same v0 in SS correlation! %d %d", v0dataTrigg.v0Id(), v0dataAssoc.v0Id());
          registry.fill(HIST("hAutoCorrelation"), 0);
          continue;
        }
        int bachIdTrigg = trigger.bachelorId();
        int bachIdAssoc = assoc.bachelorId();
        int posIdTrigg = v0dataTrigg.posTrackId();
        int negIdTrigg = v0dataTrigg.negTrackId();
        int posIdAssoc = v0dataAssoc.posTrackId();
        int negIdAssoc = v0dataAssoc.negTrackId();
        if (bachIdTrigg == bachIdAssoc) {
          // LOGF(info, "same bachelor in SS correlation! %d %d", bachIdTrigg, bachIdAssoc);
          registry.fill(HIST("hAutoCorrelation"), 1);
          continue;
        }
        // check for same tracks in v0's of cascades
        if (negIdTrigg == negIdAssoc || posIdTrigg == posIdAssoc) {
          // LOGF(info, "cascades have a v0-track in common in SS correlation!");
          registry.fill(HIST("hAutoCorrelation"), 2);
          continue;
        }
        if (trigger.sign() < 0) { // neg cascade
          if (negIdTrigg == bachIdAssoc || negIdAssoc == bachIdTrigg) {
            // LOGF(info, "bach of casc == v0-pion of other casc in neg SS correlation!");
            registry.fill(HIST("hAutoCorrelation"), 3);
            continue;
          }
        } else { // pos cascade
          if (posIdTrigg == bachIdAssoc || posIdAssoc == bachIdTrigg) {
            // LOGF(info, "bach of casc == v0-pion of other casc in pos SS correlation!");
            registry.fill(HIST("hAutoCorrelation"), 3);
            continue;
          }
        }
        registry.fill(HIST("hDeltaPhiSS"), dphi);
        registry.fill(HIST("hSparseSS"), dphi, deta, trigger.pt(), assoc.pt(), trigger.isSelected(), assoc.isSelected(), collision.posZ());
      }
    } // correlations
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascadeSelector>(cfgc),
    adaptAnalysisTask<cascadeCorrelations>(cfgc)};
}
