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
// Strangeness findable study
//
// --- deals with derived data that has been specifically
//     generated to do the findable exercise.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
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
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using recoStraCollisions = soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCents, aod::StraRawCents_003, aod::StraCollLabels>;
using reconstructedV0s = soa::Join<aod::V0CoreMCLabels, aod::V0Cores, aod::V0FoundTags, aod::V0MCCollRefs, aod::V0CollRefs, aod::V0Extras>;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

// simple helper
#define bitset(var, nbit) ((var) |= (static_cast<uint32_t>(1) << static_cast<uint32_t>(nbit)))

struct findableStudy {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master PDG code selection
  Configurable<int> pdgCode{"pdgCode", 310, "PDG code to select"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};

  void init(InitContext const&)
  {
    // Event counting
    histos.add("hCentrality", "hCentrality", kTH1D, {axisCentrality});

    // Duplicate counting
    histos.add("hNRecoV0s", "hNRecoV0s", kTH1D, {{50, -0.5, 49.5f}});
    histos.add("hNRecoV0sWrongColl", "hNRecoV0sWrongColl", kTH1D, {{50, -0.5, 49.5f}});

    // Findable versus found
    histos.add("h2dPtVsCentrality_All_Findable", "hPtVsCentrality_All_Findable", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_All_Found", "hPtVsCentrality_All_Found", kTH2D, {axisCentrality, axisPt});

    // binary feature study: encode in Y axis
    // 0 - positive and negative don't have feature
    // 1 - positive has, negative doesn't
    // 2 - negative has, positive doesn't
    // 3 - both positive and negative have feature
    //
    // all binary feature presence tests are done such that the "has feature" condition
    // is evaluated such that out of any reco V0s, at least one V0 has the desired feature.
    // this is required for generality

    const AxisSpec axisBinaryFeature{static_cast<int>(4), -0.5f, +3.5f, ""};
    histos.add("h2dPtVsCentrality_WithTPC_Findable", "hPtVsCentrality_WithTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithTPC_Found", "hPtVsCentrality_WithTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSTracker_Findable", "hPtVsCentrality_WithITSTracker_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSTracker_Found", "hPtVsCentrality_WithITSTracker_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSAfterburner_Findable", "hPtVsCentrality_WithITSAfterburner_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSAfterburner_Found", "hPtVsCentrality_WithITSAfterburner_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSTrackerTPC_Findable", "hPtVsCentrality_WithITSTrackerTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSTrackerTPC_Found", "hPtVsCentrality_WithITSTrackerTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSABTPC_Findable", "hPtVsCentrality_WithITSABTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h2dPtVsCentrality_WithITSABTPC_Found", "hPtVsCentrality_WithITSABTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
  }

  void processEvents(
    recoStraCollisions::iterator const& collision // reco collisions for collision counting
  )
  {
    histos.fill(HIST("hCentrality"), collision.centFT0C());
  }

  void processV0s(
    aod::V0MCCores::iterator const& v0,               // non-duplicated MC V0 decays
    soa::SmallGroups<reconstructedV0s> const& recv0s, // reconstructed versions of the v0
    recoStraCollisions const&,                        // reco collisions for de-reference
    aod::StraMCCollisions const&,                     // MC collisions for de-reference
    dauTracks const&                                  // daughter track extras
  )
  {
    if (v0.pdgCode() != pdgCode || v0.pdgCodePositive() != 211 || v0.pdgCodeNegative() != -211)
      return;
    if (!v0.isPhysicalPrimary())
      return;

    float rapidity = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassKaonNeutral);
    if (std::abs(rapidity) > 0.5f)
      return;

    double ptmc = std::hypot(v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC());

    // step 1: count number of times this candidate was actually reconstructed
    histos.fill(HIST("hNRecoV0s"), recv0s.size());
    bool hasWrongCollision = false;
    float centrality = 100.5f;
    bool hasBeenFound = false;

    // encode conditionals here
    uint32_t withTPC = 0, withITSTracker = 0, withITSAfterburner = 0, withITSTrackerTPC = 0, withITSABTPC = 0;

    for (auto& recv0 : recv0s) {
      if (recv0.isFound()) {
        hasBeenFound = true;
      }
      auto coll = recv0.straCollision_as<recoStraCollisions>();
      int mcCollID_fromCollision = coll.straMCCollisionId();
      int mcCollID_fromV0 = recv0.straMCCollisionId();
      if (mcCollID_fromCollision != mcCollID_fromV0) {
        hasWrongCollision = true;
      } else {
        // if this is a correctly collision-associated V0, take centrality from here
        // N.B.: this could still be an issue if collision <-> mc collision is imperfect
        centrality = coll.centFT0C();
      }

      // de-reference daughter track extras
      auto pTrack = recv0.posTrackExtra_as<dauTracks>();
      auto nTrack = recv0.negTrackExtra_as<dauTracks>();

      if (pTrack.hasTPC())
        bitset(withTPC, 0);
      if (nTrack.hasTPC())
        bitset(withTPC, 1);

      if (pTrack.hasITS() && pTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTracker, 0);
      if (nTrack.hasITS() && nTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTracker, 1);

      if (pTrack.hasITS() && pTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSAfterburner, 0);
      if (nTrack.hasITS() && nTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSAfterburner, 1);

      if (pTrack.hasTPC() && pTrack.hasITS() && pTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTrackerTPC, 0);
      if (nTrack.hasTPC() && nTrack.hasITS() && nTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTrackerTPC, 1);

      if (pTrack.hasTPC() && pTrack.hasITS() && pTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSABTPC, 0);
      if (nTrack.hasTPC() && nTrack.hasITS() && nTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSABTPC, 1);
    }
    histos.fill(HIST("h2dPtVsCentrality_All_Findable"), centrality, ptmc);
    histos.fill(HIST("h2dPtVsCentrality_WithTPC_Findable"), centrality, ptmc, withTPC);
    histos.fill(HIST("h2dPtVsCentrality_WithITSTracker_Findable"), centrality, ptmc, withITSTracker);
    histos.fill(HIST("h2dPtVsCentrality_WithITSAfterburner_Findable"), centrality, ptmc, withITSAfterburner);
    histos.fill(HIST("h2dPtVsCentrality_WithITSTrackerTPC_Findable"), centrality, ptmc, withITSTrackerTPC);
    histos.fill(HIST("h2dPtVsCentrality_WithITSABTPC_Findable"), centrality, ptmc, withITSABTPC);

    if (hasBeenFound) {
      histos.fill(HIST("h2dPtVsCentrality_All_Found"), centrality, ptmc);
      histos.fill(HIST("h2dPtVsCentrality_WithTPC_Found"), centrality, ptmc, withTPC);
      histos.fill(HIST("h2dPtVsCentrality_WithITSTracker_Found"), centrality, ptmc, withITSTracker);
      histos.fill(HIST("h2dPtVsCentrality_WithITSAfterburner_Found"), centrality, ptmc, withITSAfterburner);
      histos.fill(HIST("h2dPtVsCentrality_WithITSTrackerTPC_Found"), centrality, ptmc, withITSTrackerTPC);
      histos.fill(HIST("h2dPtVsCentrality_WithITSABTPC_Found"), centrality, ptmc, withITSABTPC);
    }

    if (hasWrongCollision) {
      histos.fill(HIST("hNRecoV0sWrongColl"), recv0s.size());
    }
  }

  PROCESS_SWITCH(findableStudy, processEvents, "process collision counters", true);
  PROCESS_SWITCH(findableStudy, processV0s, "process V0s", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<findableStudy>(cfgc)};
}
