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
#include "PWGLF/Utils/v0SelectionBits.h"
#include "PWGLF/Utils/v0SelectionGroup.h"
#include "PWGLF/Utils/v0SelectionTools.h"

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
using reconstructedV0s = soa::Join<aod::V0CoreMCLabels, aod::V0Cores, aod::V0FoundTags, aod::V0MCCollRefs, aod::V0CollRefs, aod::V0Extras, aod::V0TOFPIDs, aod::V0TOFNSigmas>;
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;

// simple checkers, but ensure 64 bit integers
#define bitset(var, nbit) ((var) |= (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint64_t>(1) << static_cast<uint64_t>(nbit)))

struct findableStudy {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master PDG code selection
  Configurable<int> pdgCode{"pdgCode", 310, "PDG code to select"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};
  ConfigurableAxis axisCentrality{"axisCentrality", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "Centrality"};

  // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
  // Full wrapper for configurables related to actual analysis
  Configurable<v0SelectionGroup> v0Selections{"v0Selections", {}, "V0 selection criteria for analysis"};
  // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+

  uint64_t maskTopological;
  uint64_t maskTrackProperties;

  uint64_t maskK0ShortSpecific;
  uint64_t maskLambdaSpecific;
  uint64_t maskAntiLambdaSpecific;

  uint64_t maskSelectionK0Short;
  uint64_t maskSelectionLambda;
  uint64_t maskSelectionAntiLambda;

  void init(InitContext const&)
  {
    v0Selections->PrintSelections(); // for the logs

    v0Selections->provideMasks(maskTopological, maskTrackProperties, maskK0ShortSpecific, maskLambdaSpecific, maskAntiLambdaSpecific);

    // Primary particle selection, central to analysis
    maskSelectionK0Short = maskTopological | maskTrackProperties | maskK0ShortSpecific;
    maskSelectionLambda = maskTopological | maskTrackProperties | maskLambdaSpecific;
    maskSelectionAntiLambda = maskTopological | maskTrackProperties | maskAntiLambdaSpecific;

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
    histos.add("h3dPtVsCentrality_WithTPC_Findable", "h3dPtVsCentrality_WithTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithTPC_Found", "h3dPtVsCentrality_WithTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSTracker_Findable", "h3dPtVsCentrality_WithITSTracker_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSTracker_Found", "h3dPtVsCentrality_WithITSTracker_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSTrackerTPC_Findable", "h3dPtVsCentrality_WithITSTrackerTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSTrackerTPC_Found", "h3dPtVsCentrality_WithITSTrackerTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSABTPC_Findable", "h3dPtVsCentrality_WithITSABTPC_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithITSABTPC_Found", "h3dPtVsCentrality_WithITSABTPC_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithSVertexerOK_Findable", "h3dPtVsCentrality_WithSVertexerOK_Findable", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});
    histos.add("h3dPtVsCentrality_WithSVertexerOK_Found", "h3dPtVsCentrality_WithSVertexerOK_Found", kTH3D, {axisCentrality, axisPt, axisBinaryFeature});

    // Pass topological criteria - N.B. this is cumulative, no need for 2 prong encoding above
    histos.add("h2dPtVsCentrality_Analysis_PassesTrackQuality", "h2dPtVsCentrality_Analysis_PassesTrackQuality", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_PassesTopological", "h2dPtVsCentrality_Analysis_PassesTopological", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_PassesThisSpecies", "h2dPtVsCentrality_Analysis_PassesThisSpecies", kTH2D, {axisCentrality, axisPt});

    // one-on-one test
    histos.add("h2dPtVsCentrality_Analysis_Topo_V0Radius", "h2dPtVsCentrality_Analysis_Topo_V0Radius", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Topo_V0RadiusMax", "h2dPtVsCentrality_Analysis_Topo_V0RadiusMax", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Topo_V0CosPA", "h2dPtVsCentrality_Analysis_Topo_V0CosPA", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Topo_DcaPosToPV", "h2dPtVsCentrality_Analysis_Topo_DcaPosToPV", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Topo_DcaNegToPV", "h2dPtVsCentrality_Analysis_Topo_DcaNegToPV", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Topo_DcaV0Dau", "h2dPtVsCentrality_Analysis_Topo_DcaV0Dau", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Track_TPCRows", "h2dPtVsCentrality_Analysis_Track_TPCRows", kTH2D, {axisCentrality, axisPt});
    histos.add("h2dPtVsCentrality_Analysis_Track_TPCPID", "h2dPtVsCentrality_Analysis_Track_TPCPID", kTH2D, {axisCentrality, axisPt});
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
    int pdgCodePositive = 211;
    int pdgCodeNegative = -211;
    if( pdgCode == 3122 ) 
      pdgCodePositive = 2212;
    if( pdgCode == -3122 ) 
      pdgCodePositive = -2212;
    if( pdgCode == 22 ) {
      pdgCodePositive = -11;
      pdgCodeNegative = -11;
    }

    if (v0.pdgCode() != pdgCode || v0.pdgCodePositive() != pdgCodePositive || v0.pdgCodeNegative() != pdgCodeNegative)
      return;
    if (!v0.isPhysicalPrimary())
      return;

    float rapidity = 2.0;
    if (pdgCode==310)
      rapidity = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassKaonNeutral);
    if (pdgCode==22)
      rapidity = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassPhoton);
    if (pdgCode==3122 || pdgCode==-3122)
      rapidity = RecoDecay::y(std::array{v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC()}, o2::constants::physics::MassLambda0);

    if (std::abs(rapidity) > 0.5f)
      return;

    double ptmc = std::hypot(v0.pxPosMC() + v0.pxNegMC(), v0.pyPosMC() + v0.pyNegMC(), v0.pzPosMC() + v0.pzNegMC());

    // step 1: count number of times this candidate was actually reconstructed
    histos.fill(HIST("hNRecoV0s"), recv0s.size());
    bool hasWrongCollision = false;
    float centrality = 100.5f;
    bool hasBeenFound = false;

    // encode conditionals here
    uint32_t withTPC = 0; // if prongs have TPC 
    uint32_t withITSTracker = 0; // if prongs have been ITS tracked
    uint32_t withITSTrackerTPC = 0; // if prongs have TPC and are ITS tracked
    uint32_t withITSABTPC = 0; // if prongs have TPC and are ITS afterburned
    uint32_t withSVertexerOK = 0; // if prongs have acceptable tracking conditions for svertexer
    
    // Broad
    bool trackPropertiesOK = false;
    bool topologyOK = false;
    bool thisSpeciesOK = false;

    // Detailed
    bool topoV0RadiusOK = false, topoV0RadiusMaxOK = false, topoV0CosPAOK = false, topoDcaPosToPVOK = false, topoDcaNegToPVOK = false, topoDcaV0DauOK = false, trackTPCRowsOK = false, trackTPCPIDOK = false;
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

      if (
        (pTrack.hasTPC() && pTrack.hasITS()) || 
        (!pTrack.hasTPC() && pTrack.itsNCls() >= 6)
        )
        bitset(withSVertexerOK, 0);
      if (
        (nTrack.hasTPC() && nTrack.hasITS()) || 
        (!nTrack.hasTPC() && nTrack.itsNCls() >= 6)
      )
        bitset(withSVertexerOK, 1);

      if (pTrack.hasITS() && pTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTracker, 0);
      if (nTrack.hasITS() && nTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTracker, 1);

      if (pTrack.hasTPC() && pTrack.hasITS() && pTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTrackerTPC, 0);
      if (nTrack.hasTPC() && nTrack.hasITS() && nTrack.itsChi2PerNcl() > -10.0f)
        bitset(withITSTrackerTPC, 1);

      if (pTrack.hasTPC() && pTrack.hasITS() && pTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSABTPC, 0);
      if (nTrack.hasTPC() && nTrack.hasITS() && nTrack.itsChi2PerNcl() < -10.0f)
        bitset(withITSABTPC, 1);

      // determine if this V0 would go to analysis or not
      if( recv0.isFound() ){ 
        uint64_t selMap = v0data::computeReconstructionBitmap(recv0, pTrack, nTrack, coll, recv0.yLambda(), recv0.yK0Short(), v0Selections);
        
        // Consider in all cases
        selMap = selMap | (uint64_t(1) << v0data::selConsiderK0Short) | (uint64_t(1) << v0data::selConsiderLambda) | (uint64_t(1) << v0data::selConsiderAntiLambda);
        
        // global selection checker
        bool validTrackProperties = v0Selections->verifyMask(selMap, maskTrackProperties);
        bool validTopology = v0Selections->verifyMask(selMap, maskTopological);

        uint64_t thisSpeciesMask = maskK0ShortSpecific; 
        if(pdgCode==3122) 
          thisSpeciesMask = maskLambdaSpecific;
        if(pdgCode==-3122) 
          thisSpeciesMask = maskAntiLambdaSpecific;
        // add other species masks as necessary

        bool validThisSpecies = v0Selections->verifyMask(selMap, thisSpeciesMask);
        if ( validTrackProperties ) trackPropertiesOK = true; // stay true even if last recv0 is false 
        if ( validTrackProperties && validTopology ) topologyOK = true; // stay true even if last recv0 is false 
        if ( validTrackProperties && validTopology && validThisSpecies ) thisSpeciesOK = true; // stay true even if last recv0 is false 

        // specific selection (not cumulative)
        topoV0RadiusOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selRadius);
        topoV0RadiusMaxOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selRadiusMax);
        topoV0CosPAOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selCosPA);
        topoDcaPosToPVOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selDCAPosToPV);
        topoDcaNegToPVOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selDCANegToPV);
        topoDcaV0DauOK |= v0Selections->verifyMask(selMap, uint64_t(1) << v0data::selDCAV0Dau);
        trackTPCRowsOK |= v0Selections->verifyMask(selMap, (uint64_t(1) << v0data::selPosGoodTPCTrack) | (uint64_t(1) << v0data::selNegGoodTPCTrack) );

        uint64_t tpcPidMask = (uint64_t(1) << v0data::selTPCPIDPositivePion) | (uint64_t(1) << v0data::selTPCPIDNegativePion); 
        if(pdgCode==3122) 
          tpcPidMask = (uint64_t(1) << v0data::selTPCPIDPositiveProton) | (uint64_t(1) << v0data::selTPCPIDNegativePion); 
        if(pdgCode==-3122) 
          tpcPidMask = (uint64_t(1) << v0data::selTPCPIDPositivePion) | (uint64_t(1) << v0data::selTPCPIDNegativeProton); 

        trackTPCPIDOK |= v0Selections->verifyMask(selMap, tpcPidMask);
      }
    }


    histos.fill(HIST("h2dPtVsCentrality_All_Findable"), centrality, ptmc);
    histos.fill(HIST("h3dPtVsCentrality_WithTPC_Findable"), centrality, ptmc, withTPC);
    histos.fill(HIST("h3dPtVsCentrality_WithITSTracker_Findable"), centrality, ptmc, withITSTracker);
    histos.fill(HIST("h3dPtVsCentrality_WithITSTrackerTPC_Findable"), centrality, ptmc, withITSTrackerTPC);
    histos.fill(HIST("h3dPtVsCentrality_WithITSABTPC_Findable"), centrality, ptmc, withITSABTPC);
    histos.fill(HIST("h3dPtVsCentrality_WithSVertexerOK_Findable"), centrality, ptmc, withSVertexerOK);

    if (hasBeenFound) {
      histos.fill(HIST("h2dPtVsCentrality_All_Found"), centrality, ptmc);
      histos.fill(HIST("h3dPtVsCentrality_WithTPC_Found"), centrality, ptmc, withTPC);
      histos.fill(HIST("h3dPtVsCentrality_WithITSTracker_Found"), centrality, ptmc, withITSTracker);
      histos.fill(HIST("h3dPtVsCentrality_WithITSTrackerTPC_Found"), centrality, ptmc, withITSTrackerTPC);
      histos.fill(HIST("h3dPtVsCentrality_WithITSABTPC_Found"), centrality, ptmc, withITSABTPC);
      histos.fill(HIST("h3dPtVsCentrality_WithSVertexerOK_Found"), centrality, ptmc, withSVertexerOK);
    }

    // broad analysis Level selections -> last axis is 0 not, 1 yes 
    if(trackPropertiesOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_PassesTrackQuality"), centrality, ptmc);
    if(topologyOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_PassesTopological"), centrality, ptmc);
    if(thisSpeciesOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_PassesThisSpecies"), centrality, ptmc);

    // specifics (could be a bit cleaner but ok)
    if(topoV0RadiusOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_V0Radius"), centrality, ptmc);
    if(topoV0RadiusMaxOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_V0RadiusMax"), centrality, ptmc);
    if(topoV0CosPAOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_V0CosPA"), centrality, ptmc);
    if(topoDcaPosToPVOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_DcaPosToPV"), centrality, ptmc);
    if(topoDcaNegToPVOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_DcaNegToPV"), centrality, ptmc);
    if(topoDcaV0DauOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Topo_DcaV0Dau"), centrality, ptmc);
    if(trackTPCRowsOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Track_TPCRows"), centrality, ptmc);
    if(trackTPCPIDOK)
      histos.fill(HIST("h2dPtVsCentrality_Analysis_Track_TPCPID"), centrality, ptmc);

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
