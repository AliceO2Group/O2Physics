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
//__________________________________________________
// this task provides general links between collisions
// and strange objects reconstructed in various ways.
// It is meant to help with providing auxiliary information
// when dealing with derived data.

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using TracksWithExtra = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe>;

struct strangederivedbuilder {
  //__________________________________________________
  // fundamental building blocks of derived data
  Produces<aod::StraCollision> strangeColl;      // characterises collisions
  Produces<aod::StraCents> strangeCents;         // characterises collisions / centrality
  Produces<aod::V0CollRefs> v0collref;           // references collisions from V0s
  Produces<aod::CascCollRefs> casccollref;       // references collisions from cascades
  Produces<aod::KFCascCollRefs> kfcasccollref;   // references collisions from KF cascades
  Produces<aod::TraCascCollRefs> tracasccollref; // references collisions from tracked cascades

  //__________________________________________________
  // track extra references
  Produces<aod::DauTrackExtras> dauTrackExtras;   // daughter track detector properties
  Produces<aod::DauTrackTPCPIDs> dauTrackTPCPIDs; // daughter track TPC PID 
  Produces<aod::V0Extras> v0Extras;               // references DauTracks from V0s
  Produces<aod::CascExtras> cascExtras;           // references DauTracks from cascades
  Produces<aod::KFCascExtras> kfcascExtras;       // references DauTracks from KF cascades
  Produces<aod::TraCascExtras> tracascExtras;     // references DauTracks from tracked cascades

  //__________________________________________________
  // cascade interlinks
  Produces<aod::CascToTraRefs> cascToTraRefs;   // cascades -> tracked
  Produces<aod::CascToKFRefs> cascToKFRefs;     // cascades -> KF
  Produces<aod::TraToCascRefs> traToCascRefs;   // tracked -> cascades
  Produces<aod::KFToCascRefs> kfToCascRefs;   // KF -> cascades
  

  //__________________________________________________
  // correlation information between cascades: standard<->KF, standard<->tracked, KF<->tracked

  Configurable<bool> fillEmptyCollisions{"fillEmptyCollisions", false, "fill collision entries without candidates"};

  // For manual sliceBy
  Preslice<aod::V0Datas> V0perCollision = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> CascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::KFCascDatas> KFCascperCollision = o2::aod::cascdata::collisionId;
  Preslice<aod::TraCascDatas> TraCascperCollision = o2::aod::cascdata::collisionId;

  void init(InitContext& context)
  {
  }

  void processCollisions(soa::Join<aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As> const& collisions, aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades)
  {
    int currentCollIdx = -1;
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisColl = V0s.sliceBy(V0perCollision, collIdx);
      auto CascTable_thisColl = Cascades.sliceBy(CascperCollision, collIdx);
      auto KFCascTable_thisColl = KFCascades.sliceBy(KFCascperCollision, collIdx);
      auto TraCascTable_thisColl = TraCascades.sliceBy(TraCascperCollision, collIdx);
      bool strange = V0Table_thisColl.size() > 0 ||
                     CascTable_thisColl.size() > 0 ||
                     KFCascTable_thisColl.size() > 0 ||
                     TraCascTable_thisColl.size() > 0;
      // casc table sliced
      if (strange || fillEmptyCollisions) {
        if (currentCollIdx != collIdx) {
          strangeColl(collision.posX(), collision.posY(), collision.posZ());
          strangeCents(collision.centFT0M(), collision.centFT0A(),
                       collision.centFT0C(), collision.centFV0A());
          currentCollIdx = collIdx;
        }
      }
      for (int i = 0; i < V0Table_thisColl.size(); i++)
        v0collref(strangeColl.lastIndex());
      for (int i = 0; i < CascTable_thisColl.size(); i++)
        casccollref(strangeColl.lastIndex());
      for (int i = 0; i < KFCascTable_thisColl.size(); i++)
        kfcasccollref(strangeColl.lastIndex());
      for (int i = 0; i < TraCascTable_thisColl.size(); i++)
        tracasccollref(strangeColl.lastIndex());
    }
  }

  void processTrackExtras(aod::V0Datas const& V0s, aod::CascDatas const& Cascades, aod::KFCascDatas const& KFCascades, aod::TraCascDatas const& TraCascades, TracksWithExtra const& tracksExtra)
  {
    std::vector<int> trackMap(tracksExtra.size(), -1); // index -1: not used

    //__________________________________________________
    // mark tracks that belong to V0s
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
    }

    //__________________________________________________
    // index tracks that belong to CascDatas
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to KFCascDatas
    for (auto const& casc : KFCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // index tracks that belong to TraCascDatas
    for (auto const& casc : TraCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      auto strangeTrack = casc.strangeTrack_as<TracksWithExtra>();
      trackMap[posTrack.globalIndex()] = 0;
      trackMap[negTrack.globalIndex()] = 0;
      trackMap[bachTrack.globalIndex()] = 0;
      trackMap[strangeTrack.globalIndex()] = 0;
    }
    //__________________________________________________
    // Figure out the numbering of the new tracks table
    // assume filling per order
    int nTracks = 0;
    for (int i = 0; i < trackMap.size(); i++) {
      if (trackMap[i] >= 0) {
        trackMap[i] = nTracks++;
      }
    }
    //__________________________________________________
    // populate track references
    for (auto const& v0 : V0s) {
      auto const& posTrack = v0.posTrack_as<TracksWithExtra>();
      auto const& negTrack = v0.negTrack_as<TracksWithExtra>();
      v0Extras(trackMap[posTrack.globalIndex()],
               trackMap[negTrack.globalIndex()]); // joinable with V0Datas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : Cascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      cascExtras(trackMap[posTrack.globalIndex()],
                 trackMap[negTrack.globalIndex()],
                 trackMap[bachTrack.globalIndex()]); // joinable with CascDatas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : KFCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      kfcascExtras(trackMap[posTrack.globalIndex()],
                   trackMap[negTrack.globalIndex()],
                   trackMap[bachTrack.globalIndex()]); // joinable with KFCascDatas
    }
    //__________________________________________________
    // populate track references
    for (auto const& casc : TraCascades) {
      auto bachTrack = casc.bachelor_as<TracksWithExtra>();
      auto v0 = casc.v0();
      auto posTrack = v0.posTrack_as<TracksWithExtra>();
      auto negTrack = v0.negTrack_as<TracksWithExtra>();
      auto strangeTrack = casc.strangeTrack_as<TracksWithExtra>();
      tracascExtras(trackMap[posTrack.globalIndex()],
                    trackMap[negTrack.globalIndex()],
                    trackMap[bachTrack.globalIndex()],
                    trackMap[strangeTrack.globalIndex()]); // joinable with TraCascDatas
    }
    //__________________________________________________
    // circle back and populate actual DauTrackExtra table
    for (auto const& tr : tracksExtra) {
      if (trackMap[tr.globalIndex()] >= 0) {
        dauTrackExtras(tr.detectorMap(), tr.itsClusterSizes(),
                       tr.tpcNClsFound(), tr.tpcNClsCrossedRows());
        dauTrackTPCPIDs(tr.tpcSignal(), tr.tpcNSigmaEl(),
                        tr.tpcNSigmaPi(), tr.tpcNSigmaKa(),
                        tr.tpcNSigmaPr(), tr.tpcNSigmaHe());
      }
    }
    // done!
  }

  using interlinkedCascades = soa::Join<aod::Cascades, aod::CascDataLink, aod::KFCascDataLink, aod::TraCascDataLink>;

  void processCascadeInterlink(interlinkedCascades const& masterCascades, aod::CascIndices const& Cascades, aod::KFCascIndices const& KFCascades, aod::TraCascIndices const& TraCascades)
  {
    //Standard to tracked
    for (auto const& c : Cascades) {
      int indexTracked = -1, indexKF = -1;
      if(c.has_cascade()){
        auto cascade = c.cascade_as<interlinkedCascades>();
        indexTracked = cascade.traCascDataId();
        indexKF = cascade.kfCascDataId();
      }
      cascToTraRefs(indexTracked);
      cascToKFRefs(indexKF);
    }
    //Tracked to standard
    for (auto const& c : TraCascades) {
      int index = -1;
      if(c.has_cascade()){
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      traToCascRefs(index);
    }
    //Tracked to KF
    for (auto const& c : KFCascades) {
      int index = -1;
      if(c.has_cascade()){
        auto cascade = c.cascade_as<interlinkedCascades>();
        index = cascade.cascDataId();
      }
      kfToCascRefs(index);
    }
  }

  PROCESS_SWITCH(strangederivedbuilder, processCollisions, "Produce collisions", true);
  PROCESS_SWITCH(strangederivedbuilder, processTrackExtras, "Produce track extra information", true);
  PROCESS_SWITCH(strangederivedbuilder, processCascadeInterlink, "Produce tables connecting cascades", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangederivedbuilder>(cfgc)};
}
