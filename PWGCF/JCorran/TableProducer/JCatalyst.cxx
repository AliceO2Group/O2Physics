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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>

// centrality
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

////TODO: remove redundant:
#include "Framework/HistogramRegistry.h"

#include "DCAFitter/DCAFitterN.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
////

#include <Math/Vector4D.h>
#include <Math/LorentzVector.h>
#include <TRandom.h>

#include "PWGCF/JCorran/DataModel/JCatalyst.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using namespace ROOT;
using namespace ROOT::Math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::aod
{
namespace jmultiplicity
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float); //! Collision centrality or multiplicity
} // namespace cfmultiplicity
DECLARE_SOA_TABLE(JMultiplicities, "AOD", "JMULTIPLICITY", jmultiplicity::Multiplicity); //! Transient multiplicity table
} // namespace o2::aod

struct JCatalyst {
 public:
  O2_DEFINE_CONFIGURABLE(trigger, int, 0, "Trigger choice: 0 = none, 7 = sel7, 8 = sel8");
  O2_DEFINE_CONFIGURABLE(zvertex, double, 8.0, "Accepted z-vertex range");
  O2_DEFINE_CONFIGURABLE(ptmin, double, 0.2, "Minimal pT for tracks");
  O2_DEFINE_CONFIGURABLE(ptmax, double, 5.0, "Maximal pT for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(charge, int, 0, "Particle charge: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(trackingMode, int, 0, "Tracking mode: 0 = global; 1 = hybrid");
  O2_DEFINE_CONFIGURABLE(collisionFlags, uint16_t, aod::collision::CollisionFlagsRun2::Run2VertexerTracks, "Request collision flags if non-zero: 0 = off, 1 = Run2VertexerTracks")

  Configurable<bool> cutOutliers{"cutOutliers", false, "Cut outlier events"};

  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "CCDB repository URL"};
  O2_DEFINE_CONFIGURABLE(mapCentFlattening, std::string, "", "CCDB path to centrality flattening map");
  //O2_DEFINE_CONFIGURABLE(mapNUACorrection, std::string, "", "CCDB path to NUA correction map");
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < zvertex;
  Filter collisionVertexTypeFilter = (collisionFlags == 0) || ((aod::collision::flags & collisionFlags) == collisionFlags);

  Filter trackFilter = (nabs(aod::track::eta) < etamax) && (aod::track::pt > ptmin) && (aod::track::pt < ptmax);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  TH1D* pcentFlatteningMap = 0;
  //TList* pnuaMapList = 0;

  Produces<aod::JTracks> particleTrack;
  Produces<aod::JCollisions> collisionData;

  void init(InitContext const& ic)
  {
    LOGF(info, "JCatalyst init()");

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    if (!mapCentFlattening.value.empty()) {
      pcentFlatteningMap = ccdb->getForTimeStamp<TH1D>(mapCentFlattening.value, nolaterthan.value);
      if (pcentFlatteningMap)
        LOGF(info, "Centrality flattening enabled. Loaded %s.", mapCentFlattening.value.c_str());
      else
        LOGF(info, "Failed to load centrality flattening histogram %s. Flattening will be disabled.", mapCentFlattening.value.c_str());
    }
    /*if (!mapNUACorrection.value.empty()) {
      pnuaMapList = ccdb->getForTimeStamp<TList>(mapNUACorrection.value, nolaterthan.value);
      if (pnuaMapList)
        LOGF(info, "NUA correction enabled. Loaded %s.", mapNUACorrection.value.c_str());
      else
        LOGF(info, "Failed to load NUA correction catalog %s. Correction will be disabled.", mapNUACorrection.value.c_str());
    }*/
  }

  //void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::BCsWithTimestamps>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection> const& tracks)
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::JMultiplicities>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>> const& tracks)
  {
    switch(trigger){
	case 7:
		if (!collision.alias_bit(kINT7) || !collision.sel7())
		  return;
		break;
	case 8:
	 	if(!collision.sel8())
		  return;
		break;
	default:
		break;
	}

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    collisionData(bc.runNumber(),collision.posZ(),collision.multiplicity());

	for(auto &track : tracks){
	  particleTrack(collisionData.lastIndex(),track.pt(),track.eta(),track.phi(),track.sign());
	}
    /*if (pcentFlatteningMap) {
      Int_t bin = pcentFlatteningMap->GetXaxis()->FindBin(cent[centEst]);
      if (gRandom->Uniform(0, 1) > pcentFlatteningMap->GetBinContent(bin))
        return;
    }*/

    /*TH1* pweightMap = pnuaMapList ? (TH1*)pnuaMapList->FindObject(Form("PhiWeights_%u_%02u", bc.runNumber(), cbin)) : 0;

    for (auto& track : tracks) {
      if (trackingMode == 0 && !track.isGlobalTrack())
        continue;
      else if (trackingMode == 1 && !track.isGlobalTrackSDD())
        continue;

      float pt = track.pt();
      if (pt < ptmin || pt > ptmax)
        continue;

      Int_t ch = track.sign();
      if (charge != 0 && charge * ch < 0)
        continue;

      float eta = track.eta();
      float phi = track.phi();

      float phiWeight = 1.0;
      if (pweightMap) {
        Int_t bin = pweightMap->FindBin(phi, eta, zvertex);
        phiWeight = pweightMap->GetBinContent(bin);
      }

      particleTrack(track.collisionId(), pt, eta, phi, phiWeight, 1.0f);
    }*/
    // LOGF(info,"event %u processed with %u tracks.",collisionId,tracks.size());
  }
};

struct JMultiplicitySelector{
	Produces<aod::JMultiplicities> output;

  O2_DEFINE_CONFIGURABLE(ptmin, float, 0.2f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(etamax, float, 0.8f, "Eta range for tracks")

  Filter trackFilter = (nabs(aod::track::eta) < etamax) && (aod::track::pt > ptmin);
  Filter trackSelection = (requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true);

  void processTracks(aod::Collision const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    output(tracks.size());
  }
  PROCESS_SWITCH(JMultiplicitySelector, processTracks, "Select track count as multiplicity", false);

  void processRun2V0M(aod::CentRun2V0Ms const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centRun2V0M());
    }
  }
  PROCESS_SWITCH(JMultiplicitySelector, processRun2V0M, "Select V0M centrality as multiplicity", false);

  void processRun2CL0(aod::CentRun2CL0s const& centralities)
  {
    for (auto& c : centralities) {
      output(c.centRun2CL0());
    }
  }
  PROCESS_SWITCH(JMultiplicitySelector, processRun2CL0, "Select CL0 centrality as multiplicity", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JCatalyst>(cfgc),
    adaptAnalysisTask<JMultiplicitySelector>(cfgc)};
}
