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

struct JCatalyst {
 public:
  O2_DEFINE_CONFIGURABLE(zvertex, double, 8.0, "Accepted z-vertex range");
  O2_DEFINE_CONFIGURABLE(ptmin, double, 0.2, "Minimal pT for tracks");
  O2_DEFINE_CONFIGURABLE(ptmax, double, 5.0, "Maximal pT for tracks");
  O2_DEFINE_CONFIGURABLE(charge, int, 0, "Particle charge: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(trackingMode, int, 0, "Tracking mode: 0 = global; 1 = hybrid");
  O2_DEFINE_CONFIGURABLE(centEst, int, 0, "Centrality estimator: 0 = V0M; 1 = CL0; 2 = CL1");
  Configurable<bool> cutOutliers{"cutOutliers", false, "Cut outlier events"};

  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "CCDB repository URL"};
  O2_DEFINE_CONFIGURABLE(mapCentFlattening, std::string, "", "CCDB path to centrality flattening map");
  O2_DEFINE_CONFIGURABLE(mapNUACorrection, std::string, "", "CCDB path to NUA correction map");
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  TH1D* pcentFlatteningMap = 0;
  TList* pnuaMapList = 0;

  Produces<aod::ParticleTrack> particleTrack;
  Produces<aod::CollisionData> collisionData;

  Int_t GetCentBin(Double_t cent)
  {
    for (Int_t i = 1, bins = sizeof(jflucCentBins) / sizeof(jflucCentBins[0]);
         i < bins; ++i) {
      if (cent < jflucCentBins[i])
        return i - 1;
    }
    return -1;
  }

  void init(InitContext const& ic)
  {
    LOGF(info, "JCatalyst init()");

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // gRandom->SetSeed(15122022); // for centrality flattening

    if (!mapCentFlattening.value.empty()) {
      pcentFlatteningMap = ccdb->getForTimeStamp<TH1D>(mapCentFlattening.value, nolaterthan.value);
      if (pcentFlatteningMap)
        LOGF(info, "Centrality flattening enabled. Loaded %s.", mapCentFlattening.value.c_str());
      else
        LOGF(info, "Failed to load centrality flattening histogram %s. Flattening will be disabled.", mapCentFlattening.value.c_str());
    }
    if (!mapNUACorrection.value.empty()) {
      pnuaMapList = ccdb->getForTimeStamp<TList>(mapNUACorrection.value, nolaterthan.value);
      if (pnuaMapList)
        LOGF(info, "NUA correction enabled. Loaded %s.", mapNUACorrection.value.c_str());
      else
        LOGF(info, "Failed to load NUA correction catalog %s. Correction will be disabled.", mapNUACorrection.value.c_str());
    }
  }

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection> const& tracks)
  {
    Double_t cent[3] = {
      collision.centRun2V0M(),
      collision.centRun2CL0(),
      collision.centRun2CL1()};
    Int_t cbin = GetCentBin(cent[centEst]);

    collisionData(collision.globalIndex(), cent[centEst], cbin);

    if (cbin < 0)
      return;
    if (!collision.alias_bit(kINT7) || !collision.sel7())
      return;
    if (std::abs(collision.posZ()) > zvertex)
      return;

    if (pcentFlatteningMap) {
      Int_t bin = pcentFlatteningMap->GetXaxis()->FindBin(cent[centEst]);
      if (gRandom->Uniform(0, 1) > pcentFlatteningMap->GetBinContent(bin))
        return;
    }

    // TODO: outlier cutting

    if (cutOutliers.value) {
      double centCL0 = collision.centRun2CL0();
      double center = 0.973488 * centCL0 + 0.0157497;
      double sigma = 0.673612 + centCL0 * (0.0290718 + centCL0 * (-0.000546728 + centCL0 * 5.82749e-06));
      if (cent[0] < center - 5.0 * sigma || cent[0] > center + 5.5 * sigma || cent[0] < 0.0 || cent[0] > 60.0)
        return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    TH1* pweightMap = pnuaMapList ? (TH1*)pnuaMapList->FindObject(Form("PhiWeights_%u_%02u", bc.runNumber(), cbin)) : 0;

    for (auto& track : tracks) {
      if (trackingMode == 0 && !track.isGlobalTrack())
        continue;
      else if (trackingMode == 1 && !track.isGlobalTrackSDD())
        continue;

      Double_t pt = track.pt();
      if (pt < ptmin || pt > ptmax)
        continue;

      Int_t ch = track.sign();
      if (charge != 0 && charge * ch < 0)
        continue;

      Double_t eta = track.eta();
      Double_t phi = track.phi();

      Double_t phiWeight = 1.0;
      if (pweightMap) {
        Int_t bin = pweightMap->FindBin(phi, eta, zvertex);
        phiWeight = pweightMap->GetBinContent(bin);
      }

      particleTrack(track.collisionId(), pt, eta, phi, phiWeight, 1.0f);
    }
    // LOGF(info,"event %u processed with %u tracks.",collisionId,tracks.size());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JCatalyst>(cfgc)};
}
