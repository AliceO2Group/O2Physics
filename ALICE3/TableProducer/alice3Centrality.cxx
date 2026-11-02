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
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \author David Dobrigkeit Chinellato <david.dobrigkeit.chinellato@cern.ch>, UNICAMP/CERN
/// \brief Centrality producer for ALICE3
/// \file alice3Centrality.cxx

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3Centrality {
  Produces<aod::CentRun2V0Ms> cent;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", 1, "latest acceptable timestamp of creation for the object"};
  Configurable<float> minEta{"minEta", -4.0f, "Minimum eta in range"};
  Configurable<float> maxEta{"maxEta", 4.0f, "Maximum eta in range"};
  Configurable<float> maxDCA{"maxDCA", 0.0025f, "Max DCAxy and DCAz for counted tracks"};
  Configurable<float> vtxZ{"vtxZ", 10.0f, "Max event vertex z position allowed"};
  Configurable<int> minNumContrib{"minNumContrib", 1, "Minimum required number of primary vertex contributors"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/ALICE3/Centrality", "path to the ccdb object"};

  ConfigurableAxis axisMult{"axisMult", {10000, 0, 10000}, "Reconstructed tracks"};
  ConfigurableAxis axisCent{"axisCent", {150, 0, 150}, "Percentile"};

  Filter trackFilter = (aod::track::eta >= minEta) && (aod::track::eta <= maxEta) && (nabs(aod::track::dcaXY) <= maxDCA) && (nabs(aod::track::dcaZ) <= maxDCA);

  bool centralityLoaded = false;
  TH1D* hCumMultALICE3 = nullptr;

  void init(InitContext&)
  {
    TString tit = Form("%.3f < #it{#eta} < %.3f", minEta.value, maxEta.value);
    histos.add("centrality/numberOfTracks", tit, kTH1D, {axisMult});
    histos.add("centrality/centralityDistribution", "Centrality test", kTH1D, {axisCent});

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(const o2::aod::Collision& collision, const soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA>>& tracks)
  {
    if (!centralityLoaded) {
      hCumMultALICE3 = ccdb->getForTimeStamp<TH1D>(ccdbPath.value, ccdbNoLaterThan.value);
      centralityLoaded = true;
      LOGF(info, "ALICE 3 centrality calibration loaded!");
    }

    if (collision.numContrib() < minNumContrib.value) {
      histos.fill(HIST("centrality/centralityDistribution"), 101);
      cent(101);
      return;
    }
    if (std::fabs(collision.posZ()) > vtxZ.value) {
      histos.fill(HIST("centrality/centralityDistribution"), 102);
      cent(102);
      return;
    }

    const auto nTracks = tracks.size();
    histos.fill(HIST("centrality/numberOfTracks"), nTracks);

    float centALICE3 = hCumMultALICE3->GetBinContent(hCumMultALICE3->FindBin(nTracks));
    histos.fill(HIST("centrality/centralityDistribution"), centALICE3);
    cent(centALICE3);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3Centrality>(cfgc)};
}
