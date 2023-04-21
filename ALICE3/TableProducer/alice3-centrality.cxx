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

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ALICE3Centrality {
  Produces<aod::CentRun2V0Ms> cent;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> MinEta{"MinEta", -4.0f, "Minimum eta in range"};
  Configurable<float> MaxEta{"MaxEta", 4.0f, "Maximum eta in range"};
  Configurable<float> MaxMult{"MaxMult", 10000.f, "Maximum multiplicity in range"};
  Configurable<float> MaxDCA{"MaxDCA", 0.0025f, "Max DCAxy and DCAz for counted tracks"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  bool centralityLoaded = false;
  TH1D* hCumMultALICE3 = nullptr;

  void init(InitContext&)
  {
    const AxisSpec axisMult{MaxMult.value > 10000.f ? 10000 : (int)MaxMult, 0, MaxMult, "Reconstructed tracks"};
    const AxisSpec axisCent{150, 0, 150, "Percentile"};
    TString tit = Form("%.3f < #it{#eta} < %.3f", MinEta.value, MaxEta.value);
    histos.add("centrality/numberOfTracks", tit, kTH1D, {axisMult});
    histos.add("centrality/centralityDistribution", "Centrality test", kTH1D, {axisCent});

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  int nevs = 0;
  void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks)
  {
    if (!centralityLoaded) {
      hCumMultALICE3 = ccdb->getForTimeStamp<TH1D>("Analysis/ALICE3/Centrality", -1);
      centralityLoaded = true;
      LOGF(info, "ALICE 3 centrality calibration loaded!");
    }

    int nTracks = 0;
    if (collision.numContrib() < 1) {
      histos.fill(HIST("centrality/centralityDistribution"), 101);
      cent(101);
      return;
    }
    if (fabs(collision.posZ()) > 10) {
      histos.fill(HIST("centrality/centralityDistribution"), 102);
      cent(102);
      return;
    }
    for (const auto& track : tracks) {
      if (track.eta() < MinEta || track.eta() > MaxEta) {
        continue;
      }
      if (abs(track.dcaXY()) > MaxDCA || abs(track.dcaZ()) > MaxDCA) {
        continue;
      }
      nTracks++;
    }
    LOG(info) << nevs++ << ") Event " << collision.globalIndex() << " has " << nTracks << " tracks";
    histos.fill(HIST("centrality/numberOfTracks"), nTracks);

    float centALICE3 = hCumMultALICE3->GetBinContent(hCumMultALICE3->FindBin(nTracks));
    histos.fill(HIST("centrality/centralityDistribution"), centALICE3);
    cent(centALICE3);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ALICE3Centrality>(cfgc, TaskName{"alice3-centrality"})};
}
