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

/// \file dptDptFilterQa.cxx
/// \brief basic checks for the behavior of the filter task
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptDptFilter.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define DPTDPTFILTERLOGCOLLISIONS debug
// #define DPTDPTFILTERLOGTRACKS debug

namespace o2::analysis::dptdptfilterqa
{
typedef enum { kRECO = 0,
               kGEN } innerdatatype;
static constexpr std::string_view Dirname[] = {"reconstructed/", "generated/"};
} // namespace o2::analysis::dptdptfilterqa

// Checking the filtered tables
struct DptDptFilterQa {
  Configurable<std::string> cfgDataType{"cfgDataType", "data", "Data type: data, MC, FastMC, OnTheFlyMC. Default data"};
  HistogramRegistry histos{"DptDptFilterQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  o2::analysis::dptdptfilter::DataType datatype;

  template <o2::analysis::dptdptfilterqa::innerdatatype dir>
  void createHistograms()
  {

    using namespace o2::analysis::dptdptfilterqa;

    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksOne").Data(), "Tracks as track one", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksTwo").Data(), "Tracks as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksOneAndTwo").Data(), "Tracks as track one and as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksNone").Data(), "Not selected tracks", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksOneUnsel").Data(), "Tracks as track one", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksTwoUnsel").Data(), "Tracks as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksOneAndTwoUnsel").Data(), "Tracks as track one and as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "TracksNoneUnsel").Data(), "Not selected tracks", kTH1F, {{1500, 0.0, 1500.0}});
    histos.add(TString::Format("%s%s", Dirname[dir].data(), "SelectedEvents").Data(), "Selected events", kTH1F, {{2, 0.0, 2.0}});
    histos.get<TH1>(HIST(Dirname[dir]) + HIST("SelectedEvents"))->GetXaxis()->SetBinLabel(1, "Not selected events");
    histos.get<TH1>(HIST(Dirname[dir]) + HIST("SelectedEvents"))->GetXaxis()->SetBinLabel(2, "Selected events");
  };

  void init(InitContext const&)
  {
    using namespace o2::analysis::dptdptfilter;
    using namespace o2::analysis::dptdptfilterqa;

    switch (getDataType(cfgDataType)) {
      case kData:
        createHistograms<kRECO>();
        break;
      case kMC:
        createHistograms<kRECO>();
        createHistograms<kGEN>();
        break;
      case kFastMC:
      case kOnTheFly:
        createHistograms<kGEN>();
        break;
      default:
        LOGF(fatal, "Data type %s not supported", (std::string)cfgDataType);
        break;
    }
  }

  template <o2::analysis::dptdptfilterqa::innerdatatype dir, typename FilteredCollision, typename FilteredTracks>
  void processQATask(FilteredCollision const& collision,
                     FilteredTracks const& tracks)
  {
    using namespace o2::analysis::dptdptfilterqa;
    static constexpr int kNoOfIclusiveParticles = 2; /* number of inclusive charged particles, aka positive and negative */

    if (collision.collisionaccepted() != uint8_t(true)) {
      histos.fill(HIST(Dirname[dir]) + HIST("SelectedEvents"), 0.5);
    } else {
      histos.fill(HIST(Dirname[dir]) + HIST("SelectedEvents"), 1.5);
    }

    int nTracksOne = 0;
    int nTracksTwo = 0;
    int nTracksOneAndTwo = 0;
    int nTracksNone = 0;
    for (auto const& track : tracks) {
      if (!(track.trackacceptedid() < 0) && !(track.trackacceptedid() < kNoOfIclusiveParticles)) {
        LOGF(fatal, "Task not prepared for identified particles");
      }
      if (track.trackacceptedid() != 0 && track.trackacceptedid() != 1) {
        nTracksNone++;
      }
      if (track.trackacceptedid() == 0) {
        nTracksOne++;
      }
      if (track.trackacceptedid() == 1) {
        nTracksTwo++;
      }
    }
    if (collision.collisionaccepted() != uint8_t(true)) {
      /* control for non selected events */
      histos.fill(HIST(Dirname[dir]) + HIST("TracksOneUnsel"), nTracksOne);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksTwoUnsel"), nTracksTwo);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksNoneUnsel"), nTracksNone);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksOneAndTwoUnsel"), nTracksOneAndTwo);
    } else {
      histos.fill(HIST(Dirname[dir]) + HIST("TracksOne"), nTracksOne);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksTwo"), nTracksTwo);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksNone"), nTracksNone);
      histos.fill(HIST(Dirname[dir]) + HIST("TracksOneAndTwo"), nTracksOneAndTwo);
    }
  }

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyacceptedtracks = (int8_t(0) <= aod::dptdptfilter::trackacceptedid);

  void processGeneratorLevel(soa::Filtered<aod::DptDptCFAcceptedTrueCollisions>::iterator const& collision, soa::Filtered<aod::ScannedTrueTracks> const& tracks)
  {
    using namespace o2::analysis::dptdptfilterqa;
    LOGF(DPTDPTFILTERLOGCOLLISIONS, "New filtered generated collision with BC id %d and with %d accepted tracks", collision.bcId(), tracks.size());
    processQATask<kGEN>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptFilterQa, processGeneratorLevel, "Process generator level filter task QA", true);

  void processDetectorLevel(soa::Filtered<aod::DptDptCFAcceptedCollisions>::iterator const& collision, soa::Filtered<aod::ScannedTracks> const& tracks)
  {
    using namespace o2::analysis::dptdptfilterqa;
    LOGF(DPTDPTFILTERLOGCOLLISIONS, "New filtered collision with BC id %d and with %d accepted tracks", collision.bcId(), tracks.size());
    processQATask<kRECO>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptFilterQa, processDetectorLevel, "Process detector level filter task QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DptDptFilterQa>(cfgc)};
  return workflow;
}
