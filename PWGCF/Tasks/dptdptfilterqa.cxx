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

#include <cmath>

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define DPTDPTFILTERLOGCOLLISIONS debug
#define DPTDPTFILTERLOGTRACKS debug

namespace o2::analysis::dptdptfilterqa
{
typedef enum { kRECO = 0,
               kGEN } innerdatatype;
static constexpr std::string_view dirname[] = {"reconstructed/", "generated/"};
} // namespace o2::analysis::dptdptfilterqa

// Checking the filtered tables
struct DptDptFilterQA {
  Configurable<std::string> cfgDataType{"datatype", "data", "Data type: data, MC, FastMC, OnTheFlyMC. Default data"};
  HistogramRegistry histos{"DptDptFilterQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  o2::analysis::dptdptfilter::DataType datatype;

  template <o2::analysis::dptdptfilterqa::innerdatatype dir>
  void createHistograms()
  {

    using namespace o2::analysis::dptdptfilterqa;

    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksOne").Data(), "Tracks as track one", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksTwo").Data(), "Tracks as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksOneAndTwo").Data(), "Tracks as track one and as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksNone").Data(), "Not selected tracks", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksOneUnsel").Data(), "Tracks as track one", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksTwoUnsel").Data(), "Tracks as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksOneAndTwoUnsel").Data(), "Tracks as track one and as track two", kTH1F, {{1500, 0.0, 1500.0, "number of tracks"}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "TracksNoneUnsel").Data(), "Not selected tracks", kTH1F, {{1500, 0.0, 1500.0}});
    histos.add(TString::Format("%s%s", dirname[dir].data(), "SelectedEvents").Data(), "Selected events", kTH1F, {{2, 0.0, 2.0}});
    histos.get<TH1>(HIST(dirname[dir]) + HIST("SelectedEvents"))->GetXaxis()->SetBinLabel(1, "Not selected events");
    histos.get<TH1>(HIST(dirname[dir]) + HIST("SelectedEvents"))->GetXaxis()->SetBinLabel(2, "Selected events");
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

    if (collision.collisionaccepted() != uint8_t(true)) {
      histos.fill(HIST(dirname[dir]) + HIST("SelectedEvents"), 0.5);
    } else {
      histos.fill(HIST(dirname[dir]) + HIST("SelectedEvents"), 1.5);
    }

    int ntracks_one = 0;
    int ntracks_two = 0;
    int ntracks_one_and_two = 0;
    int ntracks_none = 0;
    for (auto& track : tracks) {
      if (!(track.trackacceptedid() < 0) && !(track.trackacceptedid() < 2)) {
        LOGF(fatal, "Task not prepared for identified particles");
      }
      if (track.trackacceptedid() != 0 && track.trackacceptedid() != 1) {
        ntracks_none++;
      }
      if (track.trackacceptedid() == 0) {
        ntracks_one++;
      }
      if (track.trackacceptedid() == 1) {
        ntracks_two++;
      }
    }
    if (collision.collisionaccepted() != uint8_t(true)) {
      /* control for non selected events */
      histos.fill(HIST(dirname[dir]) + HIST("TracksOneUnsel"), ntracks_one);
      histos.fill(HIST(dirname[dir]) + HIST("TracksTwoUnsel"), ntracks_two);
      histos.fill(HIST(dirname[dir]) + HIST("TracksNoneUnsel"), ntracks_none);
      histos.fill(HIST(dirname[dir]) + HIST("TracksOneAndTwoUnsel"), ntracks_one_and_two);
    } else {
      histos.fill(HIST(dirname[dir]) + HIST("TracksOne"), ntracks_one);
      histos.fill(HIST(dirname[dir]) + HIST("TracksTwo"), ntracks_two);
      histos.fill(HIST(dirname[dir]) + HIST("TracksNone"), ntracks_none);
      histos.fill(HIST(dirname[dir]) + HIST("TracksOneAndTwo"), ntracks_one_and_two);
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
  PROCESS_SWITCH(DptDptFilterQA, processGeneratorLevel, "Process generator level filter task QA", true);

  void processDetectorLevel(soa::Filtered<aod::DptDptCFAcceptedCollisions>::iterator const& collision, soa::Filtered<aod::ScannedTracks> const& tracks)
  {
    using namespace o2::analysis::dptdptfilterqa;
    LOGF(DPTDPTFILTERLOGCOLLISIONS, "New filtered collision with BC id %d and with %d accepted tracks", collision.bcId(), tracks.size());
    processQATask<kRECO>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptFilterQA, processDetectorLevel, "Process detector level filter task QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DptDptFilterQA>(cfgc)};
  return workflow;
}
