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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

#include "DataModel/UDDerived.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TH1D.h>

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UPCSpectraProviderTask {

  Produces<aod::UDTracks> outputTracks;

  Filter trackFilter = (requireGlobalTrackInFilter());

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& col, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>> const& tracks)
  {
    bool checkV0 = col.selection_bit(kIsBBV0A) || col.selection_bit(kIsBBV0C) || !col.selection_bit(kNoBGV0A) || !col.selection_bit(kNoBGV0C);
    if (checkV0) {
      return;
    }
    bool checkFDD = col.selection_bit(kIsBBFDA) || col.selection_bit(kIsBBFDC) || !col.selection_bit(kNoBGFDA) || !col.selection_bit(kNoBGFDC);
    if (checkFDD) {
      return;
    }
    if (!col.alias_bit(kCUP9)) {
      return;
    }
    if (tracks.size() != 2) {
      return;
    }
    auto track1 = tracks.begin();
    auto track2 = track1 + 1;
    if (track1.sign() * track2.sign() >= 0) {
      return;
    }
    UChar_t clustermap1 = track1.itsClusterMap();
    UChar_t clustermap2 = track2.itsClusterMap();
    bool checkClusMap = TESTBIT(clustermap1, 0) && TESTBIT(clustermap1, 1) && TESTBIT(clustermap2, 0) && TESTBIT(clustermap2, 1);
    if (!checkClusMap) {
      return;
    }
    /*
    if ((p.Pt() >= 0.1) || (signalTPC1 + signalTPC2 > 140.)) {
      return;
    }
    */
    outputTracks(track1.pt(), track1.eta(), track1.phi(), track1.tpcSignal(),
                 track2.pt(), track2.eta(), track2.phi(), track2.tpcSignal());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCSpectraProviderTask>(cfgc, TaskName{"upcspectra-task-skim-provider"})};
}
