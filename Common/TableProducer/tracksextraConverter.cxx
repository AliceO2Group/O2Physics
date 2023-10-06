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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;

// Converts TracksExtra table from version 000 to 001
struct tracksextraConverter {
  Produces<aod::StoredTracksExtra_001> tracksExtra_001;
  void process(aod::TracksExtra_000 const& tracksExtra_000)
  {
    // dummy itsClusterSizes, table filled with maximum uint32_t value
    uint32_t itsClusterSizes = 0xFFFFFFFF;

    for (auto& track0 : tracksExtra_000) {
      tracksExtra_001(track0.tpcInnerParam(),
                      track0.flags(),
                      track0.itsClusterMap(),
                      itsClusterSizes,
                      track0.tpcNClsFindable(),
                      track0.tpcNClsFindableMinusFound(),
                      track0.tpcNClsFindableMinusCrossedRows(),
                      track0.tpcNClsShared(),
                      track0.trdPattern(),
                      track0.itsChi2NCl(),
                      track0.tpcChi2NCl(),
                      track0.trdChi2(),
                      track0.tofChi2(),
                      track0.tpcSignal(),
                      track0.trdSignal(),
                      track0.length(),
                      track0.tofExpMom(),
                      track0.trackEtaEmcal(),
                      track0.trackPhiEmcal(),
                      track0.trackTime(),
                      track0.trackTimeRes());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<tracksextraConverter>(cfgc),
  };
}
