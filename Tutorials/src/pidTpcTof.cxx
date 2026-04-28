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

///
/// \file   pidTpcTof.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to test the PID features and utilities
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

// O2Physics includes
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

using namespace o2;
using namespace o2::aod::pidutils;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the response table
struct pidTpcTof {
  void init(o2::framework::InitContext&)
  {
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                         aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                         aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                         aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                         aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl>;
  // Define slice per collision
  void process(Trks::iterator const& track)
  {
    // Access TPC and TOF nSigma statically for each particle mass hypothesis
    static_for<0, 8>([&](auto i) {
      tofNSigma<i>(track);
      tofExpSigma<i>(track);
      tofExpSignal<i>(track);
      tofExpSignalDiff<i>(track);

      tpcNSigma<i>(track);
      tpcExpSigma<i>(track);
      tpcExpSignal<i>(track);
      tpcExpSignalDiff<i>(track);
    });

    // Access TPC and TOF nSigma dinamically for each particle mass hypothesis
    for (int i = 0; i < 9; i++) {
      tofNSigma(i, track);
      tofExpSigma(i, track);
      tofExpSignal(i, track);
      tofExpSignalDiff(i, track);

      tpcNSigma(i, track);
      tpcExpSigma(i, track);
      tpcExpSignal(i, track);
      tpcExpSignalDiff(i, track);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<pidTpcTof>(cfgc)};
}
