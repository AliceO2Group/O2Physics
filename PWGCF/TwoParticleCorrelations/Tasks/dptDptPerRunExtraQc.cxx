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

/// \file dptDptPerRunExtraQc.cxx
/// \brief basic per run check of the per analyzed species p vs TPC IW momentum
/// \author victor.gonzalez.sebastian@gmail.com

#include <array>
#include <cmath>
#include <unordered_map>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

namespace perrunextraqc
{
std::unordered_map<int, TProfile3D*> gRunMapPvsTpcIwP;
TProfile3D* gCurrentRunPvsPtcIwP;
} // namespace perrunextraqc

struct DptDptPerRunExtraQc {
  int mRunNumber{-1};
  AxisSpec qaPAxis{150, 0.1, 5.0};

  HistogramRegistry mHistos{"PerRunExtraQcHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  void initRunNumber(aod::BCsWithTimestamps::iterator const& bc)
  {
    using namespace perrunextraqc;

    if (mRunNumber == bc.runNumber()) {
      return;
    } else {
      mRunNumber = bc.runNumber();
      if (gRunMapPvsTpcIwP.find(mRunNumber) == gRunMapPvsTpcIwP.end()) {
        gRunMapPvsTpcIwP[mRunNumber] = mHistos.add<TProfile3D>(TString::Format("Reco/%d_pVsTpcIwP", mRunNumber).Data(), ";species;p (GeV/#it{c}); p_{tpciw} (GeV/#it{c})", {HistType::kTProfile3D, {{10, -0.5, 9.5}, qaPAxis, qaPAxis}}).get();
      }
      gCurrentRunPvsPtcIwP = gRunMapPvsTpcIwP[mRunNumber];
    }
  }

  void init(InitContext&)
  {
    using namespace perrunextraqc;

    qaPAxis.makeLogarithmic();
  }

  template <typename PassedTracks>
  void processTracks(PassedTracks const& tracks)
  {
    using namespace perrunextraqc;

    for (const auto& track : tracks) {
      gCurrentRunPvsPtcIwP->Fill(track.trackacceptedid(), track.p(), track.tpcInnerParam(), track.pt());
    }
  }

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyacceptedtracks = (aod::dptdptfilter::trackacceptedid >= int8_t(0));

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision, soa::Filtered<soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo>> const& tracks, aod::BCsWithTimestamps const&)
  {
    using namespace analysis::dptdptfilter;

    if (!collision.collisionaccepted()) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initRunNumber(bc);
    processTracks(tracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DptDptPerRunExtraQc>(cfgc)};
}
