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

#include <CCDB/BasicCCDBManager.h>
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "tables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

float fReal_fTPCSignalN(float mbb0R, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.019869f * mbb0R) + (0.0012031f * a1pt) + (-0.0031766f * atgl) + (-0.0058023f * atgl * mbb0R) + (0.00087494f * a1pt * mbb0R) + (0.0020074f * side) + (-0.0010434f * a1pt * a1pt) + (0.011812f)) * occ / 1.e3f + //
         ((0.009032f * mbb0R) + (0.0011737f * a1pt) + (-0.0010241f * atgl) + (-0.0075789f * atgl * mbb0R) + (0.00029324f * a1pt * mbb0R) + (0.00052475f * side) + (-0.00045413f * a1pt * a1pt) + (0.0024879f)) * fOccTPCN + //
         ((0.004255f * mbb0R) + (0.0011954f * a1pt) + (0.0054092f * atgl) + (-0.0033655f * atgl * mbb0R) + (0.00052243f * a1pt * mbb0R) + (-0.0002969f * side) + (-0.00074909f * a1pt * a1pt) + (-0.0075754f)) * fTrackOccMeanN + //
         ((-0.07925f * mbb0R) + (-0.03737f * a1pt) + (0.0017054f * atgl) + (0.093686f * atgl * mbb0R) + (0.023925f * a1pt * mbb0R) + (-0.0083407f * side) + (0.00336f * a1pt * a1pt) + (1.0461f));
};

float clamp(float value, float lo, float hi)
{
  return value < lo ? lo : (value > hi ? hi : value);
}

struct ProduceDynamicExtensionCalib {
  Produces<aod::TracksQACorrected> tpcex;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher fetcher;

  Preslice<aod::Tracks> perColl = aod::track::collisionId;

  using BCs = soa::Join<aod::BCs, aod::Timestamps>;
  using Collisions = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;

  int runNumber{0};
  int colId{-100};
  int bcId{-100};
  int trkId{-100};
  Collisions::iterator col;
  BCs::iterator bc;
  Tracks::iterator track;

  void process(Collisions const& collisions, BCs const& bcs, /*aod::FT0s const& ft0s,*/ Tracks const& tracks, aod::TracksQAVersion const& tracksQA)
  {
    col = collisions.begin();
    bc = bcs.begin();
    runNumber = bc.runNumber();
    track = tracks.begin();
    tpcex.reserve(tracksQA.size());
    for (auto& trackqa : tracksQA) {
      if (!trackqa.has_track()) {
        tpcex(0);
        continue;
      }
      if (trackqa.trackId() != trkId) {
        track.setCursor(trackqa.trackId());
      }
      if (!track.has_collision()) {
        tpcex(0);
        continue;
      }
      if (track.collisionId() != colId) {
        colId = track.collisionId();
        col.setCursor(colId);
      }
      if (!col.has_foundBC()) {
        tpcex(0);
        continue;
      }
      if (col.foundBCId() != bcId) {
        bc.setCursor(col.foundBCId());
        if (bc.runNumber() != runNumber) {
          runNumber = bc.runNumber();
        }
      }

      float rate = fetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;
      float occ = col.trackOccupancyInTimeRange();
      float fOccTPCN = clamp(col.multTPC() / 1100.f, 0.f, 12.f);

      float correction0 = fReal_fTPCSignalN(clamp(50.f / track.tpcSignal(), 0.05f, 1.05f), std::abs(track.signed1Pt()), std::abs(track.tgl()), track.tgl() > 0 ? 1.f : 0.f, occ, fOccTPCN, rate / 5.f);
      float correction1 = fReal_fTPCSignalN(clamp(correction0 * 50.f / track.tpcSignal(), 0.05f, 1.05f), std::abs(track.signed1Pt()), std::abs(track.tgl()), track.tgl() > 0 ? 1.f : 0.f, occ, fOccTPCN, rate / 5.);

      tpcex(track.tpcSignal() / correction1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<ProduceDynamicExtensionCalib>(cfgc)};
}
