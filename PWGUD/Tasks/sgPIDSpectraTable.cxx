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
//
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/SGTables.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SGPIDSpectraTable {
  Produces<aod::SGEvents> SGevents;
  Produces<aod::SGTracks> SGtracks;
  SGSelector sgSelector;

  // configurables
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", .1, "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> GS_cut{"GS", 0., "Gap-side A=0, C=1, AC = 2, No Gap = -1, All events = 3"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  Configurable<int> occ_cut{"occ_cut", 200, "Maximum Occupancy"};
  Configurable<int> occ_bit1_cut{"occ_bit1_cut", 0, "Check NoCollInTimeRangeStandard"};
  Configurable<int> occ_bit2_cut{"occ_bit2_cut", 0, "Check NoCollInRofStandard"};
  Configurable<int> occ_bit3_cut{"occ_bit3_cut", 0, "Check NoHighMultCollInPrevRof"};
  Configurable<float> ir_cut{"ir_cut", 100, "Maximum IR"};
  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // Collision histograms
  }

  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksPIDExtra, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(coll, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    if (GS_cut != 3) {
      if (truegapSide != GS_cut)
        return;
    }
    // fill collision histograms
    // check occupancies:
    if (occ_bit1_cut && !coll.trs())
      return;
    if (occ_bit2_cut && !coll.trofs())
      return;
    if (occ_bit3_cut && !coll.hmpr())
      return;
    if (coll.occupancyInTime() > occ_cut)
      return;
    if (coll.hadronicRate() > ir_cut)
      return;
    // int truegapSide = sgSelector.trueGap(dgcand, FV0_cut, ZDC_cut);
    // select PV contributors
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // check rho0 signals
    float tpcpi, tpcka, tpcel, tpcpr, tofpi, tofka, tofpr, tofel;
    float tpcde, tofde, tpcmu, tofmu;
    TVector3 a;
    int goodtracks = 0;
    for (auto t : tracks) {
      if (trackselector(t, parameters)) {
        goodtracks++;
      }
    }
    if (!goodtracks)
      return;
    SGevents(coll.runNumber(), coll.flags(), truegapSide, coll.energyCommonZNA(), coll.energyCommonZNC(), goodtracks, coll.occupancyInTime(), coll.hadronicRate());
    // SGevents(coll.runNumber(), coll.flags());
    for (auto t : tracks) {
      if (trackselector(t, parameters)) {
        a.SetXYZ(t.px(), t.py(), t.pz());
        tpcpi = t.hasTPC() ? t.tpcNSigmaPi() : -999;
        tpcmu = t.hasTPC() ? t.tpcNSigmaMu() : -999;
        tpcka = t.hasTPC() ? t.tpcNSigmaKa() : -999;
        tpcpr = t.hasTPC() ? t.tpcNSigmaPr() : -999;
        tpcel = t.hasTPC() ? t.tpcNSigmaEl() : -999;
        tofpi = t.hasTOF() ? t.tofNSigmaPi() : -999;
        tofmu = t.hasTOF() ? t.tofNSigmaMu() : -999;
        tofka = t.hasTOF() ? t.tofNSigmaKa() : -999;
        tofpr = t.hasTOF() ? t.tofNSigmaPr() : -999;
        tofel = t.hasTOF() ? t.tofNSigmaEl() : -999;
        tpcde = t.hasTPC() ? t.tpcNSigmaDe() : -999;
        tofde = t.hasTOF() ? t.tofNSigmaDe() : -999;
        SGtracks(SGevents.lastIndex(), a.Pt(), a.Eta(), a.Phi(), t.sign(), tpcpi, tpcka, tpcpr, tpcel, tofpi, tofka, tofpr, tofel, tpcmu, tofmu, tpcde, tofde);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGPIDSpectraTable>(cfgc, TaskName{"sgpidspectratable"}),
  };
}
