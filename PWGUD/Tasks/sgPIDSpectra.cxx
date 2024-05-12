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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace pid_spectra{
  DECLARE_SOA_COLUMN(GS, gs, int);
  DECLARE_SOA_COLUMN(ZA, za, int);
  DECLARE_SOA_COLUMN(ZC, zc, int);
  DECLARE_SOA_COLUMN(SIGN, sign, std::vector<int>);
  DECLARE_SOA_COLUMN(PX, px, std::vector<float>);
  DECLARE_SOA_COLUMN(PY, py, std::vector<float>);
  DECLARE_SOA_COLUMN(PZ, pz, std::vector<float>);
  DECLARE_SOA_COLUMN(TPCPION, tpcpion, std::vector<int>);
  DECLARE_SOA_COLUMN(TPCKAON, tpckaon, std::vector<int>);
  DECLARE_SOA_COLUMN(TPCPROTON, tpcproton, std::vector<int>);
  DECLARE_SOA_COLUMN(TOFPION, tofpion, std::vector<int>);
  DECLARE_SOA_COLUMN(TOFKAON, tofkaon, std::vector<int>);
  DECLARE_SOA_COLUMN(TOFPROTON, tofproton, std::vector<int>);
}
namespace o2::aod {
  DECLARE_SOA_TABLE(PID_spectra, "AOD", "PID_SPECTRA",
                  pid_spectra::GS, pid_spectra::ZA, pid_spectra::ZC, pid_spectra::SIGN, pid_spectra::PX, pid_spectra::PY, pid_spectra::PZ, pid_spectra::TPCPION, pid_spectra::TPCKAON, pid_spectra::TPCPROTON, pid_spectra::TOFPION, pid_spectra::TOFKAON, pid_spectra::TOFPROTON);
}
struct SGPIDSpectra {
  Produces<o2::aod::PID_spectra> pid_spectra;
  SGSelector sgSelector;

  // configurables
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  //Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // Collision histograms
    registry.add("collisions/GapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("collisions/TrueGapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{4, -1.5, 2.5}}});
  }

  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    // fill collision histograms
    registry.get<TH1>(HIST("collisions/GapSide"))->Fill(coll.gapSide(), 1.);
    // int truegapSide = sgSelector.trueGap(dgcand, FV0_cut, ZDC_cut);
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(coll, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    int gs = truegapSide;
    registry.get<TH1>(HIST("collisions/TrueGapSide"))->Fill(truegapSide, 1.);
    // select PV contributors
    float zna = -1;
    float znc = -1;
    int an = 0;
    int cn = 0;
    if (coll.energyCommonZNC() > 0)
      znc = coll.energyCommonZNC();
    if (coll.energyCommonZNA() > 0)
      zna = coll.energyCommonZNA();
    if (zna > 0 && zna < 4)
      an = 1;
    else if (zna > 4 && zna < 6.8)
      an = 2;
    else if (zna > 6.8 && zna < 10)
      an = 3;
    if (znc > 0 && znc < 4)
      cn = 1;
    else if (znc > 4 && znc < 6.8)
      cn = 2;
    else if (znc > 6.8 && znc < 10)
      cn = 3;

    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // check rho0 signals
    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<int> sign;
    std::vector<int> tpcpion;
    std::vector<int> tpckaon;
    std::vector<int> tpcproton;
    std::vector<int> tofpion;
    std::vector<int> tofkaon;
    std::vector<int> tofproton;
    for (auto t : tracks) {
      if (trackselector(t, parameters)) {
        px.push_back(t.px());
        py.push_back(t.py());
        pz.push_back(t.pz());
        sign.push_back(t.sign());
        if (t.hasTOF()){
	   tofproton.push_back(t.tofNSigmaPr()); 
	   tofkaon.push_back(t.tofNSigmaKa()); 
	   tofpion.push_back(t.tofNSigmaPi()); 
	} else {
           tofproton.push_back(-999.);
           tofkaon.push_back(-999.);
           tofpion.push_back(-999.);
	}
        if (t.hasTPC()){
	   tpcproton.push_back(t.tpcNSigmaPr()); 
	   tpckaon.push_back(t.tpcNSigmaKa()); 
	   tpcpion.push_back(t.tpcNSigmaPi()); 
	} else {
           tpcproton.push_back(-999.);
           tpckaon.push_back(-999.);
           tpcpion.push_back(-999.);
	}
      }
    }
    //Fill Tables here
    pid_spectra(gs, an, cn, sign, px, py, pz, tpcpion, tpckaon, tpcproton, tofpion, tofkaon, tofproton);
 }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGPIDSpectra>(cfgc, TaskName{"sgpidspectra"}),
  };
}
