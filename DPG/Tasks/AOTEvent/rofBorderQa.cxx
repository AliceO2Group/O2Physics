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
/// \brief QA task to study ROF border effect for different event, track and particle selection parameters
/// \author Igor Altsybeev, Igor.Altsybeev@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPECSObject.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TMathBase.h>
#include <TString.h>

#include <RtypesCore.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// main class
struct RofBorderQaTask {
  // for vertex vs time:
  bool flagShowInfo = false;
  int lastRunNumber = -1;
  int nBCsPerOrbit = 3564;

  // bc position correlations
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int64_t bcSOR = -1; // global bc of the start of the first orbit
  int64_t orbitSOR = -1;
  int64_t nBCsPerTF = 1; // 128*3564; // duration of TF in bcs
  uint32_t nOrbitsPerTF = 0;

  // Configurable<int> flagPbPb{"flagPbPb", 0, "0 - pp, 1 - PbPb"};

  // ##### hist registries
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext const&)
  {
    AxisSpec axisBC{3601, -0.5, 3600.5, "bc"};

    // ##### tracks for All collisions
    histos.add("hFoundBC_nAllTracks", "hFoundBC_nAllTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nTracksPV", "hFoundBC_nTracksPV", kTH1D, {axisBC});
    histos.add("hFoundBC_nGlobalTracks", "hFoundBC_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nITStracks", "hFoundBC_nITStracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nTPCtracks", "hFoundBC_nTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nTOFtracks", "hFoundBC_nTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nTRDtracks", "hFoundBC_nTRDtracks", kTH1D, {axisBC});

    histos.add("hFoundBC_nTRDTOFtracks", "hFoundBC_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nITSTPCTRDtracks", "hFoundBC_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nITSTPCTOFtracks", "hFoundBC_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_nITSTPCTRDTOFtracks", "hFoundBC_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

    // ##### tracks for TVX collisions
    histos.add("hFoundBC_kTVX_nAllTracks", "hFoundBC_kTVX_nAllTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nTracksPV", "hFoundBC_kTVX_nTracksPV", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nGlobalTracks", "hFoundBC_kTVX_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITStracks", "hFoundBC_kTVX_nITStracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nTPCtracks", "hFoundBC_kTVX_nTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCtracks", "hFoundBC_kTVX_nITSTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nTOFtracks", "hFoundBC_kTVX_nTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nTRDtracks", "hFoundBC_kTVX_nTRDtracks", kTH1D, {axisBC});

    histos.add("hFoundBC_kTVX_nTRDTOFtracks", "hFoundBC_kTVX_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks", "hFoundBC_kTVX_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTOFtracks", "hFoundBC_kTVX_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDTOFtracks", "hFoundBC_kTVX_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

    histos.add("hFoundBC_kTVX_counter", "hFoundBC_kTVX_counter", kTH1D, {axisBC});

    // n ITS layers per track in pt bins
    AxisSpec axisPtBinsForROFstudy{{0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0}, "p_{T}"};
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt", "hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt", kTH2D, {axisBC, axisPtBinsForROFstudy});

    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pi", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pi", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pi", "hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pi", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_ka", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_ka", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_ka", "hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_ka", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pr", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pr", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pr", "hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pr", kTH2D, {axisBC, axisPtBinsForROFstudy});

    // n ITS layers per track
    histos.add("hBC_kTVX_nITSlayers_for_ITSTPCtracks", "hBC_kTVX_nITSlayers_for_ITSTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_kTVXinTRD", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_kTVXinTRD", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pi", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pi", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_ka", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_ka", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pr", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pr", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSlayers_for_GlobalTracks", "hFoundBC_kTVX_nITSlayers_for_GlobalTracks", kTH1D, {axisBC});

    // counters
    histos.add("hBC_kTVX_counter_ITSTPCtracks", "hBC_kTVX_counter_ITSTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks", "hFoundBC_kTVX_counter_ITSTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_kTVXinTRD", "hFoundBC_kTVX_counter_ITSTPCtracks_kTVXinTRD", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_pi", "hFoundBC_kTVX_counter_ITSTPCtracks_pi", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_ka", "hFoundBC_kTVX_counter_ITSTPCtracks_ka", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_pr", "hFoundBC_kTVX_counter_ITSTPCtracks_pr", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_counter_GlobalTracks", "hFoundBC_kTVX_counter_GlobalTracks", kTH1D, {axisBC});

    // mult bins
    AxisSpec axisMultBins{{
                            -0.5,
                            4.5,
                            9.5,
                            19.5,
                            39.5,
                            79.5,
                            199.5,
                          },
                          "n PV ITSTPC tracks"};
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_multBins", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_multBins", kTH2D, {axisBC, axisMultBins});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_multBins", "hFoundBC_kTVX_counter_ITSTPCtracks_multBins", kTH2D, {axisBC, axisMultBins});

    AxisSpec axisClusterSizes{15, 0.5, 15.5, "cluster size #times cos(#Lambda)"};
    histos.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_clsSizeBins", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_clsSizeBins", kTH2D, {axisBC, axisClusterSizes});
    histos.add("hFoundBC_kTVX_counter_ITSTPCtracks_clsSizeBins", "hFoundBC_kTVX_counter_ITSTPCtracks_clsSizeBins", kTH2D, {axisBC, axisClusterSizes});

    // ##### tracks for TVXinTRD collisions
    histos.add("hFoundBC_kTVXinTRD_nAllTracks", "hFoundBC_kTVXinTRD_nAllTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nTracksPV", "hFoundBC_kTVXinTRD_nTracksPV", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nGlobalTracks", "hFoundBC_kTVXinTRD_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nITStracks", "hFoundBC_kTVXinTRD_nITStracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nTPCtracks", "hFoundBC_kTVXinTRD_nTPCtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nTOFtracks", "hFoundBC_kTVXinTRD_nTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nTRDtracks", "hFoundBC_kTVXinTRD_nTRDtracks", kTH1D, {axisBC});

    histos.add("hFoundBC_kTVXinTRD_nTRDTOFtracks", "hFoundBC_kTVXinTRD_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nITSTPCTRDtracks", "hFoundBC_kTVXinTRD_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nITSTPCTOFtracks", "hFoundBC_kTVXinTRD_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks", "hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

    // to study ITS ROF dips vs multiplicity
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000", kTH1D, {axisBC});
    histos.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000", kTH1D, {axisBC});

    // to study ITS ROF dips vs pt
    AxisSpec axisITSlayers{7, -0.5, 7 - 0.5, "layer id"};
    histos.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

    histos.add("hITSlayerCounts_vs_BC_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});

    // pions
    histos.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

    histos.add("hITSlayerCounts_vs_BC_pi_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_pi_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_pi_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});

    // protons
    histos.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
    histos.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

    histos.add("hITSlayerCounts_vs_BC_p_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_p_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    histos.add("hITSlayerCounts_vs_BC_p_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
    // histos.add("hITSlayerCounts_vs_BC_p_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});

    // bcInTF
    AxisSpec axisBCinTF32orbits{32 * 3564, -0.5, 32 * 3564, "bcInTF"};
    histos.add("hBCinTF_kTVX_nTracksPV", "hBCinTF_kTVX_nTracksPV;bc in TF; n tracks", kTH1D, {axisBCinTF32orbits});
    histos.add("hBCinTF_kTVX_nITSTPCTRDtracks", "hBCinTF_kTVX_nITSTPCTRDtracks;bc in TF; n tracks", kTH1D, {axisBCinTF32orbits});
  }

  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi>;
  // using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::FT0sCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected>;
  // using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::FT0sCorrected>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>; //, aod::Run3MatchedToBCSparse>;

  void processRun3(
    Colls::iterator const& collision,
    aod::FT0s const&,
    BCsRun3 const& bcs,
    aod::Origins const& /*origins*/,
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
              aod::TrackSelectionExtension,
              aod::TracksDCA, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
              aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFEl> const& tracks)
  // aod::V0Datas const& v0s, DaughterTracks const&)
  {
    auto bc = collision.bc_as<BCsRun3>();
    auto collBC = bc.globalBC() % 3564;
    uint64_t globalFoundBC = 9999;
    if (collision.has_foundBC()) {
      auto bcFound = collision.foundBC_as<BCsRun3>();
      globalFoundBC = bcFound.globalBC() % 3564;
    } else {
      return;
    }
    // histosEvent.fill(HIST("hBC_Bef"), collBC);

    if (fabs(collision.posZ()) > 10)
      return;

    // routine for each new run
    int runNumber = bc.runNumber();

    // #### begin of the code from qaPrimVtxVsTime.cxx
    if (lastRunNumber != runNumber) {
      /// execute the code in this scope only once, i.e. when the current run is considered for the first time in this DF
      lastRunNumber = runNumber;
      int64_t tsSOR = 0;
      // int64_t tsEOR = 0;

      /// reject AO2Ds for which no CCDB access is possible
      if (runNumber < 500000) {
        LOG(warning) << ">>> run number " << runNumber << " < 500000. access to CCDB not possible. Exiting";
        return;
      }
      /// If we are here, the current run was never considered before.

      // ##### code to find TF borders, Feb 1, 2024
      int64_t ts = bcs.iteratorAt(0).timestamp();
      // access orbit-reset timestamp
      auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
      int64_t tsOrbitReset = (*ctpx)[0]; // us

      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", runNumber);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
      nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF
      tsSOR = grpecs->getTimeStart();        // ms
      // tsEOR = grpecs->getTimeEnd();                   // ms

      // duration of TF in bcs
      nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;

      // IA - try without "-1" shift
      orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF; // - 1;
      bcSOR = orbitSOR * nBCsPerOrbit;
    }

    // bc in Time Frame:
    int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;

    // ### special sub-loop over tracks for numerator of tracks / T0ampl ratios
    if (!collision.has_foundBC())
      return;
    int nAllTracks = 0;
    int nTracksPV = 0;

    int nITStracks = 0;
    int nTPCtracks = 0;
    int nITSTPCtracks = 0;
    int nTOFtracks = 0;
    int nTRDtracks = 0;

    int nTRDTOFtracks = 0;
    int nITSTPCTRDtracks = 0;
    int nITSTPCTOFtracks = 0;
    int nITSTPCTRDTOFtracks = 0;

    int nGlobalTracks = 0;

    for (auto& track : tracks) {
      nAllTracks++;
      if (!track.isPVContributor()) {
        continue;
      }
      nTracksPV++;

      if (track.isGlobalTrack())
        nGlobalTracks++;

      nITStracks += track.hasITS() && !track.hasTPC();
      nTPCtracks += track.hasTPC();
      nITSTPCtracks += track.hasITS() && track.hasTPC();
      nTOFtracks += track.hasTOF();
      nTRDtracks += track.hasTRD() && !track.hasTOF();

      nTRDTOFtracks += track.hasTRD() && track.hasTOF();
      nITSTPCTRDtracks += track.hasITS() && track.hasTPC() && track.hasTRD();
      nITSTPCTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF();
      nITSTPCTRDTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF() && track.hasTRD();

      float eta = track.eta();
      float pt = track.pt();

      bool isPion = (TMath::Abs(track.tpcNSigmaPi()) < 3 && TMath::Abs(track.tpcNSigmaKa()) > 3 && TMath::Abs(track.tpcNSigmaPr()) > 3 && TMath::Abs(track.tpcNSigmaEl()) > 1) ? true : false;
      bool isKaon = (TMath::Abs(track.tpcNSigmaKa()) < 3 && TMath::Abs(track.tpcNSigmaPi()) > 3 && TMath::Abs(track.tpcNSigmaPr()) > 3 && TMath::Abs(track.tpcNSigmaEl()) > 1) ? true : false;
      bool isProton = (TMath::Abs(track.tpcNSigmaPr()) < 3 && TMath::Abs(track.tpcNSigmaKa()) > 3 && TMath::Abs(track.tpcNSigmaPi()) > 3 && TMath::Abs(track.tpcNSigmaEl()) > 1) ? true : false;

      if (track.hasITS() && track.hasTPC() && fabs(eta) < 0.8) {
        histos.fill(HIST("hBC_kTVX_nITSlayers_for_ITSTPCtracks"), collBC, track.itsNCls());
        histos.fill(HIST("hBC_kTVX_counter_ITSTPCtracks"), collBC);
      }

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) && track.hasITS() && track.hasTPC() && fabs(eta) < 0.8) {

        histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks"), globalFoundBC, track.itsNCls());
        histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks"), globalFoundBC);

        // ### vs mult
        histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_multBins"), globalFoundBC, nTracksPV, track.itsNCls());
        histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_multBins"), globalFoundBC, nTracksPV);

        // ### vs pt
        histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt"), globalFoundBC, pt, track.itsNCls());
        histos.fill(HIST("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt"), globalFoundBC, pt);

        if (collision.alias_bit(kTVXinTRD)) {
          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_kTVXinTRD"), globalFoundBC, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_kTVXinTRD"), globalFoundBC);
        }

        // ### with pid selection
        if (isPion) {
          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pi"), globalFoundBC, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_pi"), globalFoundBC);

          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pi"), globalFoundBC, pt, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pi"), globalFoundBC, pt);
        }
        if (isKaon) {
          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_ka"), globalFoundBC, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_ka"), globalFoundBC);

          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_ka"), globalFoundBC, pt, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_ka"), globalFoundBC, pt);
        }
        if (isProton) {
          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_pr"), globalFoundBC, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_pr"), globalFoundBC);

          histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt_pr"), globalFoundBC, pt, track.itsNCls());
          histos.fill(HIST("hFoundBC_kTVX_counter_nITSTPCtracks_vs_pt_pr"), globalFoundBC, pt);
        }
      }
      // ### global tracks
      if (track.isGlobalTrack() && fabs(eta) < 0.8) {
        histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_GlobalTracks"), globalFoundBC, track.itsNCls());
        histos.fill(HIST("hFoundBC_kTVX_counter_GlobalTracks"), globalFoundBC);
      }

      // ### vs cluster size
      if (track.itsNCls() >= 5) {
        float averageClusterSize = 0.;
        for (int i = 0; i < 7; i++) { // info stored in 4 bits
          averageClusterSize += (((1 << 4) - 1) & (track.itsClusterSizes() >> 4 * i));
        }
        averageClusterSize /= track.itsNCls();
        histos.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_clsSizeBins"), globalFoundBC, averageClusterSize, track.itsNCls());
        histos.fill(HIST("hFoundBC_kTVX_counter_ITSTPCtracks_clsSizeBins"), globalFoundBC, averageClusterSize);
      }
    }

    // all collisions
    histos.fill(HIST("hFoundBC_nAllTracks"), globalFoundBC, nAllTracks);
    histos.fill(HIST("hFoundBC_nTracksPV"), globalFoundBC, nTracksPV);
    histos.fill(HIST("hFoundBC_nGlobalTracks"), globalFoundBC, nGlobalTracks);
    histos.fill(HIST("hFoundBC_nITStracks"), globalFoundBC, nITStracks);
    histos.fill(HIST("hFoundBC_nTPCtracks"), globalFoundBC, nTPCtracks);
    histos.fill(HIST("hFoundBC_nTOFtracks"), globalFoundBC, nTOFtracks);
    histos.fill(HIST("hFoundBC_nTRDtracks"), globalFoundBC, nTRDtracks);

    histos.fill(HIST("hFoundBC_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
    histos.fill(HIST("hFoundBC_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
    histos.fill(HIST("hFoundBC_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
    histos.fill(HIST("hFoundBC_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);

    // TVX collisions
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      histos.fill(HIST("hFoundBC_kTVX_nAllTracks"), globalFoundBC, nAllTracks);
      histos.fill(HIST("hFoundBC_kTVX_nTracksPV"), globalFoundBC, nTracksPV);
      histos.fill(HIST("hFoundBC_kTVX_nGlobalTracks"), globalFoundBC, nGlobalTracks);
      histos.fill(HIST("hFoundBC_kTVX_nITStracks"), globalFoundBC, nITStracks);
      histos.fill(HIST("hFoundBC_kTVX_nTPCtracks"), globalFoundBC, nTPCtracks);
      histos.fill(HIST("hFoundBC_kTVX_nITSTPCtracks"), globalFoundBC, nITSTPCtracks);
      histos.fill(HIST("hFoundBC_kTVX_nTOFtracks"), globalFoundBC, nTOFtracks);
      histos.fill(HIST("hFoundBC_kTVX_nTRDtracks"), globalFoundBC, nTRDtracks);

      histos.fill(HIST("hFoundBC_kTVX_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
      histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
      histos.fill(HIST("hFoundBC_kTVX_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
      histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);

      histos.fill(HIST("hFoundBC_kTVX_counter"), globalFoundBC);

      // to study ITS ROF dips vs multiplicity
      if (nTracksPV >= 1 && nTracksPV < 20)
        histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20"), globalFoundBC, nITSTPCTRDtracks);
      else if (nTracksPV >= 20 && nTracksPV < 100)
        histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100"), globalFoundBC, nITSTPCTRDtracks);
      else if (nTracksPV >= 100 && nTracksPV < 400)
        histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400"), globalFoundBC, nITSTPCTRDtracks);
      else if (nTracksPV >= 400 && nTracksPV < 1000)
        histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000"), globalFoundBC, nITSTPCTRDtracks);
      else if (nTracksPV >= 1000)
        histos.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000"), globalFoundBC, nITSTPCTRDtracks);

      // vs bcInTF:
      histos.fill(HIST("hBCinTF_kTVX_nTracksPV"), bcInTF, nTracksPV);
      histos.fill(HIST("hBCinTF_kTVX_nITSTPCTRDtracks"), bcInTF, nITSTPCTRDtracks);
    }

    // TVXinTRD collisions
    if (collision.alias_bit(kTVXinTRD)) {
      histos.fill(HIST("hFoundBC_kTVXinTRD_nAllTracks"), globalFoundBC, nAllTracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nTracksPV"), globalFoundBC, nTracksPV);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nGlobalTracks"), globalFoundBC, nGlobalTracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nITStracks"), globalFoundBC, nITStracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nTPCtracks"), globalFoundBC, nTPCtracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nTOFtracks"), globalFoundBC, nTOFtracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nTRDtracks"), globalFoundBC, nTRDtracks);

      histos.fill(HIST("hFoundBC_kTVXinTRD_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
      histos.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);
    }
  }

  PROCESS_SWITCH(RofBorderQaTask, processRun3, "Process RofBorderQaTask", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<RofBorderQaTask>(cfgc)};
}
