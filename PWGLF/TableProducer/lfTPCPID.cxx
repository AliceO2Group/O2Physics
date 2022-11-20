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
/// \file   lfTPCPID.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2022-11-20
/// \brief  Task to produce the PID information for the TPC for the purpose of the Light flavor PWG
///

// ROOT includes
#include "TFile.h"
#include "TSystem.h"

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/CcdbApi.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "TableHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the response table
struct lfTpcPid {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = aod::Collisions;

  // Tables to produce
  Produces<o2::aod::pidTPCEl> tablePIDEl;
  Produces<o2::aod::pidTPCMu> tablePIDMu;
  Produces<o2::aod::pidTPCPi> tablePIDPi;
  Produces<o2::aod::pidTPCKa> tablePIDKa;
  Produces<o2::aod::pidTPCPr> tablePIDPr;
  Produces<o2::aod::pidTPCDe> tablePIDDe;
  Produces<o2::aod::pidTPCTr> tablePIDTr;
  Produces<o2::aod::pidTPCHe> tablePIDHe;
  Produces<o2::aod::pidTPCAl> tablePIDAl;
  // TPC PID Response
  o2::ccdb::CcdbApi ccdbApi;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<uint64_t> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  struct bbParams {
    float bb1 = 0.03209809958934784f;
    float bb2 = 19.9768009185791f;
    float bb3 = 2.5266601063857674e-16f;
    float bb4 = 2.7212300300598145f;
    float bb5 = 6.080920219421387f;
  } bbEl, bbMu, bbPi, bbKa, bbPr, bbDe, bbTr, bbHe, bbAl;

  template <o2::track::PID::ID id, typename T>
  float BetheBlochLf(const T& track, const bbParams& params)
  {
    static constexpr float invmass = 1.f / o2::track::pid_constants::sMasses[id];
    return o2::tpc::BetheBlochAleph(track.tpcInnerParam() * invmass, params.bb1, params.bb2, params.bb3, params.bb4, params.bb5);
  }

  float BetheBlochEl(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Electron>(track, bbEl); }
  float BetheBlochMu(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Muon>(track, bbMu); }
  float BetheBlochPi(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Pion>(track, bbPi); }
  float BetheBlochKa(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Kaon>(track, bbKa); }
  float BetheBlochPr(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Proton>(track, bbPr); }
  float BetheBlochDe(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Deuteron>(track, bbDe); }
  float BetheBlochTr(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Triton>(track, bbTr); }
  float BetheBlochHe(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Helium3>(track, bbHe); }
  float BetheBlochAl(const Trks::iterator& track) { return BetheBlochLf<o2::track::PID::Alpha>(track, bbAl); }

  void init(o2::framework::InitContext& initContext)
  {
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      enableFlagIfTableRequired(initContext, "pidTPC" + particle, flag);
    };
    enableFlag("El", pidEl);
    enableFlag("Mu", pidMu);
    enableFlag("Pi", pidPi);
    enableFlag("Ka", pidKa);
    enableFlag("Pr", pidPr);
    enableFlag("De", pidDe);
    enableFlag("Tr", pidTr);
    enableFlag("He", pidHe);
    enableFlag("Al", pidAl);

    // Get the parameters
    const std::string path = ccdbPath.value;
    const auto time = ccdbTimestamp.value;
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (time != 0) {
      LOGP(info, "Initialising TPC PID response for fixed timestamp {}:", time);
      ccdb->setTimestamp(time);
    } else {
      LOGP(info, "Initialising default TPC PID response:");
    }
  }

  void process(Coll const& collisions,
               Trks const& tracks,
               aod::BCsWithTimestamps const&)
  {

    const uint64_t tracks_size = tracks.size();

    auto reserveTable = [&tracks_size](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks_size);
    };
    // Prepare memory for enabled tables
    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    // int lastCollisionId = -1; // Last collision ID analysed

    for (auto const& trk : tracks) {
      // Loop on Tracks
      // if (useCCDBParam && ccdbTimestamp.value == 0 && trk.has_collision() && trk.collisionId() != lastCollisionId) { // Updating parametrization only if the initial timestamp is 0
      //   lastCollisionId = trk.collisionId();
      //   const auto& bc = collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>();
      // }
      // Check and fill enabled tables

#define doFillTable(Particle)                                                                                       \
  if (pid##Particle.value == 1) {                                                                                   \
    aod::pidutils::packInTable<aod::pidtpc_tiny::binning>(trk.tpcSignal() - BetheBlochPi(trk), tablePID##Particle); \
  }

      doFillTable(El);
      doFillTable(Mu);
      doFillTable(Pi);
      doFillTable(Ka);
      doFillTable(Pr);
      doFillTable(De);
      doFillTable(Tr);
      doFillTable(He);
      doFillTable(Al);

#undef doFillTable
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lfTpcPid>(cfgc)};
}
