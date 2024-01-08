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
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include <fairlogger/Logger.h>
#include <Framework/AnalysisDataModel.h>

#include <vector>

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct singleTrackSelector {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  // Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  Configurable<std::vector<int>> _particlesToKeep{"particlesToKeepPDGs", std::vector<int>{2212, 1000010020}, "PDG codes of perticles for which the 'singletrackselector' tables will be created (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> keepWithinNsigmaTPC{"keepWithinNsigmaTPC", std::vector<float>{-4.0f, 4.0f}, "TPC range for preselection of particles specified with PDG"};
  Configurable<std::vector<int>> _particlesToReject{"particlesToRejectPDGs", std::vector<int>{}, "PDG codes of perticles that will be rejected with TOF (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-5.0f, 5.0f}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};

  Configurable<float> _min_P{"min_P", 0.f, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.f, "upper mometum limit"};
  Configurable<float> _eta{"eta", 100.f, "abs eta value limit"};
  Configurable<float> _dcaXY{"dcaXY", 1000.f, "Maximum dca of track in xy"};
  Configurable<float> _dcaZ{"dcaZ", 1000.f, "Maximum dca of track in xy"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidEvTimeFlags, aod::TracksDCA,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                         aod::pidTPCFullPr, aod::pidTPCFullDe,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa,
                         aod::pidTOFFullPr, aod::pidTOFFullDe,
                         aod::TrackSelection, aod::pidTOFbeta>;

  using Coll = soa::Join<aod::Collisions, aod::TPCMults, aod::EvSels, aod::CentFT0Ms>;

  Produces<o2::aod::SingleTrackSels> tableRow;
  Produces<o2::aod::SingleCollSels> tableRowColl;
  Produces<o2::aod::SingleTrkExtras> tableRowExtra;
  Produces<o2::aod::SingleTrkMCs> tableRowMC;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (aod::evsel::sel8 == true));
  Filter vertexFilter = ((o2::aod::collision::posZ < 15.f) && (o2::aod::collision::posZ > -15.f));
  Filter trackFilter = ((o2::aod::track::itsChi2NCl <= 36.f) && (o2::aod::track::itsChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl <= 4.f));

  Filter pFilter = o2::aod::track::p > _min_P&& o2::aod::track::p < _max_P;
  Filter etaFilter = nabs(o2::aod::track::eta) < _eta;
  Filter dcaFilter = (o2::aod::track::dcaXY <= _dcaXY) && (o2::aod::track::dcaZ <= _dcaZ);

  int mRunNumber = 0;
  float d_bz = 0.f;

  std::vector<int> particlesToKeep;
  std::vector<int> particlesToReject;

  void init(InitContext& context)
  {

    particlesToKeep = _particlesToKeep;
    particlesToReject = _particlesToReject;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // inspired by PWGLF/TableProducer/lambdakzerobuilder.cxx
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    d_bz = 0.f;

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    d_bz = 0.1 * d_bz;
  }

  template <bool isMC, typename Col, typename Trks>
  inline void fillTheTables(Col collision, Trks const& tracks)
  {
    bool skip_track = false; // flag used for track rejection

    tableRowColl(collision.multTPC(),
                 collision.centFT0M(),
                 collision.posZ(),
                 d_bz);

    for (auto& track : tracks) {
      if constexpr (isMC) {
        if (!track.has_mcParticle())
          continue;
      }
      skip_track = false;
      if (Configurable<bool>{"rejectNotPropagatedTrks", true, "rejects tracks that are not propagated to the primary vertex"} && track.trackType() != aod::track::Track) {
        continue;
      }
      for (auto i : particlesToReject) {
        // if satisfied, want to continue in the upper loop (over tracks) -- skip the current track
        // cannot use simple 'continue' since it will be applied to the current loop, so have to use a flag
        if (o2::aod::singletrackselector::TOFselection(track, std::make_pair(i, rejectWithinNsigmaTOF))) {
          skip_track = true;
          break;
        }
      }

      if (skip_track)
        continue;

      for (auto ii : particlesToKeep)
        if (o2::aod::singletrackselector::TPCselection(track, std::make_pair(ii, keepWithinNsigmaTPC))) {

          tableRow(tableRowColl.lastIndex(),
                   track.p(),
                   track.eta(),
                   track.phi(),
                   track.sign(),
                   track.tpcNClsFound(),
                   track.tpcNClsShared(),
                   track.itsNCls(),
                   singletrackselector::packInTable<singletrackselector::binning::dca>(track.dcaXY()),
                   singletrackselector::packInTable<singletrackselector::binning::dca>(track.dcaZ()),
                   singletrackselector::packInTable<singletrackselector::binning::chi2>(track.tpcChi2NCl()),
                   singletrackselector::packInTable<singletrackselector::binning::chi2>(track.itsChi2NCl()),
                   singletrackselector::packInTable<singletrackselector::binning::rowsOverFindable>(track.tpcCrossedRowsOverFindableCls()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tofNSigmaPi()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tpcNSigmaPi()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tofNSigmaKa()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tpcNSigmaKa()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tofNSigmaPr()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tpcNSigmaPr()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tofNSigmaDe()),
                   singletrackselector::packInTable<singletrackselector::binning::nsigma>(track.tpcNSigmaDe()));

          tableRowExtra(track.tpcInnerParam(),
                        track.tpcSignal(),
                        track.beta());

          if constexpr (isMC) {
            int origin = -1;
            if (track.mcParticle().isPhysicalPrimary())
              origin = 0;
            if (!track.mcParticle().isPhysicalPrimary() && track.mcParticle().producedByGenerator())
              origin = 1;
            if (!track.mcParticle().isPhysicalPrimary() && !track.mcParticle().producedByGenerator())
              origin = 2;

            if (origin == -1)
              LOGF(fatal, "Could not define the origin (primary/weak decay/material) of the track!!!");

            tableRowMC(track.mcParticle().pdgCode(),
                       origin,
                       track.mcParticle().p(),
                       track.mcParticle().eta(),
                       track.mcParticle().phi());
          }

          break; // break the loop with particlesToKeep after the 'if' condition is satisfied -- don't want double entries
        }
    }
  }

  void processData(soa::Filtered<Coll>::iterator const& collision, soa::Filtered<Trks> const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    fillTheTables<false>(collision, tracks);
  }
  PROCESS_SWITCH(singleTrackSelector, processData, "process data", true);

  void processMC(soa::Filtered<Coll>::iterator const& collision, soa::Filtered<soa::Join<Trks, aod::McTrackLabels>> const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    fillTheTables<true>(collision, tracks);
  }
  PROCESS_SWITCH(singleTrackSelector, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelector>(cfgc)};
}
