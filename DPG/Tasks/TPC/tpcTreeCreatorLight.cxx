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

/// \file tpcTreeCreatorLight.cxx
/// \brief Task to produce table with clean selections for TPC PID calibration
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>

#include "tpcTreeCreatorLight.h"

#include <CCDB/BasicCCDBManager.h>

#include <cmath>
/// ROOT
#include "TRandom3.h"
/// O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
/// O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;

struct TreeWriterTPCTOF {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;

  /// Tables to be produced
  Produces<o2::aod::TPCTOFTree> rowTPCTOFTree;

  /// Configurables
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> applyTrkSel{"applyTrkSel", 1, "Flag to apply track selection: 0 -> no track selection, 1 -> track selection"};
  Configurable<double> dwnSmplFactor{"dwnSmplFactor", 1., "downsampling factor, default fraction to keep is 1. == 100%"};
  /// Proton
  Configurable<float> maxMomTPCOnlyPr{"maxMomTPCOnlyPr", 0.6, "Maximum momentum for TPC only cut proton"};
  Configurable<float> nSigmaTPCOnlyPr{"nSigmaTPCOnlyPr", 4., "number of sigma for TPC only cut proton"};
  Configurable<float> nSigmaTPC_TPCTOF_Pr{"nSigmaTPC_TPCTOF_Pr", 4., "number of sigma for TPC cut for TPC and TOF combined proton"};
  Configurable<float> nSigmaTOF_TPCTOF_Pr{"nSigmaTOF_TPCTOF_Pr", 3., "number of sigma for TOF cut for TPC and TOF combined proton"};

  /// Kaon
  Configurable<float> maxMomTPCOnlyKa{"maxMomTPCOnlyKa", 0.4, "Maximum momentum for TPC only cut kaon"};
  Configurable<float> maxMomHardCutOnlyKa{"maxMomHardCutOnlyKa", 50, "Maximum TPC inner momentum for kaons"};
  Configurable<float> nSigmaTPCOnlyKa{"nSigmaTPCOnlyKa", 4., "number of sigma for TPC only cut kaon"};
  Configurable<float> nSigmaTPC_TPCTOF_Ka{"nSigmaTPC_TPCTOF_Ka", 4., "number of sigma for TPC cut for TPC and TOF combined kaon"};
  Configurable<float> nSigmaTOF_TPCTOF_Ka{"nSigmaTOF_TPCTOF_Ka", 3., "number of sigma for TOF cut for TPC and TOF combined kaon"};
  /// Pion
  Configurable<float> maxMomTPCOnlyPi{"maxMomTPCOnlyPi", 0.5, "Maximum momentum for TPC only cut pion"};
  Configurable<float> nSigmaTPCOnlyPi{"nSigmaTPCOnlyPi", 4., "number of sigma for TPC only cut pion"};
  Configurable<float> nSigmaTPC_TPCTOF_Pi{"nSigmaTPC_TPCTOF_Pi", 4., "number of sigma for TPC cut for TPC and TOF combined pion"};
  Configurable<float> nSigmaTOF_TPCTOF_Pi{"nSigmaTOF_TPCTOF_Pi", 4., "number of sigma for TOF cut for TPC and TOF combined pion"};

  /// Electron
  Configurable<float> minMomTPCOnlyEl{"minMomTPCOnlyEl", 0.3, "Minimum momentum for TPC only cut electron"};
  Configurable<float> maxMomTPCOnlyEl{"maxMomTPCOnlyEl", 0.4, "Maximum momentum for TPC only cut electron"};
  Configurable<float> mindEdxTPCOnlyEl{"mindEdxTPCOnlyEl", 70., "dE/dx min for TPC only cut electron"};
  Configurable<float> maxdEdxTPCOnlyEl{"maxdEdxTPCOnlyEl", 100., "dE/dx max for TPC only cut electron"};

  int TrackPID = 0;

  /// Function to fill trees
  template <typename T, typename C>
  void fillSkimmedTPCTOFTable(T const& track, C const& collision, int runnumber, double dwnSmplFactor)
  {

    const double ncl = track.tpcNClsFound();
    const double p = track.tpcInnerParam();
    const int multTPC = collision.multTPC();

    const double pseudoRndm = track.pt() * 1000. - (int64_t)(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      rowTPCTOFTree(track.tpcSignal(),
                    p,
                    track.tgl(),
                    track.signed1Pt(),
                    track.eta(),
                    track.phi(),
                    track.y(),
                    multTPC / 11000.,
                    std::sqrt(nClNorm / ncl),
                    TrackPID,
                    runnumber);
    }
  };

  /// Event selection
  template <typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& /*tracks*/)
  {
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }
    return true;
  };

  /// Track selection
  template <typename CollisionType, typename TrackType>
  bool isTrackSelected(const CollisionType&, const TrackType& track)
  {
    if (!track.isGlobalTrack()) { // Skipping non global tracks
      return false;
    }
    if (!track.hasITS()) { // Skipping tracks without ITS
      return false;
    }
    if (!track.hasTPC()) { // Skipping tracks without TPC
      return false;
    }
    return true;
  };

  template <typename TrackType>
  bool KeepTrack(const TrackType& trk)
  {
    bool ProtonTrack = false;
    bool PionTrack = false;
    bool KaonTrack = false;
    bool ElectronTrack = false;
    TrackPID = 0;
    if ((trk.tpcInnerParam() < maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr) || (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr)) {
      ProtonTrack = true;
      TrackPID |= kProtonTrack;
    } else {
      TrackPID &= ~kProtonTrack;
    }
    if ((trk.tpcInnerParam() < maxMomHardCutOnlyKa) && ((trk.tpcInnerParam() < maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa) || (trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka))) {
      KaonTrack = true;
      TrackPID |= kKaonTrack;
    } else {
      TrackPID &= ~kKaonTrack;
    }
    if ((trk.tpcInnerParam() < maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi) || (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi)) {
      PionTrack = true;
      TrackPID |= kPionTrack;
    } else {
      TrackPID &= ~kPionTrack;
    }
    if ((trk.tpcInnerParam() > minMomTPCOnlyEl && trk.tpcInnerParam() < maxMomTPCOnlyEl && trk.tpcSignal() > mindEdxTPCOnlyEl && trk.tpcSignal() < maxdEdxTPCOnlyEl)) {
      ElectronTrack = true;
      TrackPID |= kElectronTrack;
    } else {
      TrackPID &= ~kElectronTrack;
    }
    if (ProtonTrack || KaonTrack || PionTrack || ElectronTrack) {
      return true;
    } else {
      return false;
    }
  }

  void init(o2::framework::InitContext&)
  {
  }
  void process(Coll::iterator const& collision, Trks const& tracks, aod::BCsWithTimestamps const&)
  {
    /// Check event selection
    if (!isEventSelected(collision, tracks)) {
      return;
    }

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const int runnumber = bc.runNumber();

    rowTPCTOFTree.reserve(tracks.size());
    for (auto const& trk : tracks) {

      /// Check track selection
      if (applyTrkSel == 1 && !isTrackSelected(collision, trk)) {
        continue;
      }
      if (KeepTrack(trk)) {
        fillSkimmedTPCTOFTable(trk, collision, runnumber, dwnSmplFactor);
      }

    } /// Loop tracks
  } /// process
}; /// struct TreeWriterTPCTOF

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<TreeWriterTPCTOF>(cfgc)};
  return workflow;
}
