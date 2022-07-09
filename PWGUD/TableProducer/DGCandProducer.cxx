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
// \brief Saves relevant information of DG candidates
//
//     options:
//           DiffCuts.mNDtcoll(4)
//           DiffCuts.mMinNBCs(7)
//           DiffCuts.mMinNTracks(0)
//           DiffCuts.mMaxNTracks(10000)
//           DiffCuts.mMinNetCharge(0)
//           DiffCuts.mMaxNetCharge(0)
//           DiffCuts.mPidHypo(211)
//           DiffCuts.mMinPosz(-1000.)
//           DiffCuts.mMaxPosz(1000.)
//           DiffCuts.mMinPt(0.)
//           DiffCuts.mMaxPt(1000.)
//           DiffCuts.mMinEta(-1.)
//           DiffCuts.mMaxEta(1.)
//           DiffCuts.mMinIVM(0.)
//           DiffCuts.mMaxIVM(1000.)
//           DiffCuts.mMaxnSigmaTPC(1000.)
//           DiffCuts.mMaxnSigmaTOF(1000.)
//           DiffCutsX.mFITAmpLimits({0., 0., 0., 0., 0.})
//
//     usage: copts="--configuration json://DGCandProducerConfig.json --aod-writer-json DGCandProducerWriter.json -b"
//
//           o2-analysis-timestamp $copts |
//           o2-analysis-track-propagation $copts |
//           o2-analysis-multiplicity-table $copts |
//           o2-analysis-ft0-corrected-table $copts |
//           o2-analysis-event-selection $copts |
//           o2-analysis-trackextension $copts |
//           o2-analysis-trackselection $copts |
//           o2-analysis-pid-tpc-full $copts |
//           o2-analysis-pid-tof-base $copts |
//           o2-analysis-pid-tof-full $copts |
//           o2-analysis-ud-dgcand-producer $copts > DGCandProducer.log
//
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  20.05.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "EventFiltering/PWGUD/DGHelpers.h"
#include "PWGUD/Core/UDHelperFunctions.h"
#include "PWGUD/DataModel/DGCandidates.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandProducer {

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // DG selector
  DGSelector dgSelector;

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
  }

  Produces<aod::DGCandidates> outputCollisions;
  Produces<aod::DGTracks> outputTracks;

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  void process(CC const& collision,
               BCs const& bcs,
               TCs& tracks,
               FWs& fwdtracks,
               aod::Zdcs& zdcs,
               aod::FT0s& ft0s,
               aod::FV0As& fv0as,
               aod::FDDs& fdds)
  {
    // nominal BC
    auto bc = collision.bc_as<BCs>();

    // obtain slice of compatible BCs
    auto bcRange = compatibleBCs(collision, diffCuts.NDtcoll(), bcs, diffCuts.minNBCs());

    // apply DG selection
    auto isDGEvent = dgSelector.IsSelected(diffCuts, collision, bc, bcRange, tracks, fwdtracks);

    // save DG candidates
    if (isDGEvent == 0) {

      // update DGCandidates tables
      outputCollisions(bc.runNumber(), bc.timestamp(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       collision.numContrib(), netCharge(tracks), rPVtrwTOF(tracks, collision.numContrib()));

      // update DGTracks tables
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          outputTracks(outputCollisions.lastIndex(), track.pt(), track.eta(), track.phi(), track.sign(),
                       track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                       track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandProducer>(cfgc, TaskName{"dgcandproducer"}),
  };
}
