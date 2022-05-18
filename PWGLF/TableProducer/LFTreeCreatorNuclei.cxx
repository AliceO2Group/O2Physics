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
/// \file LFTreeCreatorNuclei.cxx
/// \brief Writer of the nuclei candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch> and Francesca Bellini <fbellini@cern.ch>
///

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFNucleiTables.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Writes the full information in an output TTree
struct LfTreeCreatorNuclei {
  Produces<o2::aod::LfCandNucleusFull> rowCandidateFull;
  Produces<o2::aod::LfCandNucleusFullEvents> rowCandidateFullEvents;

  void init(o2::framework::InitContext&)
  {
  }

  //track
  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -8.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +8.0, "Value of the Nsigma cut"};
  //events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> useEvsel{"useEvsel", true, "Use sel8 for run3 Event Selection"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (requireGlobalTrackInFilter());
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::TrackSelection,
                                                  aod::pidTOFbeta, aod::TOFSignal,
                                                  aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                  aod::pidTPCFullPr, aod::pidTOFFullPr,
                                                  aod::pidTPCFullDe, aod::pidTOFFullDe,
                                                  aod::pidTPCFullHe, aod::pidTOFFullHe>>;

  void process(soa::Filtered<soa::Join<aod::Collisions,
                                       aod::EvSels, aod::Mults>>::iterator const& collision,
               TrackCandidates const& tracks,
               aod::BCs const&)
  {
    if (useEvsel && !collision.sel8()) {
      return;
    }
    // Filling event properties
    rowCandidateFullEvents(
      collision.bcId(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      collision.multFV0M(),
      collision.sel8(),
      collision.bc().runNumber());

    // Filling candidate properties
    rowCandidateFull.reserve(tracks.size());
    for (auto& track : tracks) {
      rowCandidateFull(
        rowCandidateFullEvents.lastIndex(),
        track.dcaXY(),
        track.dcaZ(),
        track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
        track.tpcNSigmaDe(), track.tpcNSigmaHe(),
        track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
        track.tofNSigmaDe(), track.tofNSigmaHe(),
        track.hasTOF(),
        track.tpcInnerParam(),
        track.tpcSignal(),
        track.beta(),
        track.px(),
        track.py(),
        track.pz(),
        track.pt(),
        track.p(),
        track.eta(),
        track.phi(),
        track.sign(),
        track.tpcNClsCrossedRows(),
        track.tpcCrossedRowsOverFindableCls(),
        track.tpcChi2NCl(),
        track.itsChi2NCl());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorNuclei>(cfgc)};
}
