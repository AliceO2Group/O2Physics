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

/// \file LFTreeCreatorNuclei.cxx
/// \brief Writer of the nuclei candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch> and Francesca Bellini <fbellini@cern.ch>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
//included from NucleiSpectraTask.cxx
#include "Framework/AnalysisDataModel.h" 
#include "Framework/ASoAHelpers.h"  
#include "Framework/HistogramRegistry.h"

#include "Common/Core/PID/PIDResponse.h" 
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h" 
#include "Common/Core/trackUtilities.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::lf_cand_nucleus;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigTPCPi0, nsigTPCPi0, float);
DECLARE_SOA_COLUMN(NSigTPCPr0, nsigTPCPr0, float);
DECLARE_SOA_COLUMN(NSigTPCDe0, nsigTPCD0, float);
DECLARE_SOA_COLUMN(NSigTPC3He0, nsigTPC3He0, float);
DECLARE_SOA_COLUMN(NSigTOFPi0, nsigTOFPi0, float);
DECLARE_SOA_COLUMN(NSigTOFPr0, nsigTOFPr0, float);
DECLARE_SOA_COLUMN(NSigTOFDe0, nsigTOFD0, float);
DECLARE_SOA_COLUMN(NSigTOF3He0, nsigTOF3He0, float);
DECLARE_SOA_COLUMN(TOFmatch, tofMatch, bool);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(LfCandNucleusFull, "AOD", "HFCANDP3Full",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::NSigTPCPi0,
                  full::NSigTPCPr0,
                  full::NSigTPCDe0,
                  full::NSigTPC3He0,
                  full::NSigTOFPi0,
                  full::NSigTOFPr0,
                  full::NSigTOFDe0,
                  full::NSigTOF3He0,
                  full::TOFmatch,
                  full::Pt,
                  full::P,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::sign,
                  full::MCflag);

DECLARE_SOA_TABLE(LfCandNucleusFullEvents, "AOD", "LFNUCLFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct CandidateTreeWriter {
  Produces<o2::aod::LfCandNucleusFull> rowCandidateFull;
  Produces<o2::aod::LfCandNucleusFullEvents> rowCandidateFullEvents;

  void init(o2::framework::InitContext&)
  {
  }

  Configurable<float> yMin{"yMin", -0.5, "Maximum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Minimum rapidity"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> nsigmacutLow{"nsigmacutLow", -8.0, "Value of the Nsigma cut"};
  Configurable<float> nsigmacutHigh{"nsigmacutHigh", +8.0, "Value of the Nsigma cut"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::isGlobalTrack == (uint8_t) true);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCFullHe, aod::pidTOFFullHe, aod::TrackSelection, aod::TOFSignal>>;

  void process(soa::Filtered<soa::Join<aod::Collisions const& collisions,
               aod::EvSels>>::iterator const& collision, 
               TrackCandidates const& tracks)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto fillTable = [&](int CandFlag,
                           int FunctionSelection,
                           float FunctionY,
                           float FunctionE) {
        if (FunctionSelection >= 1) {
          rowCandidateFull(
            collision.bcId(),
            collision.numContrib(),
            track.tpcNSigmaPi(),
            track.tpcNSigmaPr(),
            track.tpcNSigmaDe(),
            track.tpcNSigmaHe(),
            track.tofNSigmaPi(),
            track.tofNSigmaPr(),
            track.tofNSigmaD(),
            track.tofNSigmaHe(),
            track.hasTOF(),
            track.pt(),
            candidate.p(),
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            FunctionE,
            track.sign(),
            candidate.flagMCMatchRec());
        }
      };
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<CandidateTreeWriter>(cfgc, TaskName{"lf-tree-creator-nucleus"}));
  return workflow;
}
