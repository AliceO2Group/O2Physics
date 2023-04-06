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
// V0 MC finder task
// -----------------
//
//    This task allows for the re-creation of the V0 table (not the V0Data table)
//    using pure MC information. It serves the purpose of cross-checking the
//    maximum efficiency attainable in a perfect vertexing algorithm.
//
//    Nota bene: special attention could still be dedicated to the PV
//               reconstruction efficiency.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;

using LabeledTracks = soa::Join<aod::TracksIU, aod::McTrackLabels>;

struct lambdakzeromcfinder {
  Produces<aod::V0s> v0;
  Produces<aod::McFullV0Labels> fullv0labels;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for selecting which particles to generate
  Configurable<bool> findK0Short{"findK0Short", true, "findK0Short"};
  Configurable<bool> findLambda{"findLambda", true, "findLambda"};
  Configurable<bool> findAntiLambda{"findAntiLambda", true, "findAntiLambda"};
  Configurable<bool> findHyperTriton{"findHyperTriton", false, "findHyperTriton"};
  Configurable<bool> findAntiHyperTriton{"findAntiHyperTriton", false, "findAntiHyperTriton"};

  void init(InitContext& context)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});
  }

  template <typename TmcParticle, typename TTrackList>
  void PopulateV0s(TmcParticle const& mcParticle, TTrackList const& trackList, int bestCollisionIndex)
  {
    int trackIndexPositive = -1;
    int trackIndexNegative = -1;
    if (mcParticle.has_daughters()) {
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      if (daughters.size() == 2) {
        for (auto const& daughter : daughters) { // might be better ways of doing this but ok
          for (auto const& track : trackList) {
            if (track.mcParticleId() == daughter.globalIndex()) {
              // determine which charge this particle has
              if (track.sign() > 0) {
                trackIndexPositive = track.globalIndex();
                if (trackIndexNegative >= 0)
                  break;
              } else {
                trackIndexNegative = track.globalIndex();
                if (trackIndexPositive >= 0)
                  break;
              }
            }
          }
        }
      }
    }
    if (trackIndexPositive >= 0 && trackIndexNegative >= 0) {
      v0(bestCollisionIndex, trackIndexPositive, trackIndexNegative);
      fullv0labels(mcParticle.globalIndex());
    }
  }

  void process(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, LabeledTracks const& tracks, aod::McParticles const& mcParticles)
  {
    // Resolve collision (note: this loop is only over recoed collisions)
    histos.fill(HIST("hNTimesCollRecoed"), collisions.size());
    if (collisions.size() < 1)
      return; // not recoed, can't do anything, skip
    int biggestNContribs = -1;
    int bestCollisionIndex = -1;
    for (auto& collision : collisions) {
      if (biggestNContribs < collision.numContrib()) {
        biggestNContribs = collision.numContrib();
        bestCollisionIndex = collision.globalIndex();
      }
    }

    // Iterate over MC collisions, identify particles that are desired
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == 310 && findK0Short) {
        PopulateV0s(mcParticle, tracks, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == 3122 && findLambda) {
        PopulateV0s(mcParticle, tracks, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == -3122 && findAntiLambda) {
        PopulateV0s(mcParticle, tracks, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == 1010010030 && findHyperTriton) {
        PopulateV0s(mcParticle, tracks, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == -1010010030 && findAntiHyperTriton) {
        PopulateV0s(mcParticle, tracks, bestCollisionIndex);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeromcfinder>(cfgc, TaskName{"lf-lambdakzeromcfinder"})};
}
