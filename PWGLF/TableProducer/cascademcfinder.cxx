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
//    This task allows for the re-creation of the cascade table (not the CascData table)
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
using LabeledFullV0s = soa::Join<aod::V0s, aod::McFullV0Labels>;

struct cascademcfinder {
  Produces<aod::Cascades> cascades;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for selecting which particles to generate
  Configurable<bool> findXiMinus{"findXiMinus", true, "findXiMinus"};
  Configurable<bool> findXiPlus{"findXiPlus", true, "findXiPlus"};
  Configurable<bool> findOmegaMinus{"findOmegaMinus", true, "findOmegaMinus"};
  Configurable<bool> findOmegaPlus{"findOmegaPlus", true, "findOmegaPlus"};

  void init(InitContext& context)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});
  }

  template <typename TmcParticle, typename TTrackList, typename TV0List>
  void PopulateCascade(TmcParticle const& mcParticle, TTrackList const& trackList, TV0List const& v0s, int bestCollisionIndex)
  {
    int trackIndexBachelor = -1;
    int trackIndexV0 = -1;
    if (mcParticle.has_daughters()) {
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      if (daughters.size() == 2) {
        for (auto const& daughter : daughters) { // might be better ways of doing this but ok
          // option 1: this is the lambda daughter
          if (TMath::Abs(daughter.pdgCode()) == 3122) {
            for (auto const& v0 : v0s) {
              if (v0.mcParticleId() == daughter.globalIndex()) {
                trackIndexV0 = v0.globalIndex();
                break;
              }
            }
          } // end lambda search
          if (TMath::Abs(daughter.pdgCode()) == 211 || TMath::Abs(daughter.pdgCode()) == 321) {
            for (auto const& track : trackList) {
              if (track.mcParticleId() == daughter.globalIndex()) {
                trackIndexBachelor = track.globalIndex();
                break;
              }
            }
          } // end bachelor search
        }
      }
    }
    if (trackIndexBachelor >= 0 && trackIndexV0 >= 0) {
      cascades(bestCollisionIndex, trackIndexV0, trackIndexBachelor);
    }
  }

  void process(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, LabeledTracks const& tracks, aod::McParticles const& mcParticles, LabeledFullV0s const& v0s)
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
      if (mcParticle.pdgCode() == 3312 && findXiMinus) {
        PopulateCascade(mcParticle, tracks, v0s, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == -3312 && findXiPlus) {
        PopulateCascade(mcParticle, tracks, v0s, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == 3334 && findOmegaMinus) {
        PopulateCascade(mcParticle, tracks, v0s, bestCollisionIndex);
      }
      if (mcParticle.pdgCode() == -3334 && findOmegaPlus) {
        PopulateCascade(mcParticle, tracks, v0s, bestCollisionIndex);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascademcfinder>(cfgc, TaskName{"lf-cascademcfinder"})};
}
