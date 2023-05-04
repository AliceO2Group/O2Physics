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
#include "PWGLF/DataModel/LFQATables.h"
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

using LabeledTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

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
  Configurable<bool> requireITS{"requireITS", false, "require ITS information used in tracks"};

  Preslice<aod::McParticle> perMcCollision = aod::mcparticle::mcCollisionId;

  std::vector<int> v0collisionId;
  std::vector<int> v0positiveIndex;
  std::vector<int> v0negativeIndex;
  std::vector<int> v0mcLabel;

  void init(InitContext& context)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};
    const AxisSpec axisPt{(int)100, +0.0f, +10.0f, "p_{T} (GeV/c)"};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});

    histos.add("hPtK0ShortGenerated", "hPtK0ShortGenerated", kTH1F, {axisPt});
    histos.add("hPtLambdaGenerated", "hPtLambdaGenerated", kTH1F, {axisPt});
    histos.add("hPtAntiLambdaGenerated", "hPtAntiLambdaGenerated", kTH1F, {axisPt});
    histos.add("hPtHypertritonGenerated", "hPtHypertritonGenerated", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonGenerated", "hPtAntiHypertritonGenerated", kTH1F, {axisPt});

    histos.add("hPtK0ShortReconstructed", "hPtK0ShortReconstructed", kTH1F, {axisPt});
    histos.add("hPtLambdaReconstructed", "hPtLambdaReconstructed", kTH1F, {axisPt});
    histos.add("hPtAntiLambdaReconstructed", "hPtAntiLambdaReconstructed", kTH1F, {axisPt});
    histos.add("hPtHypertritonReconstructed", "hPtHypertritonReconstructed", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonReconstructed", "hPtAntiHypertritonReconstructed", kTH1F, {axisPt});

    histos.add("hPtK0ShortGlobal", "hPtK0ShortGlobal", kTH1F, {axisPt});
    histos.add("hPtLambdaGlobal", "hPtLambdaGlobal", kTH1F, {axisPt});
    histos.add("hPtAntiLambdaGlobal", "hPtAntiLambdaGlobal", kTH1F, {axisPt});
    histos.add("hPtHypertritonGlobal", "hPtHypertritonGlobal", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonGlobal", "hPtAntiHypertritonGlobal", kTH1F, {axisPt});

    histos.add("hPtK0ShortGlobalWithPV", "hPtK0ShortGlobalWithPV", kTH1F, {axisPt});
    histos.add("hPtLambdaGlobalWithPV", "hPtLambdaGlobalWithPV", kTH1F, {axisPt});
    histos.add("hPtAntiLambdaGlobalWithPV", "hPtAntiLambdaGlobalWithPV", kTH1F, {axisPt});
    histos.add("hPtHypertritonGlobalWithPV", "hPtHypertritonGlobalWithPV", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonGlobalWithPV", "hPtAntiHypertritonGlobalWithPV", kTH1F, {axisPt});
  }

  // for sorting
  template <typename T>
  std::vector<std::size_t> sort_indices(const std::vector<T>& v)
  {
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
                     [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
    return idx;
  }

  template <typename TmcParticle, typename TTrackList>
  bool ProcessV0(TmcParticle const& mcParticle, TTrackList const& trackList, int bestCollisionIndex, bool& positiveITS, bool& negativeITS, bool& positiveTPC, bool& negativeTPC)
  {
    bool reconstructed = false;
    positiveITS = false;
    negativeITS = false;
    positiveTPC = false;
    negativeTPC = false;
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
                if (track.hasITS())
                  positiveITS = true;
                if (track.hasTPC())
                  positiveTPC = true;
                if (trackIndexNegative >= 0)
                  break;
              } else {
                trackIndexNegative = track.globalIndex();
                if (track.hasITS())
                  negativeITS = true;
                if (track.hasTPC())
                  negativeTPC = true;
                if (trackIndexPositive >= 0)
                  break;
              }
            }
          }
        }
      }
    }
    if (trackIndexPositive >= 0 && trackIndexNegative >= 0 && (!requireITS || (requireITS && positiveITS && negativeITS))) {
      reconstructed = true;
      v0collisionId.emplace_back(bestCollisionIndex);
      v0positiveIndex.emplace_back(trackIndexPositive);
      v0negativeIndex.emplace_back(trackIndexNegative);
      v0mcLabel.emplace_back(mcParticle.globalIndex());
    }
    return reconstructed;
  }

  void process(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, LabeledTracks const& tracks, aod::McParticles const& allMcParticles)
  {
    v0collisionId.clear();
    v0positiveIndex.clear();
    v0negativeIndex.clear();
    v0mcLabel.clear();

    // Step 1: sweep over all mcCollisions and find all relevant candidates
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hNTimesCollRecoed"), mcCollision.hasRecoCollision());
      int bestCollisionIndex = mcCollision.bestCollisionIndex();

      auto mcParticles = allMcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      bool positiveITS = false;
      bool negativeITS = false;
      bool positiveTPC = false;
      bool negativeTPC = false;
      bool reconstructed = false;
      for (auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() == 310 && findK0Short) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtK0ShortGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtK0ShortReconstructed"), mcParticle.pt());
            if (reconstructed && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtK0ShortGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtK0ShortGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 3122 && findLambda) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtLambdaGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtLambdaReconstructed"), mcParticle.pt());
            if (reconstructed && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtLambdaGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtLambdaGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -3122 && findAntiLambda) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiLambdaGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiLambdaReconstructed"), mcParticle.pt());
            if (reconstructed && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtAntiLambdaGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtAntiLambdaGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 1010010030 && findHyperTriton) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -1010010030 && findAntiHyperTriton) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtAntiHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && positiveTPC && negativeTPC)
              histos.fill(HIST("hPtAntiHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(v0collisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      if (v0collisionId[ic] >= 0) {
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic]);
        fullv0labels(v0mcLabel[ic]);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeromcfinder>(cfgc, TaskName{"lf-lambdakzeromcfinder"})};
}
