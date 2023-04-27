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
using LabeledFullV0s = soa::Join<aod::V0s, aod::McFullV0Labels>;

struct cascademcfinder {
  Produces<aod::Cascades> cascades;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for selecting which particles to generate
  Configurable<bool> findXiMinus{"findXiMinus", true, "findXiMinus"};
  Configurable<bool> findXiPlus{"findXiPlus", true, "findXiPlus"};
  Configurable<bool> findOmegaMinus{"findOmegaMinus", true, "findOmegaMinus"};
  Configurable<bool> findOmegaPlus{"findOmegaPlus", true, "findOmegaPlus"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS information used in tracks"};

  Preslice<aod::McParticle> perMcCollision = aod::mcparticle::mcCollisionId;

  std::vector<int> casccollisionId;
  std::vector<int> cascv0Index;
  std::vector<int> cascbachelorIndex;

  void init(InitContext& context)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};
    const AxisSpec axisPt{(int)100, +0.0f, +10.0f, "p_{T} (GeV/c)"};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});

    histos.add("hPtXiMinusGenerated", "hPtXiMinusGenerated", kTH1F, {axisPt});
    histos.add("hPtXiPlusGenerated", "hPtXiPlusGenerated", kTH1F, {axisPt});
    histos.add("hPtOmegaMinusGenerated", "hPtOmegaMinusGenerated", kTH1F, {axisPt});
    histos.add("hPtOmegaPlusGenerated", "hPtOmegaPlusGenerated", kTH1F, {axisPt});

    histos.add("hPtXiMinusReconstructed", "hPtXiMinusReconstructed", kTH1F, {axisPt});
    histos.add("hPtXiPlusReconstructed", "hPtXiPlusReconstructed", kTH1F, {axisPt});
    histos.add("hPtOmegaMinusReconstructed", "hPtOmegaMinusReconstructed", kTH1F, {axisPt});
    histos.add("hPtOmegaPlusReconstructed", "hPtOmegaPlusReconstructed", kTH1F, {axisPt});

    histos.add("hPtXiMinusReconstructedWithPV", "hPtXiMinusReconstructedWithPV", kTH1F, {axisPt});
    histos.add("hPtXiPlusReconstructedWithPV", "hPtXiPlusReconstructedWithPV", kTH1F, {axisPt});
    histos.add("hPtOmegaMinusReconstructedWithPV", "hPtOmegaMinusReconstructedWithPV", kTH1F, {axisPt});
    histos.add("hPtOmegaPlusReconstructedWithPV", "hPtOmegaPlusReconstructedWithPV", kTH1F, {axisPt});

    histos.add("hPtXiMinusReconstructedWithPV_ITS", "hPtXiMinusReconstructedWithPV_ITS", kTH1F, {axisPt});
    histos.add("hPtXiPlusReconstructedWithPV_ITS", "hPtXiPlusReconstructedWithPV_ITS", kTH1F, {axisPt});
    histos.add("hPtOmegaMinusReconstructedWithPV_ITS", "hPtOmegaMinusReconstructedWithPV_ITS", kTH1F, {axisPt});
    histos.add("hPtOmegaPlusReconstructedWithPV_ITS", "hPtOmegaPlusReconstructedWithPV_ITS", kTH1F, {axisPt});
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

  template <typename TmcParticle, typename TTrackList, typename TV0List>
  bool ProcessCascade(TmcParticle const& mcParticle, TTrackList const& trackList, TV0List const& v0s, int bestCollisionIndex, bool& positiveITS, bool& negativeITS, bool& bachelorITS)
  {
    bool reconstructed = false;
    positiveITS = false;
    negativeITS = false;
    bachelorITS = false;
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
                auto const& positiveTrack = v0.template posTrack_as<LabeledTracks>();
                auto const& negativeTrack = v0.template negTrack_as<LabeledTracks>();
                if (positiveTrack.hasITS())
                  positiveITS = true;
                if (negativeTrack.hasITS())
                  negativeITS = true;
                break;
              }
            }
          } // end lambda search
          if (TMath::Abs(daughter.pdgCode()) == 211 || TMath::Abs(daughter.pdgCode()) == 321) {
            for (auto const& track : trackList) {
              if (track.mcParticleId() == daughter.globalIndex()) {
                trackIndexBachelor = track.globalIndex();
                if (track.hasITS())
                  bachelorITS = true;
                break;
              }
            }
          } // end bachelor search
        }
      }
    }
    if (trackIndexBachelor >= 0 && trackIndexV0 >= 0 && (!requireITS || (requireITS && positiveITS && negativeITS && bachelorITS))) {
      reconstructed = true;
      casccollisionId.emplace_back(bestCollisionIndex);
      cascbachelorIndex.emplace_back(trackIndexBachelor);
      cascv0Index.emplace_back(trackIndexV0);
    }
    return reconstructed;
  }

  void process(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, LabeledTracks const& tracks, aod::McParticles const& allMcParticles, LabeledFullV0s const& v0s)
  {
    casccollisionId.clear();
    cascbachelorIndex.clear();
    cascv0Index.clear();

    // Step 1: sweep over all mcCollisions and find all relevant candidates
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hNTimesCollRecoed"), mcCollision.hasRecoCollision());
      int bestCollisionIndex = mcCollision.bestCollisionIndex();

      auto mcParticles = allMcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      bool positiveITS = false;
      bool negativeITS = false;
      bool bachelorITS = false;
      bool reconstructed = false;
      for (auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() == 3312 && findXiMinus) {
          reconstructed = ProcessCascade(mcParticle, tracks, v0s, bestCollisionIndex, positiveITS, negativeITS, bachelorITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtXiMinusGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtXiMinusReconstructed"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0)
              histos.fill(HIST("hPtXiMinusReconstructedWithPV"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && bachelorITS)
              histos.fill(HIST("hPtXiMinusReconstructedWithPV_ITS"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -3312 && findXiPlus) {
          reconstructed = ProcessCascade(mcParticle, tracks, v0s, bestCollisionIndex, positiveITS, negativeITS, bachelorITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtXiPlusGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtXiPlusReconstructed"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0)
              histos.fill(HIST("hPtXiPlusReconstructedWithPV"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && bachelorITS)
              histos.fill(HIST("hPtXiPlusReconstructedWithPV_ITS"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 3334 && findOmegaMinus) {
          reconstructed = ProcessCascade(mcParticle, tracks, v0s, bestCollisionIndex, positiveITS, negativeITS, bachelorITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtOmegaMinusGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtOmegaMinusReconstructed"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0)
              histos.fill(HIST("hPtOmegaMinusReconstructedWithPV"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && bachelorITS)
              histos.fill(HIST("hPtOmegaMinusReconstructedWithPV_ITS"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -3334 && findOmegaPlus) {
          reconstructed = ProcessCascade(mcParticle, tracks, v0s, bestCollisionIndex, positiveITS, negativeITS, bachelorITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtOmegaPlusGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtOmegaPlusReconstructed"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0)
              histos.fill(HIST("hPtOmegaPlusReconstructedWithPV"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveITS && negativeITS && bachelorITS)
              histos.fill(HIST("hPtOmegaPlusReconstructedWithPV_ITS"), mcParticle.pt());
          }
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(casccollisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      if (casccollisionId[ic] >= 0) {
        cascades(casccollisionId[ic], cascv0Index[ic], cascbachelorIndex[ic]);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascademcfinder>(cfgc, TaskName{"lf-cascademcfinder"})};
}
