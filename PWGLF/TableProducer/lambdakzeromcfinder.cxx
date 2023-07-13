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

#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"

#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

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
  Configurable<bool> findGamma{"findGamma", true, "findGamma"};
  Configurable<bool> findK0Short{"findK0Short", true, "findK0Short"};
  Configurable<bool> findLambda{"findLambda", true, "findLambda"};
  Configurable<bool> findAntiLambda{"findAntiLambda", true, "findAntiLambda"};
  Configurable<bool> findHyperTriton{"findHyperTriton", false, "findHyperTriton"};
  Configurable<bool> findAntiHyperTriton{"findAntiHyperTriton", false, "findAntiHyperTriton"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS information used in tracks"};
  Configurable<bool> doUnassociatedV0s{"doUnassociatedV0s", true, "generate also unassociated V0s (for cascades!)"};
  Configurable<bool> doQA{"doQA", true, "do qa plots"};
  Configurable<bool> doSameCollisionOnly{"doSameCollisionOnly", false, "stick to decays in which tracks are assoc to same collision"};
  Configurable<int> qaNbins{"qaNbins", 200, "qa plots: binning"};
  Configurable<float> yPreFilter{"yPreFilter", 2.5, "broad y pre-filter for speed"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Preslice<aod::McParticle> perMcCollision = aod::mcparticle::mcCollisionId;

  std::vector<int> v0collisionId;
  std::vector<int> v0positiveIndex;
  std::vector<int> v0negativeIndex;
  std::vector<int> v0mcLabel;

  std::vector<int> searchedV0PDG;
  std::vector<int> searchedV0PositivePDG;
  std::vector<int> searchedV0NegativePDG;
  std::vector<float> searchedV0PositiveMass;
  std::vector<float> searchedV0NegativeMass;

  void init(InitContext& context)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{static_cast<int>(10), -0.5f, +9.5f, ""};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});

    histos.add("hPtGammaGenerated", "hPtGammaGenerated", kTH1F, {axisPtQA});
    histos.add("hPtK0ShortGenerated", "hPtK0ShortGenerated", kTH1F, {axisPtQA});
    histos.add("hPtLambdaGenerated", "hPtLambdaGenerated", kTH1F, {axisPtQA});
    histos.add("hPtAntiLambdaGenerated", "hPtAntiLambdaGenerated", kTH1F, {axisPtQA});
    histos.add("hPtHypertritonGenerated", "hPtHypertritonGenerated", kTH1F, {axisPtQA});
    histos.add("hPtAntiHypertritonGenerated", "hPtAntiHypertritonGenerated", kTH1F, {axisPtQA});

    histos.add("hPtGammaReconstructed", "hPtGammaReconstructed", kTH1F, {axisPtQA});
    histos.add("hPtK0ShortReconstructed", "hPtK0ShortReconstructed", kTH1F, {axisPtQA});
    histos.add("hPtLambdaReconstructed", "hPtLambdaReconstructed", kTH1F, {axisPtQA});
    histos.add("hPtAntiLambdaReconstructed", "hPtAntiLambdaReconstructed", kTH1F, {axisPtQA});
    histos.add("hPtHypertritonReconstructed", "hPtHypertritonReconstructed", kTH1F, {axisPtQA});
    histos.add("hPtAntiHypertritonReconstructed", "hPtAntiHypertritonReconstructed", kTH1F, {axisPtQA});

    histos.add("hPtGammaGlobal", "hPtGammaGlobal", kTH1F, {axisPtQA});
    histos.add("hPtK0ShortGlobal", "hPtK0ShortGlobal", kTH1F, {axisPtQA});
    histos.add("hPtLambdaGlobal", "hPtLambdaGlobal", kTH1F, {axisPtQA});
    histos.add("hPtAntiLambdaGlobal", "hPtAntiLambdaGlobal", kTH1F, {axisPtQA});
    histos.add("hPtHypertritonGlobal", "hPtHypertritonGlobal", kTH1F, {axisPtQA});
    histos.add("hPtAntiHypertritonGlobal", "hPtAntiHypertritonGlobal", kTH1F, {axisPtQA});

    histos.add("hPtGammaGlobalWithPV", "hPtGammaGlobalWithPV", kTH1F, {axisPtQA});
    histos.add("hPtK0ShortGlobalWithPV", "hPtK0ShortGlobalWithPV", kTH1F, {axisPtQA});
    histos.add("hPtLambdaGlobalWithPV", "hPtLambdaGlobalWithPV", kTH1F, {axisPtQA});
    histos.add("hPtAntiLambdaGlobalWithPV", "hPtAntiLambdaGlobalWithPV", kTH1F, {axisPtQA});
    histos.add("hPtHypertritonGlobalWithPV", "hPtHypertritonGlobalWithPV", kTH1F, {axisPtQA});
    histos.add("hPtAntiHypertritonGlobalWithPV", "hPtAntiHypertritonGlobalWithPV", kTH1F, {axisPtQA});

    if (doQA) {
      histos.add("hPtGammaDaughters", "hPtGammaDaughters", kTH2F, {axisPtQA, axisPtQA});
      histos.add("hPtK0ShortDaughters", "hPtK0ShortDaughters", kTH2F, {axisPtQA, axisPtQA});
      histos.add("hPtLambdaDaughters", "hPtLambdaDaughters", kTH2F, {axisPtQA, axisPtQA});
      histos.add("hPtAntiLambdaDaughters", "hPtAntiLambdaDaughters", kTH2F, {axisPtQA, axisPtQA});
      histos.add("hPtHypertritonDaughters", "hPtHypertritonDaughters", kTH2F, {axisPtQA, axisPtQA});
      histos.add("hPtAntiHypertritonDaughters", "hPtAntiHypertritonDaughters", kTH2F, {axisPtQA, axisPtQA});
    }

    // initialise search vectors
    if (findGamma) {
      searchedV0PDG.emplace_back(22);
      searchedV0PositivePDG.emplace_back(-11);
      searchedV0NegativePDG.emplace_back(+11);
      searchedV0PositiveMass.emplace_back(0.0f);
      searchedV0NegativeMass.emplace_back(0.0f);
    }
    if (findK0Short) {
      searchedV0PDG.emplace_back(310);
      searchedV0PositivePDG.emplace_back(+211);
      searchedV0NegativePDG.emplace_back(-211);
      searchedV0PositiveMass.emplace_back(o2::constants::physics::MassPionCharged);
      searchedV0NegativeMass.emplace_back(o2::constants::physics::MassPionCharged);
    }
    if (findLambda) {
      searchedV0PDG.emplace_back(3122);
      searchedV0PositivePDG.emplace_back(+2212);
      searchedV0NegativePDG.emplace_back(-211);
      searchedV0PositiveMass.emplace_back(o2::constants::physics::MassProton);
      searchedV0NegativeMass.emplace_back(o2::constants::physics::MassPionCharged);
    }
    if (findAntiLambda) {
      searchedV0PDG.emplace_back(-3122);
      searchedV0PositivePDG.emplace_back(+211);
      searchedV0NegativePDG.emplace_back(-2212);
      searchedV0PositiveMass.emplace_back(o2::constants::physics::MassPionCharged);
      searchedV0NegativeMass.emplace_back(o2::constants::physics::MassProton);
    }
    if (findHyperTriton) {
      searchedV0PDG.emplace_back(+1010010030);
      searchedV0PositivePDG.emplace_back(+1000020030);
      searchedV0NegativePDG.emplace_back(-211);
      searchedV0PositiveMass.emplace_back(o2::constants::physics::MassHelium3);
      searchedV0NegativeMass.emplace_back(o2::constants::physics::MassPionCharged);
    }
    if (findAntiHyperTriton) {
      searchedV0PDG.emplace_back(-1010010030);
      searchedV0PositivePDG.emplace_back(+211);
      searchedV0NegativePDG.emplace_back(-1000020030);
      searchedV0PositiveMass.emplace_back(o2::constants::physics::MassPionCharged);
      searchedV0NegativeMass.emplace_back(o2::constants::physics::MassHelium3);
    }
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
  bool ProcessV0(TmcParticle const& mcParticle, TTrackList const& trackList, int bestCollisionIndex, bool& positiveITS, bool& negativeITS, bool& positiveTPC, bool& negativeTPC, bool& positiveTPCITS, bool& negativeTPCITS)
  {
    bool reconstructed = false;
    positiveITS = false;
    negativeITS = false;
    positiveTPC = false;
    negativeTPC = false;
    positiveTPCITS = false;
    negativeTPCITS = false;
    int trackIndexPositive = -1;
    int trackIndexNegative = -1;
    float posPt = -1.0f;
    float negPt = -1.0f;

    int positivePdg = 211;
    int negativePdg = -211;
    if (mcParticle.pdgCode() == 22) {
      positivePdg = -11;
      negativePdg = +11;
    }
    if (mcParticle.pdgCode() == 3122) {
      positivePdg = 2212;
    }
    if (mcParticle.pdgCode() == -3122) {
      negativePdg = -2212;
    }
    if (mcParticle.pdgCode() == 1010010030) {
      positivePdg = 1000020030;
    }
    if (mcParticle.pdgCode() == -1010010030) {
      negativePdg = -1000020030;
    }

    if (mcParticle.has_daughters()) {
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      if (daughters.size() >= 2) {
        for (auto const& daughter : daughters) { // might be better ways of doing this but ok
          if ((daughter.getProcess() != 4 && mcParticle.pdgCode() != 22) || daughter.getProcess() != 5)
            continue; // skip deltarays (if ever), stick to decay products only
          if (daughter.pdgCode() == positivePdg) {
            for (auto const& track : trackList) {
              if (track.mcParticleId() == daughter.globalIndex()) {
                if (track.hasITS())
                  positiveITS = true;
                if (track.hasTPC()) {
                  positiveTPC = true;
                  trackIndexPositive = track.globalIndex(); // assign only if TPC present
                  posPt = track.pt();
                  if (track.hasITS())
                    positiveTPCITS = true;
                }
              } // end daughter ID check
            }   // end track list loop
          }     // end positive pdg check
          if (daughter.pdgCode() == negativePdg) {
            for (auto const& track : trackList) {
              if (track.mcParticleId() == daughter.globalIndex()) {
                if (track.hasITS())
                  negativeITS = true;
                if (track.hasTPC()) {
                  negativeTPC = true;
                  trackIndexNegative = track.globalIndex(); // assign only if TPC present
                  negPt = track.pt();
                  if (track.hasITS())
                    negativeTPCITS = true;
                }
              } // end daughter ID check
            }   // end track list loop
          }     // end positive pdg check
        }
      }
    }
    if (trackIndexPositive >= 0 && trackIndexNegative >= 0 && (!requireITS || (requireITS && positiveITS && negativeITS))) {
      reconstructed = true;
      v0collisionId.emplace_back(bestCollisionIndex);
      v0positiveIndex.emplace_back(trackIndexPositive);
      v0negativeIndex.emplace_back(trackIndexNegative);
      v0mcLabel.emplace_back(mcParticle.globalIndex());
      if (doQA) {
        if (mcParticle.pdgCode() == 22)
          histos.fill(HIST("hPtGammaDaughters"), posPt, negPt);
        if (mcParticle.pdgCode() == 310)
          histos.fill(HIST("hPtK0ShortDaughters"), posPt, negPt);
        if (mcParticle.pdgCode() == 3122)
          histos.fill(HIST("hPtLambdaDaughters"), posPt, negPt);
        if (mcParticle.pdgCode() == -3122)
          histos.fill(HIST("hPtAntiLambdaDaughters"), posPt, negPt);
        if (mcParticle.pdgCode() == 1010010030)
          histos.fill(HIST("hPtHypertritonDaughters"), posPt, negPt);
        if (mcParticle.pdgCode() == -1010010030)
          histos.fill(HIST("hPtAntiHypertritonDaughters"), posPt, negPt);
      }
    }
    return reconstructed;
  }

  void processFromMcParticles(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, LabeledTracks const& tracks, aod::McParticles const& allMcParticles)
  {
    v0collisionId.clear();
    v0positiveIndex.clear();
    v0negativeIndex.clear();
    v0mcLabel.clear();

    // Step 1: sweep over all mcCollisions and find all relevant candidates
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hNTimesCollRecoed"), mcCollision.numRecoCollision());
      int bestCollisionIndex = mcCollision.bestCollisionIndex();

      auto mcParticles = allMcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      bool positiveITS = false;
      bool negativeITS = false;
      bool positiveTPC = false;
      bool negativeTPC = false;
      bool positiveTPCITS = false;
      bool negativeTPCITS = false;
      bool reconstructed = false;
      for (auto& mcParticle : mcParticles) {
        if (fabs(mcParticle.y()) > yPreFilter)
          continue; // go declarative at a later stage but pre-filter here

        if (mcParticle.pdgCode() == 22 && findGamma) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtGammaGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtGammaReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtGammaGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtGammaGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 310 && findK0Short) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtK0ShortGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtK0ShortReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtK0ShortGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtK0ShortGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 3122 && findLambda) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtLambdaGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtLambdaReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtLambdaGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtLambdaGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -3122 && findAntiLambda) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiLambdaGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiLambdaReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtAntiLambdaGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtAntiLambdaGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 1010010030 && findHyperTriton) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -1010010030 && findAntiHyperTriton) {
          reconstructed = ProcessV0(mcParticle, tracks, bestCollisionIndex, positiveITS, negativeITS, positiveTPC, negativeTPC, positiveTPCITS, negativeTPCITS);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtAntiHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && positiveTPCITS && negativeTPCITS)
              histos.fill(HIST("hPtAntiHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(v0collisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      if (v0collisionId[ic] >= 0 || doUnassociatedV0s) {
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic]);
        fullv0labels(v0mcLabel[ic]);
      }
    }
  }

  // this improved process function does not start from MC particles, as that would be too costly
  // rather, it starts from appropriately detected prongs of the correct charge.
  Partition<LabeledTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<LabeledTracks> negTracks = aod::track::signed1Pt < 0.0f;

  void processFromSingleProngs(aod::Collisions const& collisions, LabeledTracks const& tracks, soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, aod::McParticles const& allMcParticles)
  {
    v0collisionId.clear();
    v0positiveIndex.clear();
    v0negativeIndex.clear();
    v0mcLabel.clear();

    // This will take place once per TF!
    for (auto& posTrack : posTracks) { //<- no grouping, deliberately
      int v0pdgIndex = -1;
      int motherIndex = -1;
      if (!posTrack.has_mcParticle())
        continue;
      auto posParticle = posTrack.mcParticle_as<aod::McParticles>();
      if (!posParticle.has_mothers())
        continue;
      for (auto& posMotherParticle : posParticle.mothers_as<aod::McParticles>()) {
        // determine if mother particle satisfies any condition curently being searched for
        for (uint16_t ipdg = 0; ipdg < searchedV0PDG.size(); ipdg++)
          if (searchedV0PDG[ipdg] == posMotherParticle.pdgCode()) {
            v0pdgIndex = ipdg; // index mapping to desired V0 species
            motherIndex = posMotherParticle.globalIndex();
            continue;
          }
        if (v0pdgIndex < 0 || posParticle.pdgCode() != searchedV0PositivePDG[v0pdgIndex])
          continue; // not interesting, skip

        // if we got here, we need to search for the other prong
        for (auto& negTrack : negTracks) { //<- no grouping, deliberately
          if (doSameCollisionOnly && negTrack.collisionId() != posTrack.collisionId())
            continue; // skip if requested to look only at the same collision (fixme: could be better)
          if (!negTrack.has_mcParticle())
            continue;
          auto negParticle = negTrack.mcParticle_as<aod::McParticles>();
          if (!negParticle.has_mothers())
            continue;
          for (auto& negMotherParticle : negParticle.mothers_as<aod::McParticles>()) {
            if (negMotherParticle.globalIndex() == posMotherParticle.globalIndex() && negMotherParticle.pdgCode() == searchedV0NegativePDG[v0pdgIndex]) {
              // de-reference best collision
              int bestCollisionIndex = -1;
              auto mcCollision = posParticle.mcCollision_as<soa::Join<aod::McCollisions, aod::McCollsExtra>>();
              if (mcCollision.numRecoCollision())
                bestCollisionIndex = mcCollision.bestCollisionIndex();

              // place in list to be passed along, please
              v0collisionId.emplace_back(bestCollisionIndex);
              v0positiveIndex.emplace_back(posTrack.globalIndex());
              v0negativeIndex.emplace_back(negTrack.globalIndex());
              v0mcLabel.emplace_back(motherIndex);
            }
          }
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(v0collisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      if (v0collisionId[ic] >= 0 || doUnassociatedV0s) {
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic]);
        fullv0labels(v0mcLabel[ic]);
      }
    }
  }

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  /// basic building options (one of them must be chosen)
  PROCESS_SWITCH(lambdakzeromcfinder, processFromMcParticles, "Switch to generate from mc particle list (slower)", false);
  PROCESS_SWITCH(lambdakzeromcfinder, processFromSingleProngs, "Switch to generate from single prong combinations (faster)", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeromcfinder>(cfgc, TaskName{"lf-lambdakzeromcfinder"})};
}
