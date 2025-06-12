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

#include <cmath>
#include <array>
#include <cstdlib>
#include <vector>

#include "Math/Vector4D.h"
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TPDGCode.h>

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
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;

using LabeledTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;
using FullMcParticles = soa::Join<aod::McParticles, aod::ParticlesToTracks>;

struct lambdakzeromcfinder {
  Produces<aod::FindableV0s> v0;
  Produces<aod::McFullV0Labels> fullv0labels;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for selecting which particles to generate
  Configurable<bool> findGamma{"findGamma", true, "findGamma"};
  Configurable<bool> findK0Short{"findK0Short", true, "findK0Short"};
  Configurable<bool> findLambda{"findLambda", true, "findLambda"};
  Configurable<bool> findAntiLambda{"findAntiLambda", true, "findAntiLambda"};
  Configurable<bool> findHyperTriton{"findHyperTriton", false, "findHyperTriton"};
  Configurable<bool> findAntiHyperTriton{"findAntiHyperTriton", false, "findAntiHyperTriton"};
  Configurable<bool> requireTPC{"requireTPC", true, "require TPC"};
  Configurable<bool> skipTPConly{"skipTPConly", false, "skip tracks that are TPC-only"};
  Configurable<bool> storeSingleTPCOnlyProng{"storeSingleTPCOnlyProng", false, "in case a TPC-only track is found, do not allow another TPC-only for the same mcParticle. Works only in MC particle path."};
  Configurable<bool> doAssociatedV0s{"doAssociatedV0s", true, "generate collision-associated V0s (for cascades!)"};
  Configurable<bool> doUnassociatedV0s{"doUnassociatedV0s", true, "generate also unassociated V0s (for cascades, UPC)"};
  Configurable<bool> doSameCollisionOnly{"doSameCollisionOnly", false, "stick to decays in which tracks are assoc to same collision"};
  Configurable<int> qaNbins{"qaNbins", 200, "qa plots: binning"};
  Configurable<float> yPreFilter{"yPreFilter", 2.5, "broad y pre-filter for speed"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};

  Preslice<FullMcParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  std::vector<int> v0collisionId;
  std::vector<int> v0positiveIndex;
  std::vector<int> v0negativeIndex;
  std::vector<int> v0mcLabel;

  std::vector<int> searchedV0PDG;
  std::vector<int> searchedV0PositivePDG;
  std::vector<int> searchedV0NegativePDG;
  std::vector<float> searchedV0PositiveMass;
  std::vector<float> searchedV0NegativeMass;

  void init(InitContext&)
  {
    // initialize histograms
    const AxisSpec axisNTimesRecoed{static_cast<int>(10), -0.5f, +9.5f, ""};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesRecoed});

    // store number of recoed V0s and number of recoed V0s with no collision association
    histos.add("hNCollisionAssociation", "hNCollisionAssociation", kTH1F, {axisNTimesRecoed});

    // warning: this stores (composite) number of copies of tracks
    histos.add("hNTimesRecoedGamma", "hNTimesRecoedGamma", kTH2F, {axisNTimesRecoed, axisPtQA});
    histos.add("hNTimesRecoedK0Short", "hNTimesRecoedK0Short", kTH2F, {axisNTimesRecoed, axisPtQA});
    histos.add("hNTimesRecoedLambda", "hNTimesRecoedLambda", kTH2F, {axisNTimesRecoed, axisPtQA});
    histos.add("hNTimesRecoedAntiLambda", "hNTimesRecoedAntiLambda", kTH2F, {axisNTimesRecoed, axisPtQA});
    histos.add("hNTimesRecoedHypertriton", "hNTimesRecoedHypertriton", kTH2F, {axisNTimesRecoed, axisPtQA});
    histos.add("hNTimesRecoedAntiHypertriton", "hNTimesRecoedAntiHypertriton", kTH2F, {axisNTimesRecoed, axisPtQA});

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

  template <typename TmcParticle>
  int ProcessV0(TmcParticle const& mcParticle, int bestCollisionIndex)
  {
    int nPosReco = 0;
    int nNegReco = 0;
    std::vector<int> trackIndexPositive;
    std::vector<int> trackIndexNegative;

    int positivePdg = 211;
    int negativePdg = -211;
    int relevantProcess = 4; // normal search: decay
    if (mcParticle.pdgCode() == 22) {
      positivePdg = -11;
      negativePdg = +11;
      relevantProcess = 5; // look for pair production if photon
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
      auto const& daughters = mcParticle.template daughters_as<FullMcParticles>();
      if (daughters.size() >= 2) {
        for (auto const& daughter : daughters) { // might be better ways of doing this but ok
          if (daughter.getProcess() != relevantProcess)
            continue; // skip deltarays (if ever), stick to decay products only
          if (daughter.pdgCode() == positivePdg) {
            auto const& thisDaughterTracks = daughter.template tracks_as<LabeledTracks>();
            bool tpcOnlyFound = false;
            for (auto const& track : thisDaughterTracks) {
              if (track.detectorMap() == o2::aod::track::TPC) {
                if (tpcOnlyFound == true && storeSingleTPCOnlyProng)
                  continue; // in case a previous TPC-only version of this mcParticle was found + we want to store only one copy, skip
                if (skipTPConly)
                  continue;
                tpcOnlyFound = true;
              }
              if (track.sign() > 0 && (track.hasTPC() || !requireTPC)) {
                trackIndexPositive.push_back(track.globalIndex()); // assign only if TPC present
                nPosReco++;
              }
            } // end track list loop
          }   // end positive pdg check
          if (daughter.pdgCode() == negativePdg) {
            auto const& thisDaughterTracks = daughter.template tracks_as<LabeledTracks>();
            bool tpcOnlyFound = false;
            for (auto const& track : thisDaughterTracks) {
              if (track.detectorMap() == o2::aod::track::TPC) {
                if (tpcOnlyFound == true && storeSingleTPCOnlyProng)
                  continue; // in case a previous TPC-only version of this mcParticle was found + we want to store only one copy, skip
                if (skipTPConly)
                  continue;
                tpcOnlyFound = true;
              }
              if (track.sign() < 0 && (track.hasTPC() || !requireTPC)) {
                trackIndexNegative.push_back(track.globalIndex()); // assign only if TPC present
                nNegReco++;
              }
            } // end track list loop
          }   // end negative pdg check
        }
      }
    }
    // account for track duplicates
    int reconstructed = 0;
    for (int ip = 0; ip < nPosReco; ip++) {
      for (int in = 0; in < nNegReco; in++) {
        reconstructed++;
        v0collisionId.emplace_back(bestCollisionIndex);
        v0positiveIndex.emplace_back(trackIndexPositive[ip]);
        v0negativeIndex.emplace_back(trackIndexNegative[in]);
        v0mcLabel.emplace_back(mcParticle.globalIndex());
      }
    }
    return reconstructed;
  }

  void processFromMcParticles(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, LabeledTracks const& /*tracks*/, FullMcParticles const& allMcParticles)
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
      for (auto& mcParticle : mcParticles) {
        if (fabs(mcParticle.y()) > yPreFilter)
          continue; // go declarative at a later stage but pre-filter here

        if (mcParticle.pdgCode() == 22 && findGamma) {
          histos.fill(HIST("hNTimesRecoedGamma"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == 310 && findK0Short) {
          histos.fill(HIST("hNTimesRecoedK0Short"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == 3122 && findLambda) {
          histos.fill(HIST("hNTimesRecoedLambda"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == -3122 && findAntiLambda) {
          histos.fill(HIST("hNTimesRecoedAntiLambda"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == 1010010030 && findHyperTriton) {
          histos.fill(HIST("hNTimesRecoedHypertriton"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == -1010010030 && findAntiHyperTriton) {
          histos.fill(HIST("hNTimesRecoedAntiHypertriton"), ProcessV0(mcParticle, bestCollisionIndex), mcParticle.pt());
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(v0collisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      histos.fill(HIST("hNCollisionAssociation"), 0.0f); // any correctly recoed
      if (v0collisionId[ic] >= 0)
        histos.fill(HIST("hNCollisionAssociation"), 1.0f); // reconstructed with a collision associated to it

      if (v0collisionId[ic] < 0 && doUnassociatedV0s) {
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic], 1);
        fullv0labels(v0mcLabel[ic]);
      }
      if (v0collisionId[ic] >= 0 && doAssociatedV0s) {
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic], 1);
        fullv0labels(v0mcLabel[ic]);
      }
    }
  }

  // this improved process function does not start from MC particles, as that would be too costly
  // rather, it starts from appropriately detected prongs of the correct charge.
  Partition<LabeledTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<LabeledTracks> negTracks = aod::track::signed1Pt < 0.0f;

  void processFromSingleProngs(aod::Collisions const& /*collisions*/, LabeledTracks const& /*tracks*/, soa::Join<aod::McCollisions, aod::McCollsExtra> const& /*mcCollisions*/, aod::McParticles const& /*allMcParticles*/)
  {
    v0collisionId.clear();
    v0positiveIndex.clear();
    v0negativeIndex.clear();
    v0mcLabel.clear();

    // This will take place once per TF!
    for (auto& posTrack : posTracks) { //<- no grouping, deliberately
      int v0pdgIndex = -1;
      int motherIndex = -1;
      int bestCollisionIndex = -1;
      if (!posTrack.has_mcParticle())
        continue; // skip unindexed particles
      if (posTrack.detectorMap() == o2::aod::track::TPC && skipTPConly)
        continue;
      if (!posTrack.hasTPC() && requireTPC)
        continue; // skip particles without TPC
      auto posParticle = posTrack.mcParticle_as<aod::McParticles>();
      if (posParticle.getProcess() != 4)
        continue; // skip particles not coming from a decay
      if (!posParticle.has_mothers())
        continue; // skip particles without decay mothers
      for (auto& posMotherParticle : posParticle.mothers_as<aod::McParticles>()) {
        // determine if mother particle satisfies any condition curently being searched for
        for (std::size_t ipdg = 0; ipdg < searchedV0PDG.size(); ipdg++)
          if (searchedV0PDG[ipdg] == posMotherParticle.pdgCode() && fabs(posMotherParticle.y()) < yPreFilter) {
            v0pdgIndex = ipdg; // index mapping to desired V0 species
            motherIndex = posMotherParticle.globalIndex();

            // de-reference best collision
            auto mcCollision = posMotherParticle.mcCollision_as<soa::Join<aod::McCollisions, aod::McCollsExtra>>();
            if (mcCollision.numRecoCollision())
              bestCollisionIndex = mcCollision.bestCollisionIndex();
            continue;
          }
        if (v0pdgIndex < 0 || posParticle.pdgCode() != searchedV0PositivePDG[v0pdgIndex])
          continue; // not interesting, skip

        // if we got here, we need to search for the other prong
        for (auto& negTrack : negTracks) { //<- no grouping, deliberately
          if (doSameCollisionOnly && negTrack.collisionId() != posTrack.collisionId())
            continue; // skip if requested to look only at the same collision (fixme: could be better)
          if (!negTrack.has_mcParticle())
            continue; // skip unindexed particles
          if (negTrack.detectorMap() == o2::aod::track::TPC && skipTPConly)
            continue;
          if (!negTrack.hasTPC() && requireTPC)
            continue; // skip particles without TPC
          auto negParticle = negTrack.mcParticle_as<aod::McParticles>();
          if (negParticle.getProcess() != 4)
            continue; // skip particles not coming from a decay
          if (!negParticle.has_mothers())
            continue; // skip particles without decay mothers
          for (auto& negMotherParticle : negParticle.mothers_as<aod::McParticles>()) {
            if (negMotherParticle.globalIndex() == posMotherParticle.globalIndex() && negParticle.pdgCode() == searchedV0NegativePDG[v0pdgIndex]) {

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
        v0(v0collisionId[ic], v0positiveIndex[ic], v0negativeIndex[ic], 1);
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
