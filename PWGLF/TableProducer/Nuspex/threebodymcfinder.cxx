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
// 3-body MC finder task
// -----------------
//
//    This task allows for the re-creation of the Decay3Body table with
//    perfect MC information. It is meant to be used to understand
//    baseline svertexer efficiency and to allow for tuning
//    of the 3-body decay reconstruction using MC.
//
//    Nota bene: special attention could still be dedicated to the PV
//               reconstruction efficiency.
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace ROOT::Math;

using LabeledTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

struct threebodymcfinder {
  Produces<aod::Decay3Bodys> d3b;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for selecting which particles to generate
  Configurable<bool> findHyHe4{"findHyHe4", true, "findHyHe4"};
  Configurable<bool> findHypertriton{"findHypertriton", false, "findHypertriton"};

  Configurable<bool> requireITS{"requireITS", false, "require ITS information used in tracks"};

  Preslice<aod::McParticle> perMcCollision = aod::mcparticle::mcCollisionId;

  std::vector<int> d3bcollisionId;
  std::vector<int> d3bprong0Index;
  std::vector<int> d3bprong1Index;
  std::vector<int> d3bprong2Index;

  void init(InitContext&)
  {
    // initialize histograms
    const AxisSpec axisNTimesCollRecoed{(int)10, -0.5f, +9.5f, ""};
    const AxisSpec axisPt{(int)100, +0.0f, +10.0f, "p_{T} (GeV/c)"};

    histos.add("hNTimesCollRecoed", "hNTimesCollRecoed", kTH1F, {axisNTimesCollRecoed});
    histos.add("hNDaughters", "hNDaughters", kTH1F, {axisNTimesCollRecoed});

    histos.add("hPtHyHe4Generated", "hPtHyHe4Generated", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4Generated", "hPtAntiHyHe4Generated", kTH1F, {axisPt});
    histos.add("hPtHypertritonGenerated", "hPtHypertritonGenerated", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonGenerated", "hPtAntiHypertritonGenerated", kTH1F, {axisPt});

    histos.add("hPtHyHe4Reconstructed", "hPtHyHe4Reconstructed", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4Reconstructed", "hPtAntiHyHe4Reconstructed", kTH1F, {axisPt});
    histos.add("hPtHypertritonReconstructed", "hPtHypertritonReconstructed", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonReconstructed", "hPtAntiHypertritonReconstructed", kTH1F, {axisPt});

    histos.add("hPtHyHe4Global", "hPtHyHe4Global", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4Global", "hPtAntiHyHe4Global", kTH1F, {axisPt});
    histos.add("hPtHypertritonGlobal", "hPtHypertritonGlobal", kTH1F, {axisPt});
    histos.add("hPtAntiHypertritonGlobal", "hPtAntiHypertritonGlobal", kTH1F, {axisPt});

    histos.add("hPtHyHe4GlobalWithPV", "hPtHyHe4GlobalWithPV", kTH1F, {axisPt});
    histos.add("hPtAntiHyHe4GlobalWithPV", "hPtAntiHyHe4GlobalWithPV", kTH1F, {axisPt});
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
  bool ProcessThreeBody(TmcParticle const& mcParticle, TTrackList const& trackList, int bestCollisionIndex, bool& prong0ITS, bool& prong1ITS, bool& prong2ITS, bool& prong0TPC, bool& prong1TPC, bool& prong2TPC)
  {
    bool reconstructed = false;
    prong0ITS = false;
    prong1ITS = false;
    prong2ITS = false;
    prong0TPC = false;
    prong1TPC = false;
    prong2TPC = false;
    int trackIndexProng0 = -1, trackIndexProng1 = -1, trackIndexProng2 = -1;
    int prong0pdg = -1, prong1pdg = -1, prong2pdg = -1;

    // expected daughter pdgs
    if (mcParticle.pdgCode() == +1010020040) {
      prong0pdg = +1000020030;
      prong1pdg = +2212;
      prong2pdg = -211;
    }
    if (mcParticle.pdgCode() == -1010020040) {
      prong0pdg = -1000020030;
      prong1pdg = -2212;
      prong2pdg = +211;
    }
    if (mcParticle.pdgCode() == +1010010030) {
      prong0pdg = +1000010020;
      prong1pdg = +2212;
      prong2pdg = -211;
    }
    if (mcParticle.pdgCode() == -1010010030) {
      prong0pdg = -1000010020;
      prong1pdg = -2212;
      prong2pdg = +211;
    }

    if (mcParticle.has_daughters()) {
      auto const& daughters = mcParticle.template daughters_as<aod::McParticles>();
      histos.fill(HIST("hNDaughters"), daughters.size());
      if (daughters.size() >= 3) {               // consider also delta-rays
        for (auto const& daughter : daughters) { // might be better ways of doing this but ok
          for (auto const& track : trackList) {
            if (track.mcParticleId() == daughter.globalIndex()) {
              // determine which charge this particle has
              if (daughter.pdgCode() == prong0pdg && daughter.getProcess() == 4) {
                trackIndexProng0 = track.globalIndex();
                if (track.hasITS())
                  prong0ITS = true;
                if (track.hasTPC())
                  prong0TPC = true;
              }
              if (daughter.pdgCode() == prong1pdg && daughter.getProcess() == 4) {
                trackIndexProng1 = track.globalIndex();
                if (track.hasITS())
                  prong1ITS = true;
                if (track.hasTPC())
                  prong1TPC = true;
              }
              if (daughter.pdgCode() == prong2pdg && daughter.getProcess() == 4) {
                trackIndexProng2 = track.globalIndex();
                if (track.hasITS())
                  prong2ITS = true;
                if (track.hasTPC())
                  prong2TPC = true;
              }
              if (trackIndexProng0 >= 0 && trackIndexProng1 >= 0 && trackIndexProng2 >= 0)
                break; // all found
            }
          }
        }
      }
    }
    if (trackIndexProng0 >= 0 && trackIndexProng1 >= 0 && trackIndexProng2 >= 0 && (!requireITS || (requireITS && prong0ITS && prong1ITS && prong2ITS))) {
      reconstructed = true;
      d3bcollisionId.emplace_back(bestCollisionIndex);
      d3bprong0Index.emplace_back(trackIndexProng0);
      d3bprong1Index.emplace_back(trackIndexProng1);
      d3bprong2Index.emplace_back(trackIndexProng2);
    }
    return reconstructed;
  }

  void process(soa::Join<aod::McCollisions, aod::McCollsExtra> const& mcCollisions, LabeledTracks const& tracks, aod::McParticles const& allMcParticles)
  {
    d3bcollisionId.clear();
    d3bprong0Index.clear();
    d3bprong1Index.clear();
    d3bprong2Index.clear();

    // Step 1: sweep over all mcCollisions and find all relevant candidates
    for (auto const& mcCollision : mcCollisions) {
      histos.fill(HIST("hNTimesCollRecoed"), mcCollision.numRecoCollision());
      int bestCollisionIndex = mcCollision.bestCollisionIndex();

      auto mcParticles = allMcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      bool prong0ITS = false, prong1ITS = false, prong2ITS = false;
      bool prong0TPC = false, prong1TPC = false, prong2TPC = false;
      bool reconstructed = false;
      for (auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() == 1010020040 && findHyHe4) {
          reconstructed = ProcessThreeBody(mcParticle, tracks, bestCollisionIndex, prong0ITS, prong1ITS, prong2ITS, prong0TPC, prong1TPC, prong2TPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtHyHe4Generated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtHyHe4Reconstructed"), mcParticle.pt());
            if (reconstructed && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtHyHe4Global"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtHyHe4GlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -1010020040 && findHyHe4) {
          reconstructed = ProcessThreeBody(mcParticle, tracks, bestCollisionIndex, prong0ITS, prong1ITS, prong2ITS, prong0TPC, prong1TPC, prong2TPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiHyHe4Generated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiHyHe4Reconstructed"), mcParticle.pt());
            if (reconstructed && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtAntiHyHe4Global"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtAntiHyHe4GlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == 1010010030 && findHypertriton) {
          reconstructed = ProcessThreeBody(mcParticle, tracks, bestCollisionIndex, prong0ITS, prong1ITS, prong2ITS, prong0TPC, prong1TPC, prong2TPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
        if (mcParticle.pdgCode() == -1010010030 && findHypertriton) {
          reconstructed = ProcessThreeBody(mcParticle, tracks, bestCollisionIndex, prong0ITS, prong1ITS, prong2ITS, prong0TPC, prong1TPC, prong2TPC);
          if (fabs(mcParticle.y()) < 0.5) {
            histos.fill(HIST("hPtAntiHyperTritonGenerated"), mcParticle.pt());
            if (reconstructed)
              histos.fill(HIST("hPtAntiHyperTritonReconstructed"), mcParticle.pt());
            if (reconstructed && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtAntiHyperTritonGlobal"), mcParticle.pt());
            if (reconstructed && bestCollisionIndex >= 0 && prong0ITS && prong1ITS && prong2ITS && prong0TPC && prong1TPC && prong2TPC)
              histos.fill(HIST("hPtAntiHyperTritonGlobalWithPV"), mcParticle.pt());
          }
        }
      }
    }

    // sort according to collision ID
    auto sortedIndices = sort_indices(d3bcollisionId);

    // V0 list established, populate
    for (auto ic : sortedIndices) {
      if (d3bcollisionId[ic] >= 0) {
        d3b(d3bcollisionId[ic], d3bprong0Index[ic], d3bprong1Index[ic], d3bprong2Index[ic]);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<threebodymcfinder>(cfgc, TaskName{"lf-threebodymcfinder"})};
}
