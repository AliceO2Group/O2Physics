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
/// \file    alice3-dilepton.cxx
/// \author  s.scheid@cern.ch, daiki.sekihata@cern.ch
///

#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/prefilterDilepton.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>

#include "Math/GenVector/VectorUtil.h"
#include <Math/Vector4D.h>

#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3DileptonEventCentSelection {

  Service<o2::framework::O2DatabasePDG> inspdg;

  Produces<aod::DiEventCentCuts> eventCentSel;
  Configurable<float> minFITPart{"minFITPart", 0.f, "Minimum number of charged particles in the FIT acceptance"};
  Configurable<float> maxFITPart{"maxFITPart", 0.f, "Maximum number of charged particles in the FIT acceptance"};

  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    registry.add("Generated/Before/ParticlesFit", "Charged Particles in Fit acceptance per event", kTH1F, {{15000, 0, 15000}});
    registry.add("Generated/After/ParticlesFit", "Charged Particles in Fit acceptance per event", kTH1F, {{15000, 0, 15000}});
    registry.add("Generated/Before/ParticlesEta", "Charged Particles per event", kTH1F, {{1200, -6, 6}});
    registry.add("Generated/After/ParticlesEta", "Charged Particles per event", kTH1F, {{1200, -6, 6}});
  }

  using MyEvents = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;

  void processGen(MyEvents::iterator const& event, o2::aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    if (!event.has_mcCollision()) {
      eventCentSel(0);
    } else {
      auto mccollision = event.mcCollision();
      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());
      int nParticlesFIT = 0;
      for (const auto& mcParticle : mcParticles_per_coll) {
        if (mcParticle.isPhysicalPrimary()) {
          auto pdgParticle = inspdg->GetParticle(mcParticle.pdgCode());
          if (pdgParticle) {
            float charge = pdgParticle->Charge() / 3.f; // Charge in units of |e|
            if (std::abs(charge) >= 1.) {
              registry.fill(HIST("Generated/Before/ParticlesEta"), mcParticle.eta());
              if ((2.2 < mcParticle.eta() && mcParticle.eta() < 5.0) || (-3.4 < mcParticle.eta() && mcParticle.eta() < -2.3)) {
                nParticlesFIT++;
              }
            }
          }
        }
      } // end of mc particle loop
      registry.fill(HIST("Generated/Before/ParticlesFit"), nParticlesFIT);
      if (nParticlesFIT > minFITPart && nParticlesFIT < maxFITPart) {
        registry.fill(HIST("Generated/After/ParticlesFit"), nParticlesFIT);
        for (const auto& mcParticle : mcParticles_per_coll) {
          if (mcParticle.isPhysicalPrimary()) {
            auto pdgParticle = inspdg->GetParticle(mcParticle.pdgCode());
            if (pdgParticle) {
              float charge = pdgParticle->Charge() / 3.f; // Charge in units of |e|
              if (std::abs(charge) >= 1.) {
                registry.fill(HIST("Generated/After/ParticlesEta"), mcParticle.eta());
              }
            }
          }
        } // end of mc particle loop
        eventCentSel(1);
      } else {
        eventCentSel(0);
      }
    } // if mc collision loop
  } // end of process

  void processDummy(MyEvents::iterator const&)
  {
    eventCentSel(1);
  }

  PROCESS_SWITCH(Alice3DileptonEventCentSelection, processGen, "Select centrality", false);
  PROCESS_SWITCH(Alice3DileptonEventCentSelection, processDummy, "Dummy", true);
};

struct Alice3DileptonPrefilter {

  Produces<aod::DiTrackPrefilter> pfb_derived;

  SliceCache cache_mc;
  SliceCache cache_rec;

  Service<o2::framework::O2DatabasePDG> inspdg;

  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  Configurable<float> etaMin{"eta-min", -5.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 5.f, "Upper limit in eta"};
  Configurable<float> maxMass{"max-mass", 5.f, "Upper limit in mass"};
  Configurable<float> maxOp{"max-Op", 5.f, "Upper limit in opening angle"};
  Configurable<float> nSigmaEleCutOuterTOF{"nSigmaEleCutOuterTOF", 3., "Electron inclusion in outer TOF"};
  Configurable<float> nSigmaEleCutInnerTOF{"nSigmaEleCutInnerTOF", 3., "Electron inclusion in inner TOF"};
  Configurable<float> nSigmaPionCutOuterTOF{"nSigmaPionCutOuterTOF", 3., "Pion exclusion in outer TOF"};
  Configurable<float> nSigmaPionCutInnerTOF{"nSigmaPionCutInnerTOF", 3., "Pion exclusion in inner TOF"};
  Configurable<float> nSigmaElectronRich{"nSigmaElectronRich", 3., "Electron inclusion RICH"};
  Configurable<float> nSigmaPionRich{"nSigmaPionRich", 3., "Pion exclusion RICH"};

  std::unordered_map<int, uint16_t> map_pfb; // map track.globalIndex -> prefilter bit

  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using MyTracksMCs = soa::Join<aod::Tracks, aod::UpgradeTofs, aod::UpgradeTofMCs, aod::UpgradeRichs, aod::TracksAlice3>;
  using MyTracksMC = MyTracksMCs::iterator;
  using DileptonCollisions = soa::Join<aod::Collisions, aod::DiEventCentCuts>;
  using DileptonCollision = DileptonCollisions::iterator;

  Preslice<MyTracksMCs> perCollision = aod::track::collisionId;
  Partition<MyTracksMCs> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyTracksMCs> negTracks = o2::aod::track::signed1Pt < 0.f;

  void init(InitContext&)
  {
    const AxisSpec axisM{500, 0, 5, "#it{m}_{ll} (GeV/#it{c}^{2})"};
    const AxisSpec axisPt{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisSigmaEl{200, -10, 10, "n#sigma_{El}"};
    const AxisSpec axisTrackLengthOuterTOF{300, 0., 300., "Track length (cm)"};
    const AxisSpec axisEta{1000, -5, 5, "#it{#eta}"};
    const AxisSpec axisPhi{360, 0, TMath::TwoPi(), "#it{#varphi} (rad.)"};
    const AxisSpec axisProdx{2000, -100, 100, "Prod. Vertex X (cm)"};
    const AxisSpec axisPrody{2000, -100, 100, "Prod. Vertex Y (cm)"};
    const AxisSpec axisProdz{2000, -100, 100, "Prod. Vertex Z (cm)"};

    registry.add("Reconstructed/Pair/ULS/Mass_Pt", "Pair Mass vs. Pt", kTH2F, {axisM, axisPt}, true);
    registry.add("ReconstructedFiltered/Pair/ULS/Mass_Pt", "Pair Mass vs. Pt", kTH2F, {axisM, axisPt}, true);
    registry.add("Reconstructed/Track/Pt", "Track Pt", kTH1F, {axisPt});
    registry.add("Reconstructed/Track/Eta", "Track Eta", kTH1F, {axisEta});
    registry.add("Reconstructed/Track/Phi", "Track Phi", kTH1F, {axisPhi});
    registry.add("Reconstructed/Track/Eta_Pt", "Eta vs. Pt", kTH2F, {axisPt, axisEta}, true);
    registry.add("Reconstructed/Track/SigmaOTofvspt", "Track #sigma oTOF", kTH2F, {axisPt, axisSigmaEl});
    registry.add("Reconstructed/Track/SigmaITofvspt", "Track #sigma iTOF", kTH2F, {axisPt, axisSigmaEl});
    registry.add("Reconstructed/Track/SigmaRichvspt", "Track #sigma RICH", kTH2F, {axisPt, axisSigmaEl});
    registry.add("Reconstructed/Track/outerTOFTrackLength", "Track length outer TOF", kTH1F, {axisTrackLengthOuterTOF});
  }

  void processPreFilter(DileptonCollisions const& collisions, MyTracksMCs const& tracks)
  {
    for (const auto& track : tracks) {
      map_pfb[track.globalIndex()] = 0;
    } // end of track loop

    Int_t countercoll = 0;
    for (const auto& collision : collisions) {
      countercoll++;

      auto negTracks_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache_rec);
      auto posTracks_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache_rec);

      if (collision.isEventCentSelected() == 0) {
        for (const auto& pos : posTracks_coll) {
          map_pfb[pos.globalIndex()] = 0;
        }
        for (const auto& neg : negTracks_coll) {
          map_pfb[neg.globalIndex()] = 0;
        }
        continue;
      }

      for (const auto& pos : posTracks_coll) {
        if (!pos.isReconstructed()) {
          continue;
        }
        if (pos.eta() < etaMin || etaMax < pos.eta()) {
          continue;
        }
        if (pos.pt() < ptMin || ptMax < pos.pt()) {
          continue;
        }
        if ((std::abs(pos.nSigmaElectronRich()) < nSigmaElectronRich && nSigmaPionRich < std::abs(pos.nSigmaPionRich())) || (std::abs(pos.nSigmaElectronOuterTOF()) < nSigmaEleCutOuterTOF && nSigmaPionCutOuterTOF < std::abs(pos.nSigmaPionOuterTOF())) || (std::abs(pos.nSigmaElectronInnerTOF()) < nSigmaEleCutInnerTOF && nSigmaPionCutInnerTOF < std::abs(pos.nSigmaPionInnerTOF()))) {
          registry.fill(HIST("Reconstructed/Track/SigmaOTofvspt"), pos.pt(), pos.nSigmaElectronOuterTOF());
          registry.fill(HIST("Reconstructed/Track/SigmaITofvspt"), pos.pt(), pos.nSigmaElectronInnerTOF());
          registry.fill(HIST("Reconstructed/Track/SigmaRichvspt"), pos.pt(), pos.nSigmaElectronRich());
          registry.fill(HIST("Reconstructed/Track/outerTOFTrackLength"), pos.outerTOFTrackLength());
          registry.fill(HIST("Reconstructed/Track/Pt"), pos.pt());
          registry.fill(HIST("Reconstructed/Track/Eta"), pos.eta());
          registry.fill(HIST("Reconstructed/Track/Phi"), pos.phi());
          registry.fill(HIST("Reconstructed/Track/Eta_Pt"), pos.pt(), pos.eta());
        }
      }
      for (const auto& ele : negTracks_coll) {
        if (!ele.isReconstructed()) {
          continue;
        }
        if (ele.eta() < etaMin || etaMax < ele.eta()) {
          continue;
        }
        if (ele.pt() < ptMin || ptMax < ele.pt()) {
          continue;
        }
        if ((std::abs(ele.nSigmaElectronRich()) < nSigmaElectronRich && nSigmaPionRich < std::abs(ele.nSigmaPionRich())) || (std::abs(ele.nSigmaElectronOuterTOF()) < nSigmaEleCutOuterTOF && nSigmaPionCutOuterTOF < std::abs(ele.nSigmaPionOuterTOF())) || (std::abs(ele.nSigmaElectronInnerTOF()) < nSigmaEleCutInnerTOF && nSigmaPionCutInnerTOF < std::abs(ele.nSigmaPionInnerTOF()))) {
          registry.fill(HIST("Reconstructed/Track/SigmaOTofvspt"), ele.pt(), ele.nSigmaElectronOuterTOF());
          registry.fill(HIST("Reconstructed/Track/SigmaITofvspt"), ele.pt(), ele.nSigmaElectronInnerTOF());
          registry.fill(HIST("Reconstructed/Track/SigmaRichvspt"), ele.pt(), ele.nSigmaElectronRich());
          registry.fill(HIST("Reconstructed/Track/outerTOFTrackLength"), ele.outerTOFTrackLength());
          registry.fill(HIST("Reconstructed/Track/Pt"), ele.pt());
          registry.fill(HIST("Reconstructed/Track/Eta"), ele.eta());
          registry.fill(HIST("Reconstructed/Track/Phi"), ele.phi());
          registry.fill(HIST("Reconstructed/Track/Eta_Pt"), ele.pt(), ele.eta());
        }
      }

      for (const auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_coll, negTracks_coll))) { // ULS
        if (!pos.isReconstructed()) {
          continue;
        }
        if (pos.eta() < etaMin || etaMax < pos.eta()) {
          continue;
        }
        if (pos.pt() < ptMin || ptMax < pos.pt()) {
          continue;
        }
        if (!ele.isReconstructed()) {
          continue;
        }
        if (ele.eta() < etaMin || etaMax < ele.eta()) {
          continue;
        }
        if (ele.pt() < ptMin || ptMax < ele.pt()) {
          continue;
        }
        if ((std::abs(pos.nSigmaElectronRich()) < nSigmaElectronRich && nSigmaPionRich < std::abs(pos.nSigmaPionRich())) || (std::abs(pos.nSigmaElectronOuterTOF()) < nSigmaEleCutOuterTOF && nSigmaPionCutOuterTOF < std::abs(pos.nSigmaPionOuterTOF())) || (std::abs(pos.nSigmaElectronInnerTOF()) < nSigmaEleCutInnerTOF && nSigmaPionCutInnerTOF < std::abs(pos.nSigmaPionInnerTOF()))) {
          if ((std::abs(ele.nSigmaElectronRich()) < nSigmaElectronRich && nSigmaPionRich < std::abs(ele.nSigmaPionRich())) || (std::abs(ele.nSigmaElectronOuterTOF()) < nSigmaEleCutOuterTOF && nSigmaPionCutOuterTOF < std::abs(ele.nSigmaPionOuterTOF())) || (std::abs(ele.nSigmaElectronInnerTOF()) < nSigmaEleCutInnerTOF && nSigmaPionCutInnerTOF < std::abs(ele.nSigmaPionInnerTOF()))) {
            ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
            ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
            float angle = ROOT::Math::VectorUtil::Angle(v1, v2);
            o2::math_utils::bringToPMPi(angle);

            registry.fill(HIST("Reconstructed/Pair/ULS/Mass_Pt"), v12.M(), v12.Pt());

            if (v12.M() < maxMass && angle < maxOp) {
              map_pfb[pos.globalIndex()] = 1;
              map_pfb[ele.globalIndex()] = 1;
              registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Mass_Pt"), v12.M(), v12.Pt());
            }
          }
        }
      }
    } // end of collision

    Int_t counter = 0;
    for (const auto& track : tracks) {
      counter++;
      pfb_derived(map_pfb[track.globalIndex()]);
    } // end of track loop
    map_pfb.clear();
  }

  void processDummy(MyTracksMCs const& tracks)
  {
    Int_t counter = 0;
    for (const auto& track : tracks) {
      counter++;
      pfb_derived(0);
    }
  }

  PROCESS_SWITCH(Alice3DileptonPrefilter, processPreFilter, "Run prefilter", false);
  PROCESS_SWITCH(Alice3DileptonPrefilter, processDummy, "Dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3DileptonEventCentSelection>(cfgc, TaskName{"alice3-dilepton-event-cent-selection"}),
    adaptAnalysisTask<Alice3DileptonPrefilter>(cfgc, TaskName{"alice3-dilepton-prefilter"})};
}
