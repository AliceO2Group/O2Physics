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
// ========================
//
// This code produces reduced events for photon analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using TracksMC = soa::Join<aod::Tracks, aod::McTrackLabels>;
using TrackMC = TracksMC::iterator;

struct createEMReducedMCEvent {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kDalitzEE = 0x8,
    kUndef = -1,
  };
  Produces<o2::aod::EMReducedEvents> events;
  Produces<o2::aod::EMReducedMCEvents> mcevents;
  Produces<o2::aod::EMReducedMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryTrackMCLabels> emprimarytrackmclabels;
  Produces<o2::aod::V0KFEMReducedEventIds> v0kfeventid;
  Produces<o2::aod::DalitzEEEMReducedEventIds> dalitzeventid;
  Produces<o2::aod::PHOSEMReducedEventIds> phoseventid;
  Produces<o2::aod::EMCEMReducedEventIds> emceventid;

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0Photons> perCollision_pcm = aod::v0photon::collisionId;
  Preslice<aod::DalitzEEs> perCollision_dalitz = aod::dalitzee::collisionId;
  Preslice<aod::EMPrimaryTracks> perCollision_emprmtrk = aod::emprimarytrack::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  using MyCollisions = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs, aod::McCollisionLabels>;

  template <uint8_t system, typename TTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TDielectorns, typename TEMPrimaryTracks>
  void skimmingMC(MyCollisions const& collisions, aod::BCs const&, aod::McCollisions const&, aod::McParticles const& mcTracks, TTracks const& o2tracks, TPCMs const& v0photons, TPCMLegs const& v0legs, TPHOSs const& phosclusters, TEMCs const& emcclusters, TDielectorns const& dielectrons, TEMPrimaryTracks const& emprimarytracks)
  {
    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    // std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter

    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);
      int ng_pcm = 0;
      int ng_phos = 0;
      int ng_emc = 0;
      int ng_dilepton = 0;

      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2);

      auto mcCollision = collision.mcCollision();
      // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      bool is_phoscpv_readout = collision.alias_bit(kTVXinPHOS);
      bool is_emc_readout = collision.alias_bit(kTVXinEMC);

      if constexpr (static_cast<bool>(system & kPCM)) {
        auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
        ng_pcm = v0photons_coll.size();
        for (int iv0 = 0; iv0 < v0photons_coll.size(); iv0++) {
          v0kfeventid(events.lastIndex() + 1);
        }
      }
      if constexpr (static_cast<bool>(system & kPHOS)) {
        auto phos_coll = phosclusters.sliceBy(perCollision_phos, collision.globalIndex());
        ng_phos = phos_coll.size();
        for (int iphos = 0; iphos < phos_coll.size(); iphos++) {
          phoseventid(events.lastIndex() + 1);
        }
      }
      if constexpr (static_cast<bool>(system & kEMC)) {
        auto emc_coll = emcclusters.sliceBy(perCollision_emc, collision.globalIndex());
        ng_emc = emc_coll.size();
        for (int iemc = 0; iemc < emc_coll.size(); iemc++) {
          emceventid(events.lastIndex() + 1);
        }
      }
      if constexpr (static_cast<bool>(system & kDalitzEE)) {
        auto dielectrons_coll = dielectrons.sliceBy(perCollision_dalitz, collision.globalIndex());
        ng_dilepton = dielectrons_coll.size();
        for (int iee = 0; iee < dielectrons_coll.size(); iee++) {
          dalitzeventid(events.lastIndex() + 1);
        }
      }

      // store event selection decisions
      uint64_t tag = collision.selection_raw();

      events(collision.globalIndex(), tag, collision.bc().runNumber(), collision.bc().triggerMask(), collision.sel8(),
             is_phoscpv_readout, is_emc_readout,
             collision.posX(), collision.posY(), collision.posZ(),
             collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes(),
             collision.multTPC(), collision.multFV0A(), collision.multFV0C(), collision.multFT0A(), collision.multFT0C(),
             collision.multFDDA(), collision.multFDDC(), collision.multZNA(), collision.multZNC(), collision.multTracklets(), collision.multNTracksPV(), collision.multNTracksPVeta1(),
             collision.centFT0M(), collision.centFT0A(), collision.centFT0C(), collision.centNTPV(),
             ng_pcm, ng_phos, ng_emc);

      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        // mcevents(mcCollision.globalIndex(), mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter());
        mcevents(mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
      }

      mceventlabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

      // store mc particles
      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (auto& mctrack : groupedMcTracks) {
        if (mctrack.pt() < 1e-2 || abs(mctrack.y()) > 1.5 || abs(mctrack.vz()) > 250 || sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) > 260) {
          continue;
        }
        int pdg = mctrack.pdgCode();

        // Note that pi0 from weak decay gives producedByGenerator() = false
        if (
          abs(pdg) != 11      // electron
          && (abs(pdg) != 22) // photon
          // light mesons
          && (abs(pdg) != 111) // pi0
          && (abs(pdg) != 113) // rho(770)
          // && (abs(pdg) != 211) // changed pion
          && (abs(pdg) != 221) // eta
          && (abs(pdg) != 223) // omega(782)
          && (abs(pdg) != 331) // eta'(958)
          && (abs(pdg) != 333) // phi(1020)
          // strange hadrons
          && (abs(pdg) != 310)  // K0S
          && (abs(pdg) != 130)  // K0L
          && (abs(pdg) != 3122) // Lambda
        ) {
          continue;
        }

        // LOGF(info,"index = %d , mc track pdg = %d , producedByGenerator =  %d , isPhysicalPrimary = %d", mctrack.index(), mctrack.pdgCode(), mctrack.producedByGenerator(), mctrack.isPhysicalPrimary());

        // these are used as denominator for efficiency. (i.e. generated information)
        if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mctrack.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
          // fMCFlags[mctrack.globalIndex()] = mcflags;
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
          fCounters[0]++;
        }
      } // end of mc track loop

      if constexpr (static_cast<bool>(system & kPCM)) {
        auto groupedV0s = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
        for (auto& v0 : groupedV0s) {
          auto ele = v0.template negTrack_as<aod::V0Legs>();
          auto pos = v0.template posTrack_as<aod::V0Legs>();

          auto o2track_ele = ele.template track_as<TracksMC>();
          auto o2track_pos = pos.template track_as<TracksMC>();
          if (!o2track_ele.has_mcParticle() || !o2track_pos.has_mcParticle()) {
            continue; // If no MC particle is found, skip the v0
          }

          for (auto& leg : {pos, ele}) { // be carefull of order {pos, ele}!
            auto o2track = leg.template track_as<TracksMC>();
            auto mctrack = o2track.template mcParticle_as<aod::McParticles>();

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mctrack.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
              // fMCFlags[mctrack.globalIndex()] = mcflags;
              fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
              fCounters[0]++;
            }
            v0legmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());
          } // end of leg loop
        }   // end of v0 loop
      }
      if constexpr (static_cast<bool>(system & kDalitzEE)) {
        // for dalitz ee
        auto emprimarytracks_coll = emprimarytracks.sliceBy(perCollision_emprmtrk, collision.globalIndex());
        for (auto& emprimarytrack : emprimarytracks_coll) {
          auto o2track = o2tracks.iteratorAt(emprimarytrack.trackId());
          if (!o2track.has_mcParticle()) {
            continue; // If no MC particle is found, skip the dilepton
          }
          auto mctrack = o2track.template mcParticle_as<aod::McParticles>();

          // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            // fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }
          emprimarytrackmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());
        } // end of em primary track loop
      }
    } // end of collision loop

    //  Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      auto mctrack = mcTracks.iteratorAt(oldLabel);
      // uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (auto& m : mctrack.mothersIds()) {
          if (m < mcTracks.size()) { // protect against bad mother indices
            if (fNewLabels.find(m) != fNewLabels.end()) {
              mothers.push_back(fNewLabels.find(m)->second);
            }
          } else {
            std::cout << "Mother label (" << m << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
            std::cout << " Check the MC generator" << std::endl;
          }
        }
      }

      // TODO: Check that the daughter slice in the skimmed table works as expected
      //       Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcTracks.size()) { // protect against bad daughter indices
            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
            std::cout << " Check the MC generator" << std::endl;
          }
        }
      }
      int daughterRange[2] = {-1, -1};
      if (daughters.size() > 0) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      emmcparticles(fEventIdx.find(oldLabel)->second, mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
                    mothers, daughterRange,
                    mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
                    mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt());
    } // end loop over labels

    fNewLabels.clear();
    fNewLabelsReversed.clear();
    // fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
    fCounters[0] = 0;
    fCounters[1] = 0;
  } //  end of skimmingMC

  void processMC_PCM(soa::SmallGroups<MyCollisions> const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    skimmingMC<kPCM>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_DalitzEE(soa::SmallGroups<MyCollisions> const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::DalitzEEs const& dielectrons, aod::EMPrimaryTracks const& emprimarytracks)
  {
    const uint8_t sysflag = kPCM | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, nullptr, dielectrons, emprimarytracks);
  }
  void processMC_PHOS(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters)
  {
    skimmingMC<kPHOS>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::SkimEMCClusters const& emcclusters)
  {
    skimmingMC<kEMC>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_DalitzEE(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::DalitzEEs const& dielectrons, aod::EMPrimaryTracks const& emprimarytracks)
  {
    const uint8_t sysflag = kPCM | kPHOS | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, nullptr, dielectrons, emprimarytracks);
  }
  void processMC_PCM_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_EMC_DalitzEE(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::SkimEMCClusters const& emcclusters, aod::DalitzEEs const& dielectrons, aod::EMPrimaryTracks const& emprimarytracks)
  {
    const uint8_t sysflag = kPCM | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, dielectrons, emprimarytracks);
  }
  void processMC_PHOS_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC_DalitzEE(MyCollisions const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0Photons const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters, aod::DalitzEEs const& dielectrons, aod::EMPrimaryTracks const& emprimarytracks)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, dielectrons, emprimarytracks);
  }

  void processDummy(MyCollisions const& collisions) {}

  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM, "create em mc event table for PCM", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_DalitzEE, "create em mc event table for PCM, DalitzEE", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PHOS, "create em mc event table for PHOS", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_EMC, "create em mc event table for EMCal", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_PHOS, "create em mc event table for PCM, PHOS", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_PHOS_DalitzEE, "create em mc event table for PCM, PHOS, DalitzEE", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_EMC, "create em mc event table for PCM, EMCal", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_EMC_DalitzEE, "create em mc event table for PCM, EMCal, DalitzEE", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PHOS_EMC, "create em mc event table for PHOS, EMCal", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_PHOS_EMC, "create em mc event table for PCM, PHOS, EMCal", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processMC_PCM_PHOS_EMC_DalitzEE, "create em mc event table for PCM, PHOS, EMCal, DalitzEE", false);
  PROCESS_SWITCH(createEMReducedMCEvent, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createEMReducedMCEvent>(cfgc, TaskName{"create-emreduced-mc-event"})};
}
