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
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels>;
using TracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;

struct AssociateMCInfo {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kDalitzEE = 0x8,
    kDalitzMuMu = 0x10,
  };

  Produces<o2::aod::EMReducedMCEvents> mcevents;
  Produces<o2::aod::EMReducedMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMPrimaryMuonMCLabels> emprimarymuonmclabels;

  Configurable<float> max_rxy_gen{"max_rxy_gen", 100, "max rxy to store generated information"};

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <uint8_t system, typename TTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const&, aod::McParticles const& mcTracks, TTracks const& o2tracks, TPCMs const& v0photons, TPCMLegs const& v0legs, TPHOSs const& phosclusters, TEMCs const& emcclusters, TEMPrimaryElectrons const& emprimaryelectrons, TEMPrimaryMuons const& emprimarymuons)
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

      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2);

      auto mcCollision = collision.mcCollision();

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
        if (mctrack.pt() < 1e-3 || abs(mctrack.y()) > 1.5 || abs(mctrack.vz()) > 250 || sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) > max_rxy_gen) {
          continue;
        }
        int pdg = mctrack.pdgCode();

        // Note that pi0 from weak decay gives producedByGenerator() = false
        if (
          abs(pdg) != 11      // electron
          && (abs(pdg) != 13) // muon
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

          auto o2track_ele = o2tracks.iteratorAt(pos.trackId());
          auto o2track_pos = o2tracks.iteratorAt(ele.trackId());

          if (!o2track_ele.has_mcParticle() || !o2track_pos.has_mcParticle()) {
            continue; // If no MC particle is found, skip the v0
          }

          for (auto& leg : {pos, ele}) { // be carefull of order {pos, ele}!
            // auto o2track = leg.template track_as<TracksMC>();
            auto o2track = o2tracks.iteratorAt(leg.trackId());
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
        auto emprimaryelectrons_coll = emprimaryelectrons.sliceBy(perCollision_el, collision.globalIndex());
        for (auto& emprimaryelectron : emprimaryelectrons_coll) {
          auto o2track = o2tracks.iteratorAt(emprimaryelectron.trackId());
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
          emprimaryelectronmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());
        } // end of em primary track loop
      }
      if constexpr (static_cast<bool>(system & kDalitzMuMu)) {
        // for dalitz mumu
        auto emprimarymuons_coll = emprimarymuons.sliceBy(perCollision_mu, collision.globalIndex());
        for (auto& emprimarymuon : emprimarymuons_coll) {
          auto o2track = o2tracks.iteratorAt(emprimarymuon.trackId());
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
          emprimarymuonmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());
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

  void processMC_PCM(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs)
  {
    skimmingMC<kPCM>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, nullptr, emprimaryelectrons, nullptr);
  }
  void processMC_PCM_DalitzEE_DalitzMuMu(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons, aod::EMPrimaryMuons const& emprimarymuons)
  {
    const uint8_t sysflag = kPCM | kDalitzEE | kDalitzMuMu;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, nullptr, emprimaryelectrons, emprimarymuons);
  }

  void processMC_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters)
  {
    skimmingMC<kPHOS>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::SkimEMCClusters const& emcclusters)
  {
    skimmingMC<kEMC>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, nullptr, emprimaryelectrons, nullptr);
  }
  void processMC_PCM_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_EMC_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::SkimEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, emprimaryelectrons, nullptr);
  }
  void processMC_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, emprimaryelectrons, nullptr);
  }

  void processDummy(MyCollisionsMC const& collisions) {}

  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM, "create em mc event table for PCM", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_DalitzEE, "create em mc event table for PCM, DalitzEE", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_DalitzEE_DalitzMuMu, "create em mc event table for PCM, DalitzEE, DalitzMuMu", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PHOS, "create em mc event table for PHOS", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_EMC, "create em mc event table for EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS, "create em mc event table for PCM, PHOS", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_DalitzEE, "create em mc event table for PCM, PHOS, DalitzEE", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_EMC, "create em mc event table for PCM, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_EMC_DalitzEE, "create em mc event table for PCM, EMCal, DalitzEE", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PHOS_EMC, "create em mc event table for PHOS, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_EMC, "create em mc event table for PCM, PHOS, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_EMC_DalitzEE, "create em mc event table for PCM, PHOS, EMCal, DalitzEE", false);
  PROCESS_SWITCH(AssociateMCInfo, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AssociateMCInfo>(cfgc, TaskName{"associate-mc-info"})};
}
