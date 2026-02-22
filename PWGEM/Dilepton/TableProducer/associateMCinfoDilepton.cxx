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

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
// #include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/TableHelper.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <iostream>
#include <map>
#include <random>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct AssociateMCInfoDilepton {
  enum SubSystem {
    kElectron = 0x1,
    kFwdMuon = 0x2,
    // kPCM = 0x4,
  };

  using MyCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::EMEvSels, aod::EMEoIs>;
  using TracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;
  using FwdTracksMC = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
  using MFTTracksMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

  Produces<o2::aod::EMMCEvents> mcevents;
  Produces<o2::aod::EMMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::EMMCGenVectorMesons> emmcgenvms;
  // Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMPrimaryMuonMCLabels> emprimarymuonmclabels;
  Produces<o2::aod::EMMFTMCLabels> emmftmclabels;
  Produces<o2::aod::EMDummyDatas> emdummydata;

  Configurable<int> n_dummy_loop{"n_dummy_loop", 0, "for loop runs over n times"};
  Configurable<float> down_scaling_omega{"down_scaling_omega", 1.1, "down scaling factor to store omega"};
  Configurable<float> down_scaling_phi{"down_scaling_phi", 1.1, "down scaling factor to store phi"};
  Configurable<float> min_eta_gen_primary{"min_eta_gen_primary", -1.5, "min eta to store generated information"};         // smearing is applied at analysis stage. set wider value.
  Configurable<float> max_eta_gen_primary{"max_eta_gen_primary", +1.5, "max eta to store generated information"};         // smearing is applied at analysis stage. set wider value.
  Configurable<float> min_eta_gen_primary_fwd{"min_eta_gen_primary_fwd", -6.0, "min eta to store generated information"}; // smearing is applied at analysis stage. set wider value.
  Configurable<float> max_eta_gen_primary_fwd{"max_eta_gen_primary_fwd", -1.0, "max eta to store generated information"}; // smearing is applied at analysis stage. set wider value.

  HistogramRegistry registry{"EMMCEvent"};
  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isDecayDielectronInAcceptance(TMCParticle const& mctrack, TMCParticles const& mcTracks)
  {
    if (!mctrack.has_daughters()) {
      return false;
    }

    bool is_lepton = false, is_anti_lepton = false;
    for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
      auto daughter = mcTracks.iteratorAt(d);
      if (daughter.pdgCode() == 11 && (min_eta_gen_primary < daughter.eta() && daughter.eta() < max_eta_gen_primary)) {
        is_lepton = true;
      }
      if (daughter.pdgCode() == -11 && (min_eta_gen_primary < daughter.eta() && daughter.eta() < max_eta_gen_primary)) {
        is_anti_lepton = true;
      }
    }
    if (is_lepton && is_anti_lepton) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isDecayDimuonInAcceptance(TMCParticle const& mctrack, TMCParticles const& mcTracks)
  {
    if (!mctrack.has_daughters()) {
      return false;
    }

    bool is_lepton = false, is_anti_lepton = false;
    for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
      auto daughter = mcTracks.iteratorAt(d);
      if (daughter.pdgCode() == 13 && (min_eta_gen_primary_fwd < daughter.eta() && daughter.eta() < max_eta_gen_primary_fwd)) {
        is_lepton = true;
      }
      if (daughter.pdgCode() == -13 && (min_eta_gen_primary_fwd < daughter.eta() && daughter.eta() < max_eta_gen_primary_fwd)) {
        is_anti_lepton = true;
      }
    }
    if (is_lepton && is_anti_lepton) {
      return true;
    } else {
      return false;
    }
  }

  SliceCache cache;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  // Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;

  // apply rapidity cut for electrons
  Partition<aod::McParticles> mcelectrons = nabs(o2::aod::mcparticle::pdgCode) == 11 && min_eta_gen_primary < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen_primary;
  Partition<aod::McParticles> mcmuons = nabs(o2::aod::mcparticle::pdgCode) == 13 && min_eta_gen_primary_fwd < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen_primary_fwd;
  Partition<aod::McParticles> mcvectormesons = o2::aod::mcparticle::pdgCode == 223 || o2::aod::mcparticle::pdgCode == 333;

  template <uint8_t system, typename TTracks, typename TFwdTracks, typename TMFTTracks, typename TPCMs, typename TPCMLegs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const& mcCollisions, aod::McParticles const& mcTracks, TTracks const& o2tracks, TFwdTracks const& o2fwdtracks, TMFTTracks const&, TPCMs const& v0photons, TPCMLegs const&, TEMPrimaryElectrons const& emprimaryelectrons, TEMPrimaryMuons const& emprimarymuons)
  {
    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    // std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter

    // first, run loop over mc collisions to create map between aod::McCollisions and aod::EMMCEvents
    for (const auto& mcCollision : mcCollisions) {
      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        mcevents(mcCollision.globalIndex(), mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.impactParameter(), mcCollision.eventPlaneAngle());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
      }
    } // end of mc collision loop

    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }

      if (!collision.isSelected()) {
        continue;
      }

      if (!collision.isEoI()) { // events with at least 1 lepton for data reduction.
        continue;
      }

      registry.fill(HIST("hEventCounter"), 2);
      auto mcCollision = collision.mcCollision();
      mceventlabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

    } // end of reconstructed collision loop

    for (const auto& mcCollision : mcCollisions) {
      // store MC true information
      auto mcelectrons_per_mccollision = mcelectrons.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto mcmuons_per_mccollision = mcmuons.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto mcvectormesons_per_mccollision = mcvectormesons.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (const auto& mctrack : mcelectrons_per_mccollision) { // store necessary information for denominator of efficiency
        if (!mctrack.isPhysicalPrimary() && !mctrack.producedByGenerator()) {
          continue;
        }
        // auto mcCollision = mcCollisions.iteratorAt(mctrack.mcCollisionId());

        // only for temporary protection, as of 15.July.2024 (by Daiki Sekihata)
        int motherid_tmp = -999; // first mother index tmp
        if (mctrack.has_mothers()) {
          motherid_tmp = mctrack.mothersIds()[0]; // first mother index
        }
        auto mp_tmp = mcTracks.iteratorAt(motherid_tmp);
        int ndau_tmp = mp_tmp.daughtersIds()[1] - mp_tmp.daughtersIds()[0] + 1;
        if (ndau_tmp < 10) {

          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            // fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

          int motherid = -999; // first mother index
          if (mctrack.has_mothers()) {
            motherid = mctrack.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcTracks.iteratorAt(motherid);

              // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
              if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
                fNewLabels[mp.globalIndex()] = fCounters[0];
                fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
                // fMCFlags[mp.globalIndex()] = mcflags;
                fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
                fCounters[0]++;
              }

              if (mp.has_mothers()) {
                motherid = mp.mothersIds()[0]; // first mother index
              } else {
                motherid = -999;
              }
            } else {
              motherid = -999;
            }
          } // end of mother chain loop
        } // end of ndau protection
      } // end of mc electron loop

      for (const auto& mctrack : mcmuons_per_mccollision) { // store necessary information for denominator of efficiency
        if (!mctrack.isPhysicalPrimary() && !mctrack.producedByGenerator()) {
          continue;
        }
        // auto mcCollision = mcCollisions.iteratorAt(mctrack.mcCollisionId());

        // only for temporary protection, as of 15.July.2024 (by Daiki Sekihata)
        int motherid_tmp = -999; // first mother index tmp
        if (mctrack.has_mothers()) {
          motherid_tmp = mctrack.mothersIds()[0]; // first mother index
        }
        auto mp_tmp = mcTracks.iteratorAt(motherid_tmp);
        int ndau_tmp = mp_tmp.daughtersIds()[1] - mp_tmp.daughtersIds()[0] + 1;
        if (ndau_tmp < 10) {

          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            // fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

          int motherid = -999; // first mother index
          if (mctrack.has_mothers()) {
            motherid = mctrack.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcTracks.iteratorAt(motherid);

              // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
              if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
                fNewLabels[mp.globalIndex()] = fCounters[0];
                fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
                // fMCFlags[mp.globalIndex()] = mcflags;
                fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
                fCounters[0]++;
              }

              if (mp.has_mothers()) {
                motherid = mp.mothersIds()[0]; // first mother index
              } else {
                motherid = -999;
              }
            } else {
              motherid = -999;
            }
          } // end of mother chain loop
        } // end of ndau protection
      } // end of mc muon loop

      for (const auto& mctrack : mcvectormesons_per_mccollision) { // store necessary information for denominator of efficiency
        // Be careful!! dilepton rapidity is different from meson rapidity! No acceptance cut here.

        if (!mctrack.isPhysicalPrimary() && !mctrack.producedByGenerator()) {
          continue;
        }

        if (!isDecayDielectronInAcceptance(mctrack, mcTracks) && !isDecayDimuonInAcceptance(mctrack, mcTracks)) { // acceptance cut to decay dileptons
          continue;
        }

        // auto mcCollision = mcCollisions.iteratorAt(mctrack.mcCollisionId());

        int ndau = mctrack.daughtersIds()[1] - mctrack.daughtersIds()[0] + 1;
        if (ndau < 10) {

          if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mctrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
            // fMCFlags[mctrack.globalIndex()] = mcflags;
            fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }

          // store daughter of vector mesons
          if (mctrack.has_daughters()) {
            bool is_lepton_involved = false;
            for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
              // TODO: remove this check as soon as issues with MC production are fixed
              if (d < mcTracks.size()) { // protect against bad daughter indices
                auto daughter = mcTracks.iteratorAt(d);
                if (std::abs(daughter.pdgCode()) == 11 || std::abs(daughter.pdgCode()) == 13) {
                  is_lepton_involved = true;
                  break;
                }
              } else {
                std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
                std::cout << " Check the MC generator" << std::endl;
              }
            }

            if (is_lepton_involved) {
              // LOGF(info, "daughter range in original MC stack pdg = %d | %d - %d , n dau = %d", mctrack.pdgCode(), mctrack.daughtersIds()[0], mctrack.daughtersIds()[1], mctrack.daughtersIds()[1] -mctrack.daughtersIds()[0] +1);
              for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
                // TODO: remove this check as soon as issues with MC production are fixed
                if (d < mcTracks.size()) { // protect against bad daughter indices
                  auto daughter = mcTracks.iteratorAt(d);
                  // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
                  if (!(fNewLabels.find(daughter.globalIndex()) != fNewLabels.end())) {
                    fNewLabels[daughter.globalIndex()] = fCounters[0];
                    fNewLabelsReversed[fCounters[0]] = daughter.globalIndex();
                    // fMCFlags[daughter.globalIndex()] = mcflags;
                    fEventIdx[daughter.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
                    fCounters[0]++;
                  }
                } else {
                  std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
                  std::cout << " Check the MC generator" << std::endl;
                }
              } // end of daughter loop
            }
          }
        } // end of ndau protection
      } // end of generated vector mesons loop

    } // end of mc collision loop

    // if constexpr (static_cast<bool>(system & kPCM)) {
    //   for (const auto& v0 : v0photons) {
    //     auto collision_from_v0 = collisions.iteratorAt(v0.collisionId());
    //     if (!collision_from_v0.has_mcCollision()) {
    //       continue;
    //     }

    //     auto ele = v0.template negTrack_as<TPCMLegs>();
    //     auto pos = v0.template posTrack_as<TPCMLegs>();

    //     auto o2track_ele = o2tracks.iteratorAt(ele.trackId());
    //     auto o2track_pos = o2tracks.iteratorAt(pos.trackId());

    //     if (!o2track_ele.has_mcParticle() || !o2track_pos.has_mcParticle()) {
    //       continue; // If no MC particle is found, skip the v0
    //     }

    //     for (const auto& leg : {pos, ele}) { // be carefull of order {pos, ele}!
    //       auto o2track = o2tracks.iteratorAt(leg.trackId());
    //       auto mctrack = o2track.template mcParticle_as<aod::McParticles>();
    //       // LOGF(info, "mctrack.globalIndex() = %d, mctrack.index() = %d", mctrack.globalIndex(), mctrack.index()); // these are exactly the same.

    //       // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
    //       if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
    //         fNewLabels[mctrack.globalIndex()] = fCounters[0];
    //         fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
    //         // fMCFlags[mctrack.globalIndex()] = mcflags;
    //         fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mctrack.mcCollisionId())->second;
    //         fCounters[0]++;
    //       }
    //       v0legmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());

    //       // Next, store mother-chain of this reconstructed track.
    //       int motherid = -999; // first mother index
    //       if (mctrack.has_mothers()) {
    //         motherid = mctrack.mothersIds()[0]; // first mother index
    //       }
    //       while (motherid > -1) {
    //         if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
    //           auto mp = mcTracks.iteratorAt(motherid);

    //           // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
    //           if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
    //             fNewLabels[mp.globalIndex()] = fCounters[0];
    //             fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
    //             // fMCFlags[mp.globalIndex()] = mcflags;
    //             fEventIdx[mp.globalIndex()] = fEventLabels.find(mp.mcCollisionId())->second;
    //             fCounters[0]++;
    //           }

    //           if (mp.has_mothers()) {
    //             motherid = mp.mothersIds()[0]; // first mother index
    //           } else {
    //             motherid = -999;
    //           }
    //         } else {
    //           motherid = -999;
    //         }
    //       } // end of mother chain loop
    //     } // end of leg loop
    //   } // end of v0 loop
    // }

    if constexpr (static_cast<bool>(system & kElectron)) {
      // auto emprimaryelectrons_coll = emprimaryelectrons.sliceBy(perCollision_el, collision.globalIndex());
      for (const auto& emprimaryelectron : emprimaryelectrons) {
        auto collision_from_el = collisions.iteratorAt(emprimaryelectron.collisionId());
        if (!collision_from_el.has_mcCollision()) {
          continue;
        }

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
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mctrack.mcCollisionId())->second;
          fCounters[0]++;
        }
        emprimaryelectronmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());

        // Next, store mother-chain of this reconstructed track.
        int motherid = -999; // first mother index
        if (mctrack.has_mothers()) {
          motherid = mctrack.mothersIds()[0]; // first mother index
        }
        while (motherid > -1) {
          if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mcTracks.iteratorAt(motherid);

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mp.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
              // fMCFlags[mp.globalIndex()] = mcflags;
              fEventIdx[mp.globalIndex()] = fEventLabels.find(mp.mcCollisionId())->second;
              fCounters[0]++;
            }

            if (mp.has_mothers()) {
              motherid = mp.mothersIds()[0]; // first mother index
            } else {
              motherid = -999;
            }
          } else {
            motherid = -999;
          }
        } // end of mother chain loop

      } // end of em primary electron loop
    }

    if constexpr (static_cast<bool>(system & kFwdMuon)) {
      // auto emprimarymuons_coll = emprimarymuons.sliceBy(perCollision_mu, collision.globalIndex());
      for (const auto& emprimarymuon : emprimarymuons) {
        auto collision_from_mu = collisions.iteratorAt(emprimarymuon.collisionId());
        if (!collision_from_mu.has_mcCollision()) {
          continue;
        }

        auto o2track = o2fwdtracks.iteratorAt(emprimarymuon.fwdtrackId());
        if (!o2track.has_mcParticle()) {
          continue; // If no MC particle is found, skip the dilepton
        }
        auto mctrack = o2track.template mcParticle_as<aod::McParticles>();

        // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
        if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mctrack.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
          // fMCFlags[mctrack.globalIndex()] = mcflags;
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mctrack.mcCollisionId())->second;
          fCounters[0]++;
        }
        emprimarymuonmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());

        // Next, store mother-chain of this reconstructed track.
        int motherid = -999; // first mother index
        if (mctrack.has_mothers()) {
          motherid = mctrack.mothersIds()[0]; // first mother index
        }
        while (motherid > -1) {
          if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mcTracks.iteratorAt(motherid);

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mp.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
              // fMCFlags[mp.globalIndex()] = mcflags;
              fEventIdx[mp.globalIndex()] = fEventLabels.find(mp.mcCollisionId())->second;
              fCounters[0]++;
            }

            if (mp.has_mothers()) {
              motherid = mp.mothersIds()[0]; // first mother index
            } else {
              motherid = -999;
            }
          } else {
            motherid = -999;
          }
        } // end of mother chain loop

        // mc label for tracks registered in MFT in global muons
        if (o2track.matchMFTTrackId() > -1) {
          auto o2mfttrack = o2track.template matchMFTTrack_as<TMFTTracks>();
          if (!o2mfttrack.has_mcParticle()) {
            emmftmclabels(-1, 0);
            break;
          }

          auto mco2mfttrack = o2mfttrack.template mcParticle_as<aod::McParticles>();
          if (!(fNewLabels.find(mco2mfttrack.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mco2mfttrack.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mco2mfttrack.globalIndex();
            // fMCFlags[mco2mfttrack.globalIndex()] = mcflags;
            fEventIdx[mco2mfttrack.globalIndex()] = fEventLabels.find(mco2mfttrack.mcCollisionId())->second;
            fCounters[0]++;
          }
          emmftmclabels(fNewLabels.find(mco2mfttrack.index())->second, o2track.mcMask());

          // Next, store mother-chain of this reconstructed track.
          int motherid = -999; // first mother index
          if (mctrack.has_mothers()) {
            motherid = mctrack.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcTracks.iteratorAt(motherid);

              // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
              if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
                fNewLabels[mp.globalIndex()] = fCounters[0];
                fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
                // fMCFlags[mp.globalIndex()] = mcflags;
                fEventIdx[mp.globalIndex()] = fEventLabels.find(mp.mcCollisionId())->second;
                fCounters[0]++;
              }

              if (mp.has_mothers()) {
                motherid = mp.mothersIds()[0]; // first mother index
              } else {
                motherid = -999;
              }
            } else {
              motherid = -999;
            }
          } // end of mother chain loop
        } else {
          emmftmclabels(-1, 0);
        }
      } // end of em primary muon loop
    }

    //  Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      auto mctrack = mcTracks.iteratorAt(oldLabel);
      // uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (const auto& m : mctrack.mothersIds()) {
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

      // Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        // int ndau = mctrack.daughtersIds()[1] - mctrack.daughtersIds()[0] + 1;
        // LOGF(info, "daughter range in original MC stack pdg = %d | %d - %d , n dau = %d", mctrack.pdgCode(), mctrack.daughtersIds()[0], mctrack.daughtersIds()[1], mctrack.daughtersIds()[1] -mctrack.daughtersIds()[0] +1);
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcTracks.size()) { // protect against bad daughter indices
            // auto dau_tmp = mcTracks.iteratorAt(d);
            // // LOGF(info, "daughter pdg = %d", dau_tmp.pdgCode());
            // if ((mctrack.pdgCode() == 223 || mctrack.pdgCode() == 333) && (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
            //   if (fNewLabels.find(d) == fNewLabels.end() && (abs(dau_tmp.pdgCode()) == 11 || abs(dau_tmp.pdgCode()) == 13)) {
            //     LOGF(info, "daughter lepton is not found mctrack.globalIndex() = %d, mctrack.producedByGenerator() == %d, ndau = %d | dau_tmp.globalIndex() = %d, dau_tmp.pdgCode() = %d, dau_tmp.producedByGenerator() = %d, dau_tmp.pt() = %f, dau_tmp.eta() = %f, dau_tmp.phi() = %f", mctrack.globalIndex(), mctrack.producedByGenerator(), ndau, dau_tmp.globalIndex(), dau_tmp.pdgCode(), dau_tmp.producedByGenerator(), dau_tmp.pt(), dau_tmp.eta(), dau_tmp.phi());
            //   }
            // }

            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
            std::cout << " Check the MC generator" << std::endl;
          }
        }
      }

      emmcparticles(fEventIdx.find(oldLabel)->second, mctrack.pdgCode(), mctrack.flags(), mctrack.statusCode(),
                    mothers, daughters,
                    mctrack.px(), mctrack.py(), mctrack.pz(), mctrack.e(),
                    mctrack.vx(), mctrack.vy(), mctrack.vz());

      mothers.clear();
      mothers.shrink_to_fit();
      daughters.clear();
      daughters.shrink_to_fit();
    } // end loop over labels

    // only for omega, phi mesons
    for (const auto& mcCollision : mcCollisions) {
      auto mcvectormesons_per_mccollision = mcvectormesons.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (const auto& mctrack : mcvectormesons_per_mccollision) { // store necessary information for denominator of efficiency
        if (!mctrack.isPhysicalPrimary() && !mctrack.producedByGenerator()) {
          continue;
        }

        if (mctrack.pdgCode() == 223) {
          if (dist01(engine) < down_scaling_omega) {
            emmcgenvms(fEventLabels[mcCollision.globalIndex()], mctrack.pdgCode(), mctrack.flags(), mctrack.px(), mctrack.py(), mctrack.pz(), mctrack.e(), down_scaling_omega.value);
          }
        } else if (mctrack.pdgCode() == 333) {
          if (dist01(engine) < down_scaling_phi) {
            emmcgenvms(fEventLabels[mcCollision.globalIndex()], mctrack.pdgCode(), mctrack.flags(), mctrack.px(), mctrack.py(), mctrack.pz(), mctrack.e(), down_scaling_phi.value);
          }
        }
      } // end of generated vector meson loop
    } // end of reconstructed collision loop

    fNewLabels.clear();
    fNewLabelsReversed.clear();
    // fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
    fCounters[0] = 0;
    fCounters[1] = 0;
  } //  end of skimmingMC

  void processMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, nullptr, nullptr, nullptr, emprimaryelectrons, nullptr);
  }

  void processMC_FwdMuon(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, FwdTracksMC const& o2fwdtracks, MFTTracksMC const& o2mfttracks, aod::EMPrimaryMuons const& emprimarymuons)
  {
    const uint8_t sysflag = kFwdMuon;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, o2fwdtracks, o2mfttracks, nullptr, nullptr, nullptr, emprimarymuons);
  }

  void processMC_Electron_FwdMuon(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, FwdTracksMC const& o2fwdtracks, MFTTracksMC const& o2mfttracks, aod::EMPrimaryElectrons const& emprimaryelectrons, aod::EMPrimaryMuons const& emprimarymuons)
  {
    const uint8_t sysflag = kElectron | kFwdMuon;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, o2fwdtracks, o2mfttracks, nullptr, nullptr, emprimaryelectrons, emprimarymuons);
  }

  // void processMC_Electron_FwdMuon_PCM(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, FwdTracksMC const& o2fwdtracks, MFTTracksMC const& o2mfttracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons, aod::EMPrimaryMuons const& emprimarymuons)
  // {
  //   const uint8_t sysflag = kPCM | kElectron | kFwdMuon;
  //   skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, o2fwdtracks, o2mfttracks, v0photons, v0legs, emprimaryelectrons, emprimarymuons);
  // }

  // void processMC_Electron_PCM(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons)
  // {
  //   const uint8_t sysflag = kPCM | kElectron;
  //   skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, nullptr, v0photons, v0legs, emprimaryelectrons, nullptr);
  // }

  // void processMC_PCM(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs)
  // {
  //   skimmingMC<kPCM>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, nullptr, v0photons, v0legs, nullptr, nullptr);
  // }

  void processGenDummy(MyCollisionsMC const&)
  {
    for (int i = 0; i < n_dummy_loop; i++) {
      emdummydata(
        0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
        0.f, 0.f, 0.f, 0.f, 0.f,
        0.f);
    }
  }

  void processDummy(MyCollisionsMC const&) {}

  PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_Electron, "create em mc event table for Electron", false);
  PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_FwdMuon, "create em mc event table for Forward Muon", false);
  PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_Electron_FwdMuon, "create em mc event table for Electron, FwdMuon", false);
  // PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_Electron_FwdMuon_PCM, "create em mc event table for PCM, Electron, FwdMuon", false);
  // PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_Electron_PCM, "create em mc event table for PCM, Electron", false);
  // PROCESS_SWITCH(AssociateMCInfoDilepton, processMC_PCM, "create em mc event table for PCM", false);
  PROCESS_SWITCH(AssociateMCInfoDilepton, processGenDummy, "produce dummy data", false);
  PROCESS_SWITCH(AssociateMCInfoDilepton, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AssociateMCInfoDilepton>(cfgc, TaskName{"associate-mc-info-dilepton"})};
}
