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
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::mcutil;

using MyCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels>;
using TracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;
using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCClusterMCLabels>;

struct AssociateMCInfo {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kDalitzEE = 0x8,
    kDalitzMuMu = 0x10,
  };

  Produces<o2::aod::EMMCEvents> mcevents;
  Produces<o2::aod::EMMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMPrimaryMuonMCLabels> emprimarymuonmclabels;
  Produces<o2::aod::EMEMCClusterMCLabels> ememcclustermclabels;

  Configurable<float> max_rxy_gen{"max_rxy_gen", 100, "max rxy to store generated information"};
  Configurable<float> max_Y_gen{"max_Y_gen", 0.9, "max rapidity Y to store generated information"};
  Configurable<float> margin_z_gen{"margin_z_gen", 15.f, "margin for Z of true photon conversion point to store generated information"};

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");
    registry.add<TH2>("PCM/hXY", "hRxy;X (cm);Y (cm)", kTH2F, {{400, -100, +100}, {400, -100, +100}});
    registry.add<TH2>("PCM/hRZ", "hRxy;R (cm);Z (cm)", kTH2F, {{400, -100, +100}, {200, 0, +100}});
    registry.add<TH1>("Generated/hPt_Pi0", "pT distribution of #pi^{0};p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH1>("Generated/hPt_Eta", "pT distribution of #eta;p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH1>("Generated/hPt_ChargedPion", "pT distribution of #pi^{#pm};p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH1>("Generated/hPt_ChargedKaon", "pT distribution of K^{#pm};p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH1>("Generated/hPt_K0S", "pT distribution of K0S;p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
    registry.add<TH1>("Generated/hPt_Lambda", "pT distribution of #Lambda(#bar{#Lambda});p_{T} (GeV/c)", kTH1F, {{200, 0, 20}});
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<MyEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <uint8_t system, typename TTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const&, aod::McParticles const& mcTracks, TTracks const& o2tracks, TPCMs const& v0photons, TPCMLegs const& /*v0legs*/, TPHOSs const& /*phosclusters*/, TEMCs const& emcclusters, TEMPrimaryElectrons const& emprimaryelectrons, TEMPrimaryMuons const& emprimarymuons)
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
        mcevents(mcCollision.globalIndex(), mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.t(), mcCollision.impactParameter());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
      }

      mceventlabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

      // store mc particles
      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (auto& mctrack : groupedMcTracks) { // store necessary information for denominator of efficiency
        if (mctrack.pt() < 1e-3 || abs(mctrack.vz()) > 250 || sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) > max_rxy_gen) {
          continue;
        }
        int pdg = mctrack.pdgCode();
        if (abs(pdg) > 1e+9) {
          continue;
        }

        // fill basic histograms
        if ((mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) && abs(mctrack.y()) < 0.5) {
          switch (abs(pdg)) {
            case 111:
              registry.fill(HIST("Generated/hPt_Pi0"), mctrack.pt());
              break;
            case 211:
              registry.fill(HIST("Generated/hPt_ChargedPion"), mctrack.pt());
              break;
            case 221:
              registry.fill(HIST("Generated/hPt_Eta"), mctrack.pt());
              break;
            case 310:
              registry.fill(HIST("Generated/hPt_K0S"), mctrack.pt());
              break;
            case 321:
              registry.fill(HIST("Generated/hPt_ChargedKaon"), mctrack.pt());
              break;
            case 3122:
              registry.fill(HIST("Generated/hPt_Lambda"), mctrack.pt());
              break;
            default:
              break;
          }
        }

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
          // charmonia
          && (abs(pdg) != 443)    // J/psi
          && (abs(pdg) != 100443) // psi(2S)
          // bottomonia
          && (abs(pdg) != 553)    // Upsilon(1S)
          && (abs(pdg) != 100553) // Upsilon(2S)
          && (abs(pdg) != 200553) // Upsilon(3S)

          // heavy flavor hadrons
          && (std::to_string(abs(pdg))[std::to_string(abs(pdg)).length() - 3] != '4') // charmed mesons
          && (std::to_string(abs(pdg))[std::to_string(abs(pdg)).length() - 3] != '5') // beauty mesons
          && (std::to_string(abs(pdg))[std::to_string(abs(pdg)).length() - 4] != '4') // charmed baryons
          && (std::to_string(abs(pdg))[std::to_string(abs(pdg)).length() - 4] != '5') // beauty baryons

          // strange hadrons
          // && (abs(pdg) != 310)  // K0S
          // && (abs(pdg) != 130)  // K0L
          // && (abs(pdg) != 3122) // Lambda
        ) {
          continue;
        }
        // LOGF(info,"index = %d , mc track pdg = %d , producedByGenerator =  %d , isPhysicalPrimary = %d", mctrack.index(), mctrack.pdgCode(), mctrack.producedByGenerator(), mctrack.isPhysicalPrimary());

        if (!(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) { // neither physical primary nor producedByGenerator
          if (abs(pdg) == 11) {                                                // one more check for secondary electrons. i.e. gamma->ee

            if (sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) < abs(mctrack.vz()) * std::tan(2 * std::atan(std::exp(-max_Y_gen))) - margin_z_gen) {
              continue;
            }

            if (mctrack.has_mothers()) {
              auto mp = mctrack.template mothers_first_as<aod::McParticles>(); // mother particle of electron
              int pdg_mother = mp.pdgCode();
              if (pdg_mother != 22 || !(mp.isPhysicalPrimary() || mp.producedByGenerator())) { // mother of electron is not photon, or not physical primary, or not producedByGenerator
                continue;
              }
            }
          } else { // not physical primary, not producedByGenerator, not electrons
            continue;
          }
        } else if (abs(mctrack.y()) > max_Y_gen) { // physical primary or producedByGenerator, but outside of acceptance.
          continue;
        }

        if (abs(pdg) == 11 && !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) { // only for quick check, only secondary electrons should appear.
          registry.fill(HIST("PCM/hXY"), mctrack.vx(), mctrack.vy());
          registry.fill(HIST("PCM/hRZ"), mctrack.vz(), sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)));
        }

        // these are used as denominator for efficiency. (i.e. generated information)
        if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mctrack.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
          // fMCFlags[mctrack.globalIndex()] = mcflags;
          fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
          fCounters[0]++;
        }

        bool is_used_for_gen = false;
        if ((mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) && (abs(pdg) == 11 || abs(pdg) == 13) && mctrack.has_mothers()) {
          auto mp = mctrack.template mothers_first_as<aod::McParticles>(); // mother particle of electron
          int pdg_mother = abs(mp.pdgCode());

          bool is_from_sm = pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333 || pdg_mother == 443 || pdg_mother == 100443 || pdg_mother == 553 || pdg_mother == 100553 || pdg_mother == 200553;
          bool is_from_c_hadron = std::to_string(pdg_mother)[std::to_string(pdg_mother).length() - 3] == '4' || std::to_string(pdg_mother)[std::to_string(pdg_mother).length() - 4] == '4';
          bool is_from_b_hadron = std::to_string(pdg_mother)[std::to_string(pdg_mother).length() - 3] == '5' || std::to_string(pdg_mother)[std::to_string(pdg_mother).length() - 4] == '5';

          if ((is_from_sm || is_from_c_hadron || is_from_b_hadron)) {
            is_used_for_gen = true;
          }
        }

        if (is_used_for_gen) {
          // Next, store mother-chain for only HF->l, because HF->ll analysis requires correlation between HF hadrons or quarks. PLEASE DON'T do this for other partices.
          int motherid = -999; // first mother index
          if (mctrack.has_mothers()) {
            motherid = mctrack.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcTracks.iteratorAt(motherid);

              if (abs(mp.pdgCode()) < 100) { // don't store quark/gluon informaiton, because data size explodes.
                break;
              }

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
            // LOGF(info, "mctrack.globalIndex() = %d, mctrack.index() = %d", mctrack.globalIndex(), mctrack.index()); // these are exactly the same.

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mctrack.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
              // fMCFlags[mctrack.globalIndex()] = mcflags;
              fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
              fCounters[0]++;
            }
            v0legmclabels(fNewLabels.find(mctrack.index())->second, o2track.mcMask());

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
          }   // end of leg loop
        }     // end of v0 loop
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

        } // end of em primary track loop
      }
      if constexpr (static_cast<bool>(system & kEMC)) {
        // for emc photons
        auto ememcclusters_coll = emcclusters.sliceBy(perCollision_emc, collision.globalIndex());
        for (auto& ememccluster : ememcclusters_coll) {
          auto mcphoton = mcTracks.iteratorAt(ememccluster.emmcparticleIds()[0]);

          // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
          if (!(fNewLabels.find(mcphoton.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mcphoton.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mcphoton.globalIndex();
            fEventIdx[mcphoton.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
            fCounters[0]++;
          }
          ememcclustermclabels(fNewLabels.find(mcphoton.index())->second);

          // Next, store mother-chain of this reconstructed track.
          int motherid = -999; // first mother index
          if (mcphoton.has_mothers()) {
            motherid = mcphoton.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcTracks.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcTracks.iteratorAt(motherid);

              // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
              if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
                fNewLabels[mp.globalIndex()] = fCounters[0];
                fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
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

        } // end of em emc cluster loop
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

      // Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        // LOGF(info, "daughter range in original MC stack pdg = %d | %d - %d , n dau = %d", mctrack.pdgCode(), mctrack.daughtersIds()[0], mctrack.daughtersIds()[1], mctrack.daughtersIds()[1] -mctrack.daughtersIds()[0] +1);
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcTracks.size()) { // protect against bad daughter indices
            // auto dau_tmp = mcTracks.iteratorAt(d);
            // LOGF(info, "daughter pdg = %d", dau_tmp.pdgCode());
            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcTracks.size() << ")" << std::endl;
            std::cout << " Check the MC generator" << std::endl;
          }
        }
      }

      emmcparticles(fEventIdx.find(oldLabel)->second, mctrack.pdgCode(), mctrack.flags(),
                    mothers, daughters,
                    mctrack.px(), mctrack.py(), mctrack.pz(), mctrack.e(),
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
  void processMC_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, nullptr, nullptr, nullptr, emprimaryelectrons, nullptr);
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
  void processMC_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, MyEMCClusters const& emcclusters)
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
  void processMC_PCM_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_EMC_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, nullptr, emcclusters, emprimaryelectrons, nullptr);
  }
  void processMC_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC_DalitzEE(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC | kDalitzEE;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, v0photons, v0legs, phosclusters, emcclusters, emprimaryelectrons, nullptr);
  }

  void processDummy(MyCollisionsMC const&) {}

  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM, "create em mc event table for PCM", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_DalitzEE, "create em mc event table for PCM, DalitzEE", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_DalitzEE, "create em mc event table for DalitzEE", false);
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
