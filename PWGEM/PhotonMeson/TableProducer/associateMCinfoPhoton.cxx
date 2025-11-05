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
/// \file associateMCinfoPhoton.cxx
///
/// \brief This code produces reduced events for photon analyses
///
/// \author Daiki Sekihata (daiki.sekihata@cern.ch)
///

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;

using MyCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::EMEvSels>;
using TracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;
using FwdTracksMC = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCClusterMCLabels>;

struct AssociateMCInfoPhoton {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kElectron = 0x8,
  };

  Produces<o2::aod::EMMCEvents> mcevents;
  Produces<o2::aod::EMMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMEMCClusterMCLabels> ememcclustermclabels;

  Produces<o2::aod::BinnedGenPts> binnedGenPt;

  Configurable<float> max_eta_gen_secondary{"max_eta_gen_secondary", 0.9, "max eta to store generated information"};
  Configurable<float> margin_z_gen{"margin_z_gen", 15.f, "margin for Z of true photon conversion point to store generated information"};
  Configurable<float> max_rxy_gen{"max_rxy_gen", 100, "max rxy to store generated information"};
  Configurable<bool> requireGammaGammaDecay{"requireGammaGammaDecay", false, "require gamma gamma decay for generated pi0 and eta meson"};

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");

    // !!Don't change pt,eta,y binning. These binnings have to be consistent with binned data at analysis.!!
    std::vector<double> ptbins;
    for (int i = 0; i < 2; i++) {
      ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.05 GeV/c, every 0.05 GeV/c
    }
    for (int i = 2; i < 51; i++) {
      ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 4.9 GeV/c, every 0.1 GeV/c
    }
    for (int i = 51; i < 61; i++) {
      ptbins.emplace_back(0.5 * (i - 51) + 5.0); // from 5 to 9.5 GeV/c, every 0.5 GeV/c
    }
    for (int i = 61; i < 72; i++) {
      ptbins.emplace_back(1.0 * (i - 61) + 10.0); // from 10 to 20 GeV/c, every 1 GeV/c
    }
    const AxisSpec axisPt{ptbins, "p_{T} (GeV/c)"};
    const AxisSpec axisRapidity{{0.0, +0.8, +0.9}, "rapidity |y|"};

    static constexpr std::string_view ParticleNames[9] = {
      "Gamma", "Pi0", "Eta", "Omega", "Phi",
      "ChargedPion", "ChargedKaon", "K0S", "Lambda"};

    for (int i = 0; i < 9; i++) {
      registry.add<TH2>(Form("Generated/h2PtY_%s", ParticleNames[i].data()), Form("Generated %s", ParticleNames[i].data()), kTH2F, {axisPt, axisRapidity}, true);
    }

    // reserve space for generated vectors if that process enabled
    auto hBinFinder = registry.get<TH2>(HIST("Generated/h2PtY_Pi0"));
    LOGF(info, "Binned generated processing enabled. Initialising with %i elements...", hBinFinder->GetNcells());
    genGamma.resize(hBinFinder->GetNcells(), 0);
    genPi0.resize(hBinFinder->GetNcells(), 0);
    genEta.resize(hBinFinder->GetNcells(), 0);
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0PhotonsKF> perCollisionPCM = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectronsFromDalitz> perCollisionEl = aod::emprimaryelectron::collisionId;
  Preslice<aod::PHOSClusters> perCollisionPHOS = aod::skimmedcluster::collisionId;
  Preslice<MyEMCClusters> perCollisionEMC = aod::skimmedcluster::collisionId;

  std::vector<uint16_t> genGamma; // primary, pt, y
  std::vector<uint16_t> genPi0;   // primary, pt, y
  std::vector<uint16_t> genEta;   // primary, pt, y

  template <uint8_t system, typename TTracks, typename TFwdTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TEMPrimaryElectrons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const&, aod::McParticles const& mcParticles, TTracks const& o2tracks, TFwdTracks const&, TPCMs const& v0photons, TPCMLegs const&, TPHOSs const&, TEMCs const& emcclusters, TEMPrimaryElectrons const& emprimaryelectrons)
  {
    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    // std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter
    auto hBinFinder = registry.get<TH2>(HIST("Generated/h2PtY_Gamma"));

    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }

      if (!collision.isSelected()) {
        continue;
      }

      registry.fill(HIST("hEventCounter"), 2);

      std::fill(genGamma.begin(), genGamma.end(), 0);
      std::fill(genPi0.begin(), genPi0.end(), 0);
      std::fill(genEta.begin(), genEta.end(), 0);

      auto mcCollision = collision.mcCollision();
      // store mc particles
      auto groupedMcParticles = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (const auto& mcParticle : groupedMcParticles) { // store necessary information for denominator of efficiency
        if ((mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator()) && std::fabs(mcParticle.y()) < 0.9f && mcParticle.pt() < 20.f) {
          auto binNumber = hBinFinder->FindBin(mcParticle.pt(), std::fabs(mcParticle.y())); // caution: pack
          switch (std::abs(mcParticle.pdgCode())) {
            case 22:
              registry.fill(HIST("Generated/h2PtY_Gamma"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genGamma[binNumber]++;
              break;
            case 111:
              if (requireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles))
                continue;
              registry.fill(HIST("Generated/h2PtY_Pi0"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genPi0[binNumber]++;
              break;
            case 221:
              if (requireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles))
                continue;
              registry.fill(HIST("Generated/h2PtY_Eta"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genEta[binNumber]++;
              break;
            default:
              break;
          }
        }
      } // end of mc track loop

      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        mcevents(mcCollision.globalIndex(), mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.impactParameter(), mcCollision.eventPlaneAngle());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
        binnedGenPt(genGamma, genPi0, genEta);
      }

      // LOGF(info, "collision.globalIndex() = %d , mceventlabels.lastIndex() = %d", collision.globalIndex(), mceventlabels.lastIndex());
      mceventlabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

      for (const auto& mcParticle : groupedMcParticles) { // store necessary information for denominator of efficiency
        if (mcParticle.pt() < 1e-3 || std::fabs(mcParticle.vz()) > 250 || std::sqrt(std::pow(mcParticle.vx(), 2) + std::pow(mcParticle.vy(), 2)) > max_rxy_gen) {
          continue;
        }
        int pdg = mcParticle.pdgCode();
        if (std::abs(pdg) > 1e+9) {
          continue;
        }

        // Note that pi0 from weak decay gives producedByGenerator() = false
        // LOGF(info,"index = %d , mc track pdg = %d , producedByGenerator =  %d , isPhysicalPrimary = %d", mcParticle.index(), mcParticle.pdgCode(), mcParticle.producedByGenerator(), mcParticle.isPhysicalPrimary());

        if (std::abs(pdg) == 11 && mcParticle.has_mothers() && !(mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator())) { // secondary electrons. i.e. ele/pos from photon conversions.
          int motherid = mcParticle.mothersIds()[0];                                                                                    // first mother index
          auto mp = mcParticles.iteratorAt(motherid);

          if (std::sqrt(std::pow(mcParticle.vx(), 2) + std::pow(mcParticle.vy(), 2)) < std::fabs(mcParticle.vz()) * std::tan(2 * std::atan(std::exp(-max_eta_gen_secondary))) - margin_z_gen) {
            continue;
          }

          if (mp.pdgCode() == 22 && (mp.isPhysicalPrimary() || mp.producedByGenerator()) && std::fabs(mp.eta()) < max_eta_gen_secondary) {
            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mcParticle.globalIndex()) != fNewLabels.end())) { // store electron information. !!Not photon!!
              fNewLabels[mcParticle.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mcParticle.globalIndex();
              // fMCFlags[mcParticle.globalIndex()] = mcflags;
              fEventIdx[mcParticle.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
              fCounters[0]++;
            }

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) { // store conversion photon
              fNewLabels[mp.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
              // fMCFlags[mp.globalIndex()] = mcflags;
              fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
              fCounters[0]++;
            }
          }
        }
      } // end of mc track loop

    } // end of rec. collision loop

    if constexpr (static_cast<bool>(system & kPCM)) {
      for (const auto& v0 : v0photons) {
        auto collisionFromV0 = collisions.iteratorAt(v0.collisionId());
        if (!collisionFromV0.has_mcCollision()) {
          continue;
        }
        auto mcCollisionFromV0 = collisionFromV0.mcCollision();

        auto ele = v0.template negTrack_as<aod::V0Legs>();
        auto pos = v0.template posTrack_as<aod::V0Legs>();

        auto o2TrackEle = o2tracks.iteratorAt(pos.trackId());
        auto o2TrackPos = o2tracks.iteratorAt(ele.trackId());

        if (!o2TrackEle.has_mcParticle() || !o2TrackPos.has_mcParticle()) {
          continue; // If no MC particle is found, skip the v0
        }

        for (const auto& leg : {pos, ele}) { // be carefull of order {pos, ele}!
          auto o2track = o2tracks.iteratorAt(leg.trackId());
          auto mcParticle = o2track.template mcParticle_as<aod::McParticles>();
          // LOGF(info, "mcParticle.globalIndex() = %d, mcParticle.index() = %d", mcParticle.globalIndex(), mcParticle.index()); // these are exactly the same.

          // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
          if (!(fNewLabels.find(mcParticle.globalIndex()) != fNewLabels.end())) {
            fNewLabels[mcParticle.globalIndex()] = fCounters[0];
            fNewLabelsReversed[fCounters[0]] = mcParticle.globalIndex();
            // fMCFlags[mcParticle.globalIndex()] = mcflags;
            fEventIdx[mcParticle.globalIndex()] = fEventLabels.find(mcCollisionFromV0.globalIndex())->second;
            fCounters[0]++;
          }
          v0legmclabels(fNewLabels.find(mcParticle.index())->second, o2track.mcMask());

          // Next, store mother-chain of this reconstructed track.
          int motherid = -999; // first mother index
          if (mcParticle.has_mothers()) {
            motherid = mcParticle.mothersIds()[0]; // first mother index
          }
          while (motherid > -1) {
            if (motherid < mcParticles.size()) { // protect against bad mother indices. why is this needed?
              auto mp = mcParticles.iteratorAt(motherid);

              // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
              if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
                fNewLabels[mp.globalIndex()] = fCounters[0];
                fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
                // fMCFlags[mp.globalIndex()] = mcflags;
                fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollisionFromV0.globalIndex())->second;
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
        } // end of leg loop
      } // end of v0 loop
    }

    if constexpr (static_cast<bool>(system & kElectron)) {
      // auto emprimaryelectrons_coll = emprimaryelectrons.sliceBy(perCollisionEl, collision.globalIndex());
      for (const auto& emprimaryelectron : emprimaryelectrons) {
        auto collisionFromEl = collisions.iteratorAt(emprimaryelectron.collisionId());
        if (!collisionFromEl.has_mcCollision()) {
          continue;
        }
        auto mcCollisionFromEl = collisionFromEl.mcCollision();

        auto o2track = o2tracks.iteratorAt(emprimaryelectron.trackId());
        if (!o2track.has_mcParticle()) {
          continue; // If no MC particle is found, skip the dilepton
        }
        auto mcParticle = o2track.template mcParticle_as<aod::McParticles>();

        // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
        if (!(fNewLabels.find(mcParticle.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mcParticle.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mcParticle.globalIndex();
          // fMCFlags[mcParticle.globalIndex()] = mcflags;
          fEventIdx[mcParticle.globalIndex()] = fEventLabels.find(mcCollisionFromEl.globalIndex())->second;
          fCounters[0]++;
        }
        emprimaryelectronmclabels(fNewLabels.find(mcParticle.index())->second, o2track.mcMask());

        // Next, store mother-chain of this reconstructed track.
        int motherid = -999; // first mother index
        if (mcParticle.has_mothers()) {
          motherid = mcParticle.mothersIds()[0]; // first mother index
        }
        while (motherid > -1) {
          if (motherid < mcParticles.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mcParticles.iteratorAt(motherid);

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mp.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
              // fMCFlags[mp.globalIndex()] = mcflags;
              fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollisionFromEl.globalIndex())->second;
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

    if constexpr (static_cast<bool>(system & kEMC)) { // for emc photons
      // auto ememcclusters_coll = emcclusters.sliceBy(perCollisionEMC, collision.globalIndex());
      for (const auto& emccluster : emcclusters) {
        auto collisionFromEMC = collisions.iteratorAt(emccluster.collisionId());
        if (!collisionFromEMC.has_mcCollision()) {
          continue;
        }
        auto mcCollisionFromEMC = collisionFromEMC.mcCollision();

        auto mcphoton = mcParticles.iteratorAt(emccluster.emmcparticleIds()[0]);

        // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
        if (!(fNewLabels.find(mcphoton.globalIndex()) != fNewLabels.end())) {
          fNewLabels[mcphoton.globalIndex()] = fCounters[0];
          fNewLabelsReversed[fCounters[0]] = mcphoton.globalIndex();
          fEventIdx[mcphoton.globalIndex()] = fEventLabels.find(mcCollisionFromEMC.globalIndex())->second;
          fCounters[0]++;
        }
        ememcclustermclabels(fNewLabels.find(mcphoton.index())->second);

        // Next, store mother-chain of this reconstructed track.
        int motherid = -999; // first mother index
        if (mcphoton.has_mothers()) {
          motherid = mcphoton.mothersIds()[0]; // first mother index
        }
        while (motherid > -1) {
          if (motherid < mcParticles.size()) { // protect against bad mother indices. why is this needed?
            auto mp = mcParticles.iteratorAt(motherid);

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mp.globalIndex()) != fNewLabels.end())) {
              fNewLabels[mp.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mp.globalIndex();
              fEventIdx[mp.globalIndex()] = fEventLabels.find(mcCollisionFromEMC.globalIndex())->second;
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

    //  Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      auto mcParticle = mcParticles.iteratorAt(oldLabel);
      // uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mcParticle.has_mothers()) {
        for (const auto& m : mcParticle.mothersIds()) {
          if (m < mcParticles.size()) { // protect against bad mother indices
            if (fNewLabels.find(m) != fNewLabels.end()) {
              mothers.push_back(fNewLabels.find(m)->second);
            }
          } else {
            LOG(info) << "Mother label (" << m << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(info) << " Check the MC generator";
          }
        }
      }

      // Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mcParticle.has_daughters()) {
        // LOGF(info, "daughter range in original MC stack pdg = %d | %d - %d , n dau = %d", mcParticle.pdgCode(), mcParticle.daughtersIds()[0], mcParticle.daughtersIds()[1], mcParticle.daughtersIds()[1] -mcParticle.daughtersIds()[0] +1);
        for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcParticles.size()) { // protect against bad daughter indices
            // auto dau_tmp = mcParticles.iteratorAt(d);
            // LOGF(info, "daughter pdg = %d", dau_tmp.pdgCode());
            if (fNewLabels.find(d) != fNewLabels.end()) {
              daughters.push_back(fNewLabels.find(d)->second);
            }
          } else {
            LOG(info) << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(info) << " Check the MC generator";
          }
        }
      }

      emmcparticles(fEventIdx.find(oldLabel)->second, mcParticle.pdgCode(), mcParticle.flags(), mcParticle.statusCode(),
                    mothers, daughters,
                    mcParticle.px(), mcParticle.py(), mcParticle.pz(), mcParticle.e(),
                    mcParticle.vx(), mcParticle.vy(), mcParticle.vz());
    } // end loop over labels

    fNewLabels.clear();
    fNewLabelsReversed.clear();
    // fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
    fCounters[0] = 0;
    fCounters[1] = 0;
  } // end of skimmingMC

  void processMC_PCM(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs)
  {
    skimmingMC<kPCM>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, nullptr, nullptr, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM, "create em mc event table for PCM", false);

  void processMC_PCM_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectronsFromDalitz const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, nullptr, nullptr, emprimaryelectrons);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_Electron, "create em mc event table for PCM, Electron", false);

  void processMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::EMPrimaryElectronsFromDalitz const& emprimaryelectrons)
  {
    const uint8_t sysflag = kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, nullptr, nullptr, nullptr, nullptr, emprimaryelectrons);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_Electron, "create em mc event table for Electron", false);

  void processMC_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, aod::PHOSClusters const& phosclusters)
  {
    skimmingMC<kPHOS>(collisions, bcs, mccollisions, mcParticles, nullptr, nullptr, nullptr, nullptr, phosclusters, nullptr, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PHOS, "create em mc event table for PHOS", false);

  void processMC_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, MyEMCClusters const& emcclusters)
  {
    skimmingMC<kEMC>(collisions, bcs, mccollisions, mcParticles, nullptr, nullptr, nullptr, nullptr, nullptr, emcclusters, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_EMC, "create em mc event table for EMCal", false);

  void processMC_PCM_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, phosclusters, nullptr, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_PHOS, "create em mc event table for PCM, PHOS", false);

  void processMC_PCM_PHOS_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::EMPrimaryElectronsFromDalitz const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, phosclusters, nullptr, emprimaryelectrons);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_PHOS_Electron, "create em mc event table for PCM, PHOS, Electron", false);

  void processMC_PCM_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, nullptr, emcclusters, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_EMC, "create em mc event table for PCM, EMCal", false);

  void processMC_PCM_EMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters, aod::EMPrimaryElectronsFromDalitz const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kEMC | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, nullptr, emcclusters, emprimaryelectrons);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_EMC_Electron, "create em mc event table for PCM, EMCal, Electron", false);

  void processMC_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, nullptr, nullptr, nullptr, nullptr, phosclusters, emcclusters, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PHOS_EMC, "create em mc event table for PHOS, EMCal", false);

  void processMC_PCM_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, phosclusters, emcclusters, nullptr);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_PHOS_EMC, "create em mc event table for PCM, PHOS, EMCal", false);

  void processMC_PCM_PHOS_EMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters, aod::EMPrimaryElectronsFromDalitz const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcParticles, o2tracks, nullptr, v0photons, v0legs, phosclusters, emcclusters, emprimaryelectrons);
  }
  PROCESS_SWITCH(AssociateMCInfoPhoton, processMC_PCM_PHOS_EMC_Electron, "create em mc event table for PCM, PHOS, EMCal, Electron", false);

  void processDummy(MyCollisionsMC const&) {}
  PROCESS_SWITCH(AssociateMCInfoPhoton, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AssociateMCInfoPhoton>(cfgc, TaskName{"associate-mc-info-photon"})};
}
