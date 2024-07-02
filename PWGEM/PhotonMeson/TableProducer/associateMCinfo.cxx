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
using FwdTracksMC = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCClusterMCLabels>;

struct AssociateMCInfo {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kElectron = 0x8,
    kFwdMuon = 0x10,
  };

  Produces<o2::aod::EMMCEvents> mcevents;
  Produces<o2::aod::EMMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMPrimaryMuonMCLabels> emprimarymuonmclabels;
  Produces<o2::aod::EMEMCClusterMCLabels> ememcclustermclabels;

  Produces<o2::aod::BinnedGenPts> binned_gen_pt;
  // Produces<o2::aod::BinnedGenPtAccs> binned_gen_pt_acc;

  Configurable<float> max_rxy_gen{"max_rxy_gen", 100, "max rxy to store generated information"};
  Configurable<float> max_eta_gen_primary{"max_eta_gen_primary", 1.2, "max rapidity Y to store generated information"};   // smearing might be applied at analysis stage. set wider value.
  Configurable<float> min_eta_gen_primary_fwd{"min_eta_gen_primary_fwd", -4.5, "min eta to store generated information"}; // smearing might be applied at analysis stage. set wider value.
  Configurable<float> max_eta_gen_primary_fwd{"max_eta_gen_primary_fwd", -2.0, "max eta to store generated information"}; // smearing might be applied at analysis stage. set wider value.
  Configurable<float> max_eta_gen_secondary{"max_eta_gen_secondary", 0.9, "max eta to store generated information"};
  Configurable<float> margin_z_gen{"margin_z_gen", 15.f, "margin for Z of true photon conversion point to store generated information"};

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1I, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");
    registry.add<TH2>("PCM/hXY", "hRxy;X (cm);Y (cm)", kTH2F, {{1000, -100, +100}, {1000, -100, +100}});
    registry.add<TH2>("PCM/hRZ", "hRxy;R (cm);Z (cm)", kTH2F, {{1000, -100, +100}, {1000, 0, +100}});

    // !!Don't change pt,eta,y binning. These binnings have to be consistent with binned data at analysis.!!
    std::vector<double> ptbins;
    for (int i = 0; i < 2; i++) {
      ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.1 GeV/c, every 0.05 GeV/c
    }
    for (int i = 2; i < 52; i++) {
      ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 5 GeV/c, every 0.1 GeV/c
    }
    for (int i = 52; i < 62; i++) {
      ptbins.emplace_back(0.5 * (i - 52) + 5.0); // from 5 to 10 GeV/c, evety 0.5 GeV/c
    }
    for (int i = 62; i < 73; i++) {
      ptbins.emplace_back(1.0 * (i - 62) + 10.0); // from 10 to 20 GeV/c, evety 1 GeV/c
    }
    const AxisSpec axis_pt{ptbins, "p_{T} (GeV/c)"};
    const AxisSpec axis_rapidity{{0.0, +0.8, +0.9}, "rapidity |y|"};

    static constexpr std::string_view parnames[9] = {
      "Gamma", "Pi0", "Eta", "Omega", "Phi",
      "ChargedPion", "ChargedKaon", "K0S", "Lambda"};

    for (int i = 0; i < 9; i++) {
      registry.add<TH2>(Form("Generated/h2PtY_%s", parnames[i].data()), Form("Generated %s", parnames[i].data()), kTH2F, {axis_pt, axis_rapidity}, true);
    }

    // reserve space for generated vectors if that process enabled
    auto hBinFinder = registry.get<TH2>(HIST("Generated/h2PtY_Pi0"));
    LOGF(info, "Binned generated processing enabled. Initialising with %i elements...", hBinFinder->GetNcells());
    genGamma.resize(hBinFinder->GetNcells(), 0);
    genPi0.resize(hBinFinder->GetNcells(), 0);
    genEta.resize(hBinFinder->GetNcells(), 0);
    genOmega.resize(hBinFinder->GetNcells(), 0);
    genPhi.resize(hBinFinder->GetNcells(), 0);
    genChargedPion.resize(hBinFinder->GetNcells(), 0);
    genChargedKaon.resize(hBinFinder->GetNcells(), 0);
    genK0S.resize(hBinFinder->GetNcells(), 0);
    genLambda.resize(hBinFinder->GetNcells(), 0);

    static constexpr std::string_view parnames_acc[8] = {
      "Pi0_acc_gg", "Pi0_acc_eeg",
      "Eta_acc_gg", "Eta_acc_eeg", "Eta_acc_mumug", "Eta_acc_pipig",
      "Omega_acc_ee",
      "Phi_acc_ee"};
    const AxisSpec axis_eta{{0.0, +0.8, +0.9}, "pseudo-rapidity |#eta|"};

    for (int i = 0; i < 8; i++) {
      registry.add<THnSparse>(Form("Generated/hs_%s", parnames_acc[i].data()), Form("Generated %s", parnames_acc[i].data()), kTHnSparseF, {axis_pt, axis_rapidity, axis_eta, axis_eta, axis_eta}, true);
    }
    auto hBinFinder_acc = registry.get<THnSparse>(HIST("Generated/hs_Pi0_acc_gg"));
    int nbins_acc = 1;
    for (int idim = 0; idim < hBinFinder_acc->GetNdimensions(); idim++) {
      nbins_acc *= hBinFinder_acc->GetAxis(idim)->GetNbins(); // without over/underflow
    }

    LOGF(info, "Binned generated processing enabled. Initialising with %i elements... for particles in acceptance", nbins_acc);
    genPi0_acc_gg.resize(nbins_acc, 0);
    genPi0_acc_eeg.resize(nbins_acc, 0);
    genEta_acc_gg.resize(nbins_acc, 0);
    genEta_acc_eeg.resize(nbins_acc, 0);
    genEta_acc_mumug.resize(nbins_acc, 0);
    genEta_acc_pipig.resize(nbins_acc, 0);
    genOmega_acc_ee.resize(nbins_acc, 0);
    genPhi_acc_ee.resize(nbins_acc, 0);
  }

  // template <typename TMCParticle>
  // bool isBeam(TMCParticle const& p)
  // {
  //   if ((abs(p.pdgCode()) == 2212 || abs(p.pdgCode()) > 1e+9) && p.pt() < 1e-4 && p.pz() > 440.f && p.globalIndex() < 2 && !p.has_mothers()) {
  //     return true;
  //   } else {
  //     return false;
  //   }
  // }

  // template <typename TMCParticle>
  // bool isQuarkOrGluon(TMCParticle const& p)
  // {
  //   if ((1 <= abs(p.pdgCode()) && abs(p.pdgCode()) <= 6) || p.pdgCode() == 21) {
  //     return true;
  //   } else {
  //     return false;
  //   }
  // }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::EMPrimaryElectrons> perCollision_el = aod::emprimaryelectron::collisionId;
  Preslice<aod::EMPrimaryMuons> perCollision_mu = aod::emprimarymuon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<MyEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  std::vector<uint16_t> genGamma;         // primary, pt, y
  std::vector<uint16_t> genPi0;           // primary, pt, y
  std::vector<uint16_t> genPi0_acc_gg;    // primary, pt, y
  std::vector<uint16_t> genPi0_acc_eeg;   // primary, pt, y
  std::vector<uint16_t> genEta;           // primary, pt, y
  std::vector<uint16_t> genEta_acc_gg;    // primary, pt, y
  std::vector<uint16_t> genEta_acc_eeg;   // primary, pt, y
  std::vector<uint16_t> genEta_acc_mumug; // primary, pt, y
  std::vector<uint16_t> genEta_acc_pipig; // primary, pt, y
  std::vector<uint16_t> genOmega;         // primary, pt, y
  std::vector<uint16_t> genOmega_acc_ee;  // primary, pt, y
  std::vector<uint16_t> genPhi;           // primary, pt, y
  std::vector<uint16_t> genPhi_acc_ee;    // primary, pt, y
  std::vector<uint16_t> genChargedPion;   // primary, pt, y
  std::vector<uint16_t> genChargedKaon;   // primary, pt, y
  std::vector<uint16_t> genK0S;           // primary, pt, y
  std::vector<uint16_t> genLambda;        // primary, pt, y

  template <uint8_t system, typename TTracks, typename TFwdTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const&, aod::McParticles const& mcTracks, TTracks const& o2tracks, TFwdTracks const& o2fwdtracks, TPCMs const& v0photons, TPCMLegs const& /*v0legs*/, TPHOSs const& /*phosclusters*/, TEMCs const& emcclusters, TEMPrimaryElectrons const& emprimaryelectrons, TEMPrimaryMuons const& emprimarymuons)
  {
    // temporary variables used for the indexing of the skimmed MC stack
    std::map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    // std::map<uint64_t, uint16_t> fMCFlags;
    std::map<uint64_t, int> fEventIdx;
    std::map<uint64_t, int> fEventLabels;
    int fCounters[2] = {0, 0}; //! [0] - particle counter, [1] - event counter
    auto hBinFinder = registry.get<TH2>(HIST("Generated/h2PtY_Gamma"));

    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      // TODO: investigate the collisions without corresponding mcCollision
      if (!collision.has_mcCollision()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2);

      std::fill(genGamma.begin(), genGamma.end(), 0);
      std::fill(genPi0.begin(), genPi0.end(), 0);
      std::fill(genEta.begin(), genEta.end(), 0);
      std::fill(genOmega.begin(), genOmega.end(), 0);
      std::fill(genPhi.begin(), genPhi.end(), 0);
      std::fill(genChargedPion.begin(), genChargedPion.end(), 0);
      std::fill(genChargedKaon.begin(), genChargedKaon.end(), 0);
      std::fill(genK0S.begin(), genK0S.end(), 0);
      std::fill(genLambda.begin(), genLambda.end(), 0);

      std::fill(genPi0_acc_gg.begin(), genPi0_acc_gg.end(), 0);
      std::fill(genPi0_acc_eeg.begin(), genPi0_acc_eeg.end(), 0);
      std::fill(genEta_acc_gg.begin(), genEta_acc_gg.end(), 0);
      std::fill(genEta_acc_eeg.begin(), genEta_acc_eeg.end(), 0);
      std::fill(genEta_acc_mumug.begin(), genEta_acc_mumug.end(), 0);
      std::fill(genEta_acc_pipig.begin(), genEta_acc_pipig.end(), 0);
      std::fill(genOmega_acc_ee.begin(), genOmega_acc_ee.end(), 0);
      std::fill(genPhi_acc_ee.begin(), genPhi_acc_ee.end(), 0);

      auto mcCollision = collision.mcCollision();
      // store mc particles
      auto groupedMcTracks = mcTracks.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (auto& mctrack : groupedMcTracks) { // store necessary information for denominator of efficiency
        if ((mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) && abs(mctrack.y()) < 0.9f && mctrack.pt() < 20.f) {
          auto binNumber = hBinFinder->FindBin(mctrack.pt(), mctrack.y()); // caution: pack
          switch (abs(mctrack.pdgCode())) {
            case 22:
              registry.fill(HIST("Generated/h2PtY_Gamma"), mctrack.pt(), mctrack.y());
              genGamma[binNumber]++;
              break;
            case 111:
              registry.fill(HIST("Generated/h2PtY_Pi0"), mctrack.pt(), mctrack.y());
              genPi0[binNumber]++;
              if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{22, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genPi0_acc_gg[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                registry.fill(HIST("Generated/hs_Pi0_acc_gg"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), 0.0);
              } else if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-11, 11, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genPi0_acc_eeg[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                auto dau3 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 2);
                registry.fill(HIST("Generated/hs_Pi0_acc_eeg"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), abs(dau3.eta()));
              }
              break;
            case 221:
              registry.fill(HIST("Generated/h2PtY_Eta"), mctrack.pt(), mctrack.y());
              genEta[binNumber]++;
              if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{22, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genEta_acc_gg[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                registry.fill(HIST("Generated/hs_Eta_acc_gg"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), 0.0);
              } else if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-11, 11, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genEta_acc_eeg[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                auto dau3 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 2);
                registry.fill(HIST("Generated/hs_Eta_acc_eeg"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), abs(dau3.eta()));
              } else if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-13, 13, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genEta_acc_mumug[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                auto dau3 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 2);
                registry.fill(HIST("Generated/hs_Eta_acc_mumug"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), abs(dau3.eta()));
              } else if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-211, 211, 22}, -0.9, +0.9, 0, 2 * M_PI)) {
                genEta_acc_pipig[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                auto dau3 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 2);
                registry.fill(HIST("Generated/hs_Eta_acc_pipig"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), abs(dau3.eta()));
              }
              break;
            case 223:
              registry.fill(HIST("Generated/h2PtY_Omega"), mctrack.pt(), mctrack.y());
              genOmega[binNumber]++;
              if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-11, 11}, -0.9, +0.9, 0, 2 * M_PI)) {
                genOmega_acc_ee[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                registry.fill(HIST("Generated/hs_Omega_acc_ee"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), 0.f);
              }
              break;
            case 333:
              registry.fill(HIST("Generated/h2PtY_Phi"), mctrack.pt(), mctrack.y());
              genPhi[binNumber]++;
              if (o2::aod::pwgem::mcutil::IsInAcceptanceNonDerived(mctrack, mcTracks, std::vector<int>{-11, 11}, -0.9, +0.9, 0, 2 * M_PI)) {
                genPhi_acc_ee[binNumber]++;
                auto dau1 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 0);
                auto dau2 = mcTracks.iteratorAt(mctrack.daughtersIds()[0] + 1);
                registry.fill(HIST("Generated/hs_Phi_acc_ee"), mctrack.pt(), abs(mctrack.y()), abs(dau1.eta()), abs(dau2.eta()), 0.f);
              }
              break;
            case 211:
              registry.fill(HIST("Generated/h2PtY_ChargedPion"), mctrack.pt(), mctrack.y());
              genChargedPion[binNumber]++;
              break;
            case 321:
              registry.fill(HIST("Generated/h2PtY_ChargedKaon"), mctrack.pt(), mctrack.y());
              genChargedKaon[binNumber]++;
              break;
            case 310:
              registry.fill(HIST("Generated/h2PtY_K0S"), mctrack.pt(), mctrack.y());
              genK0S[binNumber]++;
              break;
            case 3122:
              registry.fill(HIST("Generated/h2PtY_Lambda"), mctrack.pt(), mctrack.y());
              genLambda[binNumber]++;
              break;
            default:
              break;
          }
        }
      } // end of mc track loop

      // make an entry for this MC event only if it was not already added to the table
      if (!(fEventLabels.find(mcCollision.globalIndex()) != fEventLabels.end())) {
        mcevents(mcCollision.globalIndex(), mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(), mcCollision.t(), mcCollision.impactParameter());
        fEventLabels[mcCollision.globalIndex()] = fCounters[1];
        fCounters[1]++;
        binned_gen_pt(genGamma, genPi0, genEta, genOmega, genPhi, genChargedPion, genChargedKaon, genK0S, genLambda);
        // binned_gen_pt_acc(
        //   genPi0_acc_gg, genPi0_acc_eeg,
        //   genEta_acc_gg, genEta_acc_eeg, genEta_acc_mumug, genEta_acc_pipig,
        //   genOmega_acc_ee,
        //   genPhi_acc_ee);
      }

      mceventlabels(fEventLabels.find(mcCollision.globalIndex())->second, collision.mcMask());

      for (auto& mctrack : groupedMcTracks) { // store necessary information for denominator of efficiency
        if (mctrack.pt() < 1e-3 || abs(mctrack.vz()) > 250 || sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) > max_rxy_gen) {
          continue;
        }
        int pdg = mctrack.pdgCode();
        if (abs(pdg) > 1e+9) {
          continue;
        }

        // Note that pi0 from weak decay gives producedByGenerator() = false
        if (
          abs(pdg) != 11      // electron
          && (abs(pdg) != 13) // muon
        ) {
          continue;
        }
        // LOGF(info,"index = %d , mc track pdg = %d , producedByGenerator =  %d , isPhysicalPrimary = %d", mctrack.index(), mctrack.pdgCode(), mctrack.producedByGenerator(), mctrack.isPhysicalPrimary());

        if (abs(pdg) == 11 && !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) { // only for quick check, only secondary electrons should appear.
          registry.fill(HIST("PCM/hXY"), mctrack.vx(), mctrack.vy());
          registry.fill(HIST("PCM/hRZ"), mctrack.vz(), sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)));
        }

        if (mctrack.isPhysicalPrimary() || mctrack.producedByGenerator()) {                                                                         // primary leptons
          if ((abs(mctrack.eta()) < max_eta_gen_primary) || (min_eta_gen_primary_fwd < mctrack.eta() && mctrack.eta() < max_eta_gen_primary_fwd)) { // primary leptons
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
          }
        } else if (abs(pdg) == 11 && mctrack.has_mothers()) { // secondary electrons. i.e. ele/pos from photon conversions.
          int motherid = mctrack.mothersIds()[0];             // first mother index
          auto mp = mcTracks.iteratorAt(motherid);

          if (sqrt(pow(mctrack.vx(), 2) + pow(mctrack.vy(), 2)) < abs(mctrack.vz()) * std::tan(2 * std::atan(std::exp(-max_eta_gen_secondary))) - margin_z_gen) {
            continue;
          }

          if (mp.pdgCode() == 22 && (mp.isPhysicalPrimary() || mp.producedByGenerator()) && abs(mp.eta()) < max_eta_gen_secondary) {
            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!(fNewLabels.find(mctrack.globalIndex()) != fNewLabels.end())) { // store electron information. !!Not photon!!
              fNewLabels[mctrack.globalIndex()] = fCounters[0];
              fNewLabelsReversed[fCounters[0]] = mctrack.globalIndex();
              // fMCFlags[mctrack.globalIndex()] = mcflags;
              fEventIdx[mctrack.globalIndex()] = fEventLabels.find(mcCollision.globalIndex())->second;
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
      if constexpr (static_cast<bool>(system & kElectron)) {
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

        } // end of em primary electron loop
      }
      if constexpr (static_cast<bool>(system & kFwdMuon)) {
        // for dalitz mumu
        auto emprimarymuons_coll = emprimarymuons.sliceBy(perCollision_mu, collision.globalIndex());
        for (auto& emprimarymuon : emprimarymuons_coll) {
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

        } // end of em primary muon loop
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
    skimmingMC<kPCM>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, nullptr, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, nullptr, nullptr, emprimaryelectrons, nullptr);
  }
  void processMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, nullptr, nullptr, nullptr, nullptr, emprimaryelectrons, nullptr);
  }
  void processMC_FwdMuon(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, FwdTracksMC const& o2fwdtracks, aod::EMPrimaryMuons const& emprimarymuons)
  {
    const uint8_t sysflag = kFwdMuon;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, o2fwdtracks, nullptr, nullptr, nullptr, nullptr, nullptr, emprimarymuons);
  }
  void processMC_PCM_Electron_FwdMuon(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, FwdTracksMC const& o2fwdtracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::EMPrimaryElectrons const& emprimaryelectrons, aod::EMPrimaryMuons const& emprimarymuons)
  {
    const uint8_t sysflag = kPCM | kElectron | kFwdMuon;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, o2fwdtracks, v0photons, v0legs, nullptr, nullptr, emprimaryelectrons, emprimarymuons);
  }
  void processMC_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters)
  {
    skimmingMC<kPHOS>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, nullptr, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, MyEMCClusters const& emcclusters)
  {
    skimmingMC<kEMC>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, nullptr, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, phosclusters, nullptr, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, phosclusters, nullptr, emprimaryelectrons, nullptr);
  }
  void processMC_PCM_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, nullptr, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_EMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, MyEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kEMC | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, nullptr, emcclusters, emprimaryelectrons, nullptr);
  }
  void processMC_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, nullptr, nullptr, nullptr, nullptr, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, phosclusters, emcclusters, nullptr, nullptr);
  }
  void processMC_PCM_PHOS_EMC_Electron(MyCollisionsMC const& collisions, aod::BCs const& bcs, aod::McCollisions const& mccollisions, aod::McParticles const& mcTracks, TracksMC const& o2tracks, aod::V0PhotonsKF const& v0photons, aod::V0Legs const& v0legs, aod::PHOSClusters const& phosclusters, MyEMCClusters const& emcclusters, aod::EMPrimaryElectrons const& emprimaryelectrons)
  {
    const uint8_t sysflag = kPCM | kPHOS | kEMC | kElectron;
    skimmingMC<sysflag>(collisions, bcs, mccollisions, mcTracks, o2tracks, nullptr, v0photons, v0legs, phosclusters, emcclusters, emprimaryelectrons, nullptr);
  }

  void processDummy(MyCollisionsMC const&) {}

  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM, "create em mc event table for PCM", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_Electron, "create em mc event table for PCM, Electron", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_Electron, "create em mc event table for Electron", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_FwdMuon, "create em mc event table for Forward Muon", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_Electron_FwdMuon, "create em mc event table for PCM, Electron, FwdMuon", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PHOS, "create em mc event table for PHOS", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_EMC, "create em mc event table for EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS, "create em mc event table for PCM, PHOS", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_Electron, "create em mc event table for PCM, PHOS, Electron", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_EMC, "create em mc event table for PCM, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_EMC_Electron, "create em mc event table for PCM, EMCal, Electron", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PHOS_EMC, "create em mc event table for PHOS, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_EMC, "create em mc event table for PCM, PHOS, EMCal", false);
  PROCESS_SWITCH(AssociateMCInfo, processMC_PCM_PHOS_EMC_Electron, "create em mc event table for PCM, PHOS, EMCal, Electron", false);
  PROCESS_SWITCH(AssociateMCInfo, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AssociateMCInfo>(cfgc, TaskName{"associate-mc-info"})};
}
