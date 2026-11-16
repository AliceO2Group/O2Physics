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

/// \file associateMCinfoPhoton.cxx
/// \brief This code produces EmMc tables were the McParticleIds get reshuffled due to not storing all McParticles
/// \author Daiki Sekihata (daiki.sekihata@cern.ch), Marvin Hemmer (marvin.hemmer@cern.ch), Nicolas Strangmann (nicolas.strangmann@cern.ch)

#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/GammaTablesRedux.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TPDGCode.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;
using namespace o2::constants::physics;

using MyCollisionsMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PMEvSels, aod::EMEoIs>;
using TracksMC = soa::Join<aod::TracksIU, aod::McTrackLabels>;
using FwdTracksMC = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
using MyEMCClusters = soa::Join<aod::MinClusters, aod::EMCClusterMCLabels_001>;

struct Counter {
  int particles{0};
  int events{0};
};

struct AssociateMCInfoPhoton {
  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kElectron = 0x8,
  };

  static constexpr uint16_t kNoMcParticleMask = (1 << 14); // isNoise bit: no genuine MC particle for this leg/track

  Produces<o2::aod::EMMCEvents> mcevents;
  Produces<o2::aod::EMMCEventLabels> mceventlabels;
  Produces<o2::aod::EMMCParticles> emmcparticles;
  Produces<o2::aod::V0LegMCLabels> v0legmclabels;
  Produces<o2::aod::EMPrimaryElectronMCLabels> emprimaryelectronmclabels;
  Produces<o2::aod::EMEMCClusterMCLabels_001> ememcclustermclabels;

  Produces<o2::aod::BinnedGenPts> binnedGenPt;

  Configurable<float> max_eta_gen_secondary{"max_eta_gen_secondary", 0.9, "max eta to store generated information"};
  Configurable<float> margin_z_gen{"margin_z_gen", 15.f, "margin for Z of true photon conversion point to store generated information"};
  Configurable<float> max_rxy_gen{"max_rxy_gen", 100, "max rxy to store generated information"};
  Configurable<bool> requireGammaGammaDecay{"requireGammaGammaDecay", false, "require gamma gamma decay for generated pi0 and eta meson"};
  Configurable<float> cfg_max_eta_photon{"cfg_max_eta_gamma", 0.8, "max eta for photons at PV"};
  Configurable<float> maxPt{"maxPt", 20.f, "max pT for BinnedGenPts table"};
  Configurable<float> maxY{"maxY", 0.9f, "max |rapidity| for BinnedGenPts table"};
  Configurable<bool> doStoreAllDaughters{"doStoreAllDaughters", false, "flag to enable storing of all photon, pi0, eta, eta' and omega daughters. This will increase the dervied data size!"};
  Configurable<bool> useIsSelected{"useIsSelected", true, "flag to enable the usage of isSelected from our EM event selection. Disabling this will increase the derived data size!"};
  Configurable<bool> useIsEOI{"useIsEOI", true, "flag to enable the usage of the EOI table"};

  Configurable<bool> fillBremsChain{"fillBremsChain", false, "MC truth histograms of the bremsstrahlung chain: photon conversion, radiation by a leg, conversion of the radiated photon (BremsChain/)"};
  Configurable<float> bremsChainMaxR{"bremsChainMaxR", 90.f, "max conversion radius (cm) for the bremsstrahlung chain histograms"};
  Configurable<float> bremsChainMinX{"bremsChainMinX", 0.05f, "min energy fraction E(radiated photon) / E(parent photon) for the bremsstrahlung chain histograms"};

  HistogramRegistry registry{"EMMCEvent"};

  void init(o2::framework::InitContext&)
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "has mc collision");

    // !!Don't change pt,eta,y binning. These binnings have to be consistent with binned data at analysis.!!
    std::vector<double> ptbins;
    ptbins.reserve(72);
    for (int i = 0; i < 2; i++) {                // o2-linter: disable=magic-number (just numbers for binning)
      ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.05 GeV/c, every 0.05 GeV/c
    }
    for (int i = 2; i < 51; i++) {              // o2-linter: disable=magic-number (just numbers for binning)
      ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 4.9 GeV/c, every 0.1 GeV/c
    }
    for (int i = 51; i < 61; i++) {              // o2-linter: disable=magic-number (just numbers for binning)
      ptbins.emplace_back(0.5 * (i - 51) + 5.0); // from 5 to 9.5 GeV/c, every 0.5 GeV/c
    }
    for (int i = 61; i < 72; i++) {               // o2-linter: disable=magic-number (just numbers for binning)
      ptbins.emplace_back(1.0 * (i - 61) + 10.0); // from 10 to 20 GeV/c, every 1 GeV/c
    }
    const AxisSpec axisPt{ptbins, "p_{T} (GeV/c)"};
    const AxisSpec axisRapidity{{0.0, +0.8, +0.9}, "rapidity |y|"};

    static constexpr uint NParticleNames = 9;
    static constexpr std::array<std::string_view, NParticleNames> ParticleNames = {
      "Gamma", "Pi0", "Eta", "Omega", "Phi",
      "ChargedPion", "ChargedKaon", "K0S", "Lambda"};

    for (uint i = 0; i < NParticleNames; i++) {
      registry.add<TH2>(Form("Generated/h2PtY_%s", ParticleNames[i].data()), Form("Generated %s", ParticleNames[i].data()), kTH2F, {axisPt, axisRapidity}, true);
      registry.add<TH2>(Form("Generated/h2PtY_Accepted_%s", ParticleNames[i].data()), Form("Generated %s", ParticleNames[i].data()), kTH2F, {axisPt, axisRapidity}, true);
    }

    if (fillBremsChain) {
      const AxisSpec axR{45, 0.f, 90.f, "#it{R}_{conv}(A) (cm)"};
      const AxisSpec axPhi{360, 0.f, 360.f, "#varphi_{conv}(A) (deg)"};
      const AxisSpec axPtA{{0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 10.0}, "#it{p}_{T}(A) (GeV/#it{c})"};
      const AxisSpec axX{{0.05, 0.1, 0.3, 0.5, 1.0}, "#it{E}(B) / #it{E}(A)"};
      registry.add("BremsChain/hGenPhiPt", "generated primary photons, |#eta| < cfg_max_eta_gamma;#varphi of the momentum (deg);#it{p}_{T} (GeV/#it{c})", kTH2F, {{360, 0.f, 360.f}, axPtA}, false);
      registry.add("BremsChain/hConvA", "true conversions of primary photons", kTH3F, {axR, axPhi, axPtA}, false);
      registry.add("BremsChain/hConvARec", "true conversions of primary photons, photon reconstructed", kTH3F, {axR, axPhi, axPtA}, false);
      registry.add("BremsChain/hChain", "photons radiated by a leg of A, cumulative stages", kTHnSparseF, {axR, axPhi, axPtA, axX, {4, -0.5f, 3.5f, "stage: 0 radiated (in the MC stack), 1 B converts, 2 B reconstructed, 3 B and A reconstructed"}}, false);
      registry.add("BremsChain/hGeometry", "radiated photons that convert: where, relative to A", kTHnSparseF, {{18, 0.f, 90.f, "#it{R}_{conv}(A) (cm)"}, {40, -5.f, 95.f, "#it{R}_{conv}(B) #minus #it{R}_{conv}(A) (cm)"}, {40, -20.f, 20.f, "#varphi_{conv}(B) #minus #varphi_{conv}(A) (deg)"}, axX, {3, 0.5f, 3.5f, "stage: 1 B converts, 2 B reconstructed, 3 B and A reconstructed"}}, false);
      registry.add("BremsChain/hEmission", "radiated photons with #it{E}(B) / #it{E}(A) > 0.1: emission point;#it{R}_{conv}(A) (cm);#it{R}_{emission} #minus #it{R}_{conv}(A) (cm)", kTH2F, {axR, {50, 0.f, 100.f}}, false);
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

  template <o2::soa::is_iterator TMCParticle, o2::soa::is_iterator TMCDaughter>
  void selectDaughtersToStore(TMCParticle& particle, int64_t mcParticleSize, TMCDaughter& daughterIter, int eventIdx, std::unordered_map<uint64_t, int>& fNewLabels, std::map<uint64_t, int>& fNewLabelsReversed, std::unordered_map<uint64_t, int>& fEventIdx, Counter& fCounter)
  {
    if (!particle.has_daughters()) {
      return;
    }
    for (int daughterId = particle.daughtersIds()[0]; daughterId <= particle.daughtersIds()[1]; ++daughterId) {
      if (daughterId < mcParticleSize) {
        daughterIter.setCursor(daughterId);
        if (!fNewLabels.contains(daughterIter.globalIndex())) {
          fNewLabels[daughterIter.globalIndex()] = fCounter.particles;
          fNewLabelsReversed[fCounter.particles] = daughterIter.globalIndex();
          fEventIdx[daughterIter.globalIndex()] = eventIdx;
          fCounter.particles++;
        }
      }
    }
  }

  template <soa::is_iterator TCollision>
  bool passesCollisionGate(TCollision const& collisionIter)
  {
    // event selection same as in createEMEventPhoton
    if (!collisionIter.has_mcCollision()) {
      return false;
    }
    if (useIsSelected.value && !collisionIter.isSelected()) {
      return false;
    }
    if (useIsEOI.value && !collisionIter.isEoI()) {
      return false;
    }
    return true;
  }

  template <o2::soa::is_iterator TMCParticle, o2::soa::is_iterator TMCCollision>
  void selectMothersToStore(int motherId, int64_t mcParticleSize, TMCParticle& motherParticle, TMCParticle& daughterIter, TMCCollision const& mcCollisionIter, std::unordered_map<uint64_t, int>& fNewLabels, std::map<uint64_t, int>& fNewLabelsReversed, std::unordered_map<uint64_t, int>& fEventIdx, std::unordered_map<uint64_t, int>& fEventLabels, Counter& fCounter)
  {
    while (motherId > -1) {
      if (motherId < mcParticleSize) {
        motherParticle.setCursor(motherId);

        auto evIt = fEventLabels.find(mcCollisionIter.globalIndex());
        if (evIt == fEventLabels.end()) {
          // Can't attribute this ancestor to a known event; stop walking rather than risk UB.
          LOG(warning) << "selectMothersToStore: mcCollision " << mcCollisionIter.globalIndex()
                       << " not found in fEventLabels, aborting mother-chain walk.";
          break;
        }
        int eventIdx = evIt->second;

        if (!fNewLabels.contains(motherParticle.globalIndex())) {
          fNewLabels[motherParticle.globalIndex()] = fCounter.particles;
          fNewLabelsReversed[fCounter.particles] = motherParticle.globalIndex();
          fEventIdx[motherParticle.globalIndex()] = eventIdx;
          fCounter.particles++;
        }
        // If this ancestor is a photon, pi0, eta, omega or eta prime store all its daughters
        if ((doStoreAllDaughters) && (std::abs(motherParticle.pdgCode()) == PDG_t::kGamma ||
                                      std::abs(motherParticle.pdgCode()) == PDG_t::kPi0 ||
                                      std::abs(motherParticle.pdgCode()) == o2::constants::physics::Pdg::kEta ||
                                      std::abs(motherParticle.pdgCode()) == o2::constants::physics::Pdg::kEtaPrime ||
                                      std::abs(motherParticle.pdgCode()) == o2::constants::physics::Pdg::kOmega)) {
          selectDaughtersToStore(motherParticle, mcParticleSize, daughterIter,
                                 eventIdx, fNewLabels, fNewLabelsReversed,
                                 fEventIdx, fCounter);
        }

        if (motherParticle.has_mothers()) {
          motherId = motherParticle.mothersIds()[0]; // first mother index
        } else {
          motherId = -999;
        }
      } else {
        motherId = -999;
      }
    } // end of mother chain loop
  }

  // conversion point of a photon: production vertex of its e+ e- daughters; false if the photon did not convert into e+ e- (in the MC stack)
  template <typename TMCParticle, typename TMCParticles>
  bool conversionPoint(TMCParticle const& photon, TMCParticles const& mcParticles, float& vx, float& vy, int64_t& idEle, int64_t& idPos)
  {
    idEle = -1;
    idPos = -1;
    if (!photon.has_daughters()) {
      return false;
    }
    for (int id = photon.daughtersIds()[0]; id <= photon.daughtersIds()[1]; ++id) {
      if (id < 0 || id >= mcParticles.size()) {
        continue;
      }
      const auto d = mcParticles.rawIteratorAt(id);
      if (d.pdgCode() == PDG_t::kElectron && idEle < 0) {
        idEle = id;
        vx = d.vx();
        vy = d.vy();
      } else if (d.pdgCode() == PDG_t::kPositron && idPos < 0) {
        idPos = id;
      }
    }
    return idEle >= 0 && idPos >= 0;
  }

  // MC truth of the bremsstrahlung chain for one MC collision (see the histogram definitions in init)
  template <typename TGroupedMCParticles, typename TMCParticles>
  void fillBremsChainTruth(TGroupedMCParticles const& groupedMcParticles, TMCParticles const& mcParticles, std::unordered_set<int64_t> const& recoPhotons)
  {
    static constexpr float RadToDeg = 57.29578f;
    static constexpr float FullTurn = 360.f, HalfTurn = 180.f;
    auto phiDeg = [](float y, float x) {
      const float phi = std::atan2(y, x) * RadToDeg;
      return phi < 0.f ? phi + FullTurn : phi;
    };
    for (const auto& photon : groupedMcParticles) {
      if (photon.pdgCode() != PDG_t::kGamma || !(photon.isPhysicalPrimary() || photon.producedByGenerator()) || std::fabs(photon.eta()) > cfg_max_eta_photon) {
        continue;
      }
      registry.fill(HIST("BremsChain/hGenPhiPt"), phiDeg(photon.py(), photon.px()), photon.pt());
      float vxA = 0.f, vyA = 0.f;
      int64_t idEle = -1, idPos = -1;
      if (!conversionPoint(photon, mcParticles, vxA, vyA, idEle, idPos)) {
        continue;
      }
      const float rA = std::hypot(vxA, vyA);
      if (rA > bremsChainMaxR) {
        continue;
      }
      const float phiA = phiDeg(vyA, vxA);
      const bool recA = recoPhotons.contains(photon.globalIndex());
      registry.fill(HIST("BremsChain/hConvA"), rA, phiA, photon.pt());
      if (recA) {
        registry.fill(HIST("BremsChain/hConvARec"), rA, phiA, photon.pt());
      }
      for (const auto& idLeg : {idEle, idPos}) {
        const auto leg = mcParticles.rawIteratorAt(idLeg);
        if (!leg.has_daughters()) {
          continue;
        }
        for (int id = leg.daughtersIds()[0]; id <= leg.daughtersIds()[1]; ++id) {
          if (id < 0 || id >= mcParticles.size()) {
            continue;
          }
          const auto radiated = mcParticles.rawIteratorAt(id);
          if (radiated.pdgCode() != PDG_t::kGamma) {
            continue;
          }
          const float x = std::min(radiated.e() / std::max(photon.e(), 1.e-6f), 0.9999f);
          if (x < bremsChainMinX) {
            continue;
          }
          registry.fill(HIST("BremsChain/hChain"), rA, phiA, photon.pt(), x, 0.f);
          if (x > 0.1f) { // o2-linter: disable=magic-number (energy fraction of the calculation)
            registry.fill(HIST("BremsChain/hEmission"), rA, std::hypot(radiated.vx(), radiated.vy()) - rA);
          }
          float vxB = 0.f, vyB = 0.f;
          int64_t idEleB = -1, idPosB = -1;
          if (!conversionPoint(radiated, mcParticles, vxB, vyB, idEleB, idPosB)) {
            continue;
          }
          const float rB = std::hypot(vxB, vyB);
          if (rB > bremsChainMaxR) {
            continue;
          }
          float dPhi = phiDeg(vyB, vxB) - phiA;
          if (dPhi > HalfTurn) {
            dPhi -= FullTurn;
          } else if (dPhi < -HalfTurn) {
            dPhi += FullTurn;
          }
          const bool recB = recoPhotons.contains(radiated.globalIndex());
          const int lastStage = recB ? (recA ? 3 : 2) : 1;
          for (int stage = 1; stage <= lastStage; ++stage) {
            registry.fill(HIST("BremsChain/hChain"), rA, phiA, photon.pt(), x, static_cast<float>(stage));
            registry.fill(HIST("BremsChain/hGeometry"), rA, std::clamp(rB - rA, -4.99f, 94.99f), std::clamp(dPhi, -19.99f, 19.99f), x, static_cast<float>(stage));
          }
        }
      }
    }
  }

  template <uint8_t system, typename TTracks, typename TFwdTracks, typename TPCMs, typename TPCMLegs, typename TPHOSs, typename TEMCs, typename TEMPrimaryElectrons>
  void skimmingMC(MyCollisionsMC const& collisions, aod::BCs const&, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, TTracks const& o2tracks, TFwdTracks const&, TPCMs const& v0photons, TPCMLegs const& legs, TPHOSs const&, TEMCs const& emcclusters, TEMPrimaryElectrons const& emprimaryelectrons)
  {
    // temporary variables used for the indexing of the skimmed MC stack
    std::unordered_map<uint64_t, int> fNewLabels;
    std::map<uint64_t, int> fNewLabelsReversed;
    // std::map<uint64_t, uint16_t> fMCFlags;
    std::unordered_map<uint64_t, int> fEventIdx;
    std::unordered_map<uint64_t, int> fEventLabels;
    Counter fCounter{.particles = 0, .events = 0}; //! [0] - particle counter, [1] - event counter
    auto hBinFinder = registry.get<TH2>(HIST("Generated/h2PtY_Gamma"));

    // collision iterator from EMCal cluster
    auto collisionIter = collisions.begin();

    // mc particle iterator for photons
    auto mcPhoton = mcParticles.begin();

    // mc particles iterator for mother
    auto motherParticle = mcParticles.begin();
    auto daughterIter = mcParticles.begin();

    // mc particles iterator for mother and other mc particles
    auto mcParticleIter = mcParticles.begin();

    // mc collision iter
    auto mcCollisionIter = mcCollisions.begin();

    // MC photons reconstructed as a stored V0 photon (both legs point to the same photon), for the bremsstrahlung chain histograms
    std::unordered_set<int64_t> recoPhotons;
    if constexpr (static_cast<bool>(system & kPCM) && !std::is_same_v<TPCMLegs, std::nullptr_t>) {
      if (fillBremsChain) {
        auto legIter = legs.begin();
        auto trackIter = o2tracks.begin();
        for (const auto& v0 : v0photons) {
          std::array<int64_t, 2> mother{-1, -2};
          int i = 0;
          for (const auto& legId : {v0.posTrackId(), v0.negTrackId()}) {
            legIter.setCursor(legId);
            trackIter.setCursor(legIter.trackId());
            if (trackIter.has_mcParticle()) {
              mcParticleIter.setCursor(trackIter.mcParticleId());
              if (mcParticleIter.has_mothers()) {
                mother[i] = mcParticleIter.mothersIds()[0];
              }
            }
            ++i;
          }
          if (mother[0] == mother[1] && mother[0] >= 0) {
            recoPhotons.insert(mother[0]);
          }
        }
      }
    }

    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      if (!passesCollisionGate(collision)) {
        continue;
      }

      registry.fill(HIST("hEventCounter"), 2);

      std::fill(genGamma.begin(), genGamma.end(), 0);
      std::fill(genPi0.begin(), genPi0.end(), 0);
      std::fill(genEta.begin(), genEta.end(), 0);

      mcCollisionIter.setCursor(collision.mcCollisionId());
      // store mc particles
      auto groupedMcParticles = mcParticles.sliceBy(perMcCollision, mcCollisionIter.globalIndex());

      for (const auto& mcParticle : groupedMcParticles) { // store necessary information for denominator of efficiency
        if ((mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator()) && std::fabs(mcParticle.y()) < maxY && mcParticle.pt() < maxPt) {
          auto binNumber = hBinFinder->FindBin(mcParticle.pt(), std::fabs(mcParticle.y())); // caution: pack

          bool isMesonAccepted = false;
          if (mcParticle.has_daughters()) {
            auto lDaughters = mcParticle.daughters_as<aod::McParticles>();

            if (areTwoPhotonDaughtersAccepted(lDaughters, cfg_max_eta_photon)) {
              isMesonAccepted = true;
            }
          }
          switch (std::abs(mcParticle.pdgCode())) {
            case PDG_t::kGamma:
              registry.fill(HIST("Generated/h2PtY_Gamma"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genGamma[binNumber]++;
              break;
            case PDG_t::kPi0:
              if (requireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles)) {
                continue;
              }
              registry.fill(HIST("Generated/h2PtY_Pi0"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genPi0[binNumber]++;
              if (isMesonAccepted) {
                registry.fill(HIST("Generated/h2PtY_Accepted_Pi0"), mcParticle.pt(), std::fabs(mcParticle.y()));
              }
              break;
            case Pdg::kEta:
              if (requireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles)) {
                continue;
              }
              registry.fill(HIST("Generated/h2PtY_Eta"), mcParticle.pt(), std::fabs(mcParticle.y()));
              genEta[binNumber]++;
              if (isMesonAccepted) {
                registry.fill(HIST("Generated/h2PtY_Accepted_Eta"), mcParticle.pt(), std::fabs(mcParticle.y()));
              }
              break;
            default:
              break;
          }
        }
      } // end of mc track loop

      // make an entry for this MC event only if it was not already added to the table
      if (!fEventLabels.contains(mcCollisionIter.globalIndex())) {
        if (fillBremsChain) {
          fillBremsChainTruth(groupedMcParticles, mcParticles, recoPhotons);
        }
        mcevents(mcCollisionIter.globalIndex(), mcCollisionIter.generatorsID(), mcCollisionIter.posX(), mcCollisionIter.posY(), mcCollisionIter.posZ(), mcCollisionIter.impactParameter(), mcCollisionIter.eventPlaneAngle());
        fEventLabels[mcCollisionIter.globalIndex()] = fCounter.events;
        fCounter.events++;
        binnedGenPt(genGamma, genPi0, genEta);
      }

      // LOGF(info, "collision.globalIndex() = %d , mceventlabels.lastIndex() = %d", collision.globalIndex(), mceventlabels.lastIndex());
      mceventlabels(fEventLabels.find(mcCollisionIter.globalIndex())->second, collision.mcMask());

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

        if (std::abs(pdg) == PDG_t::kElectron && mcParticle.has_mothers() && !(mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator())) { // secondary electrons. i.e. ele/pos from photon conversions.
          int motherid = mcParticle.mothersIds()[0];                                                                                                  // first mother index
          motherParticle.setCursor(motherid);

          if (std::sqrt(std::pow(mcParticle.vx(), 2) + std::pow(mcParticle.vy(), 2)) < std::fabs(mcParticle.vz()) * std::tan(2 * std::atan(std::exp(-max_eta_gen_secondary))) - margin_z_gen) {
            continue;
          }

          if (motherParticle.pdgCode() == PDG_t::kGamma && (motherParticle.isPhysicalPrimary() || motherParticle.producedByGenerator()) && std::fabs(motherParticle.eta()) < max_eta_gen_secondary) {
            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!fNewLabels.contains(mcParticle.globalIndex())) { // store electron information. !!Not photon!!
              fNewLabels[mcParticle.globalIndex()] = fCounter.particles;
              fNewLabelsReversed[fCounter.particles] = mcParticle.globalIndex();
              // fMCFlags[mcParticle.globalIndex()] = mcflags;
              fEventIdx[mcParticle.globalIndex()] = fEventLabels.find(mcCollisionIter.globalIndex())->second;
              fCounter.particles++;
            }

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            if (!fNewLabels.contains(motherParticle.globalIndex())) { // store conversion photon
              fNewLabels[motherParticle.globalIndex()] = fCounter.particles;
              fNewLabelsReversed[fCounter.particles] = motherParticle.globalIndex();
              // fMCFlags[motherParticle.globalIndex()] = mcflags;
              fEventIdx[motherParticle.globalIndex()] = fEventLabels.find(mcCollisionIter.globalIndex())->second;
              fCounter.particles++;
            }
          }
        }
      } // end of mc track loop

    } // end of rec. collision loop

    if constexpr (static_cast<bool>(system & kPCM) && !std::is_same_v<TPCMLegs, std::nullptr_t>) {
      // electron and positron iterators from the TPCMLegs table as well as the o2Track table
      auto ele = legs.begin();
      auto pos = legs.begin();

      auto o2TrackEle = o2tracks.begin();
      auto o2TrackPos = o2tracks.begin();
      auto o2TrackIter = o2tracks.begin();
      for (const auto& v0 : v0photons) {
        collisionIter.setCursor(v0.collisionId());

        // check if the collision is valid AND the legs have valid MC particles
        bool validMc = false;
        if (passesCollisionGate(collisionIter)) {
          mcCollisionIter.setCursor(collisionIter.mcCollisionId());

          ele.setCursor(v0.negTrackId());
          pos.setCursor(v0.posTrackId());
          o2TrackEle.setCursor(ele.trackId());
          o2TrackPos.setCursor(pos.trackId());

          validMc = o2TrackEle.has_mcParticle() && o2TrackPos.has_mcParticle();
        }

        for (const auto& leg : {pos, ele}) { // careful: order {pos, ele}, must stay fixed either way
          if (validMc) {
            o2TrackIter.setCursor(leg.trackId());
            mcParticleIter.setCursor(o2TrackIter.mcParticleId());

            auto [iter, isNew] = fNewLabels.try_emplace(mcParticleIter.globalIndex(), fCounter.particles);
            if (isNew) {
              fNewLabelsReversed[fCounter.particles] = mcParticleIter.globalIndex();
              auto evIt = fEventLabels.find(mcCollisionIter.globalIndex());
              if (evIt != fEventLabels.end()) {
                fEventIdx[mcParticleIter.globalIndex()] = evIt->second;
              }
              fCounter.particles++;
            }
            v0legmclabels(iter->second, o2TrackIter.mcMask());

            int motherid = -999;
            if (mcParticleIter.has_mothers()) {
              motherid = mcParticleIter.mothersIds()[0];
            }
            selectMothersToStore(motherid, mcParticles.size(), motherParticle, daughterIter, mcCollisionIter,
                                 fNewLabels, fNewLabelsReversed, fEventIdx, fEventLabels, fCounter);
          } else {
            v0legmclabels(-1, kNoMcParticleMask); // no MC match for this leg, but the row must still exist
          }
        } // end of loop over legs {pos, ele}
      } // end of v0 loop
    }

    if constexpr (static_cast<bool>(system & kElectron)) {
      auto o2TrackIter = o2tracks.begin();
      for (const auto& emprimaryelectron : emprimaryelectrons) {
        collisionIter.setCursor(emprimaryelectron.collisionId());

        int32_t mcLabel = -1;
        uint16_t mcMask = kNoMcParticleMask;

        if (passesCollisionGate(collisionIter)) {
          mcCollisionIter.setCursor(collisionIter.mcCollisionId());
          o2TrackIter.setCursor(emprimaryelectron.trackId());

          if (o2TrackIter.has_mcParticle()) {
            mcParticleIter.setCursor(o2TrackIter.mcParticleId());

            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            auto [iter, isNew] = fNewLabels.try_emplace(mcParticleIter.globalIndex(), fCounter.particles);
            if (isNew) {
              fNewLabelsReversed[fCounter.particles] = mcParticleIter.globalIndex();
              auto evIt = fEventLabels.find(mcCollisionIter.globalIndex());
              if (evIt != fEventLabels.end()) {
                fEventIdx[mcParticleIter.globalIndex()] = evIt->second;
              }
              fCounter.particles++;
            }
            mcLabel = iter->second;
            mcMask = o2TrackIter.mcMask();

            int motherid = -999;
            if (mcParticleIter.has_mothers()) {
              motherid = mcParticleIter.mothersIds()[0];
            }
            selectMothersToStore(motherid, mcParticles.size(), motherParticle, daughterIter, mcCollisionIter,
                                 fNewLabels, fNewLabelsReversed, fEventIdx, fEventLabels, fCounter);
          } // if (o2TrackIter.has_mcParticle())
        }
        emprimaryelectronmclabels(mcLabel, mcMask); // always exactly one row per emprimaryelectron
      } // end of em primary electron loop
    }

    if constexpr (static_cast<bool>(system & kEMC)) { // for emc photons
      for (const auto& emccluster : emcclusters) {
        collisionIter.setCursor(emccluster.collisionId());

        std::vector<int32_t> vEmcMcParticleIds;
        std::vector<float> vAmplitudes;

        if (passesCollisionGate(collisionIter)) {
          mcCollisionIter.setCursor(collisionIter.mcCollisionId());

          if (emccluster.mcParticleIds().size() != emccluster.amplitude().size()) {
            LOG(warning) << "Cluster " << emccluster.globalIndex()
                         << ": mcParticleIds/amplitude size mismatch ("
                         << emccluster.mcParticleIds().size() << " vs " << emccluster.amplitude().size() << ")";
          }
          size_t nCont = std::min(emccluster.mcParticleIds().size(), emccluster.amplitude().size());

          vEmcMcParticleIds.reserve(nCont);
          vAmplitudes.reserve(nCont);

          for (size_t iCont = 0; iCont < nCont; iCont++) {
            mcPhoton.setCursor(emccluster.mcParticleIds()[iCont]);
            // if the MC truth particle corresponding to this reconstructed track which is not already written, add it to the skimmed MC stack
            auto [iter, isNew] = fNewLabels.try_emplace(mcPhoton.globalIndex(), fCounter.particles);
            if (isNew) {
              fNewLabelsReversed[fCounter.particles] = mcPhoton.globalIndex();
              auto evIt = fEventLabels.find(mcCollisionIter.globalIndex());
              if (evIt != fEventLabels.end()) { // keep your earlier guard here too
                fEventIdx[mcPhoton.globalIndex()] = evIt->second;
              }
              fCounter.particles++;
            }
            vEmcMcParticleIds.emplace_back(iter->second);
            vAmplitudes.emplace_back(emccluster.amplitude()[iCont]);

            int motherid = -999;
            if (mcPhoton.has_mothers()) {
              motherid = mcPhoton.mothersIds()[0];
            }
            selectMothersToStore(motherid, mcParticles.size(), motherParticle, daughterIter, mcCollisionIter,
                                 fNewLabels, fNewLabelsReversed, fEventIdx, fEventLabels, fCounter);
          }
        } // end of loop over mc particles of the current emc cluster
        // whether or not the collision was okay, exactly one row is always written here so later the tables are ensured to have same size!
        ememcclustermclabels(vEmcMcParticleIds, vAmplitudes);
      }
    }

    //  Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fNewLabelsReversed) {
      mcParticleIter.setCursor(oldLabel);
      // uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mcParticleIter.has_mothers()) {
        for (const auto& m : mcParticleIter.mothersIds()) {
          if (m < mcParticles.size()) { // protect against bad mother indices
            auto iter = fNewLabels.find(m);
            if (iter != fNewLabels.end()) {
              mothers.push_back(iter->second);
            }
          } else {
            LOG(info) << "Mother label (" << m << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(info) << " Check the MC generator";
          }
        }
      }

      // Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mcParticleIter.has_daughters()) {
        // LOGF(info, "daughter range in original MC stack pdg = %d | %d - %d , n dau = %d", mcParticleIter.pdgCode(), mcParticleIter.daughtersIds()[0], mcParticleIter.daughtersIds()[1], mcParticleIter.daughtersIds()[1] -mcParticleIter.daughtersIds()[0] +1);
        for (int d = mcParticleIter.daughtersIds()[0]; d <= mcParticleIter.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcParticles.size()) { // protect against bad daughter indices
            // auto dau_tmp = mcParticles.iteratorAt(d);
            // LOGF(info, "daughter pdg = %d", dau_tmp.pdgCode());
            auto iter = fNewLabels.find(d);
            if (iter != fNewLabels.end()) {
              daughters.push_back(iter->second);
            }
          } else {
            LOG(error) << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(error) << " Check the MC generator";
          }
        }
      }

      auto evIt = fEventIdx.find(oldLabel);
      if (evIt == fEventIdx.end()) {
        LOG(warning) << "emmcparticles: no event index for particle label " << oldLabel << ", skipping.";
        continue;
      }
      emmcparticles(evIt->second, mcParticleIter.pdgCode(), mcParticleIter.flags(), mcParticleIter.statusCode(),
                    mothers, daughters,
                    mcParticleIter.px(), mcParticleIter.py(), mcParticleIter.pz(), mcParticleIter.e(),
                    mcParticleIter.vx(), mcParticleIter.vy(), mcParticleIter.vz());
    } // end loop over labels

    fNewLabels.clear();
    fNewLabelsReversed.clear();
    // fMCFlags.clear();
    fEventIdx.clear();
    fEventLabels.clear();
    fCounter.particles = 0;
    fCounter.events = 0;
  } // end of skimmingMC

  template <o2::soa::is_table Daughters>
  inline bool areTwoPhotonDaughtersAccepted(const Daughters& lDaughters, float maxEta)
  {
    if (lDaughters.size() != 2) {
      return false;
    }

    auto it0 = lDaughters.begin();
    auto it1 = it0;
    ++it1;

    return (std::fabs(it0.eta()) < maxEta &&
            std::fabs(it1.eta()) < maxEta);
  }

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

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<AssociateMCInfoPhoton>(context, TaskName{"associate-mc-info-photon"})};
}
