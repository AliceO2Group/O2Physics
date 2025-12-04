// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file gammaJetTreeProducer.cxx
/// \brief Task to produce a tree for gamma-jet analysis, including photons (and information of isolation) and charged jets
/// \author Florian Jonas <florian.jonas@cern.ch>, UC Berkeley/LBNL
/// \since 02.08.2024

// C++ system headers first
#include <TPDGCode.h>

#include <string>
#include <unordered_map>
#include <vector>

// Framework and other headers after
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/GammaJetAnalysisTree.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/filterTables.h"

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "TVector2.h"

// \struct GammaJetTreeProducer
/// \brief Task to produce a tree for gamma-jet analysis, including photons (and information of isolation) and charged and full jets
/// \author Florian Jonas <florian.jonas@cern.ch>, UC Berkeley/LBNL
/// \since 02.08.2024
///
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using emcClusters = o2::soa::Join<o2::aod::JClusters, o2::aod::JClusterTracks>;
using emcMCClusters = o2::soa::Join<o2::aod::JMcClusterLbs, o2::aod::JClusters, o2::aod::JClusterTracks>;

#include "Framework/runDataProcessing.h"

struct GammaJetTreeProducer {
  // analysis tree
  // charged jets
  // photon candidates
  Produces<aod::GjChargedJets> chargedJetsTable;   // detector level jets
  Produces<aod::GjEvents> eventsTable;             // rec events
  Produces<aod::GjGammas> gammasTable;             // detector level clusters
  Produces<aod::GjMCEvents> mcEventsTable;         // mc collisions information
  Produces<aod::GjMCParticles> mcParticlesTable;   // gen level particles (photons and pi0)
  Produces<aod::GjGammaMCInfos> gammaMCInfosTable; // detector level clusters MC information
  Produces<aod::GjChJetMCInfos> chJetMCInfosTable; // detector level charged jets MC information
  Produces<aod::GjMCJets> mcJetsTable;             // gen level jets

  HistogramRegistry mHistograms{"GammaJetTreeProducerHisto"};

  Service<o2::framework::O2DatabasePDG> pdg;

  // ---------------
  // Configureables
  // ---------------

  // event cuts
  Configurable<double> mVertexCut{"vertexCut", 10.0, "apply z-vertex cut with value in cm"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<std::string>
    trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackMinPt{"trackMinPt", 0.15, "minimum track pT cut"};
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> isoR{"isoR", 0.4, "isolation cone radius"};
  Configurable<float> perpConeJetR{"perpConeJetR", 0.4, "perpendicular cone radius used to calculate perp cone rho for jet"};
  Configurable<float> trackMatchingEoverP{"trackMatchingEoverP", 2.0, "closest track is required to have E/p < value"};
  Configurable<float> minClusterETrigger{"minClusterETrigger", 0.0, "minimum cluster energy to trigger"};
  Configurable<float> minMCGenPt{"minMCGenPt", 0.0, "minimum pt of mc gen particles to store"};

  int mRunNumber = 0;
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  const int kMaxRecursionDepth = 100;

  std::unordered_map<int32_t, int32_t> collisionMapping;
  std::unordered_map<int32_t, int> mcJetIndexMapping; // maps the global index to the index in the mc jets table (per event). This is because later we want to later construct all trees on a per event level, and we need to know at what position in the table per event this is stored
  std::vector<int> triggerMaskBits;
  std::vector<int32_t> mcCollisionsMultiRecCollisions; // used for MC. global index of MC collisions that have multiple matched rec collisions

  // kd tree for tracks and mc particles (used for fast isolation calculation)
  std::vector<float> trackEta;
  std::vector<float> trackPhi;
  std::vector<float> trackPt;
  std::vector<float> mcParticleEta;
  std::vector<float> mcParticlePhi;
  std::vector<float> mcParticlePt;
  TKDTree<int, float>* trackTree = nullptr;
  TKDTree<int, float>* mcParticleTree = nullptr;

  void init(InitContext const&)
  {
    using o2HistType = HistType;
    using o2Axis = AxisSpec;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // create histograms
    LOG(info) << "Creating histograms";

    const o2Axis ptAxis{100, 0, 200, "p_{T} (GeV/c)"};
    const o2Axis ptRecAxis{100, 0, 200, "p_{T}^{rec} (GeV/c)"};
    const o2Axis ptGenAxis{100, 0, 200, "p_{T}^{gen} (GeV/c)"};
    const o2Axis energyAxis{100, 0, 100, "E (GeV)"};
    const o2Axis m02Axis{100, 0, 3, "m02"};
    const o2Axis etaAxis{100, -1, 1, "#eta"};
    const o2Axis phiAxis{100, 0, o2::constants::math::TwoPI, "#phi"};
    const o2Axis dRAxis{100, 0, 1, "dR"};
    const o2Axis occupancyAxis{300, 0, 30000, "occupancy"};
    const o2Axis nCollisionsAxis{10, -0.5, 9.5, "nCollisions"};
    mHistograms.add("clusterE", "Energy of cluster", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("trackPt", "pT of track", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetPt", "pT of charged jet", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("chjetPtEtaPhi", "pT of charged jet", o2HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
    mHistograms.add("chjetpt_vs_constpt", "pT of charged jet vs pT of constituents", o2HistType::kTH2F, {ptRecAxis, ptGenAxis});

    // track QA THnSparse
    mHistograms.add("trackPtEtaPhi", "Track QA", o2HistType::kTHnSparseF, {ptAxis, etaAxis, phiAxis});
    mHistograms.add("trackPtEtaOccupancy", "Track QA vs occupancy", o2HistType::kTHnSparseF, {ptAxis, etaAxis, occupancyAxis});

    // QA for MC collisions to rec collision matching
    // number of reconstructed and matched collisions for each MC collision vs mc gen photon energy
    mHistograms.add("numberRecCollisionsVsPhotonPt", "Number of rec collisions vs photon energy", o2HistType::kTH2F, {nCollisionsAxis, energyAxis});

    // Cluster MC histograms
    mHistograms.add("clusterMC_E_All", "Cluster energy for photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_Photon", "Cluster energy for photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_PromptPhoton", "Cluster energy for prompt photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_DirectPromptPhoton", "Cluster energy for direct prompt photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_FragmentationPhoton", "Cluster energy for fragmentation photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_DecayPhoton", "Cluster energy for decay photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_DecayPhotonPi0", "Cluster energy for decay photons from pi0", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_DecayPhotonEta", "Cluster energy for decay photons from eta", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_MergedPi0", "Cluster energy for merged pi0s", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_MergedEta", "Cluster energy for merged etas", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_E_ConvertedPhoton", "Cluster energy for converted photons", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("clusterMC_m02_Photon", "M02 for photons", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_PromptPhoton", "M02 for prompt photons", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_DirectPromptPhoton", "M02 for direct prompt photons", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_FragmentationPhoton", "M02 for fragmentation photons", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_DecayPhoton", "M02 for decay photons", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_DecayPhotonPi0", "M02 for decay photons from pi0", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_DecayPhotonEta", "M02 for decay photons from eta", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_MergedPi0", "M02 for merged pi0s", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_MergedEta", "M02 for merged etas", o2HistType::kTH1F, {m02Axis});
    mHistograms.add("clusterMC_m02_ConvertedPhoton", "M02 for converted photons", o2HistType::kTH1F, {m02Axis});

    // MC Gen trigger particle histograms
    mHistograms.add("mcGenTrigger_Eta", "eta of mc gen trigger particle", o2HistType::kTH1F, {etaAxis});
    mHistograms.add("mcGenTrigger_Phi", "phi of mc gen trigger particle", o2HistType::kTH1F, {phiAxis});
    mHistograms.add("mcGenTrigger_Pt", "pT of mc gen trigger particle", o2HistType::kTH1F, {ptAxis});
    mHistograms.add("mcGenTrigger_E", "E of mc gen trigger particle", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_PromptPhoton", "E of mc gen trigger prompt photon", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_DirectPromptPhoton", "E of mc gen trigger direct prompt photon", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_FragmentationPhoton", "E of mc gen trigger fragmentation photon", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_DecayPhoton", "E of mc gen trigger decay photon", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_DecayPhotonPi0", "E of mc gen trigger decay photon from pi0", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_DecayPhotonEta", "E of mc gen trigger decay photon from eta", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_DecayPhotonOther", "E of mc gen trigger decay photon from other", o2HistType::kTH1F, {energyAxis});
    mHistograms.add("mcGenTrigger_E_Pi0", "E of mc gen trigger pi0", o2HistType::kTH1F, {energyAxis});

    // MC Particle level jet histograms
    mHistograms.add("mcpJetPt", "pT of mc particle level jet", o2HistType::kTH1F, {ptAxis});

    // MC Detector level jet matching jet histograms
    mHistograms.add("mcdJetPtVsTrueJetPtMatchingGeo", "pT rec (x-axis) of detector level jets vs pT true (y-axis) of mc particle level jet (geo matching)", o2HistType::kTH2F, {ptRecAxis, ptGenAxis});
    mHistograms.add("mcdJetPtVsTrueJetPtMatchingPt", "pT rec (x-axis) of detector level jets vs pT true (y-axis) of mc particle level jet (pt matching)", o2HistType::kTH2F, {ptRecAxis, ptGenAxis});

    // Event QA histogram
    const int nEventBins = 8;
    const TString eventLabels[nEventBins] = {"All", "AfterVertexCut", "AfterCollisionSelection", "AfterTriggerSelection", "AfterEMCALSelection", "AfterClusterESelection", "Has MC collision", "is not MB Gap"};
    mHistograms.add("eventQA", "Event QA", o2HistType::kTH1F, {{nEventBins, -0.5, 7.5}});
    for (int iBin = 0; iBin < nEventBins; iBin++) {
      mHistograms.get<TH1>(HIST("eventQA"))->GetXaxis()->SetBinLabel(iBin + 1, eventLabels[iBin]);
    }

    // MC collisions QA histograms)
    const int nRecCollisionBins = 4;
    const TString recCollisionLabels[nRecCollisionBins] = {"All", "1 Rec collision", "More than 1 rec collisions", "No rec collisions"};
    mHistograms.add("mcCollisionsWithRecCollisions", "MC collisions with rec collisions", o2HistType::kTH1F, {{nRecCollisionBins, -0.5, 3.5}});
    for (int iBin = 0; iBin < nRecCollisionBins; iBin++) {
      mHistograms.get<TH1>(HIST("mcCollisionsWithRecCollisions"))->GetXaxis()->SetBinLabel(iBin + 1, recCollisionLabels[iBin]);
    }
  }

  // ---------------------
  // Helper functions
  // ---------------------

  /// \brief Builds the kd tree for the tracks or mc particles (used for fast isolation calculation)
  /// \param objects The objects to build the kd tree for (tracks or mc particles)
  template <typename T>
  void buildKdTree(const T& objects)
  {
    trackEta.clear();
    trackPhi.clear();
    trackPt.clear();
    mcParticleEta.clear();
    mcParticlePhi.clear();
    mcParticlePt.clear();

    // if the track type is aod::JetTracks, we need to build the kd tree for the tracks
    if constexpr (std::is_same_v<typename std::decay_t<T>, aod::JetTracks>) {
      for (const auto& track : objects) {
        if (!isTrackSelected(track)) {
          continue;
        }
        trackEta.push_back(track.eta());
        trackPhi.push_back(track.phi());
        trackPt.push_back(track.pt());
      }
      if (trackEta.size() > 0) {
        delete trackTree;
        trackTree = new TKDTree<int, float>(trackEta.size(), 2, 1);
        trackTree->SetData(0, trackEta.data());
        trackTree->SetData(1, trackPhi.data());
        trackTree->Build();
      }
    }
    // if the track type is aod::JetParticles, we need to build the kd tree for the mc particles
    if constexpr (std::is_same_v<typename std::decay_t<T>, aod::JetParticles>) {
      for (const auto& particle : objects) {
        if (!particle.isPhysicalPrimary()) {
          continue;
        }
        if (!isCharged(particle)) {
          continue;
        }
        if (particle.pt() < trackMinPt) {
          continue;
        }
        mcParticleEta.push_back(particle.eta());
        mcParticlePhi.push_back(particle.phi());
        mcParticlePt.push_back(particle.pt());
      }
      if (mcParticleEta.size() > 0) {
        delete mcParticleTree;
        mcParticleTree = new TKDTree<int, float>(mcParticleEta.size(), 2, 1);
        mcParticleTree->SetData(0, mcParticleEta.data());
        mcParticleTree->SetData(1, mcParticlePhi.data());
        mcParticleTree->Build();
      }
    }
  }
  /// \brief Checks if a track passes the selection criteria
  /// \param track The track to be checked
  /// \return true if track passes all selection criteria, false otherwise
  bool isTrackSelected(const auto& track)
  {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
      return false;
    }
    if (track.pt() < trackMinPt) {
      return false;
    }

    return true;
  }

  /// \brief Gets the stored collision index from the collision mapping
  /// \param collision The collision to look up
  /// \return The stored collision index, or -1 if not found
  int getStoredColIndex(const auto& collision)
  {
    int32_t storedColIndex = -1;
    if (auto foundCol = collisionMapping.find(collision.globalIndex()); foundCol != collisionMapping.end()) {
      storedColIndex = foundCol->second;
    }
    return storedColIndex;
  }

  /// \brief Checks if an event passes all selection criteria
  /// \param collision The collision to check
  /// \param clusters The EMCAL clusters in the event
  /// \return true if event passes all selection criteria, false otherwise
  bool isEventAccepted(const auto& collision, const auto& clusters)
  {
    mHistograms.fill(HIST("eventQA"), 0);

    if (collision.posZ() > mVertexCut) {
      return false;
    }
    mHistograms.fill(HIST("eventQA"), 1);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return false;
    }
    mHistograms.fill(HIST("eventQA"), 2);
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return false;
    }
    mHistograms.fill(HIST("eventQA"), 3);
    if (!jetderiveddatautilities::eventEMCAL(collision)) {
      return false;
    }
    mHistograms.fill(HIST("eventQA"), 4);

    // Check if event contains a cluster with energy > minClusterETrigger
    for (const auto& cluster : clusters) {
      if (cluster.energy() > minClusterETrigger) {
        mHistograms.fill(HIST("eventQA"), 5);
        return true;
      }
    }
    return false;
  }

  /// \brief Checks if a particle is charged
  /// \param particle The MC particle to check
  /// \return true if particle has non-zero charge, false otherwise
  bool isCharged(const auto& particle)
  {
    return std::abs(pdg->GetParticle(particle.pdgCode())->Charge()) >= 1.;
  }

  /// \brief Calculates the charged particle isolation in a cone of given size using a pre-built kd tree
  /// \param particle The particle to calculate the isolation for
  /// \param radius The cone radius
  /// \param mcGenIso Whether to use the mc gen particle tree (if false, use the track tree)
  /// \return The charged particle isolation
  template <typename T>
  double ch_iso_in_cone(const T& particle, float radius = 0.4, bool mcGenIso = false)
  {
    double iso = 0;
    float point[2] = {particle.eta(), particle.phi()};
    std::vector<int> indices;

    if (!mcGenIso) {
      if (trackTree) {
        trackTree->FindInRange(point, radius, indices);
        for (const auto& index : indices) {
          iso += trackPt[index];
        }
      } else {
        LOG(error) << "Track tree not found";
        return 0;
      }
    } else {
      if (mcParticleTree) {
        mcParticleTree->FindInRange(point, radius, indices);
        for (const auto& index : indices) {
          iso += mcParticlePt[index];
        }
      } else {
        LOG(error) << "MC particle tree not found";
        return 0;
      }
    }
    return iso;
  }

  /// \brief Calculates the charged particle density in perpendicular cones
  /// \param object The reference object (cluster or jet)
  /// \param tracks The tracks to check
  /// \param radius The cone radius for density calculation
  /// \return The average charged particle density in the perpendicular cones
  template <typename T>
  double ch_perp_cone_rho(const T& object, float radius = 0.4, bool mcGenIso = false)
  {
    double ptSumLeft = 0;
    double ptSumRight = 0;

    double cPhi = TVector2::Phi_0_2pi(object.phi());

    // rotate cone left by 90 degrees
    float cPhiLeft = cPhi - o2::constants::math::PIHalf;
    float cPhiRight = cPhi + o2::constants::math::PIHalf;

    float pointLeft[2] = {object.eta(), cPhiLeft};
    float pointRight[2] = {object.eta(), cPhiRight};

    std::vector<int> indicesLeft;
    std::vector<int> indicesRight;

    if (!mcGenIso) {
      if (trackTree) {
        trackTree->FindInRange(pointLeft, radius, indicesLeft);
        trackTree->FindInRange(pointRight, radius, indicesRight);
      } else {
        LOG(error) << "Track tree not found";
        return 0;
      }

      for (const auto& index : indicesLeft) {
        ptSumLeft += trackPt[index];
      }
      for (const auto& index : indicesRight) {
        ptSumRight += trackPt[index];
      }
    } else {
      if (mcParticleTree) {
        mcParticleTree->FindInRange(pointLeft, radius, indicesLeft);
        mcParticleTree->FindInRange(pointRight, radius, indicesRight);
      } else {
        LOG(error) << "MC particle tree not found";
        return 0;
      }
      for (const auto& index : indicesLeft) {
        ptSumLeft += mcParticlePt[index];
      }
      for (const auto& index : indicesRight) {
        ptSumRight += mcParticlePt[index];
      }
    }

    float rho = (ptSumLeft + ptSumRight) / (o2::constants::math::TwoPI * radius * radius);
    return rho;
  }

  /// \brief Fills track QA histograms for a given collision
  /// \param collision The collision containing the tracks
  /// \param tracks The tracks to analyze
  void runTrackQA(const auto& collision, aod::JetTracks const& tracks)
  {
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      mHistograms.fill(HIST("trackPt"), track.pt());
      mHistograms.fill(HIST("trackPtEtaPhi"), track.pt(), track.eta(), track.phi());
      mHistograms.fill(HIST("trackPtEtaOccupancy"), track.pt(), track.eta(), collision.trackOccupancyInTimeRange());
    }
  }

  /// \brief Finds the top-most copy of a particle in the decay chain (following carbon copies)
  /// \param particle The particle to start from
  /// \return The top-most copy of the particle
  template <typename T>
  T iTopCopy(const T& particle) const
  {
    int iUp = particle.globalIndex();
    T currentParticle = particle;
    int pdgCode = particle.pdgCode();
    auto mothers = particle.template mothers_as<aod::JMcParticles>();
    while (iUp > 0 && mothers.size() == 1 && mothers[0].globalIndex() > 0 && mothers[0].pdgCode() == pdgCode) {
      iUp = mothers[0].globalIndex();
      currentParticle = mothers[0];
      mothers = currentParticle.template mothers_as<aod::JMcParticles>();
    }
    return currentParticle;
  }

  /// \brief Checks if a particle is a prompt photon
  /// \param particle The MC particle to check
  /// \return true if particle is a prompt photon, false otherwise
  bool isPromptPhoton(const auto& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) < 90) {
      return true;
    }
    return false;
  }
  /// \brief Checks if a particle is a direct prompt photon
  /// \param particle The particle to check
  /// \return true if particle is a direct prompt photon, false otherwise
  bool isDirectPromptPhoton(const auto& particle)
  {
    // check if particle isa prompt photon
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) < 90) {
      // find the top carbon copy
      auto topCopy = iTopCopy(particle);
      if (topCopy.pdgCode() == PDG_t::kGamma && std::abs(topCopy.getGenStatusCode()) < 40) { // < 40 is particle directly produced in hard scattering
        return true;
      }
    }
    return false;
  }
  /// \brief Checks if a particle is a fragmentation photon
  /// \param particle The particle to check
  /// \return true if particle is a fragmentation photon, false otherwise
  bool isFragmentationPhoton(const auto& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) < 90) {
      // find the top carbon copy
      auto topCopy = iTopCopy(particle);
      if (topCopy.pdgCode() == PDG_t::kGamma && std::abs(topCopy.getGenStatusCode()) >= 40) { // frag photon
        return true;
      }
    }
    return false;
  }
  /// \brief Checks if a particle is a decay photon
  /// \param particle The particle to check
  /// \return true if particle is a decay photon, false otherwise
  bool isDecayPhoton(const auto& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) >= 90) {
      return true;
    }
    return false;
  }
  /// \brief Checks if a particle is a decay photon from pi0
  /// \param particle The particle to check
  /// \return true if particle is a decay photon from pi0, false otherwise
  template <typename T>
  bool isDecayPhotonPi0(const T& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) >= 90) {
      // check if it has mothers that are pi0s
      const auto& mothers = particle.template mothers_as<aod::JMcParticles>();
      for (const auto& mother : mothers) {
        if (mother.pdgCode() == PDG_t::kPi0) {
          return true;
        }
      }
    }
    return false;
  }
  /// \brief Checks if a particle is a decay photon from eta
  /// \param particle The particle to check
  /// \return true if particle is a decay photon from eta, false otherwise
  template <typename T>
  bool isDecayPhotonEta(const T& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) >= 90) {
      // check if it has mothers that are etas
      const auto& mothers = particle.template mothers_as<aod::JMcParticles>();
      for (const auto& mother : mothers) {
        if (mother.pdgCode() == 221) {
          return true;
        }
      }
    }
    return false;
  }
  /// \brief Checks if a particle is a decay photon from other sources
  /// \param particle The particle to check
  /// \return true if particle is a decay photon from other sources, false otherwise
  template <typename T>
  bool isDecayPhotonOther(const T& particle)
  {
    if (particle.pdgCode() == PDG_t::kGamma && particle.isPhysicalPrimary() && std::abs(particle.getGenStatusCode()) >= 90) {
      // check if you find a pi0 mother or a eta mother
      const auto& mothers = particle.template mothers_as<aod::JMcParticles>();
      for (const auto& mother : mothers) {
        if (mother.pdgCode() == PDG_t::kPi0 || mother.pdgCode() == 221) {
          return false;
        }
      }
      return true;
    }
    return false;
  }
  /// \brief Checks if a particle is a pi0
  /// \param particle The particle to check
  /// \return true if particle is a pi0, false otherwise
  bool isPi0(const auto& particle)
  {
    if (particle.pdgCode() == PDG_t::kPi0) {
      return true;
    }
    return false;
  }

  /// \brief Gets the  bitmap for a MC particle that indicated what type of particle it is
  /// \param particle The particle to check
  /// \return A bitmap indicating the particle's origin
  uint16_t getMCParticleOrigin(const auto& particle)
  {
    uint16_t origin = 0;
    if (isPromptPhoton(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kPromptPhoton));
    }
    if (isDirectPromptPhoton(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDirectPromptPhoton));
    }
    if (isFragmentationPhoton(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kFragmentationPhoton));
    }
    if (isDecayPhoton(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhoton));
    }
    if (isDecayPhotonPi0(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonPi0));
    }
    if (isDecayPhotonEta(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonEta));
    }
    if (isDecayPhotonOther(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonOther));
    }
    if (isPi0(particle)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ParticleOrigin::kPi0));
    }
    return origin;
  }

  /// \brief Gets the index of a mother particle with specific PDG code in the decay chain (upwards)
  /// \param particle The particle to start from
  /// \param mcParticles The MC particles collection
  /// \param pdgCode The PDG code to search for
  /// \return The index of the mother particle, or -1 if not found
  template <typename T>
  int getIndexMotherChain(const T& particle, aod::JMcParticles const& mcParticles, int pdgCode, int depth = 0)
  {
    // Limit recursion depth to avoid infinite loops
    if (depth > kMaxRecursionDepth) { // 100 generations should be more than enough
      return -1;
    }
    const auto& mothers = particle.template mothers_as<aod::JMcParticles>();
    for (const auto& mother : mothers) {
      if (mother.pdgCode() == pdgCode) {
        return mother.globalIndex();
      } else {
        return getIndexMotherChain(mother, mcParticles, pdgCode, depth + 1);
      }
    }
    return -1;
  }
  // return recursive list of all daughter IDs
  /// \brief Gets all daughter particle IDs in the decay chain
  /// \param particle The particle to start from
  /// \return Vector of daughter particle IDs
  template <typename T>
  void getDaughtersInChain(const T& particle, std::vector<int>& daughters, int depth = 0)
  {
    // Limit recursion depth to avoid infinite loops
    if (depth > kMaxRecursionDepth) { // 100 generations should be more than enough
      return;
    }

    if (!particle.has_daughters()) {
      return;
    }

    const auto& daughterParticles = particle.template daughters_as<aod::JMcParticles>();
    for (const auto& daughter : daughterParticles) {
      daughters.push_back(daughter.globalIndex());
      getDaughtersInChain(daughter, daughters, depth + 1);
    }
    return;
  }
  /// \brief Finds the first physical primary particle in the decay chain (upwards)
  /// \param particle The particle to start from
  /// \return The index of the first physical primary particle, or -1 if not found
  template <typename T>
  int findPhysicalPrimaryInChain(const T& particle, int depth = 0)
  {
    // Limit recursion depth to avoid infinite loops
    if (depth > kMaxRecursionDepth) { // 100 generations should be more than enough
      return -1;
    }

    // first check if current particle is physical primary
    if (particle.isPhysicalPrimary()) {
      return particle.globalIndex();
    }

    // check if the particle has mothers
    if (!particle.has_mothers())
      return -1;

    // now get mothers
    const auto mothers = particle.template mothers_as<aod::JMcParticles>();
    if (mothers.size() == 0)
      return -1;

    // get first mother
    for (const auto& mother : mothers) {
      int primaryIndex = findPhysicalPrimaryInChain(mother, depth + 1);
      if (primaryIndex >= 0) {
        return primaryIndex;
      }
      break; // only check first mother
    }

    return -1;
  }

  /// \brief Checks if a cluster is merged from particles of a specific PDG decay to two gammas. A cluster is considered merged if the leading and subleading contribution to a cluster come from two photons that are part of a pi0 decay
  /// \param cluster The cluster to check
  /// \param mcParticles The MC particles collection
  /// \param pdgCode The PDG code to check for
  /// \return true if cluster is merged from the specified decay, false otherwise
  template <typename T, typename U>
  bool isMergedFromPDGDecay(const T& cluster, U const& mcParticles, int pdgCode)
  {
    auto inducerIDs = cluster.mcParticlesIds();
    if (inducerIDs.size() < 2) { // it can not me "merged" if it has less than 2 inducers
      return false;
    }

    bool isMerged = false;
    int motherIndex = getIndexMotherChain(mcParticles.iteratorAt(inducerIDs[0]), mcParticles, pdgCode);
    if (motherIndex != -1) {
      const auto& mother = mcParticles.iteratorAt(motherIndex);

      // get daughters of pi0 mother
      auto daughtersMother = mother.template daughters_as<aod::JMcParticles>();
      // check if there are two daughters that are both photons
      if (daughtersMother.size() == 2) {
        const auto& daughter1 = daughtersMother.iteratorAt(0);
        const auto& daughter2 = daughtersMother.iteratorAt(1);
        if (daughter1.pdgCode() == PDG_t::kGamma && daughter2.pdgCode() == PDG_t::kGamma) {
          // get the full stack of particles that these daughters create
          std::vector<int> fullDecayChain1;
          std::vector<int> fullDecayChain2;
          getDaughtersInChain(daughter1, fullDecayChain1);
          getDaughtersInChain(daughter2, fullDecayChain2);
          bool photon1Found = false;
          bool photon2Found = false;

          // check if any of the particles in the fullDecayChain are leading or subleading in the cluster
          for (const auto& particleID : fullDecayChain1) {
            if (particleID == inducerIDs[0] || particleID == inducerIDs[1]) {
              photon1Found = true;
            }
          }
          for (const auto& particleID : fullDecayChain2) {
            if (particleID == inducerIDs[0] || particleID == inducerIDs[1]) {
              photon2Found = true;
            }
          }
          if (photon1Found && photon2Found) {
            isMerged = true;
          }
        }
      }
    }
    return isMerged;
  }

  // determine cluster origin
  /// \brief Gets the origin bitmap for a cluster
  /// \param cluster The cluster to check
  /// \param mcParticles The MC particles collection
  /// \return A bitmap indicating the cluster's origin
  template <typename T, typename U>
  uint16_t getClusterOrigin(const T& cluster, U const& mcParticles)
  {
    uint16_t origin = 0;
    auto inducerIDs = cluster.mcParticlesIds();
    if (inducerIDs.size() == 0) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kUnknown));
      return origin;
    }

    // loop over all inducers and print their energy
    LOG(debug) << "Cluster with energy: " << cluster.energy() << " and nInducers: " << inducerIDs.size();
    LOG(debug) << "Number of stored amplitudes: " << cluster.amplitudeA().size();
    int aCounter = 0;
    for (const auto& inducerID : inducerIDs) {
      const auto& inducer = mcParticles.iteratorAt(inducerID);
      int motherPDG = -1;
      if (inducer.has_mothers()) {
        motherPDG = inducer.template mothers_as<aod::JMcParticles>()[0].pdgCode();
      }
      LOG(debug) << "Inducer energy: " << inducer.energy() << " amplitude: " << cluster.amplitudeA()[aCounter] << " and PDG: " << inducer.pdgCode() << " isPhysicalPrimary: " << inducer.isPhysicalPrimary() << " motherPDG: " << motherPDG;
      aCounter++;
    }

    // check if leading energy contribution is from a photon
    const auto& leadingParticle = mcParticles.iteratorAt(inducerIDs[0]);
    LOG(debug) << "Leading particle: PDG" << leadingParticle.pdgCode();
    // leading particle primary ID
    int leadingParticlePrimaryID = findPhysicalPrimaryInChain(leadingParticle);
    LOG(debug) << "Leading particle primary ID: " << leadingParticlePrimaryID;
    if (leadingParticlePrimaryID == -1) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kUnknown));
      return origin;
    }
    const auto& leadingParticlePrimary = mcParticles.iteratorAt(leadingParticlePrimaryID);
    LOG(debug) << "Leading particle primary PDG: " << leadingParticlePrimary.pdgCode();
    if (leadingParticlePrimary.pdgCode() == PDG_t::kGamma) {
      LOG(debug) << "Leading particle primary is a photon";
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kPhoton));
    }
    if (isPromptPhoton(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kPromptPhoton));
      LOG(debug) << "Leading particle primary is a prompt photon";
    }
    if (isDirectPromptPhoton(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDirectPromptPhoton));
      LOG(debug) << "Leading particle primary is a direct prompt photon";
    }
    if (isFragmentationPhoton(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kFragmentationPhoton));
      LOG(debug) << "Leading particle primary is a fragmentation photon";
    }
    if (isDecayPhoton(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhoton));
      LOG(debug) << "Leading particle primary is a decay photon";
    }
    if (isDecayPhotonPi0(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhotonPi0));
      LOG(debug) << "Leading particle primary is a decay photon from pi0";
    }
    if (isDecayPhotonEta(leadingParticlePrimary)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhotonEta));
      LOG(debug) << "Leading particle primary is a decay photon from eta";
    }

    // Do checks if a cluster is a merged pi0 decay
    // we classify a cluster as merged pi0 if the leading and subleading contribution to a cluster come from two photons that are part of a pi0 decay
    if (isMergedFromPDGDecay(cluster, mcParticles, PDG_t::kPi0)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kMergedPi0));
      LOG(debug) << "Cluster is a merged pi0";
    }
    if (isMergedFromPDGDecay(cluster, mcParticles, 221)) {
      SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kMergedEta));
      LOG(debug) << "Cluster is a merged eta";
    }

    // check if photon conversion
    // check that leading contribution is an electron or positron
    LOG(debug) << "Checking if cluster is a converted photon";
    if (std::abs(leadingParticle.pdgCode()) == PDG_t::kElectron) {
      // make sure this electron is not a physicsl primary and has mothers
      if (!leadingParticle.isPhysicalPrimary() && leadingParticle.has_mothers()) {
        const auto mothers = leadingParticle.template mothers_as<aod::JMcParticles>();
        if (mothers.size() > 0) {
          LOG(debug) << "Got the mother";
          const auto& mother = mothers[0];
          if (mother.pdgCode() == PDG_t::kGamma && mother.has_daughters()) {
            LOG(debug) << "Got the mother with PDG 22 and daughters";
            const auto& daughters = mother.template daughters_as<aod::JMcParticles>();
            // check that mother has exactly two daughters which are e+ and e-
            if (daughters.size() == 2) {
              LOG(debug) << "Got the daughters";
              if ((daughters.iteratorAt(0).pdgCode() == PDG_t::kElectron && daughters.iteratorAt(1).pdgCode() == PDG_t::kPositron) || (daughters.iteratorAt(0).pdgCode() == -PDG_t::kPositron && daughters.iteratorAt(1).pdgCode() == PDG_t::kElectron)) {
                SETBIT(origin, static_cast<uint16_t>(gjanalysis::ClusterOrigin::kConvertedPhoton));
                LOG(debug) << "Cluster is a converted photon";
              }
            }
          }
        }
      }
    }
    // display bit origin
    LOG(debug) << "Origin bits: " << std::bitset<16>(origin);
    return origin;
  }

  // ---------------------
  // Processing functions
  // ---------------------
  // WARNING: This function always has to run first in the processing chain
  /// \brief Clears collision mapping at the start of each dataframe
  /// \param collisions The collisions collection
  void processClearMaps(aod::JetCollisions const&)
  {
    collisionMapping.clear();
    mcCollisionsMultiRecCollisions.clear();
    mcJetIndexMapping.clear();
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClearMaps, "process function that clears all the maps in each dataframe", true);

  // WARNING: This function always has to run second in the processing chain
  /// \brief Processes MC event matching QA
  /// \param mcCollision The MC collision to process
  /// \param collisions The rec collisions collection
  /// \param mcgenparticles The MC particles collection
  void processMCCollisionsMatching(aod::JetMcCollision const& mcCollision, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, aod::JetParticles const& mcgenparticles)
  {
    if (mcCollision.weight() == 0) {
      return;
    }

    // determine number of rec collisions
    int nRecCollisions = 0;
    mHistograms.fill(HIST("mcCollisionsWithRecCollisions"), 0);
    for (auto const& collision : collisions) {
      if (collision.posZ() > mVertexCut) {
        continue;
      }
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        continue;
      }
      if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
        continue;
      }
      if (!jetderiveddatautilities::eventEMCAL(collision)) {
        continue;
      }
      nRecCollisions++;
    }
    if (nRecCollisions == 0) {
      mHistograms.fill(HIST("mcCollisionsWithRecCollisions"), 3);
    }
    if (nRecCollisions == 1) {
      mHistograms.fill(HIST("mcCollisionsWithRecCollisions"), 1);
    }
    if (nRecCollisions > 1) {
      mHistograms.fill(HIST("mcCollisionsWithRecCollisions"), 2);
      mcCollisionsMultiRecCollisions.push_back(mcCollision.globalIndex());
    }

    // loop over mcgenparticles
    for (auto const& particle : mcgenparticles) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      if (particle.pdgCode() != PDG_t::kGamma) {
        continue;
      }
      mHistograms.fill(HIST("numberRecCollisionsVsPhotonPt"), nRecCollisions, particle.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processMCCollisionsMatching, "Process MC event matching QA", false);

  /// \brief Processes data events in data fill event table
  /// \param collision The collision to process
  /// \param clusters The EMCAL clusters in the event
  void processEventData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, emcClusters const& clusters)
  {
    if (!isEventAccepted(collision, clusters)) {
      return;
    }

    eventsTable(collision.multFT0M(), collision.centFT0M(), collision.rho(), collision.eventSel(), collision.trackOccupancyInTimeRange(), collision.alias_raw());
    collisionMapping[collision.globalIndex()] = eventsTable.lastIndex();
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processEventData, "Process event data", true);

  using MCCol = o2::soa::Join<aod::JMcCollisions, aod::BkgChargedMcRhos>;

  /// \brief Processes MC events and fills rec and MC event tables (disable processEventData)
  /// \param collision The collision to process
  /// \param clusters The EMCAL clusters in the event
  /// \param mcCollisions The MC collisions collection
  void processEventMC(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs, JMcCollisionLbs>::iterator const& collision, emcClusters const& clusters, MCCol const&)
  {
    if (!isEventAccepted(collision, clusters)) {
      return;
    }

    // check that this event has a MC collision
    if (!collision.has_mcCollision()) {
      return;
    }
    mHistograms.fill(HIST("eventQA"), 6);

    // check if this event is not MB gap event
    if (collision.subGeneratorId() == jetderiveddatautilities::JCollisionSubGeneratorId::mbGap) {
      return;
    }
    mHistograms.fill(HIST("eventQA"), 7);

    // fill rec collision table
    eventsTable(collision.multFT0M(), collision.centFT0M(), collision.rho(), collision.eventSel(), collision.trackOccupancyInTimeRange(), collision.alias_raw());

    // fill collision mapping
    collisionMapping[collision.globalIndex()] = eventsTable.lastIndex();

    auto mcCollision = collision.mcCollision_as<MCCol>();

    bool isMultipleAssigned = false;
    // check if we are dealing with a rec collision matched to a MC collision that was matched to multiple rec collisions
    if (std::find(mcCollisionsMultiRecCollisions.begin(), mcCollisionsMultiRecCollisions.end(), mcCollision.globalIndex()) != mcCollisionsMultiRecCollisions.end()) {
      isMultipleAssigned = true;
    }
    mcEventsTable(eventsTable.lastIndex(), mcCollision.weight(), mcCollision.rho(), isMultipleAssigned);
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processEventMC, "Process MC event MC", false);

  // ---------------------
  // Processing functions can be safely added below this line
  // ---------------------

  // define cluster filter. It selects only those clusters which are of the type
  // sadly passing of the string at runtime is not possible for technical region so cluster definition is
  // an integer instead
  PresliceUnsorted<aod::JEMCTracks> EMCTrackPerTrack = aod::jemctrack::trackId;
  // Process clusters
  /// \brief Processes clusters and fills cluster table
  /// \param collision The collision to process
  /// \param clusters The EMCAL clusters to process
  /// \param tracks The tracks collection
  /// \param emctracks The EMCAL tracks collection from track matching
  void processClusters(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, emcClusters const& clusters, aod::JetTracks const& tracks, aod::JEMCTracks const& emctracks)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;

    // loop over tracks one time for QA
    runTrackQA(collision, tracks);

    // build kd tree for tracks and mc particles
    buildKdTree(tracks);

    // loop over clusters
    for (const auto& cluster : clusters) {

      // fill histograms
      mHistograms.fill(HIST("clusterE"), cluster.energy());

      double isoraw = ch_iso_in_cone(cluster, isoR, false);
      double perpconerho = ch_perp_cone_rho(cluster, isoR, false);

      // find closest matched track
      double dEta = 0;
      double dPhi = 0;
      // double dRMin = 100;
      double p = -1;

      // do track matching
      auto tracksofcluster = cluster.matchedTracks_as<aod::JetTracks>();
      for (const auto& track : tracksofcluster) {
        if (!isTrackSelected(track)) {
          continue;
        }
        auto emcTracksPerTrack = emctracks.sliceBy(EMCTrackPerTrack, track.globalIndex());
        auto emcTrack = emcTracksPerTrack.iteratorAt(0);
        // find closest track that still has E/p < trackMatchingEoverP
        if (cluster.energy() / track.p() > trackMatchingEoverP) {
          continue;
        } else {
          dEta = cluster.eta() - emcTrack.etaEmcal();
          dPhi = RecoDecay::constrainAngle(RecoDecay::constrainAngle(emcTrack.phiEmcal(), -o2::constants::math::PI) - RecoDecay::constrainAngle(cluster.phi(), -o2::constants::math::PI), -o2::constants::math::PI);
          p = track.p();
          break;
        }
      }
      gammasTable(storedColIndex, cluster.energy(), cluster.definition(), cluster.eta(), cluster.phi(), cluster.m02(), cluster.m20(), cluster.nCells(), cluster.time(), cluster.isExotic(), cluster.distanceToBadChannel(), cluster.nlm(), isoraw, perpconerho, dPhi, dEta, p);
    }

    // dummy loop over tracks
    for (const auto& track : tracks) {
      mHistograms.fill(HIST("trackPt"), track.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClusters, "Process EMCal clusters", true);

  /// \brief Processes MC cluster information (rec level)
  /// \param collision The collision to process
  /// \param mcClusters The MC clusters to process
  /// \param mcParticles The MC particles collection
  void processClustersMCInfo(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, emcMCClusters const& mcClusters, aod::JMcParticles const& mcParticles)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;
    // loop over mcClusters
    // TODO: add weights
    for (const auto& mcCluster : mcClusters) {
      mHistograms.fill(HIST("clusterMC_E_All"), mcCluster.energy());
      uint16_t origin = getClusterOrigin(mcCluster, mcParticles);
      float leadingEnergyFraction = mcCluster.amplitudeA()[0] / mcCluster.energy();
      // Fill MC origin QA histograms
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_Photon"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_Photon"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kPromptPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_PromptPhoton"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_PromptPhoton"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDirectPromptPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_DirectPromptPhoton"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_DirectPromptPhoton"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kFragmentationPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_FragmentationPhoton"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_FragmentationPhoton"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_DecayPhoton"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_DecayPhoton"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhotonPi0))) {
        mHistograms.fill(HIST("clusterMC_E_DecayPhotonPi0"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_DecayPhotonPi0"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kDecayPhotonEta))) {
        mHistograms.fill(HIST("clusterMC_E_DecayPhotonEta"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_DecayPhotonEta"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kMergedPi0))) {
        mHistograms.fill(HIST("clusterMC_E_MergedPi0"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_MergedPi0"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kMergedEta))) {
        mHistograms.fill(HIST("clusterMC_E_MergedEta"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_MergedEta"), mcCluster.m02());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ClusterOrigin::kConvertedPhoton))) {
        mHistograms.fill(HIST("clusterMC_E_ConvertedPhoton"), mcCluster.energy());
        mHistograms.fill(HIST("clusterMC_m02_ConvertedPhoton"), mcCluster.m02());
      }
      // fill table
      gammaMCInfosTable(storedColIndex, origin, leadingEnergyFraction);
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processClustersMCInfo, "Process MC cluster information", false);

  /// \brief Fills the charged jet table with jet information and calculates jet properties
  /// \param storedColIndex The stored collision index
  /// \param jet The jet to process
  /// \param tracks The tracks collection
  template <typename T, typename U>
  void fillChargedJetTable(int32_t storedColIndex, T const& jet, U const& /*tracks*/)
  {
    if (jet.pt() < jetPtMin) {
      return;
    }
    ushort nconst = 0;
    float leadingTrackPt = 0;
    for (const auto& constituent : jet.template tracks_as<aod::JetTracks>()) {
      mHistograms.fill(HIST("chjetpt_vs_constpt"), jet.pt(), constituent.pt());
      nconst++;
      if (constituent.pt() > leadingTrackPt) {
        leadingTrackPt = constituent.pt();
      }
    }
    double perpconerho = ch_perp_cone_rho(jet, perpConeJetR, false);
    chargedJetsTable(storedColIndex, jet.pt(), jet.eta(), jet.phi(), jet.r(), jet.energy(), jet.mass(), jet.area(), leadingTrackPt, perpconerho, nconst);
    mHistograms.fill(HIST("chjetPtEtaPhi"), jet.pt(), jet.eta(), jet.phi());
    mHistograms.fill(HIST("chjetPt"), jet.pt());
  }

  Filter jetCuts = aod::jet::pt > jetPtMin;
  /// \brief Processes charged jets and fills jet table
  /// \param collision The collision to process
  /// \param chargedJets The charged jets to process
  /// \param tracks The tracks collection
  void processChargedJetsData(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& chargedJets, aod::JetTracks const& tracks)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;
    // loop over charged jets
    for (const auto& jet : chargedJets) {
      fillChargedJetTable(storedColIndex, jet, tracks);
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processChargedJetsData, "Process charged jets", true);

  Preslice<aod::JetParticles> ParticlesPerMCCollisions = aod::jmcparticle::mcCollisionId;
  /// \brief Processes MC particles and fills MC particle table
  /// \param collision The collision to process
  /// \param mcgenparticles The MC particles to process
  void processMCParticles(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs, JMcCollisionLbs>::iterator const& collision, aod::JetParticles const& mcgenparticles, MCCol const&)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;

    if (!collision.has_mcCollision()) {
      return;
    }

    // only storing MC particles if we found a reconstructed collision
    auto particlesPerMcCollision = mcgenparticles.sliceBy(ParticlesPerMCCollisions, collision.mcCollisionId());

    // build kd tree for mc particles
    buildKdTree(particlesPerMcCollision);

    // Now we want to store every pi0 and every prompt photon that we find on generator level
    for (const auto& particle : particlesPerMcCollision) {
      // only store particles above a given threshold
      if (particle.pt() < minMCGenPt) {
        continue;
      }
      // Test if a particle is a physical primary according to the following definition:
      // Particles produced in the collision including products of strong and
      // electromagnetic decay and excluding feed-down from weak decays of strange
      // particles.
      if (!(particle.isPhysicalPrimary() || particle.pdgCode() == PDG_t::kPi0)) {
        continue;
      }

      // only store photons and pi0s in mcgen stack
      if (particle.pdgCode() != PDG_t::kPi0 && particle.pdgCode() != PDG_t::kGamma) {
        continue;
      }
      // check the origin of the particle
      uint16_t origin = getMCParticleOrigin(particle);
      double mcIsolation = ch_iso_in_cone(particle, isoR, true);
      mcParticlesTable(storedColIndex, particle.energy(), particle.eta(), particle.phi(), particle.pt(), particle.pdgCode(), mcIsolation, origin);

      // fill mc gen trigger particle histograms
      mHistograms.fill(HIST("mcGenTrigger_E"), particle.energy());
      mHistograms.fill(HIST("mcGenTrigger_Eta"), particle.eta());
      mHistograms.fill(HIST("mcGenTrigger_Phi"), particle.phi());
      mHistograms.fill(HIST("mcGenTrigger_Pt"), particle.pt());
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kPromptPhoton))) {
        mHistograms.fill(HIST("mcGenTrigger_E_PromptPhoton"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDirectPromptPhoton))) {
        mHistograms.fill(HIST("mcGenTrigger_E_DirectPromptPhoton"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kFragmentationPhoton))) {
        mHistograms.fill(HIST("mcGenTrigger_E_FragmentationPhoton"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhoton))) {
        mHistograms.fill(HIST("mcGenTrigger_E_DecayPhoton"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonPi0))) {
        mHistograms.fill(HIST("mcGenTrigger_E_DecayPhotonPi0"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonEta))) {
        mHistograms.fill(HIST("mcGenTrigger_E_DecayPhotonEta"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kDecayPhotonOther))) {
        mHistograms.fill(HIST("mcGenTrigger_E_DecayPhotonOther"), particle.energy());
      }
      if (origin & (1 << static_cast<uint16_t>(gjanalysis::ParticleOrigin::kPi0))) {
        mHistograms.fill(HIST("mcGenTrigger_E_Pi0"), particle.energy());
      }
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processMCParticles, "Process MC particles", false);

  // NOTE: It is important that this function runs after the processMCParticles function (where the isolation tree is built )
  Preslice<aod::ChargedMCParticleLevelJets> PJetsPerMCCollisions = aod::jmcparticle::mcCollisionId;
  /// \brief Processes MC particle level charged jets and fills MC jet table
  /// \param collision The collision to process
  /// \param chargedJets The MC particle level charged jets to process
  /// \param mcCollisions The MC collisions collection
  void processChargedJetsMCP(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs, JMcCollisionLbs>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>> const& chargedJets, MCCol const&)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;
    // loop over charged jets
    if (!collision.has_mcCollision()) {
      return;
    }
    int localIndex = 0;
    auto pjetsPerMcCollision = chargedJets.sliceBy(PJetsPerMCCollisions, collision.mcCollisionId());
    for (const auto& pjet : pjetsPerMcCollision) {
      // fill MC particle level jet table
      float perpconerho = ch_perp_cone_rho(pjet, perpConeJetR, true);
      mcJetsTable(storedColIndex, pjet.pt(), pjet.eta(), pjet.phi(), pjet.r(), pjet.energy(), pjet.mass(), pjet.area(), perpconerho);
      mcJetIndexMapping[pjet.globalIndex()] = localIndex;
      localIndex++;
      mHistograms.fill(HIST("mcpJetPt"), pjet.pt());
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processChargedJetsMCP, "Process MC particle level jets", false);

  // NOTE: It is important that this function runs after the processChargedJetsMCP function (where the mc jet index mapping is built)
  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  Filter jetCutsMCD = aod::jet::pt > jetPtMin;
  /// \brief Processes MC detector level charged jets and fills jet matching information
  /// \param collision The collision to process
  /// \param chargedJets The MC detector level charged jets to process
  /// \param tracks The tracks collection
  /// \param pjets The MC particle level jets collection (just loaded to have subscription to the table)
  void processChargedJetsMCD(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::JCollisionBCs>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& chargedJets, aod::JetTracks const& tracks, JetMCPTable const& /*pjets*/)
  {
    // event selection
    int32_t storedColIndex = getStoredColIndex(collision);
    if (storedColIndex == -1)
      return;
    // loop over charged jets
    for (const auto& jet : chargedJets) {
      fillChargedJetTable(storedColIndex, jet, tracks);

      // Fill Matching information
      int iLocalIndexGeo = -1;
      int iLocalIndexPt = -1;
      // We will always store the information for both in our tree
      if (jet.has_matchedJetGeo()) {
        const auto& pjet = jet.template matchedJetGeo_first_as<JetMCPTable>();
        iLocalIndexGeo = mcJetIndexMapping[pjet.globalIndex()];
        mHistograms.fill(HIST("mcdJetPtVsTrueJetPtMatchingGeo"), jet.pt(), pjet.pt());
      }
      if (jet.has_matchedJetPt()) {
        const auto& pjet = jet.template matchedJetPt_first_as<JetMCPTable>();
        iLocalIndexPt = mcJetIndexMapping[pjet.globalIndex()];
        mHistograms.fill(HIST("mcdJetPtVsTrueJetPtMatchingPt"), jet.pt(), pjet.pt());
      }
      chJetMCInfosTable(storedColIndex, iLocalIndexGeo, iLocalIndexPt);
    }
  }
  PROCESS_SWITCH(GammaJetTreeProducer, processChargedJetsMCD, "Process MC detector level jets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<GammaJetTreeProducer>(cfgc, TaskName{"gamma-jet-tree-producer"})};
  return workflow;
}
