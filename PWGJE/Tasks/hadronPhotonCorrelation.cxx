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
/// \file   hadronPhotonCorrelation.cxx
/// \author Peter Stratmann <peter.stratmann@cern.ch>
/// \brief This code loops over JetTracks to extract pt and angular information
///        for hadrons and photons to compute angular correlations
///

#include <map>
#include <string>
#include <concepts>
#include <vector>

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Expressions.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramSpec.h"
#include "Framework/Configurable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/MathConstants.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/Core/JetUtilities.h"

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;

struct HadronPhotonCorrelation {

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<float> etaMax{"etaMax", 0.8, "maximum eta cut"};

  AxisSpec axisPhi = {72, 0., TwoPI, "#phi"};                       // Axis for phi distribution
  AxisSpec axisDeltaPhi = {72, -PIHalf, 3 * PIHalf, "#Delta #phi"}; // Axis for Delta phi in correlations
  AxisSpec axisEta = {40, -.8, .8, "#eta"};                         // Axis for eta distribution
  AxisSpec axisDeltaEta = {80, -1.6, 1.6, "#Delta #eta"};           // Axis for Delta eta in correlations
  ConfigurableAxis axisPtTrig = {"axisPtTrig",
                                 {VARIABLE_WIDTH,
                                  3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f, 20.0f},
                                 "p_{T, trig} [GeV]"}; // Axis for trigger particle pt distribution
  ConfigurableAxis axisPtAssoc = {"axisPtAssoc",
                                  {VARIABLE_WIDTH,
                                   0.2, 0.5, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f},
                                  "p_{T, assoc} [GeV]"}; // Axis for associated particle pt distribution
  AxisSpec axisDeltaPt = {200, 0., 1.2, "#Delta p_T"};   // Axis for pt ratio between neutral hadrons and decay photons
  AxisSpec axisPid = {9, -3.5, 5.5, "pid"};              // Axis for PID of neutral hadrons
  AxisSpec axisMult = {100, 0., 99., "N_{ch}"};          // Axis for mutplipicity
  AxisSpec axisAlpha = {100, 0., 1., "alpha"};           // Axis for decay photon pt assymetry

  AxisSpec axisDeltaRDecay = {100, 0., 0.8, "#Delta R"}; // Axis for Delta R = sqrt(Delta eta^2 + Delta phi^2) between neutral hadrons and decay photons

  float ptMinTrig;
  float ptMaxTrig;
  float ptMinAssoc;
  float ptMaxAssoc;

  HistogramRegistry registry{"histogram registry"};

  // Particle ids for storing neutral hadrons
  std::map<std::string, int> pidCodes = {
    {"pi0", 1},    // pi0
    {"eta", 2},    // eta
    {"eta'", 3},   // eta'
    {"phi", 4},    // phi
    {"omega", 5}}; // omega

  Service<o2::framework::O2DatabasePDG> pdg;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    ptMinTrig = axisPtTrig->at(1);
    ptMaxTrig = axisPtTrig->back();
    ptMinAssoc = axisPtAssoc->at(1);
    ptMaxAssoc = axisPtAssoc->back();

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // Generated histograms
    //  Triggers
    registry.add("generated/triggers/hTrigMultGen", "Generated Trigger Multiplicity", kTH1F, {axisMult});
    registry.add("generated/triggers/hTrigSpectrumGen", "Generated Trigger Spectrum", kTHnSparseF, {axisPtTrig, axisEta, axisPhi});

    // Hadrons
    registry.add("generated/hadrons/hHadronCorrelGen", "Generated Trigger-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.add("generated/hadrons/hHadronMultGen", "Generated Hadron Multiplicity", kTH1F, {axisMult});
    registry.add("generated/hadrons/hHadronSpectrumGen", "Generated Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi});

    // Photons
    registry.add("generated/photons/hPhotonCorrelGen", "Generated Trigger-Photon Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.add("generated/photons/hPhotonMultGen", "Generated Photon Multiplicity", kTH1F, {axisMult});
    registry.add("generated/photons/hPhotonSpectrumGen", "Generated Photon Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi, axisPid});

    // Charged pions
    registry.add("generated/charged/hPionCorrelGen", "Generated Trigger-Pion Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.add("generated/charged/hPionMultGen", "Generated Pion Multiplicity", kTH1F, {axisMult});
    registry.add("generated/charged/hPionSpectrumGen", "Generated Pion Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi});

    ////Neutral particles
    registry.add("generated/neutral/hNeutralCorrelGen", "Generated Trigger-Neutral Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.add("generated/neutral/hNeutralMultGen", "Generated Neutral Hadron Multiplicity", kTH1F, {axisMult});
    registry.add("generated/neutral/hNeutralSpectrumGen", "Generated Neutral Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi, axisPid});                                     // Particle ID of neutral hadrons
    registry.add("generated/neutral/hNeutralDecayGen", "Generated Neutral Hadron-Decay Photon Correlation", kTHnSparseF, {axisPtAssoc, axisDeltaPt, axisDeltaRDecay, axisAlpha, axisPid}); // Correlation with decay photons

    // Reconstructed histograms
    // Triggers
    registry.add("reconstructed/triggers/hTrigMultReco", "Reconstructed Trigger Multiplicity", kTH1F, {axisMult});
    registry.add("reconstructed/triggers/hTrigSpectrumReco", "Reconstructed Trigger Spectrum", kTHnSparseF, {axisPtTrig, axisEta, axisPhi});

    // Hadrons
    registry.add("reconstructed/hadrons/hHadronCorrelReco", "Reconstructed Trigger-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.add("reconstructed/hadrons/hHadronMultReco", "Reconstructed Hadron Multiplicity", kTH1F, {axisMult});
    registry.add("reconstructed/hadrons/hHadronSpectrumReco", "Reconstructed Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi});
    registry.add("reconstructed/hadrons/hHadronPtPrimReco", "Reconstructed Primaries Spectrum", kTH1F, {axisPtAssoc});  // Primary hadron spectrum
    registry.add("reconstructed/hadrons/hHadronPtSecReco", "Reconstructed Secondaries Spectrum", kTH1F, {axisPtAssoc}); // Secondary hadron spectrum

    // Photons
    registry.add("reconstructed/photons/hPhotonCorrelReco", "Reconstructed Trigger-Photon Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.add("reconstructed/photons/hPhotonMultReco", "Reconstructed Photon Multiplicity", kTH1F, {axisMult});
    registry.add("reconstructed/photons/hPhotonSpectrumReco", "Reconstructed Photon Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi});

    // Charged Pions
    registry.add("reconstructed/charged/hPionCorrelReco", "Reconstructed Trigger-Pion Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.add("reconstructed/charged/hPionMultReco", "Reconstructed Pion Multiplicity", kTH1F, {axisMult});
    registry.add("reconstructed/charged/hPionSpectrumReco", "Reconstructed Pion Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi});
  }

  // To check if object has has_mcParticle() (i.e. is MC Track or data track)
  template <typename T>
  struct HasHasMcParticle {
    template <typename U>
    static auto test(U* ptr) -> decltype(ptr->has_mcParticle(), std::true_type{});
    static std::false_type test(...);
    static constexpr bool Value = decltype(test(std::declval<T*>()))::value;
  };

  // Initialize track
  template <typename T>
  bool initTrack(const T& track)
  {

    // if constexpr (HasHasMcParticle<T>) {
    if constexpr (HasHasMcParticle<T>::Value) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (std::abs(track.eta()) > etaMax) {
      return false;
    }

    if (track.pt() < ptMinAssoc || track.pt() > ptMaxAssoc) {
      return false;
    }

    return true;
  }

  // Initialize particle
  template <typename T>
  bool initParticle(const T& particle, bool checkIsPrimary = true)
  {

    if (checkIsPrimary && !particle.isPhysicalPrimary()) {
      return false;
    }

    if (std::abs(particle.eta()) > etaMax) {
      return false;
    }

    if (particle.pt() < ptMinAssoc || particle.pt() > ptMaxAssoc) {
      return false;
    }

    return true;
  }

  // Initialize trigger tracks
  template <typename T>
  bool initTrig(const T& track)
  {
    // if constexpr (HasHasMcParticle<T>) {
    if constexpr (HasHasMcParticle<T>::Value) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (std::abs(track.eta()) > etaMax) {
      return false;
    }

    if (track.pt() < ptMinTrig || track.pt() > ptMaxTrig) {
      return false;
    }

    return true;
  }

  // Initialize trigger particles (charged)
  template <typename T>
  bool initTrigParticle(const T& particle)
  {
    if (!particle.isPhysicalPrimary()) {
      return false;
    }

    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle || pdgParticle->Charge() == 0.) {
      return false;
    }

    if (std::abs(particle.eta()) > etaMax) {
      return false;
    }

    if (particle.pt() < ptMinTrig || particle.pt() > ptMaxTrig) {
      return false;
    }

    return true;
  }

  // Initialize V0s
  template <typename TTrack, typename TV0>
  bool initV0(TV0 const& v0)
  {
    auto pos = v0.template posTrack_as<TTrack>();
    auto neg = v0.template negTrack_as<TTrack>();
    if (!initV0leg(pos) || !initV0leg(neg)) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool initV0leg(TTrack const& track)
  {
    if (!initTrack(track)) {
      return false;
    }

    if (!track.hasTPC()) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }

    return true;
  }

  /****************************************************************************************************
  ************************************************ TRIGGER ********************************************
  ****************************************************************************************************/

  /********************************************** DATA ***********************************************/

  void processTrigsReco(JetCollision const& collision,
                        JetTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    int nTrigs = 0;
    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      if (!initTrig(track)) {
        continue;
      }

      registry.fill(HIST("reconstructed/triggers/hTrigSpectrumReco"), track.pt(), track.eta(), track.phi());
      nTrigs++;
    }
    registry.fill(HIST("reconstructed/triggers/hTrigMultReco"), nTrigs);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processTrigsReco, "trigger particle properties", true);

  /*********************************************** MC ************************************************/

  void processTrigsMCReco(JetCollision const& collision,
                          Join<JetTracks, JMcTrackLbs> const& tracks,
                          JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    int nTrigs = 0;
    for (const auto& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }

      if (!initTrig(track)) {
        continue;
      }
      registry.fill(HIST("reconstructed/triggers/hTrigSpectrumReco"), track.pt(), track.eta(), track.phi());
      nTrigs++;
    }
    registry.fill(HIST("reconstructed/triggers/hTrigMultReco"), nTrigs);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processTrigsMCReco, "trigger particle mc properties", true);

  void processTrigsMCGen(JetMcCollision const&,
                         JetParticles const& particles)
  {

    int nTrigs = 0;
    for (const auto& particle : particles) {
      if (!initTrigParticle(particle)) {
        continue;
      }
      registry.fill(HIST("generated/triggers/hTrigSpectrumGen"), particle.pt(), particle.eta(), particle.phi());
      nTrigs++;
    }
    registry.fill(HIST("generated/triggers/hTrigMultGen"), nTrigs);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processTrigsMCGen, "trigger particle mc properties", true);

  /****************************************************************************************************
  ********************************************** PHOTONS **********************************************
  ****************************************************************************************************/

  /********************************************** DATA ***********************************************/
  using MyTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  Preslice<aod::V0Datas> perCol = aod::v0::collisionId;
  void processPhotonCorrelations(JetCollision const& collision,
                                 JetTracks const& tracks,
                                 MyTracks const&,
                                 V0Datas const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    int nPhotons = 0;
    auto v0PerCollision = v0s.sliceBy(perCol, collision.globalIndex());
    for (const auto& v0 : v0PerCollision) {
      if (!initV0<MyTracks>(v0)) {
        continue;
      }
      registry.fill(HIST("reconstructed/photons/hPhotonSpectrumReco"), v0.pt(), v0.eta(), v0.phi());
      nPhotons++;

      for (const auto& track : tracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }

        if (!initTrig(track)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track.phi() - v0.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), track.eta() - v0.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelations, "hadron-photon correlation", true);

  /*********************************************** MC ************************************************/

  using MyTracksMC = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  void processPhotonCorrelationsMCReco(Join<JetCollisions, JCollisionPIs, JMcCollisionLbs>::iterator const& collision_reco,
                                       JetMcCollisions const&,
                                       JetTracks const& tracks_reco,
                                       JetParticles const&,
                                       MyTracksMC const&,
                                       V0Datas const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits)) {
      return;
    }

    int nPhotons = 0;

    auto v0PerCollision = v0s.sliceBy(perCol, collision_reco.globalIndex());
    for (const auto& v0 : v0PerCollision) {
      if (!initV0<MyTracksMC>(v0)) {
        continue;
      }
      registry.fill(HIST("reconstructed/photons/hPhotonSpectrumReco"), v0.pt(), v0.eta(), v0.phi());
      nPhotons++;

      for (const auto& track : tracks_reco) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }

        if (!initTrig(track)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track.phi() - v0.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), track.eta() - v0.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelationsMCReco, "hadron-photon correlation", true);

  void processPhotonCorrelationsMCGen(JetMcCollision const&,
                                      JetParticles const& tracks_true)
  {
    int nPhotons = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc, false)) {
        continue;
      }
      if ((PDG_t)std::abs(track_assoc.pdgCode()) != kGamma) {
        continue;
      }
      if (!track_assoc.isPhysicalPrimary() && track_assoc.getGenStatusCode() < 0) {
        continue;
      }

      // Iterate through mother particles until original mother is reached
      auto origPhoton = track_assoc;
      auto mother = track_assoc.mothers_as<JetParticles>().at(0);
      while ((PDG_t)mother.pdgCode() == kGamma) {
        origPhoton = mother;
        mother = mother.mothers_as<JetParticles>().at(0);
      }

      auto pdgMother = pdg->GetParticle(mother.pdgCode());
      if (!pdgMother)
        continue;
      int photonGeneration;
      switch (std::abs(origPhoton.getGenStatusCode())) {
        case 23: // prompt direct photons
        case 33:
          photonGeneration = -2;
          break;
        case 43: // shower photons
        case 44:
        case 51:
        case 52:
          photonGeneration = -1;
          break;
        case 91: // decay photons
          photonGeneration = pidCodes[pdgMother->GetName()];
          break;
        default:
          photonGeneration = -3;
          break;
      }

      registry.fill(HIST("generated/photons/hPhotonSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), photonGeneration);

      nPhotons++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (!initParticle(track_assoc)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);
        registry.fill(HIST("generated/photons/hPhotonCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi, photonGeneration);
      }
    }

    registry.fill(HIST("generated/photons/hPhotonMultGen"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelationsMCGen, "mc hadron-photon correlation", true);

  /****************************************************************************************************
  ***************************************** HADRONS ***************************************************
  ****************************************************************************************************/

  /********************************************** DATA ***********************************************/
  void processHadronCorrelations(JetCollision const& collision,
                                 Join<JetTracks, pidTPCEl, pidTPCMu> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    int nHadrons = 0;
    for (const auto& track_assoc : tracks) {
      if (!jetderiveddatautilities::selectTrack(track_assoc, trackSelection)) {
        continue;
      }

      if (!initTrack(track_assoc)) {
        continue;
      }
      registry.fill(HIST("reconstructed/hadrons/hHadronSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nHadrons++;

      for (const auto& track_trig : tracks) {
        if (!jetderiveddatautilities::selectTrack(track_trig, trackSelection)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/hadrons/hadrons/hHadronCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/hadrons/hHadronMultReco"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelations, "hadron-hadron correlation", true);

  /*********************************************** MC ************************************************/

  void processHadronCorrelationsMCGen(JetMcCollision const&,
                                      JetParticles const& tracks_true)
  {
    int nHadrons = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc)) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(track_assoc.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.) {
        continue;
      }
      if (std::abs(track_assoc.pdgCode()) < 100) {
        continue;
      }

      registry.fill(HIST("generated/hadrons/hHadronSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nHadrons++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);

        registry.fill(HIST("generated/hadrons/hHadronCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("generated/hadrons/hHadronMultGen"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelationsMCGen, "mc hadron-hadron correlation", true);

  void processHadronCorrelationsMCReco(Join<JetCollisions, JMcCollisionLbs>::iterator const& collision_reco,
                                       JetMcCollisions const&,
                                       Join<JetTracks, JMcTrackLbs> const& tracks_reco,
                                       JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits)) {
      return;
    }

    int nHadrons = 0;
    for (const auto& track_assoc : tracks_reco) {
      if (!jetderiveddatautilities::selectTrack(track_assoc, trackSelection)) {
        continue;
      }

      if (!initTrack(track_assoc)) {
        continue;
      }
      auto particle = track_assoc.mcParticle();
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.) {
        continue;
      }
      if (std::abs(particle.pdgCode()) < 100) {
        continue;
      }

      registry.fill(HIST("reconstructed/hadrons/hHadronSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("reconstructed/hadrons/hHadronPtPrimReco"), track_assoc.pt());
      } else {
        registry.fill(HIST("reconstructed/hadrons/hHadronPtSecReco"), track_assoc.pt());
      }
      nHadrons++;

      for (const auto& track_trig : tracks_reco) {
        if (!jetderiveddatautilities::selectTrack(track_trig, trackSelection)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);

        registry.fill(HIST("reconstructed/hadrons/hHadronCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/hadrons/hHadronMultReco"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelationsMCReco, "mc hadron-hadron correlation", true);

  /****************************************************************************************************
  *************************************** CHARGED PIONS ***********************************************
  ****************************************************************************************************/

  /********************************************** DATA ***********************************************/
  void processPionCorrelations(JetCollision const& collision,
                               Join<JetTracks, pidTPCPi> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    int nPions = 0;
    for (const auto& track_assoc : tracks) {
      if (!jetderiveddatautilities::selectTrack(track_assoc, trackSelection)) {
        continue;
      }

      if (!initTrack(track_assoc)) {
        continue;
      }
      if (std::abs(track_assoc.tpcNSigmaPi()) > 2) {
        continue;
      } // remove non-pions
      registry.fill(HIST("reconstructed/charged/hPionSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nPions++;

      for (const auto& track_trig : tracks) {
        if (!jetderiveddatautilities::selectTrack(track_trig, trackSelection)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/charged/hPionCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/charged/hPionMultReco"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelations, "hadron-pion correlation", true);

  /*********************************************** MC ************************************************/

  void processPionCorrelationsMCGen(JetMcCollision const&,
                                    JetParticles const& tracks_true)
  {
    int nPions = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc)) {
        continue;
      }
      if ((PDG_t)std::abs(track_assoc.pdgCode()) != kPiPlus) {
        continue;
      }

      registry.fill(HIST("generated/charged/hPionSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nPions++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);

        registry.fill(HIST("generated/charged/hPionCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("generated/charged/hPionMultGen"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelationsMCGen, "mc hadron-pion correlation", true);

  void processPionCorrelationsMCReco(Join<JetCollisions, JMcCollisionLbs>::iterator const& collision_reco,
                                     JetMcCollisions const&,
                                     Join<JetTracks, JMcTrackLbs> const& tracks_reco,
                                     JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits)) {
      return;
    }

    int nPions = 0;
    for (const auto& track_assoc : tracks_reco) {
      if (!jetderiveddatautilities::selectTrack(track_assoc, trackSelection)) {
        continue;
      }

      if (!initTrack(track_assoc)) {
        continue;
      }
      if ((PDG_t)std::abs(track_assoc.mcParticle().pdgCode()) != kPiPlus) {
        continue;
      }

      registry.fill(HIST("reconstructed/charged/hPionSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nPions++;

      for (const auto& track_trig : tracks_reco) {
        if (!jetderiveddatautilities::selectTrack(track_trig, trackSelection)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }

        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_trig.phi() - track_assoc.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/charged/hPionCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/charged/hPionMultReco"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelationsMCReco, "mc hadron-pion correlation", true);

  /****************************************************************************************************
  ****************************************** NEUTRALS *************************************************
  ****************************************************************************************************/

  /*********************************************** MC ************************************************/

  void processNeutralCorrelationsMCGen(JetMcCollision const&,
                                       JetParticles const& tracks_true)
  {
    int nNeutrals = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc, false)) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(track_assoc.pdgCode());
      if (!pdgParticle) {
        continue;
      } // remove unknown particles
      if (pdgParticle->Charge() != 0.) {
        continue;
      } // remove charged particles
      if (track_assoc.pdgCode() < 100 || (PDG_t)track_assoc.pdgCode() == kNeutron) {
        continue;
      } // remove non-hadrons and neutrons
      registry.fill(HIST("generated/neutral/hNeutralSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pidCodes[pdgParticle->GetName()]);
      nNeutrals++;

      // Get correlations between neutral hadrons and their respective decay photons
      auto daughters = track_assoc.daughters_as<aod::JMcParticles>();
      double alpha = -1;
      if (daughters.size() == 2) {
        auto daughter = daughters.begin();
        double pt1 = daughter.pt();
        ++daughter;
        double pt2 = daughter.pt();
        alpha = std::abs((pt1 - pt2) / (pt1 + pt2));
      }

      for (const auto& daughter : daughters) {
        if ((PDG_t)std::abs(daughter.pdgCode()) != kGamma)
          continue;
        if (!initParticle(daughter, false))
          continue;
        if (!daughter.isPhysicalPrimary() && daughter.getGenStatusCode() == -1)
          continue;
        double deltaPt = daughter.pt() / track_assoc.pt();
        double deltaEta = daughter.eta() - track_assoc.eta();
        double deltaPhi = RecoDecay::constrainAngle(daughter.phi() - track_assoc.phi(), -PIHalf);
        double deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

        registry.fill(HIST("generated/neutral/hNeutralDecayGen"), track_assoc.pt(), deltaPt, deltaR, alpha, pidCodes[pdgParticle->GetName()]);
      }

      // Get correlations between triggers and neutral hadrons
      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("generated/neutral/hNeutralCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, pidCodes[pdgParticle->GetName()]);
      }
    }
    registry.fill(HIST("generated/neutral/hNeutralMultGen"), nNeutrals);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processNeutralCorrelationsMCGen, "mc hadron-pion correlation", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HadronPhotonCorrelation>(cfgc)};
};
