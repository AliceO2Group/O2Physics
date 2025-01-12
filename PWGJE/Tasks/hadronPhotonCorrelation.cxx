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
  AxisSpec axisDeltaPhiDecay = {100, -.5, .5, "#Delta #phi"};       // Axis for Delta phi between neutral hadrons and decay photons
  AxisSpec axisEta = {40, -.8, .8, "#eta"};                         // Axis for eta distribution
  AxisSpec axisDeltaEta = {80, -1.6, 1.6, "#Delta #eta"};           // Axis for Delta eta in correlations
  AxisSpec axisDeltaEtaDecay = {100, -.6, .6, "#Delta #eta"};       // Axis for Delta eta between neutral hadrons and decay photons
  ConfigurableAxis axisPtTrig = {"axisPtTrig",
                                 {VARIABLE_WIDTH,
                                  3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f, 20.0f},
                                 "p_{T, trig} [GeV]"}; // Axis for trigger particle pt distribution
  ConfigurableAxis axisPtAssoc = {"axisPtAssoc",
                                  {VARIABLE_WIDTH,
                                   0.2, 0.5, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f},
                                  "p_{T, assoc} [GeV]"}; // Axis for associated particle pt distribution
  AxisSpec axisDeltaPt = {200, 0., 1.2, "#Delta p_T"};   // Axis for pt ratio between neutral hadrons and decay photons
  AxisSpec axisPid = {7, -1.5, 5.5, "pid"};              // Axis for PID of neutral hadrons
  AxisSpec axisMult = {100, 0., 99., "N_{ch}"};          // Axis for mutplipicity

  float ptMinTrig;
  float ptMaxTrig;
  float ptMinAssoc;
  float ptMaxAssoc;

  HistogramRegistry registry{"histogram registry"};

  // Particle ids for storing neutral hadrons
  std::map<int, int> pidCodes = {
    {2212, -1},
    {1, -1},
    {2, -1},
    {3, -1},
    {4, -1},
    {5, -1},
    {6, -1},
    {21, -1},  // Protons, quarks, gluons (direct)
    {111, 1},  // pi0
    {221, 2},  // eta
    {223, 3},  // eta'
    {331, 4},  // phi
    {333, 5}}; // omega

  Service<o2::framework::O2DatabasePDG> pdg;

  // Calculate difference between two azimuthal angles, projecting into range [lowerPhi, lowerPhi + pi/2]
  float calculateDelta(float phi1, float phi2, float lowerPhi)
  {
    float dphi = fmod((phi1 - phi2) + 3 * PI, TwoPI) - PI;
    dphi = RecoDecay::constrainAngle(dphi, lowerPhi);
    return dphi;
  }

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    ptMinTrig = axisPtTrig->at(1);
    ptMaxTrig = axisPtTrig->back();
    ptMinAssoc = axisPtAssoc->at(1);
    ptMaxAssoc = axisPtAssoc->back();

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
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
    registry.add("generated/neutral/hNeutralSpectrumGen", "Generated Neutral Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEta, axisPhi, axisPid});                                               // Particle ID of neutral hadrons
    registry.add("generated/neutral/hNeutralDecayGen", "Generated Neutral Hadron-Decay Photon Correlation", kTHnSparseF, {axisPtAssoc, axisDeltaPt, axisDeltaEtaDecay, axisDeltaPhiDecay, axisPid}); // Correlation with decay photons

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

  void processTrigsReco(JCollision const& collision,
                        JTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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

  void processTrigsMCReco(JCollision const& collision,
                          Join<JTracks, JMcTrackLbs> const& tracks,
                          JMcParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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

  void processTrigsMCGen(JMcCollision const&,
                         JMcParticles const& particles)
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
  void processPhotonCorrelations(JCollision const& collision,
                                 JTracks const& tracks,
                                 MyTracks const&,
                                 V0Datas const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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
        float dphi = calculateDelta(track.phi(), v0.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), track.eta() - v0.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelations, "hadron-photon correlation", true);

  /*********************************************** MC ************************************************/

  using MyTracksMC = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  void processPhotonCorrelationsMCReco(Join<JCollisions, JCollisionPIs, JMcCollisionLbs>::iterator const& collision_reco,
                                       JMcCollisions const&,
                                       JTracks const& tracks_reco,
                                       JMcParticles const&,
                                       MyTracksMC const&,
                                       V0Datas const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelection)) {
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
        float dphi = calculateDelta(track.phi(), v0.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), track.eta() - v0.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelationsMCReco, "hadron-photon correlation", true);

  void processPhotonCorrelationsMCGen(JMcCollision const&,
                                      JMcParticles const& tracks_true)
  {
    int nPhotons = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc, false)) {
        continue;
      }
      if (std::abs(track_assoc.pdgCode()) != 22) {
        continue;
      }
      if (!track_assoc.isPhysicalPrimary() && track_assoc.getGenStatusCode() == -1) {
        continue;
      }

      // Iterate through mother particles until original mother is reached
      auto mother = track_assoc.mothers_as<aod::JMcParticles>().at(0);
      while (mother.pdgCode() == 22) {
        mother = mother.mothers_as<aod::JMcParticles>().at(0);
      }

      registry.fill(HIST("generated/photons/hPhotonSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pidCodes[mother.pdgCode()]);

      nPhotons++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (!initParticle(track_assoc)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);
        registry.fill(HIST("generated/photons/hPhotonCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi, pidCodes[mother.pdgCode()]);
      }
    }

    registry.fill(HIST("generated/photons/hPhotonMultGen"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelationsMCGen, "mc hadron-photon correlation", true);

  /****************************************************************************************************
  ***************************************** HADRONS ***************************************************
  ****************************************************************************************************/

  /********************************************** DATA ***********************************************/
  void processHadronCorrelations(JCollision const& collision,
                                 Join<JTracks, pidTPCEl, pidTPCMu> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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

        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/hadrons/hadrons/hHadronCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/hadrons/hHadronMultReco"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelations, "hadron-hadron correlation", true);

  /*********************************************** MC ************************************************/

  void processHadronCorrelationsMCGen(JMcCollision const&,
                                      JMcParticles const& tracks_true)
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
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);

        registry.fill(HIST("generated/hadrons/hHadronCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("generated/hadrons/hHadronMultGen"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelationsMCGen, "mc hadron-hadron correlation", true);

  void processHadronCorrelationsMCReco(Join<JCollisions, JMcCollisionLbs>::iterator const& collision_reco,
                                       JMcCollisions const&,
                                       Join<JTracks, JMcTrackLbs> const& tracks_reco,
                                       JMcParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelection)) {
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
      // if(particle.isPhysicalPrimary()){continue;}
      // if(std::abs(particle.pdgCode())>1e9){continue;}

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

        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);

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
  void processPionCorrelations(JCollision const& collision,
                               Join<JTracks, pidTPCPi> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
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

        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/charged/hPionCorrelReco"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/charged/hPionMultReco"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelations, "hadron-pion correlation", true);

  /*********************************************** MC ************************************************/

  void processPionCorrelationsMCGen(JMcCollision const&,
                                    JMcParticles const& tracks_true)
  {
    int nPions = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc)) {
        continue;
      }
      if (std::abs(track_assoc.pdgCode()) != 211) {
        continue;
      }

      registry.fill(HIST("generated/charged/hPionSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nPions++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);

        registry.fill(HIST("generated/charged/hPionCorrelGen"), track_trig.pt(), track_assoc.pt(), track_trig.eta() - track_assoc.eta(), dphi);
      }
    }
    registry.fill(HIST("generated/charged/hPionMultGen"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelationsMCGen, "mc hadron-pion correlation", true);

  void processPionCorrelationsMCReco(Join<JCollisions, JMcCollisionLbs>::iterator const& collision_reco,
                                     JMcCollisions const&,
                                     Join<JTracks, JMcTrackLbs> const& tracks_reco,
                                     JMcParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelection)) {
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
      if (std::abs(track_assoc.mcParticle().pdgCode()) != 211) {
        continue;
      }

      registry.fill(HIST("reconstructed/charged/hPionSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi());
      nPions++;

      for (const auto& track_trig : tracks_reco) {
        if (!jetderiveddatautilities::selectTrack(track_trig, trackSelection)) {
          continue;
        }

        if (!initTrig(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_trig.phi(), track_assoc.phi(), -PIHalf);
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

  void processNeutralCorrelationsMCGen(JMcCollision const&,
                                       JMcParticles const& tracks_true)
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
      if (track_assoc.pdgCode() < 100 || track_assoc.pdgCode() == 2112) {
        continue;
      } // remove non-hadrons and protons

      registry.fill(HIST("generated/neutral/hNeutralSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pidCodes[track_assoc.pdgCode()]);
      nNeutrals++;

      // Get correlations between neutral hadrons and their respective decay photons
      for (const auto& daughter : track_assoc.daughters_as<aod::JMcParticles>()) {
        if (daughter.pdgCode() != 22)
          continue;
        if (!initParticle(daughter, false))
          continue;
        if (!daughter.isPhysicalPrimary() && daughter.getGenStatusCode() == -1)
          continue;
        registry.fill(HIST("generated/neutral/hNeutralDecayGen"), track_assoc.pt(), daughter.pt() / track_assoc.pt(), daughter.eta() - track_assoc.eta(), calculateDelta(daughter.phi(), track_assoc.phi(), -PIHalf), pidCodes[track_assoc.pdgCode()]);
      }

      // Get correlations between triggers and neutral hadrons
      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        float dphi = calculateDelta(track_assoc.phi(), track_trig.phi(), -PIHalf);
        registry.fill(HIST("generated/neutral/hNeutralCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, pidCodes[track_assoc.pdgCode()]);
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
