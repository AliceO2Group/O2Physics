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

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "Framework/InitContext.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TClonesArray.h>
#include <TLorentzVector.h> // o2-linter: disable= root/lorentz-vector (TLorentzVector is needed for TPythia8Decayer)
#include <TPDGCode.h>
#include <TParticle.h>
#include <TPythia8Decayer.h>

#include <cmath>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

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
  Configurable<int> subGeneratorIdSelections{"subGeneratorIdSelections", -1, "set sub generator id"};
  Configurable<bool> hardQCDgg2gg{"hardQCDgg2gg", true, "include gg2gg process in hardQCD"};

  Configurable<int> tpcNClsCrossedRows{"tpcNClsCrossedRows", 70, "tpcNClsCrossedRows"};
  Configurable<double> tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0.8, "tpcCrossedRowsOverFindableCls"};
  Configurable<double> tpcNSigmaPi{"tpcNSigmaPi", 2., "tpcNSigmaPi"};

  Configurable<double> pTHatMin{"pTHatMin", 5, "minimum pTHat cut"};
  Configurable<double> pTHatMax{"pTHatMax", 600, "minimum pTHat cut"};

  const int pidCodeHadronCut = 100;
  const int pythiaCodeIncomingHard = 21;
  const int pythiaCodeOutgoingHard = 23;
  const int pythiaCodeOutgoingDiff = 15;

  Configurable<float> etaMaxTrig{"etaMaxTrig", 0.8, "maximum eta cut for triggers"};
  Configurable<float> etaMaxAssoc{"etaMaxAssoc", 0.8, "maximum eta cut for associateds"};
  Configurable<int> etaBinsTrig{"etaBinsTrig", 40, "number of eta bins for triggers"};
  Configurable<int> etaBinsAssoc{"etaBinsAssoc", 40, "number of eta bins for associateds"};
  Configurable<int> phiBins{"phiBins", 72, "number of phi bins"};

  // remove this comment
  AxisSpec axisEventStats = {6, -.5, 5.5, "Stats"};
  ConfigurableAxis axisPtHat = {"axisPtHat",
                                {VARIABLE_WIDTH, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                 1., 2., 3., 4., 5., 6., 7., 8., 9.,
                                 10., 20., 30., 40., 50., 60., 70., 80., 90.,
                                 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.},
                                "p_{T}^{#hat}"};
  ConfigurableAxis axisPtTrig = {"axisPtTrig",
                                 {VARIABLE_WIDTH,
                                  3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f, 20.0f},
                                 "p_{T, trig} [GeV]"}; // Axis for trigger particle pt distribution
  ConfigurableAxis axisPtAssoc = {"axisPtAssoc",
                                  {VARIABLE_WIDTH,
                                   0.2, 0.5, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 9.0f, 11.0f, 15.0f},
                                  "p_{T, assoc} [GeV]"};                        // Axis for associated particle pt distribution
  ConfigurableAxis axisDeltaPt = {"axisDeltaPt", {200, 0., 1.2}, "#Delta p_T"}; // Axis for pt ratio between neutral hadrons and decay photons
  int nHadrons = 6;
  AxisSpec axisPid = {10, -3.5, nHadrons + 0.5, "pid"};                               // Axis for PID of neutral hadrons
  ConfigurableAxis axisMult = {"axisMult", {100, 0., 99.}, "N_{ch}"};                 // Axis for mutplipicity
  AxisSpec axisAlpha = {100, 0., 1., "alpha"};                                        // Axis for decay photon pt assymetry
  ConfigurableAxis axisDeltaRDecay = {"axisDeltaRDecay", {400, 0., 3.2}, "#Delta R"}; // Axis for Delta R = sqrt(Delta eta^2 + Delta phi^2) between neutral hadrons and decay photons

  float ptMinTrig;
  float ptMaxTrig;
  float ptMinAssoc;
  float ptMaxAssoc;

  TPythia8Decayer* decayer = new TPythia8Decayer;
  TLorentzVector* motherLV = new TLorentzVector(); // o2-linter: disable= root/lorentz-vector (TLorentzVector is needed for TPythia8Decayer)
  TClonesArray* decayParticles = new TClonesArray("TParticle", 10);

  HistogramRegistry registry{"histogram registry"};

  // Particle ids for storing neutral hadrons
  std::map<std::string, int> pidCodes = {
    {"pi0", 1},    // pi0
    {"eta", 2},    // eta
    {"eta'", 3},   // eta'
    {"phi", 4},    // phi
    {"omega", 5},  // omega
    {"Sigma0", 6}, // Sigma
    {"Sigma0_bar", 6}};

  std::map<int, int> pidToPdg = {
    {1, 111},   // pi0
    {2, 221},   // eta
    {3, 331},   // eta'
    {4, 333},   // phi
    {5, 223},   // omega
    {6, 3212}}; // Sigma

  Service<o2::framework::O2DatabasePDG> pdg;

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    ptMinTrig = axisPtTrig->at(1);
    ptMaxTrig = axisPtTrig->back();
    ptMinAssoc = axisPtAssoc->at(1);
    ptMaxAssoc = axisPtAssoc->back();

    AxisSpec axisPhi = {phiBins, 0., TwoPI, "#phi"};                                                                            // Axis for phi distribution
    AxisSpec axisDeltaPhi = {phiBins, -PIHalf, 3 * PIHalf, "#Delta #phi"};                                                      // Axis for Delta phi in correlations
    AxisSpec axisEtaTrig = {etaBinsTrig, -etaMaxTrig, etaMaxTrig, "#eta"};                                                      // Axis for eta distribution
    AxisSpec axisEtaAssoc = {etaBinsAssoc, -etaMaxAssoc, etaMaxAssoc, "#eta"};                                                  // Axis for eta distribution
    AxisSpec axisDeltaEta = {etaBinsTrig + etaBinsAssoc, -(etaMaxTrig + etaMaxAssoc), etaMaxTrig + etaMaxAssoc, "#Delta #eta"}; // Axis for Delta eta in correlations

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    // Generated histograms
    registry.add("generated/events/hEventStats", "Event statistics", kTH1F, {axisEventStats});
    registry.get<TH1>(HIST("generated/events/hEventStats"))->Sumw2();
    registry.add("generated/events/hPtHat", "pT of hard collision", kTH1F, {axisPtHat});
    registry.get<TH1>(HIST("generated/events/hPtHat"))->Sumw2();

    //  Triggers
    registry.add("generated/triggers/hTrigMultGen", "Generated Trigger Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("generated/triggers/hTrigMultGen"))->Sumw2();
    registry.add("generated/triggers/hTrigSpectrumGen", "Generated Trigger Spectrum", kTHnSparseF, {axisPtTrig, axisEtaTrig, axisPhi});
    registry.get<THnSparse>(HIST("generated/triggers/hTrigSpectrumGen"))->Sumw2();

    // Hadrons
    registry.add("generated/hadrons/hHadronCorrelGen", "Generated Trigger-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.get<THnSparse>(HIST("generated/hadrons/hHadronCorrelGen"))->Sumw2();
    registry.add("generated/hadrons/hHadronMultGen", "Generated Hadron Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("generated/hadrons/hHadronMultGen"))->Sumw2();
    registry.add("generated/hadrons/hHadronSpectrumGen", "Generated Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi});
    registry.get<THnSparse>(HIST("generated/hadrons/hHadronSpectrumGen"))->Sumw2();

    // Photons
    registry.add("generated/photons/hPhotonCorrelGen", "Generated Trigger-Photon Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/photons/hPhotonCorrelGen"))->Sumw2();
    registry.add("generated/photons/hPhotonMultGen", "Generated Photon Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("generated/photons/hPhotonMultGen"))->Sumw2();
    registry.add("generated/photons/hPhotonSpectrumGen", "Generated Photon Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/photons/hPhotonSpectrumGen"))->Sumw2();

    // Charged pions
    registry.add("generated/charged/hPionCorrelGen", "Generated Trigger-Pion Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.get<THnSparse>(HIST("generated/charged/hPionCorrelGen"))->Sumw2();
    registry.add("generated/charged/hPionMultGen", "Generated Pion Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("generated/charged/hPionMultGen"))->Sumw2();
    registry.add("generated/charged/hPionSpectrumGen", "Generated Pion Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi});
    registry.get<THnSparse>(HIST("generated/charged/hPionSpectrumGen"))->Sumw2();
    registry.add("generated/charged/hCocktailPhotonCorrelGen", "Cocktail Photon from Pion-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/charged/hCocktailPhotonCorrelGen"))->Sumw2();
    registry.add("generated/charged/hCocktailPhotonSpectrumGen", "Cocktail Photon from Pion Spectrum", kTHnSparseF, {axisPtAssoc, axisPtAssoc, axisEtaAssoc, axisPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/charged/hCocktailPhotonSpectrumGen"))->Sumw2();

    ////Neutral particles
    registry.add("generated/neutral/hNeutralCorrelGen", "Generated Trigger-Neutral Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/neutral/hNeutralCorrelGen"))->Sumw2();
    registry.add("generated/neutral/hNeutralMultGen", "Generated Neutral Hadron Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("generated/neutral/hNeutralMultGen"))->Sumw2();
    registry.add("generated/neutral/hNeutralSpectrumGen", "Generated Neutral Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/neutral/hNeutralSpectrumGen"))->Sumw2();
    registry.add("generated/neutral/hNeutralDecayGen", "Generated Neutral Hadron-Decay Photon Correlation", kTHnSparseF, {axisPtAssoc, axisDeltaPt, axisDeltaRDecay, axisAlpha, axisPid}); // Correlation with decay photons
    registry.get<THnSparse>(HIST("generated/neutral/hNeutralDecayGen"))->Sumw2();
    registry.add("generated/neutral/hCocktailPhotonCorrelGen", "Cocktail Photon from Hadron-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisPtAssoc, axisDeltaEta, axisDeltaPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/neutral/hCocktailPhotonCorrelGen"))->Sumw2();
    registry.add("generated/neutral/hCocktailPhotonSpectrumGen", "Cocktail Photon from Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisPtAssoc, axisEtaAssoc, axisPhi, axisPid});
    registry.get<THnSparse>(HIST("generated/neutral/hCocktailPhotonSpectrumGen"))->Sumw2();

    // Reconstructed histograms
    registry.add("reconstructed/events/hEventStats", "Event statistics", kTH1F, {axisEventStats});
    registry.get<TH1>(HIST("reconstructed/events/hEventStats"))->Sumw2();

    // Triggers
    registry.add("reconstructed/triggers/hTrigMultReco", "Reconstructed Trigger Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("reconstructed/triggers/hTrigMultReco"))->Sumw2();
    registry.add("reconstructed/triggers/hTrigSpectrumReco", "Reconstructed Trigger Spectrum", kTHnSparseF, {axisPtTrig, axisEtaTrig, axisPhi});
    registry.get<THnSparse>(HIST("reconstructed/triggers/hTrigSpectrumReco"))->Sumw2();

    // Hadrons
    registry.add("reconstructed/hadrons/hHadronCorrelReco", "Reconstructed Trigger-Hadron Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.get<THnSparse>(HIST("reconstructed/hadrons/hHadronCorrelReco"))->Sumw2();
    registry.add("reconstructed/hadrons/hHadronMultReco", "Reconstructed Hadron Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("reconstructed/hadrons/hHadronMultReco"))->Sumw2();
    registry.add("reconstructed/hadrons/hHadronSpectrumReco", "Reconstructed Hadron Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi});
    registry.get<THnSparse>(HIST("reconstructed/hadrons/hHadronSpectrumReco"))->Sumw2();
    registry.add("reconstructed/hadrons/hHadronPtPrimReco", "Reconstructed Primaries Spectrum", kTH1F, {axisPtAssoc}); // Primary hadron spectrum
    registry.get<TH1>(HIST("reconstructed/hadrons/hHadronPtPrimReco"))->Sumw2();
    registry.add("reconstructed/hadrons/hHadronPtSecReco", "Reconstructed Secondaries Spectrum", kTH1F, {axisPtAssoc}); // Secondary hadron spectrum
    registry.get<TH1>(HIST("reconstructed/hadrons/hHadronPtSecReco"))->Sumw2();

    // Photons
    registry.add("reconstructed/photons/hPhotonCorrelReco", "Reconstructed Trigger-Photon Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.get<THnSparse>(HIST("reconstructed/photons/hPhotonCorrelReco"))->Sumw2();
    registry.add("reconstructed/photons/hPhotonMultReco", "Reconstructed Photon Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("reconstructed/photons/hPhotonMultReco"))->Sumw2();
    registry.add("reconstructed/photons/hPhotonSpectrumReco", "Reconstructed Photon Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi});
    registry.get<THnSparse>(HIST("reconstructed/photons/hPhotonSpectrumReco"))->Sumw2();

    // Charged Pions
    registry.add("reconstructed/charged/hPionCorrelReco", "Reconstructed Trigger-Pion Correlation", kTHnSparseF, {axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi});
    registry.get<THnSparse>(HIST("reconstructed/charged/hPionCorrelReco"))->Sumw2();
    registry.add("reconstructed/charged/hPionMultReco", "Reconstructed Pion Multiplicity", kTH1F, {axisMult});
    registry.get<TH1>(HIST("reconstructed/charged/hPionMultReco"))->Sumw2();
    registry.add("reconstructed/charged/hPionSpectrumReco", "Reconstructed Pion Spectrum", kTHnSparseF, {axisPtAssoc, axisEtaAssoc, axisPhi});
    registry.get<THnSparse>(HIST("reconstructed/charged/hPionSpectrumReco"))->Sumw2();
  }

  // Get pTHat of the event
  template <typename P>
  double getCollisionPtHat(const P& particles)
  {
    double pTHat = 0;
    bool isDiffractive = false;

    for (const auto& particle : particles) {

      if (std::abs(particle.getGenStatusCode()) == pythiaCodeOutgoingHard && !isDiffractive) {
        if (particle.pt() > pTHat) {
          pTHat = particle.pt();
        }
      }
      if (std::abs(particle.getGenStatusCode()) == pythiaCodeOutgoingDiff) {
        if (isDiffractive) {
          if (particle.pt() > pTHat) {
            pTHat = particle.pt();
          }
        } else {
          pTHat = particle.pt();
          isDiffractive = true;
        }
      }
    }

    return pTHat;
  }

  template <typename P>
  bool rejectOutliers(const P& particles, double pTHat)
  {
    if (subGeneratorIdSelections == 0) {
      return false;
    }

    for (const auto& particle : particles) {
      if (particle.getGenStatusCode() > 0 && particle.pt() > pTHat) {
        return true;
      }
    }

    return false;
  }

  // If event is a HardQCD gg->gg event (which is not implemented in POWHEG directphoton)
  template <typename P>
  bool isHardQCDgg2gg(const P& particles)
  {
    for (const auto& particle : particles) {
      if (std::abs(particle.getGenStatusCode()) == pythiaCodeIncomingHard && (PDG_t)std::abs(particle.pdgCode()) != kGluon) {
        return false;
      }
      if (std::abs(particle.getGenStatusCode()) == pythiaCodeOutgoingHard && (PDG_t)std::abs(particle.pdgCode()) != kGluon) {
        return false;
      }
    }
    return true;
  }

  // Initialize collision
  template <typename C, typename P>
  bool initCollisionMC(const C& collision, const P& particles)
  {
    if (collision.getSubGeneratorId() != subGeneratorIdSelections) {
      return false;
    }

    double pTHat = getCollisionPtHat(particles);
    if (pTHat < pTHatMin || pTHat > pTHatMax) {
      return false;
    }

    if (rejectOutliers(particles, pTHat)) {
      return false;
    }

    if (!hardQCDgg2gg && isHardQCDgg2gg(particles)) {
      return false;
    }

    return true;
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

    if (std::abs(track.eta()) > etaMaxAssoc) {
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

    if (particle.getGenStatusCode() == -1) {
      return false;
    }

    if (std::abs(particle.eta()) > etaMaxAssoc) {
      return false;
    }

    if (particle.pt() < ptMinAssoc || particle.pt() > ptMaxAssoc) {
      return false;
    }

    return true;
  }

  // Initialize Pythia8 decay particle
  bool initDecayParticle(const TParticle* particle)
  {

    if (particle->GetMother(0) != 0) {
      return false;
    }

    if ((PDG_t)std::abs(particle->GetPdgCode()) != kGamma) {
      return false;
    }

    if (std::abs(particle->Eta()) > etaMaxAssoc) {
      return false;
    }

    if (particle->Pt() < ptMinAssoc || particle->Pt() > ptMaxAssoc) {
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

    if (std::abs(track.eta()) > etaMaxTrig) {
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

    if (std::abs(particle.eta()) > etaMaxTrig) {
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

    if (track.tpcNClsCrossedRows() < tpcNClsCrossedRows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < tpcCrossedRowsOverFindableCls) {
      return false;
    }

    return true;
  }

  /****************************************************************************************************
  ************************************************ EVENTS ********************************************
  ****************************************************************************************************/
  void processEventsMCGen(JetMcCollision const& collision,
                          JetParticles const& particles)
  {
    if (collision.getSubGeneratorId() == subGeneratorIdSelections) {
      registry.fill(HIST("generated/events/hEventStats"), 0);
      registry.fill(HIST("generated/events/hEventStats"), 2, collision.weight());
      registry.fill(HIST("generated/events/hEventStats"), 4, collision.xsectGen());
    }

    if (!initCollisionMC(collision, particles)) {
      return;
    }

    registry.fill(HIST("generated/events/hEventStats"), 1);
    registry.fill(HIST("generated/events/hEventStats"), 3, collision.weight());
    registry.fill(HIST("generated/events/hEventStats"), 5, collision.xsectGen());

    registry.fill(HIST("generated/events/hPtHat"), getCollisionPtHat(particles));
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processEventsMCGen, "event stats MC gen", true);

  void processEventsMCReco(JetCollisionMCD const& collision,
                           JetMcCollisions const&,
                           JetParticles const&)
  {
    if (collision.getSubGeneratorId() != subGeneratorIdSelections) {
      return;
    }
    registry.fill(HIST("reconstructed/events/hEventStats"), 0);
    registry.fill(HIST("reconstructed/events/hEventStats"), 2, collision.mcCollision().weight());
    registry.fill(HIST("reconstructed/events/hEventStats"), 4, collision.mcCollision().xsectGen());

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false)) {
      return;
    }

    registry.fill(HIST("reconstructed/events/hEventStats"), 1);
    registry.fill(HIST("reconstructed/events/hEventStats"), 3, collision.mcCollision().weight());
    registry.fill(HIST("reconstructed/events/hEventStats"), 5, collision.mcCollision().xsectGen());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processEventsMCReco, "event stats MC reco", true);

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

  void processTrigsMCReco(JetCollisionMCD const& collision,
                          JetMcCollisions const&,
                          Join<JetTracks, JMcTrackLbs> const& tracks,
                          JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, false)) {
      return;
    }
    if (collision.getSubGeneratorId() != subGeneratorIdSelections) {
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
      registry.fill(HIST("reconstructed/triggers/hTrigSpectrumReco"), track.pt(), track.eta(), track.phi(), collision.mcCollision().weight());
      nTrigs++;
    }
    registry.fill(HIST("reconstructed/triggers/hTrigMultReco"), nTrigs, collision.mcCollision().weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processTrigsMCReco, "trigger particle mc properties", true);

  void processTrigsMCGen(JetMcCollision const& collision,
                         JetParticles const& particles)
  {
    if (!initCollisionMC(collision, particles)) {
      return;
    }

    int nTrigs = 0;
    for (const auto& particle : particles) {
      if (!initTrigParticle(particle)) {
        continue;
      }
      registry.fill(HIST("generated/triggers/hTrigSpectrumGen"), particle.pt(), particle.eta(), particle.phi(), collision.weight());
      nTrigs++;
    }
    registry.fill(HIST("generated/triggers/hTrigMultGen"), nTrigs, collision.weight());
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
        float dphi = RecoDecay::constrainAngle(v0.phi() - track.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), v0.eta() - track.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelations, "hadron-photon correlation", true);

  /*********************************************** MC ************************************************/

  using MyTracksMC = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  void processPhotonCorrelationsMCReco(Join<JetCollisionsMCD, JCollisionPIs>::iterator const& collision_reco,
                                       JetMcCollisions const&,
                                       JetTracks const& tracks_reco,
                                       JetParticles const&,
                                       MyTracksMC const&,
                                       V0Datas const& v0s)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits, false)) {
      return;
    }
    if (collision_reco.mcCollision().getSubGeneratorId() != subGeneratorIdSelections) {
      return;
    }

    int nPhotons = 0;

    auto v0PerCollision = v0s.sliceBy(perCol, collision_reco.globalIndex());
    for (const auto& v0 : v0PerCollision) {
      if (!initV0<MyTracksMC>(v0)) {
        continue;
      }
      registry.fill(HIST("reconstructed/photons/hPhotonSpectrumReco"), v0.pt(), v0.eta(), v0.phi(), collision_reco.mcCollision().weight());
      nPhotons++;

      for (const auto& track : tracks_reco) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
          continue;
        }

        if (!initTrig(track)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(v0.phi() - track.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/photons/hPhotonCorrelReco"), track.pt(), v0.pt(), v0.eta() - track.eta(), dphi, collision_reco.mcCollision().weight());
      }
    }
    registry.fill(HIST("reconstructed/photons/hPhotonMultReco"), nPhotons, collision_reco.mcCollision().weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonCorrelationsMCReco, "hadron-photon correlation", true);

  void processPhotonCorrelationsMCGen(JetMcCollision const& collision,
                                      JetParticles const& tracks_true)
  {
    if (!initCollisionMC(collision, tracks_true)) {
      return;
    }

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
      if (!pdgMother) {
        continue;
      }
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

      registry.fill(HIST("generated/photons/hPhotonSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), photonGeneration, collision.weight());

      nPhotons++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("generated/photons/hPhotonCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, photonGeneration, collision.weight());
      }
    }

    registry.fill(HIST("generated/photons/hPhotonMultGen"), nPhotons, collision.weight());
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
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/hadrons/hHadronCorrelReco"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/hadrons/hHadronMultReco"), nHadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelations, "hadron-hadron correlation", true);

  /*********************************************** MC ************************************************/

  void processHadronCorrelationsMCGen(JetMcCollision const& collision,
                                      JetParticles const& tracks_true)
  {
    if (!initCollisionMC(collision, tracks_true)) {
      return;
    }

    int nHadrons = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc)) {
        continue;
      }
      auto pdgParticle = pdg->GetParticle(track_assoc.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.) {
        continue;
      }
      if (std::abs(track_assoc.pdgCode()) < pidCodeHadronCut) {
        continue;
      }

      registry.fill(HIST("generated/hadrons/hHadronSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), collision.weight());
      nHadrons++;

      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);

        registry.fill(HIST("generated/hadrons/hHadronCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, collision.weight());
      }
    }
    registry.fill(HIST("generated/hadrons/hHadronMultGen"), nHadrons, collision.weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processHadronCorrelationsMCGen, "mc hadron-hadron correlation", true);

  void processHadronCorrelationsMCReco(JetCollisionMCD const& collision_reco,
                                       JetMcCollisions const&,
                                       JetTracksMCD const& tracks_reco,
                                       JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits, false)) {
      return;
    }
    if (collision_reco.mcCollision().getSubGeneratorId() != subGeneratorIdSelections) {
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
      if (std::abs(particle.pdgCode()) < pidCodeHadronCut) {
        continue;
      }

      registry.fill(HIST("reconstructed/hadrons/hHadronSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), collision_reco.mcCollision().weight());

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("reconstructed/hadrons/hHadronPtPrimReco"), track_assoc.pt(), collision_reco.mcCollision().weight());
      } else {
        registry.fill(HIST("reconstructed/hadrons/hHadronPtSecReco"), track_assoc.pt(), collision_reco.mcCollision().weight());
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
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);

        registry.fill(HIST("reconstructed/hadrons/hHadronCorrelReco"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, collision_reco.mcCollision().weight());
      }
    }
    registry.fill(HIST("reconstructed/hadrons/hHadronMultReco"), nHadrons, collision_reco.mcCollision().weight());
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
      if (std::abs(track_assoc.tpcNSigmaPi()) > tpcNSigmaPi) {
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
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/charged/hPionCorrelReco"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi);
      }
    }
    registry.fill(HIST("reconstructed/charged/hPionMultReco"), nPions);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelations, "hadron-pion correlation", true);

  /*********************************************** MC ************************************************/

  void processPionCorrelationsMCGen(JetMcCollision const& collision,
                                    JetParticles const& tracks_true)
  {
    if (!initCollisionMC(collision, tracks_true)) {
      return;
    }

    int nPions = 0;
    for (const auto& track_assoc : tracks_true) {
      if (!initParticle(track_assoc)) {
        continue;
      }
      if ((PDG_t)std::abs(track_assoc.pdgCode()) != kPiPlus) {
        continue;
      }

      registry.fill(HIST("generated/charged/hPionSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), collision.weight());
      nPions++;

      // Get correlations
      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        if (track_trig == track_assoc) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);

        registry.fill(HIST("generated/charged/hPionCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, collision.weight());
      }

      // Use PYTHIA to simulate decay
      decayer->Init();
      for (int pid = 1; pid <= nHadrons; pid++) {
        TParticlePDG* pdgParticle = nullptr;

        auto it = pidToPdg.find(pid);
        if (it != pidToPdg.end()) {
          pdgParticle = pdg->GetParticle(it->second);
        } else {
          continue;
        }

        motherLV->SetPtEtaPhiM(track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pdgParticle->Mass());
        decayer->Decay(pdgParticle->PdgCode(), motherLV);
        decayer->ImportParticles(decayParticles);
        for (int i = 0; i < decayParticles->GetEntriesFast(); ++i) {
          TParticle* daughter = static_cast<TParticle*>(decayParticles->At(i));

          if (!initDecayParticle(daughter)) {
            continue;
          }

          registry.fill(HIST("generated/charged/hCocktailPhotonSpectrumGen"), track_assoc.pt(), daughter->Pt(), daughter->Eta(), daughter->Phi(), pidCodes[pdgParticle->GetName()], collision.weight());

          for (const auto& track_trig : tracks_true) {
            if (!initTrigParticle(track_trig)) {
              continue;
            }
            if (track_trig == track_assoc) {
              continue;
            }
            float dphi = RecoDecay::constrainAngle(daughter->Phi() - track_trig.phi(), -PIHalf);
            registry.fill(HIST("generated/charged/hCocktailPhotonCorrelGen"), track_trig.pt(), track_assoc.pt(), daughter->Pt(), daughter->Eta() - track_trig.eta(), dphi, pidCodes[pdgParticle->GetName()], collision.weight());
          }
        }
      }
    }
    registry.fill(HIST("generated/charged/hPionMultGen"), nPions, collision.weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelationsMCGen, "mc hadron-pion correlation", true);

  void processPionCorrelationsMCReco(JetCollisionMCD const& collision_reco,
                                     JetMcCollisions const&,
                                     JetTracksMCD const& tracks_reco,
                                     JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision_reco, eventSelectionBits, false)) {
      return;
    }
    if (collision_reco.mcCollision().getSubGeneratorId() != subGeneratorIdSelections) {
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

      registry.fill(HIST("reconstructed/charged/hPionSpectrumReco"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), collision_reco.mcCollision().weight());
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
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("reconstructed/charged/hPionCorrelReco"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, collision_reco.mcCollision().weight());
      }
    }
    registry.fill(HIST("reconstructed/charged/hPionMultReco"), nPions, collision_reco.mcCollision().weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPionCorrelationsMCReco, "mc hadron-pion correlation", true);

  /****************************************************************************************************
  ****************************************** NEUTRALS *************************************************
  ****************************************************************************************************/

  /*********************************************** MC ************************************************/

  void processNeutralCorrelationsMCGen(JetMcCollision const& collision,
                                       JetParticles const& tracks_true)
  {
    if (!initCollisionMC(collision, tracks_true)) {
      return;
    }

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
      if (std::abs(track_assoc.pdgCode()) < pidCodeHadronCut || (PDG_t)track_assoc.pdgCode() == kNeutron) {
        continue;
      } // remove non-hadrons and neutrons

      registry.fill(HIST("generated/neutral/hNeutralSpectrumGen"), track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pidCodes[pdgParticle->GetName()], collision.weight());
      nNeutrals++;

      // Get correlations between neutral hadrons and their respective decay photons
      auto daughters = track_assoc.daughters_as<aod::JMcParticles>();
      double alpha = -1;
      int nPhotonsPionDecay = 2;
      if (daughters.size() == nPhotonsPionDecay) {
        auto daughter = daughters.begin();
        double pt1 = daughter.pt();
        ++daughter;
        double pt2 = daughter.pt();
        alpha = std::abs((pt1 - pt2) / (pt1 + pt2));
      }

      for (const auto& daughter : daughters) {
        if ((PDG_t)std::abs(daughter.pdgCode()) != kGamma) {
          continue;
        }
        if (!daughter.isPhysicalPrimary() && daughter.getGenStatusCode() == -1) {
          continue;
        }
        double deltaPt = daughter.pt() / track_assoc.pt();
        double deltaEta = daughter.eta() - track_assoc.eta();
        double deltaPhi = RecoDecay::constrainAngle(daughter.phi() - track_assoc.phi(), -PIHalf);
        double deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

        registry.fill(HIST("generated/neutral/hNeutralDecayGen"), track_assoc.pt(), deltaPt, deltaR, alpha, pidCodes[pdgParticle->GetName()], collision.weight());
      }

      // Get correlations between triggers and neutral hadrons
      for (const auto& track_trig : tracks_true) {
        if (!initTrigParticle(track_trig)) {
          continue;
        }
        float dphi = RecoDecay::constrainAngle(track_assoc.phi() - track_trig.phi(), -PIHalf);
        registry.fill(HIST("generated/neutral/hNeutralCorrelGen"), track_trig.pt(), track_assoc.pt(), track_assoc.eta() - track_trig.eta(), dphi, pidCodes[pdgParticle->GetName()], collision.weight());
      }

      // Use PYTHIA to simulate decay
      decayer->Init();
      motherLV->SetPtEtaPhiM(track_assoc.pt(), track_assoc.eta(), track_assoc.phi(), pdgParticle->Mass());
      decayer->Decay(track_assoc.pdgCode(), motherLV);
      decayer->ImportParticles(decayParticles);
      for (int i = 0; i < decayParticles->GetEntriesFast(); ++i) {
        TParticle* daughter = static_cast<TParticle*>(decayParticles->At(i));

        if (!initDecayParticle(daughter)) {
          continue;
        }

        registry.fill(HIST("generated/neutral/hCocktailPhotonSpectrumGen"), track_assoc.pt(), daughter->Pt(), daughter->Eta(), daughter->Phi(), pidCodes[pdgParticle->GetName()], collision.weight());

        for (const auto& track_trig : tracks_true) {
          if (!initTrigParticle(track_trig)) {
            continue;
          }
          float dphi = RecoDecay::constrainAngle(daughter->Phi() - track_trig.phi(), -PIHalf);
          registry.fill(HIST("generated/neutral/hCocktailPhotonCorrelGen"), track_trig.pt(), track_assoc.pt(), daughter->Pt(), daughter->Eta() - track_trig.eta(), dphi, pidCodes[pdgParticle->GetName()], collision.weight());
        }
      }
    }
    registry.fill(HIST("generated/neutral/hNeutralMultGen"), nNeutrals, collision.weight());
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processNeutralCorrelationsMCGen, "mc hadron-pion correlation", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HadronPhotonCorrelation>(cfgc)};
};
