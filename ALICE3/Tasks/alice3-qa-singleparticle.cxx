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
/// \file    alice3-qa-singleparticle.cxx
/// \author  Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief   Task to monitor the single particle QA, at the particle and track level, showing the tracked and the origin of particles
///

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TMCProcess.h>
#include <TString.h>

#include <RtypesCore.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3SingleParticle {
  Service<o2::framework::O2DatabasePDG> pdg{};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> pdgCode{"pdgCode", 2212, "PDG code of the particle of interest"};

  struct : ConfigurableGroup {
    std::string prefix = "settings";
    Configurable<int> isStable{"isStable", 0, "Flag to check stable particles"};
    Configurable<float> charge{"charge", 1.f, "Particle charge to scale the reconstructed momentum"};
    Configurable<float> nchEtaMax{"nchEtaMax", 0.5f, "|eta| acceptance of the Nch estimator (use the Ditto tune mActivityEtaMax)"};
    Configurable<bool> doPrint{"doPrint", false, "Flag to print debug messages"};
  } settings;

  // The pT, eta, y and production-vertex-z axes also define the selection of the particles of interest
  struct : ConfigurableGroup {
    std::string prefix = "settingsAxes";
    ConfigurableAxis axisMomentum{"axisMomentum", {500, 0.f, 10.f}, "Binning for p and pT"};
    ConfigurableAxis axisMomentumComponent{"axisMomentumComponent", {500, -10.f, 10.f}, "Binning for px, py and pz"};
    ConfigurableAxis axisEnergy{"axisEnergy", {500, 0.f, 1000.f}, "Binning for the energy"};
    ConfigurableAxis axisEta{"axisEta", {500, -3.f, 3.f}, "Binning for eta"};
    ConfigurableAxis axisY{"axisY", {200, -3.f, 3.f}, "Binning for rapidity"};
    ConfigurableAxis axisVtxXY{"axisVtxXY", {100, -1.f, 1.f}, "Binning for the event vertex in X and Y"};
    ConfigurableAxis axisVtxZ{"axisVtxZ", {100, -20.f, 20.f}, "Binning for the event vertex in Z"};
    ConfigurableAxis axisProdXY{"axisProdXY", {100, -1.f, 1.f}, "Binning for the production vertex in X and Y"};
    ConfigurableAxis axisProdZ{"axisProdZ", {100, -10.f, 10.f}, "Binning for the production vertex in Z"};
    ConfigurableAxis axisProdRadius{"axisProdRadius", {100, 0.f, 2.f}, "Binning for the production vertex radius"};
    ConfigurableAxis axisNch{"axisNch", {150, 0.f, 150.f}, "Charged-particle multiplicity"};
    ConfigurableAxis axisNSignal{"axisNSignal", {20, -0.5f, 19.5f}, "Number of selected particles of interest per collision"};
    ConfigurableAxis distanceAxis{"distanceAxis", {1000, 0, 1e-10}, "Distance (cm)"};
  } settingsAxes;

  std::vector<int64_t> particlesOfInterest; // global indices of the MC particles selected in genericProcessParticles, per DF

  // Selection range [min, max) of the particles of interest, taken from an axis so that every selected particle is inside its histograms
  struct Range {
    Range() = default;
    Range(double min_, double max_, int64_t bins_) : min(min_), max(max_), bins(bins_) {}
    double min = 0.;
    double max = 0.;
    int64_t bins = 0;
    int getBins() const { return static_cast<int>(bins); }
    [[nodiscard]] bool contains(double value) const { return value >= min && value < max; }
  };
  Range ptRange, etaRange, yRange, prodXYrange, prodZRange;

  static Range getRange(const ConfigurableAxis& axis)
  {
    const AxisSpec spec{axis};
    return Range(spec.binEdges.front(), spec.binEdges.back(), spec.getNbins());
  }

  void init(InitContext&)
  {
    int nEnabled = 0;
    if (doprocessStandard) {
      nEnabled++;
    }
    if (doprocessParticleOnly) {
      nEnabled++;
    }
    if (doprocessNonIU) {
      nEnabled++;
    }
    if (nEnabled == 0 || nEnabled > 1) {
      LOG(fatal) << "You can't process with " << nEnabled << " process functions, pick one";
    }

    pdg->AddParticle("deuteron", "deuteron", 1.8756134, kTRUE, 0.0, 3, "Nucleus", 1000010020);
    pdg->AddAntiParticle("anti-deuteron", -1000010020);

    pdg->AddParticle("triton", "triton", 2.8089218, kTRUE, 0.0, 3, "Nucleus", 1000010030);
    pdg->AddAntiParticle("anti-triton", -1000010030);

    pdg->AddParticle("helium3", "helium3", o2::constants::physics::MassHelium3, kTRUE, 0.0, 6, "Nucleus", 1000020030);
    pdg->AddAntiParticle("anti-helium3", -1000020030);

    pdg->AddParticle("helium4", "helium4", o2::constants::physics::MassAlpha, kTRUE, 0.0, 6, "Nucleus", 1000020040);
    pdg->AddAntiParticle("anti-helium4", -1000020040);

    ptRange = getRange(settingsAxes.axisMomentum);
    etaRange = getRange(settingsAxes.axisEta);
    yRange = getRange(settingsAxes.axisY);
    prodXYrange = getRange(settingsAxes.axisProdXY);
    prodZRange = getRange(settingsAxes.axisProdZ);
    LOG(info) << "Selecting particles of interest in pT [" << ptRange.min << ", " << ptRange.max << "), eta [" << etaRange.min << ", " << etaRange.max
              << "), y [" << yRange.min << ", " << yRange.max << "), production z [" << prodZRange.min << ", " << prodZRange.max << ")";

    const TString titlePdgCode = Form("%i", pdgCode.value);
    const AxisSpec axisPDGs{100, 0.f, 100.f, "PDG code"};
    const AxisSpec axisCharge{21, -10.5f, 10.5f, "Charge (e/3)"};
    const AxisSpec axisP{settingsAxes.axisMomentum, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt{settingsAxes.axisMomentum, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPtMC{settingsAxes.axisMomentum, "#it{p}_{T}^{MC} (GeV/#it{c})"};
    const AxisSpec axisPx{settingsAxes.axisMomentumComponent, "#it{p}_{x} (GeV/#it{c})"};
    const AxisSpec axisPy{settingsAxes.axisMomentumComponent, "#it{p}_{y} (GeV/#it{c})"};
    const AxisSpec axisPz{settingsAxes.axisMomentumComponent, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec axisEta{settingsAxes.axisEta, "#it{#eta}"};
    const AxisSpec axisY{settingsAxes.axisY, "#it{y}"};
    const AxisSpec axisE{etaRange.getBins(), 0, 1000, "E"};
    const AxisSpec axisProdx{settingsAxes.axisProdXY, "Prod. Vertex X (cm)"};
    const AxisSpec axisPrody{settingsAxes.axisProdXY, "Prod. Vertex Y (cm)"};
    const AxisSpec axisProdz{settingsAxes.axisProdZ, "Prod. Vertex Z (cm)"};
    const AxisSpec axisProdRadius{prodXYrange.getBins(), 0., 2. * prodXYrange.max, "Prod. Vertex Radius (cm)"};
    const AxisSpec axisDistx{settingsAxes.distanceAxis, "Distance X (fm)"};
    const AxisSpec axisDisty{settingsAxes.distanceAxis, "Distance Y (fm)"};
    const AxisSpec axisDistz{settingsAxes.distanceAxis, "Distance Z (fm)"};
    const AxisSpec axisProdDistance{settingsAxes.distanceAxis, "Distance Radius (fm)"};
    const AxisSpec axisVx{settingsAxes.axisVtxXY, "Vertex X (cm)"};
    const AxisSpec axisVy{settingsAxes.axisVtxXY, "Vertex Y (cm)"};
    const AxisSpec axisVz{settingsAxes.axisVtxZ, "Vertex Z (cm)"};

    histos.add("AnalysedParticlePdgCode", "Analysed particle PDG code", kTH1D, {{1, -0.5, 0.5}});
    histos.get<TH1>(HIST("AnalysedParticlePdgCode"))->SetBinContent(1, pdgCode);

    const auto* particle = pdg->GetParticle(pdgCode);
    histos.add("AnalysedParticleMass", "Analysed particle mass", kTH1D, {{1, -0.5, 0.5}});
    histos.get<TH1>(HIST("AnalysedParticleMass"))->SetBinContent(1, particle ? particle->Mass() : -1.);

    histos.add("event/VtxX", "Vertex X", kTH1D, {axisVx});
    histos.add("event/VtxY", "Vertex Y", kTH1D, {axisVy});
    histos.add("event/VtxZ", "Vertex Z", kTH1D, {axisVz});
    histos.add("event/Nch", "Events;N_{ch};Events", kTH1D, {settingsAxes.axisNch});
    histos.add("particle/NSignalVsNch", "Particles of interest per event " + titlePdgCode + ";N_{ch};N_{signal}", kTH2D, {settingsAxes.axisNch, settingsAxes.axisNSignal});

    histos.add("particle/PDGs", "Particle PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/PDGsPrimaries", "Particle PDGs of Primaries", kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/PDGsSecondaries", "Particle PDGs of Secondaries", kTH2D, {axisPDGs, axisCharge});
    auto h = histos.add<TH1>("particle/producedProcess", "Particle process of production " + titlePdgCode, kTH1D, {{kMaxMCProcess, -0.5, -0.5 + kMaxMCProcess, "getProcess"}});
    for (int i = 0; i < kMaxMCProcess; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, TMCProcessName[i]);
    }
    histos.add("particle/Pt", "Particle Pt " + titlePdgCode, kTH1D, {axisPt});
    histos.add("particle/P", "Particle P " + titlePdgCode, kTH1D, {axisP});
    histos.add("particle/primariesPt", "Particle Pt (primary) " + titlePdgCode, kTH1D, {axisPt});
    histos.add("particle/primariesP", "Particle P (primary) " + titlePdgCode, kTH1D, {axisP});
    histos.add("particle/secondariesPt", "Particle Pt (secondary) " + titlePdgCode, kTH1D, {axisPt});
    histos.add("particle/secondariesP", "Particle P (secondary) " + titlePdgCode, kTH1D, {axisP});
    histos.add("particle/prodVx", "Particle Prod. Vertex X " + titlePdgCode, kTH1D, {axisProdx});
    histos.add("particle/prodVy", "Particle Prod. Vertex Y " + titlePdgCode, kTH1D, {axisPrody});
    histos.add("particle/prodVz", "Particle Prod. Vertex Z " + titlePdgCode, kTH1D, {axisProdz});
    histos.add("particle/prodRadius", "Particle Prod. Vertex Radius XY " + titlePdgCode, kTH1D, {axisProdRadius});
    histos.add("particle/prodRadius3D", "Particle Prod. Vertex Radius XYZ " + titlePdgCode, kTH1D, {axisProdRadius});
    histos.add("particle/prodDistance", "Particle Prod. Vertex Distance XY " + titlePdgCode, kTH1D, {axisProdDistance});
    histos.add("particle/prodDistance3D", "Particle Prod. Vertex Distance XYZ " + titlePdgCode, kTH1D, {axisProdDistance});
    histos.add("particle/prodVxVy", "Particle Prod. Vertex X vs Vertex Y " + titlePdgCode, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyPrm", "Particle Prod. Vertex X vs Vertex Y Primary " + titlePdgCode, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyStr", "Particle Prod. Vertex X vs Vertex Y Decay " + titlePdgCode, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVyMat", "Particle Prod. Vertex X vs Vertex Y Material " + titlePdgCode, kTH2D, {axisProdx, axisPrody});
    histos.add("particle/prodVxVsPt", "Particle Prod. Vertex X " + titlePdgCode, kTH2D, {axisPt, axisProdx});
    histos.add("particle/prodVyVsPt", "Particle Prod. Vertex Y " + titlePdgCode, kTH2D, {axisPt, axisPrody});
    histos.add("particle/prodVzVsPt", "Particle Prod. Vertex Z " + titlePdgCode, kTH2D, {axisPt, axisProdz});
    histos.add("particle/prodRadiusVsPt", "Particle Prod. Vertex Radius XY " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/prodRadiusVsEta", "Particle Prod. Vertex Radius XY " + titlePdgCode, kTH2D, {axisEta, axisProdRadius});
    histos.add("particle/prodRadius3DVsPt", "Particle Prod. Vertex Radius XYZ " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/Eta", "Particle Eta " + titlePdgCode, kTH1D, {axisEta});
    histos.add("particle/primariesEta", "Particle Eta (primary) " + titlePdgCode, kTH1D, {axisEta});
    histos.add("particle/secondariesEta", "Particle Eta (secondary) " + titlePdgCode, kTH1D, {axisEta});
    histos.add("particle/Y", "Particle Y " + titlePdgCode, kTH1D, {axisY});
    histos.add("particle/primariesY", "Particle Y (primary)" + titlePdgCode, kTH1D, {axisY});
    histos.add("particle/secondariesY", "Particle Y (secondary)" + titlePdgCode, kTH1D, {axisY});
    histos.add("particle/EvsPz", "Particle E vs Pz " + titlePdgCode, kTH2D, {axisE, axisPz});
    histos.add("particle/YvzPz", "Particle Y vs Pz " + titlePdgCode, kTH2D, {axisY, axisPz});
    histos.add("particle/EtavzPz", "Particle Eta vs Pz " + titlePdgCode, kTH2D, {axisEta, axisPz});
    histos.add("particle/PtvzPz", "Particle Pt vs Pz " + titlePdgCode, kTH2D, {axisPt, axisPz});
    histos.add("particle/PvzPz", "Particle P vs Pz " + titlePdgCode, kTH2D, {axisP, axisPz});
    histos.add("particle/Px", "Particle Px " + titlePdgCode, kTH1D, {axisPx});
    histos.add("particle/Py", "Particle Py " + titlePdgCode, kTH1D, {axisPy});
    histos.add("particle/Pz", "Particle Pz " + titlePdgCode, kTH1D, {axisPz});
    histos.add("particle/decayDistanceX", "Particle decay distance " + titlePdgCode, kTH1D, {axisProdx});
    histos.add("particle/decayDistanceY", "Particle decay distance " + titlePdgCode, kTH1D, {axisPrody});
    histos.add("particle/decayDistanceZ", "Particle decay distance " + titlePdgCode, kTH1D, {axisProdz});
    histos.add("particle/decayDistanceVsPt", "Particle decay distance " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/decayCtau", "Particle decay distance " + titlePdgCode, kTH1D, {{1000, 0, 10000, "c#tau (#mum)"}});

    histos.add("particle/daughters/Number", "Number of Daughters " + titlePdgCode, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/daughters/PDGs", "Daughters PDGs " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/PDGsPrimaries", "Daughters PDGs Primaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/PDGsSecondaries", "Daughters PDGs Secondaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/daughters/prodVx", "Daughters Prod. Vertex X " + titlePdgCode, kTH1D, {axisProdx});
    histos.add("particle/daughters/prodVy", "Daughters Prod. Vertex Y " + titlePdgCode, kTH1D, {axisPrody});
    histos.add("particle/daughters/prodVz", "Daughters Prod. Vertex Z " + titlePdgCode, kTH1D, {axisProdz});
    histos.add("particle/daughters/prodRadiusVsPt", "Daughters Prod. Vertex Radius " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/daughters/prodRadius3DVsPt", "Daughters Prod. Vertex Radius XYZ " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});

    histos.add("particle/mothers/Number", "Number of Mothers " + titlePdgCode, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/mothers/PDGs", "Mothers PDGs " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/PDGsPrimaries", "Mothers PDGs Primaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/PDGsSecondaries", "Mothers PDGs Secondaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/prodVx", "Mothers Prod. Vertex X " + titlePdgCode, kTH1D, {axisProdx});
    histos.add("particle/mothers/prodVy", "Mothers Prod. Vertex Y " + titlePdgCode, kTH1D, {axisPrody});
    histos.add("particle/mothers/prodVz", "Mothers Prod. Vertex Z " + titlePdgCode, kTH1D, {axisProdz});
    histos.add("particle/mothers/prodRadiusVsPt", "Mothers Prod. Vertex Radius " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/mothers/prodRadius3DVsPt", "Mothers Prod. Vertex Radius XYZ " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});

    // Go up one generation
    histos.add("particle/mothers/mothers/Number", "Number of Mothers mothers " + titlePdgCode, kTH1D, {{20, -0.5, 19.5}});
    histos.add("particle/mothers/mothers/PDGs", "Mothers mothers PDGs " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/PDGsPrimaries", "Mothers mothers PDGs Primaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/PDGsSecondaries", "Mothers mothers PDGs Secondaries of " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("particle/mothers/mothers/prodVx", "Mothers mothers Prod. Vertex X " + titlePdgCode, kTH1D, {axisProdx});
    histos.add("particle/mothers/mothers/prodVy", "Mothers mothers Prod. Vertex Y " + titlePdgCode, kTH1D, {axisPrody});
    histos.add("particle/mothers/mothers/prodVz", "Mothers mothers Prod. Vertex Z " + titlePdgCode, kTH1D, {axisProdz});
    histos.add("particle/mothers/mothers/prodRadiusVsPt", "Mothers mothers Prod. Vertex Radius " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});
    histos.add("particle/mothers/mothers/prodRadius3DVsPt", "Mothers mothers Prod. Vertex Radius XYZ " + titlePdgCode, kTH2D, {axisPt, axisProdRadius});

    if (doprocessParticleOnly) {
      return;
    }
    histos.add("track/PDGs", "Track PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("track/withoutParticle", "Tracks without particles", kTH1D, {axisPt});
    histos.add("track/tofPDGs", "Track wTOF PDGs", kTH2D, {axisPDGs, axisCharge});
    histos.add("track/Pt", "Track Pt " + titlePdgCode, kTH1D, {axisPt});
    histos.add("track/primariesPt", "Track Pt (primary) " + titlePdgCode, kTH1D, {axisPt});
    histos.add("track/secondariesPt", "Track Pt (secondary) " + titlePdgCode, kTH1D, {axisPt});
    histos.add("track/PtvsMCPt", "Track Pt vs MC Pt " + titlePdgCode, kTH2D, {axisPt, axisPt});
    histos.add("track/Eta", "Track Eta " + titlePdgCode, kTH1D, {axisEta});
    histos.add("track/primariesEta", "Track Eta (primary) " + titlePdgCode, kTH1D, {axisEta});
    histos.add("track/secondariesEta", "Track Eta (secondary) " + titlePdgCode, kTH1D, {axisEta});
    histos.add("track/primaries", "Source for primaries " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
    histos.add("track/secondaries", "Source for secondaries " + titlePdgCode, kTH2D, {axisPDGs, axisCharge});
  }

  template <typename ParticleType>
  double getCharge(ParticleType const& particle)
  {
    auto pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return 10.f;
    }
    return pdgParticle->Charge();
  }

  template <typename ParticleType>
  const char* getPdgCodeString(ParticleType const& particle)
  {
    return Form("%i", particle.pdgCode());
  }

  // Nch estimator: charged physical primaries in |eta| < settings.nchEtaMax
  template <typename ParticleType>
  bool isNchParticle(ParticleType const& particle)
  {
    if (!particle.isPhysicalPrimary() || std::abs(particle.eta()) >= settings.nchEtaMax) {
      return false;
    }
    const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
    return pdgParticle && pdgParticle->Charge() != 0;
  }

  void fillSignalMultiplicity(const std::vector<int>& nchPerCollision, const std::vector<int>& nSignalPerCollision)
  {
    for (std::size_t iCollision = 0; iCollision < nchPerCollision.size(); ++iCollision) {
      histos.fill(HIST("event/Nch"), nchPerCollision[iCollision]);
      histos.fill(HIST("particle/NSignalVsNch"), nchPerCollision[iCollision], nSignalPerCollision[iCollision]);
    }
  }
  // Track histograms. Must be called after genericProcessParticles, which fills particlesOfInterest.
  template <typename TTracks>
  void genericProcessTracks(const TTracks& tracks)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        histos.fill(HIST("track/withoutParticle"), track.pt());
        continue;
      }
      const auto mcParticle = track.template mcParticle_as<aod::McParticles>();
      histos.get<TH2>(HIST("track/PDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      if (track.hasTOF()) {
        histos.get<TH2>(HIST("track/tofPDGs"))->Fill(getPdgCodeString(mcParticle), getCharge(mcParticle), 1.f);
      }
      if (settings.isStable.value == 0) {
        if (settings.doPrint) {
          LOG(info) << mcParticle.pdgCode() << " asked for " << pdgCode.value;
        }
        if (!mcParticle.has_mothers()) {
          continue;
        }
        // auto mothers = mcParticle.mothers();
        auto mothers = mcParticle.template mothers_as<aod::McParticles>();
        const auto particleIsInteresting = std::find(particlesOfInterest.begin(), particlesOfInterest.end(), mothers[0].globalIndex()) != particlesOfInterest.end();
        if (!particleIsInteresting) {
          continue;
        }
        if (settings.doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " comes from a " << mothers[0].pdgCode() << " and is a " << mcParticle.pdgCode();
        }
      } else {
        if (pdgCode != 0 && mcParticle.pdgCode() != pdgCode) {
          continue;
        }
        histos.fill(HIST("track/Pt"), track.pt() * settings.charge);
        histos.fill(HIST("track/PtvsMCPt"), track.pt() * settings.charge, mcParticle.pt());
        histos.fill(HIST("track/Eta"), track.eta());
        if (mcParticle.isPhysicalPrimary()) {
          histos.fill(HIST("track/primariesPt"), track.pt() * settings.charge);
          histos.fill(HIST("track/primariesEta"), track.eta());
        } else {
          histos.fill(HIST("track/secondariesPt"), track.pt() * settings.charge);
          histos.fill(HIST("track/secondariesEta"), track.eta());
        }

        if (!mcParticle.has_mothers()) {
          if (settings.doPrint) {
            LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode();
          }
          continue;
        }
        // auto mothers = mcParticle.mothers();
        auto mothers = mcParticle.template mothers_as<aod::McParticles>();
        if (mcParticle.isPhysicalPrimary()) {
          histos.get<TH2>(HIST("track/primaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        } else {
          histos.get<TH2>(HIST("track/secondaries"))->Fill(getPdgCodeString(mothers[0]), getCharge(mothers[0]), 1.f);
        }
        if (settings.doPrint) {
          LOG(info) << "Track " << track.globalIndex() << " is a " << mcParticle.pdgCode() << " and comes from a " << mothers[0].pdgCode() << " and is " << (mcParticle.isPhysicalPrimary() ? "" : "not") << " a primary";
        }
      }
    }
  }

  // Event and MC particle histograms. Stores the selected MC particles in particlesOfInterest.
  void genericProcessParticles(const o2::aod::McCollisions& colls,
                               const aod::McParticles& mcParticles)
  {
    for (const auto& col : colls) {
      histos.fill(HIST("event/VtxX"), col.posX());
      histos.fill(HIST("event/VtxY"), col.posY());
      histos.fill(HIST("event/VtxZ"), col.posZ());
    }
    std::vector<int> nchPerCollision(colls.size(), 0);
    std::vector<int> nSignalPerCollision(colls.size(), 0);
    particlesOfInterest.clear();
    for (const auto& mcParticle : mcParticles) {
      const auto& pdgString = getPdgCodeString(mcParticle);
      const auto& pdgCharge = getCharge(mcParticle);
      histos.get<TH2>(HIST("particle/PDGs"))->Fill(pdgString, pdgCharge, 1.f);
      if (mcParticle.isPhysicalPrimary()) {
        histos.get<TH2>(HIST("particle/PDGsPrimaries"))->Fill(pdgString, pdgCharge, 1.f);
      } else {
        histos.get<TH2>(HIST("particle/PDGsSecondaries"))->Fill(pdgString, pdgCharge, 1.f);
      }
      if (isNchParticle(mcParticle)) {
        ++nchPerCollision.at(mcParticle.mcCollisionId());
      }
      if (pdgCode != 0 && mcParticle.pdgCode() != pdgCode) {
        continue;
      }
      if (!yRange.contains(mcParticle.y())) {
        continue;
      }
      if (!etaRange.contains(mcParticle.eta())) {
        continue;
      }
      if (!ptRange.contains(mcParticle.pt())) {
        continue;
      }
      if (!prodZRange.contains(mcParticle.vz())) {
        continue;
      }
      ++nSignalPerCollision.at(mcParticle.mcCollisionId());
      histos.fill(HIST("particle/producedProcess"), mcParticle.getProcess());
      histos.fill(HIST("particle/Pt"), mcParticle.pt());
      histos.fill(HIST("particle/P"), mcParticle.p());
      histos.fill(HIST("particle/Eta"), mcParticle.eta());
      if (mcParticle.isPhysicalPrimary()) {
        histos.fill(HIST("particle/primariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/primariesP"), mcParticle.p());
        histos.fill(HIST("particle/primariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/primariesY"), mcParticle.y());
      } else {
        histos.fill(HIST("particle/secondariesPt"), mcParticle.pt());
        histos.fill(HIST("particle/secondariesP"), mcParticle.p());
        histos.fill(HIST("particle/secondariesEta"), mcParticle.eta());
        histos.fill(HIST("particle/secondariesY"), mcParticle.y());
      }
      histos.fill(HIST("particle/EvsPz"), mcParticle.e(), mcParticle.pz());
      histos.fill(HIST("particle/Y"), mcParticle.y());
      histos.fill(HIST("particle/YvzPz"), mcParticle.y(), mcParticle.pz());
      histos.fill(HIST("particle/EtavzPz"), mcParticle.eta(), mcParticle.pz());
      histos.fill(HIST("particle/PvzPz"), mcParticle.p(), mcParticle.pz());
      histos.fill(HIST("particle/PtvzPz"), mcParticle.pt(), mcParticle.pz());
      histos.fill(HIST("particle/Px"), mcParticle.px());
      histos.fill(HIST("particle/Py"), mcParticle.py());
      histos.fill(HIST("particle/Pz"), mcParticle.pz());
      histos.fill(HIST("particle/prodVx"), mcParticle.vx());
      histos.fill(HIST("particle/prodVy"), mcParticle.vy());
      histos.fill(HIST("particle/prodVz"), mcParticle.vz());
      if (mcParticle.has_daughters()) {
        const auto& daughters = mcParticle.daughters_as<aod::McParticles>();
        histos.fill(HIST("particle/daughters/Number"), daughters.size());
        bool filledDistance = false;
        for (const auto& daughter : daughters) {
          const auto& pdgStringDau = getPdgCodeString(daughter);
          const auto& pdgChargeDau = getCharge(daughter);

          histos.get<TH2>(HIST("particle/daughters/PDGs"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/daughters/PDGsPrimaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/daughters/PDGsSecondaries"))->Fill(pdgStringDau, pdgChargeDau, 1.f);
          }

          histos.fill(HIST("particle/daughters/prodVx"), daughter.vx());
          histos.fill(HIST("particle/daughters/prodVy"), daughter.vy());
          histos.fill(HIST("particle/daughters/prodVz"), daughter.vz());
          if (filledDistance) {
            continue;
          }
          filledDistance = true;
          histos.fill(HIST("particle/daughters/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(daughter.vx() * daughter.vx() + daughter.vy() * daughter.vy()));
          histos.fill(HIST("particle/daughters/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(daughter.vx() * daughter.vx() + daughter.vy() * daughter.vy() + daughter.vz() * daughter.vz()));

          const float dx = daughter.vx() - mcParticle.vx();
          const float dy = daughter.vy() - mcParticle.vy();
          const float dz = daughter.vz() - mcParticle.vz();
          histos.fill(HIST("particle/decayDistanceX"), dx);
          histos.fill(HIST("particle/decayDistanceY"), dy);
          histos.fill(HIST("particle/decayDistanceZ"), dz);
          histos.fill(HIST("particle/decayDistanceVsPt"), mcParticle.pt(), std::sqrt((dx) * (dx) + (dy) * (dy) + (dz) * (dz)));
          if (pdg->GetParticle(mcParticle.pdgCode())) {
            const float mass = pdg->GetParticle(mcParticle.pdgCode())->Mass();
            histos.fill(HIST("particle/decayCtau"), std::sqrt((dx) * (dx) + (dy) * (dy) + (dz) * (dz)) / mcParticle.pt() * mass * 10 * 1000);
          }
        }
      } else {
        histos.fill(HIST("particle/daughters/Number"), 0.f);
      }
      if (mcParticle.has_mothers()) {
        const auto& mothers = mcParticle.mothers_as<aod::McParticles>();
        for (const auto& mother : mothers) {
          const auto& pdgStringMot = getPdgCodeString(mother);
          const auto& pdgChargeMot = getCharge(mother);

          histos.get<TH2>(HIST("particle/mothers/PDGs"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          if (mcParticle.isPhysicalPrimary()) {
            histos.get<TH2>(HIST("particle/mothers/PDGsPrimaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          } else {
            histos.get<TH2>(HIST("particle/mothers/PDGsSecondaries"))->Fill(pdgStringMot, pdgChargeMot, 1.f);
          }

          histos.fill(HIST("particle/mothers/prodVx"), mother.vx());
          histos.fill(HIST("particle/mothers/prodVy"), mother.vy());
          histos.fill(HIST("particle/mothers/prodVz"), mother.vz());
          histos.fill(HIST("particle/mothers/prodRadiusVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy()));
          histos.fill(HIST("particle/mothers/prodRadius3DVsPt"), mother.pt(), std::sqrt(mother.vx() * mother.vx() + mother.vy() * mother.vy() + mother.vz() * mother.vz()));
          if (mother.has_mothers()) {
            const auto& mothers2 = mother.mothers_as<aod::McParticles>();
            for (const auto& mother2 : mothers2) {
              const auto& pdgStringMot2 = getPdgCodeString(mother2);
              const auto& pdgChargeMot2 = getCharge(mother2);

              histos.get<TH2>(HIST("particle/mothers/mothers/PDGs"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              if (mcParticle.isPhysicalPrimary()) {
                histos.get<TH2>(HIST("particle/mothers/mothers/PDGsPrimaries"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              } else {
                histos.get<TH2>(HIST("particle/mothers/mothers/PDGsSecondaries"))->Fill(pdgStringMot2, pdgChargeMot2, 1.f);
              }

              histos.fill(HIST("particle/mothers/mothers/prodVx"), mother2.vx());
              histos.fill(HIST("particle/mothers/mothers/prodVy"), mother2.vy());
              histos.fill(HIST("particle/mothers/mothers/prodVz"), mother2.vz());
              histos.fill(HIST("particle/mothers/mothers/prodRadiusVsPt"), mother2.pt(), std::sqrt(mother2.vx() * mother2.vx() + mother2.vy() * mother2.vy()));
              histos.fill(HIST("particle/mothers/mothers/prodRadius3DVsPt"), mother2.pt(), std::sqrt(mother2.vx() * mother2.vx() + mother2.vy() * mother2.vy() + mother2.vz() * mother2.vz()));
            }
          }
        }
      }

      histos.fill(HIST("particle/prodRadius"), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodVxVy"), mcParticle.vx(), mcParticle.vy());
      if (!mcParticle.isPhysicalPrimary()) {
        if (mcParticle.getProcess() == TMCProcess::kPDecay) {
          histos.fill(HIST("particle/prodVxVyStr"), mcParticle.vx(), mcParticle.vy());
        } else {
          histos.fill(HIST("particle/prodVxVyMat"), mcParticle.vx(), mcParticle.vy());
        }
      } else {
        histos.fill(HIST("particle/prodVxVyPrm"), mcParticle.vx(), mcParticle.vy());
      }
      histos.fill(HIST("particle/prodVxVsPt"), mcParticle.pt(), mcParticle.vx());
      histos.fill(HIST("particle/prodVyVsPt"), mcParticle.pt(), mcParticle.vy());
      histos.fill(HIST("particle/prodVzVsPt"), mcParticle.pt(), mcParticle.vz());
      histos.fill(HIST("particle/prodRadiusVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadiusVsEta"), mcParticle.eta(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy()));
      histos.fill(HIST("particle/prodRadius3DVsPt"), mcParticle.pt(), std::sqrt(mcParticle.vx() * mcParticle.vx() + mcParticle.vy() * mcParticle.vy() + mcParticle.vz() * mcParticle.vz()));
      particlesOfInterest.push_back(mcParticle.globalIndex());
    }
    fillSignalMultiplicity(nchPerCollision, nSignalPerCollision);
  }

  void processStandard(const o2::aod::McCollisions& colls,
                       const soa::Join<o2::aod::TracksIU, o2::aod::McTrackLabels, o2::aod::TracksExtra>& tracks,
                       const aod::McParticles& mcParticles)
  {
    genericProcessParticles(colls, mcParticles);
    genericProcessTracks(tracks);
  }

  void processParticleOnly(const o2::aod::McCollisions& colls,
                           const aod::McParticles& mcParticles)
  {
    genericProcessParticles(colls, mcParticles);
  }

  void processNonIU(const o2::aod::McCollisions& colls,
                    const soa::Join<o2::aod::Tracks, o2::aod::McTrackLabels, o2::aod::TracksExtra>& tracks,
                    const aod::McParticles& mcParticles)
  {
    genericProcessParticles(colls, mcParticles);
    genericProcessTracks(tracks);
  }

  PROCESS_SWITCH(Alice3SingleParticle, processStandard, "Process IU tracks", false);
  PROCESS_SWITCH(Alice3SingleParticle, processParticleOnly, "Process Particle only", true);
  PROCESS_SWITCH(Alice3SingleParticle, processNonIU, "Process non IU tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3SingleParticle>(cfgc)};
}
