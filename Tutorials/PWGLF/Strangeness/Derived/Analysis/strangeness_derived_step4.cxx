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
/// \brief Step4 of the Strangeness tutorial
/// \author Romain Schotter
/// based on the original codes from:
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/EventSelection.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all cascades and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of Cascades
// STEP 2
// Apply TPC PID selections on cascade daughter tracks
// STEP 3
// Apply TOF PID selections on cascade daugther tracks (if info is available)
// STEP 4
// Check the MC information of the cascades

struct strangeness_derived_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rXi{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rOmega{"omega", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rGenParticles{"genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable parameters for cascade selection
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.98, "Casc CosPA"};
  Configurable<double> cascadesetting_v0cospa{"cascadesetting_v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "DCA cascade daughters"};
  Configurable<float> cascadesetting_dcav0dau{"cascadesetting_dcav0dau", 1.0, "DCA v0 daughters"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.06, "DCA bachelor to PV"};
  Configurable<float> cascadesetting_dcapostopv{"cascadesetting_dcapostopv", 0.06, "DCA positive to PV"};
  Configurable<float> cascadesetting_dcanegtopv{"cascadesetting_dcanegtopv", 0.06, "DCA negative to PV"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "minimum V0 DCA to PV"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.5, "cascradius"};
  Configurable<float> cascadesetting_v0radius{"cascadesetting_v0radius", 1.2, "v0radius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "v0 mass window"};
  Configurable<float> cascadesetting_competingmassrej{"cascadesetting_competingmassrej", 0.008, "Competing mass rejection"};

  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCKaon{"NSigmaTPCKaon", 4, "NSigmaTPCKaon"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  // Configurable parameters for TOF PID selection
  Configurable<float> NSigmaTOFPion{"NSigmaTOFPion", 3, "NSigmaTOFPion"};
  Configurable<float> NSigmaTOFKaon{"NSigmaTOFKaon", 3, "NSigmaTOFKaon"};
  Configurable<float> NSigmaTOFProton{"NSigmaTOFProton", 3, "NSigmaTOFProton"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec XiMassAxis = {100, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec OmegaMassAxis = {100, 1.63f, 1.7f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});

    // Xi/Omega reconstruction
    rXi.add("hMassXi", "hMassXi", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiSelectedWithTOF", "hMassXiSelectedWithTOF", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hMassXiTrueRec", "hMassXiTrueRec", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hPtXiTrueRec", "hPtXiTrueRec", {HistType::kTH1F, {ptAxis}});
    rXi.add("hMassXiTrueRecWithTOF", "hMassXiTrueRecWithTOF", {HistType::kTH1F, {XiMassAxis}});
    rXi.add("hPtXiTrueRecWithTOF", "hPtXiTrueRecWithTOF", {HistType::kTH1F, {ptAxis}});

    rOmega.add("hMassOmega", "hMassOmega", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hMassOmegaSelected", "hMassOmegaSelected", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hMassOmegaSelectedWithTOF", "hMassOmegaSelectedWithTOF", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hMassOmegaTrueRec", "hMassOmegaTrueRec", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hPtOmegaTrueRec", "hPtOmegaTrueRec", {HistType::kTH1F, {ptAxis}});
    rOmega.add("hMassOmegaTrueRecWithTOF", "hMassOmegaTrueRecWithTOF", {HistType::kTH1F, {OmegaMassAxis}});
    rOmega.add("hPtOmegaTrueRecWithTOF", "hPtOmegaTrueRecWithTOF", {HistType::kTH1F, {ptAxis}});

    // Xi/Omega topological cuts
    rXi.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rXi.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});

    rOmega.add("hCascDCAV0Daughters", "hCascDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rOmega.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.95f, 1.f}}});

    // Generated level histograms
    rEventSelection.add("hVertexZGen", "hVertexZGen", {HistType::kTH1F, {vertexZAxis}});
    rGenParticles.add("hPtXiGen", "hPtXiGen", {HistType::kTH1F, {{ptAxis}}});
    rGenParticles.add("hPtOmegaGen", "hPtOmegaGen", {HistType::kTH1F, {{ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  // Filters on Cascades
  // Cannot filter on dynamic columns
  Filter preFilterCascades = (aod::cascdata::dcaV0daughters < cascadesetting_dcav0dau &&
                              nabs(aod::cascdata::dcapostopv) > cascadesetting_dcapostopv &&
                              nabs(aod::cascdata::dcanegtopv) > cascadesetting_dcanegtopv &&
                              nabs(aod::cascdata::dcabachtopv) > cascadesetting_dcabachtopv &&
                              aod::cascdata::dcacascdaughters < cascadesetting_dcacascdau);

  // Defining the type of the daughter tracks
  using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs, aod::DauTrackMCIds>;

  void processRecMC(soa::Filtered<soa::Join<aod::StraCollisions, aod::StraEvSels>>::iterator const& collision,
                    soa::Filtered<soa::Join<aod::CascCores, aod::CascExtras, aod::CascTOFNSigmas, aod::CascCoreMCLabels>> const& Cascades,
                    dauTracks const&,
                    aod::CascMCCores const& /*cascmccores*/)
  {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // Cascades
    for (const auto& casc : Cascades) {
      const auto& bachDaughterTrackCasc = casc.bachTrackExtra_as<dauTracks>();
      const auto& posDaughterTrackCasc = casc.posTrackExtra_as<dauTracks>();
      const auto& negDaughterTrackCasc = casc.negTrackExtra_as<dauTracks>();

      rXi.fill(HIST("hMassXi"), casc.mXi());
      rOmega.fill(HIST("hMassOmega"), casc.mOmega());

      // Cut on dynamic columns
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_cospa)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_v0cospa)
        continue;
      if (std::abs(casc.mLambda() - o2::constants::physics::MassLambda) > cascadesetting_v0masswindow)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascadesetting_mindcav0topv)
        continue;
      if (casc.cascradius() < cascadesetting_cascradius)
        continue;
      if (casc.v0radius() < cascadesetting_v0radius)
        continue;

      // PID selection
      if (casc.sign() < 0) {
        if (std::abs(posDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (std::abs(negDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      } else {
        if (std::abs(negDaughterTrackCasc.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (std::abs(posDaughterTrackCasc.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      }

      // TOF PID check
      bool xiPassTOFSelection = true;
      bool omegaPassTOFSelection = true;
      if (casc.sign() < 0) {
        if (posDaughterTrackCasc.hasTOF()) {
          if (std::abs(casc.tofNSigmaXiLaPr()) > NSigmaTOFProton) {
            xiPassTOFSelection &= false;
          }
          if (std::abs(casc.tofNSigmaOmLaPr()) > NSigmaTOFProton) {
            omegaPassTOFSelection &= false;
          }
        }
        if (negDaughterTrackCasc.hasTOF()) {
          if (std::abs(casc.tofNSigmaXiLaPi()) > NSigmaTOFPion) {
            xiPassTOFSelection &= false;
          }
          if (std::abs(casc.tofNSigmaOmLaPi()) > NSigmaTOFPion) {
            omegaPassTOFSelection &= false;
          }
        }
      } else {
        if (posDaughterTrackCasc.hasTOF()) {
          if (std::abs(casc.tofNSigmaXiLaPi()) > NSigmaTOFPion) {
            xiPassTOFSelection &= false;
          }
          if (std::abs(casc.tofNSigmaOmLaPi()) > NSigmaTOFPion) {
            omegaPassTOFSelection &= false;
          }
        }
        if (negDaughterTrackCasc.hasTOF()) {
          if (std::abs(casc.tofNSigmaXiLaPr()) > NSigmaTOFProton) {
            xiPassTOFSelection &= false;
          }
          if (std::abs(casc.tofNSigmaOmLaPr()) > NSigmaTOFProton) {
            omegaPassTOFSelection &= false;
          }
        }
      }

      if (bachDaughterTrackCasc.hasTOF()) {
        if (std::abs(casc.tofNSigmaXiPi()) > NSigmaTOFPion) {
          xiPassTOFSelection &= false;
        }
        if (std::abs(casc.tofNSigmaOmKa()) > NSigmaTOFKaon) {
          omegaPassTOFSelection &= false;
        }
      }

      if (bachDaughterTrackCasc.hasTOF()) {
        if (std::abs(casc.tofNSigmaXiPi()) > NSigmaTOFPion) {
          xiPassTOFSelection &= false;
        }
        if (std::abs(casc.tofNSigmaOmKa()) > NSigmaTOFKaon) {
          omegaPassTOFSelection &= false;
        }
      }

      // Fill histograms! (if possible)
      if (std::abs(bachDaughterTrackCasc.tpcNSigmaPi()) < NSigmaTPCPion) { // Xi case
        rXi.fill(HIST("hMassXiSelected"), casc.mXi());
        if (xiPassTOFSelection) {
          rXi.fill(HIST("hMassXiSelectedWithTOF"), casc.mXi());
        }

        rXi.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
        rXi.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      }
      if (std::abs(bachDaughterTrackCasc.tpcNSigmaKa()) < NSigmaTPCKaon) {                                  // Omega case
        if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascadesetting_competingmassrej) { // competing mass rejection, only in case of Omega
          rOmega.fill(HIST("hMassOmegaSelected"), casc.mOmega());
          if (omegaPassTOFSelection) {
            rOmega.fill(HIST("hMassOmegaSelectedWithTOF"), casc.mOmega());
          }

          rOmega.fill(HIST("hCascDCAV0Daughters"), casc.dcaV0daughters());
          rOmega.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        }
      }

      // MC truth info
      if (!casc.has_cascMCCore()) {
        continue;
      }
      auto cascmccore = casc.cascMCCore_as<aod::CascMCCores>();

      // Checking that the cascade is a true Xi
      if (std::abs(cascmccore.pdgCode()) == PDG_t::kXiMinus) {
        if (std::abs(bachDaughterTrackCasc.tpcNSigmaPi()) < NSigmaTPCPion) { // Xi case
          rXi.fill(HIST("hMassXiTrueRec"), casc.mXi());
          rXi.fill(HIST("hPtXiTrueRec"), casc.pt());
          if (xiPassTOFSelection) {
            rXi.fill(HIST("hMassXiTrueRecWithTOF"), casc.mXi());
            rXi.fill(HIST("hPtXiTrueRecWithTOF"), casc.pt());
          }
        }
      }
      if (std::abs(cascmccore.pdgCode()) == PDG_t::kOmegaMinus) {
        if (std::abs(bachDaughterTrackCasc.tpcNSigmaKa()) < NSigmaTPCKaon) {                                  // Omega case
          if (std::abs(casc.mXi() - o2::constants::physics::MassXiMinus) > cascadesetting_competingmassrej) { // competing mass rejection, only in case of Omega
            rOmega.fill(HIST("hMassOmegaTrueRec"), casc.mOmega());
            rOmega.fill(HIST("hPtOmegaTrueRec"), casc.pt());
            if (omegaPassTOFSelection) {
              rOmega.fill(HIST("hMassOmegaTrueRecWithTOF"), casc.mOmega());
              rOmega.fill(HIST("hPtOmegaTrueRecWithTOF"), casc.pt());
            }
          }
        }
      }
    }
  }

  void processGenMC(soa::Filtered<aod::StraMCCollisions>::iterator const& mcCollision,
                    const soa::SmallGroups<soa::Join<aod::StraCollisions, o2::aod::StraEvSels, aod::StraCollLabels>>& collisions,
                    const soa::SmallGroups<soa::Join<aod::CascMCCores, aod::CascMCCollRefs>>& cascMC)
  {
    if (collisions.size() < 1) // to process generated collisions that've been reconstructed at least once
      return;
    rEventSelection.fill(HIST("hVertexZGen"), mcCollision.posZ());

    for (const auto& cascmc : cascMC) {
      if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus) {
        rGenParticles.fill(HIST("hPtXiGen"), cascmc.ptMC());
      }
      if (std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus) {
        rGenParticles.fill(HIST("hPtOmegaGen"), cascmc.ptMC());
      }
    }
  }

  PROCESS_SWITCH(strangeness_derived_tutorial, processRecMC, "Process Run 3 mc, reconstructed", true);
  PROCESS_SWITCH(strangeness_derived_tutorial, processGenMC, "Process Run 3 mc, generated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeness_derived_tutorial>(cfgc)};
}
