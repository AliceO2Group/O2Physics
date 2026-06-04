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
/// \file coalescenceTreeProducer.cxx
/// \brief Tree producer for deuteron, triton, helium-3, and hypertriton coalescence studies.
///
/// The task loops over generated MC particles and builds bound-state
/// candidates according to the selected species. For each candidate,
/// the constituent particles are boosted to the candidate rest frame,
/// propagated to a common time, and used to compute the Jacobi relative
/// momenta.
///
/// Candidate triplets/pairs are stored in a reduced tree if:
/// - p_rho < pRhoMax
/// - p_lambda < pLambdaMax (for three-body coalescence only)
///
/// The output tree stores the kinematic and space-time coordinates of
/// the constituent baryons for offline coalescence studies.
///
/// \author Alberto Calivà <alberto.caliva@cern.ch>

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
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TTree.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

// Lightweight particle container
struct ReducedParticle {
  int64_t idPart = 0;
  int pdgPart = 0;
  float xPart = 0.f;
  float yPart = 0.f;
  float zPart = 0.f;
  float tPart = 0.f;
  float pxPart = 0.f;
  float pyPart = 0.f;
  float pzPart = 0.f;
};

struct CoalescenceTreeProducer {

  // Histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable analysis parameters
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum |z vertex| in cm"};
  Configurable<float> etaMax{"etaMax", 1.5f, "Maximum |eta| for generated particles"};
  Configurable<float> pRhoMax{"pRhoMax", 1.0f, "Maximum Jacobi p_rho in GeV/c"};
  Configurable<float> pLambdaMax{"pLambdaMax", 1.0f, "Maximum Jacobi p_lambda in GeV/c"};

  // Output tree storing bound-state candidates
  OutputObj<TTree> treeBoundState{"treeBoundState", OutputObjHandlingPolicy::AnalysisObject};

  // Variables to be stored in the tree
  Long64_t eventID = 0;                       // Event identifier
  Long64_t nMatterCandidatesPerEvent = 0;     // Number of matter candidates per event
  Long64_t nAntimatterCandidatesPerEvent = 0; // Number of antimatter candidates per event
  Long64_t idB1 = 0, idB2 = 0, idB3 = 0;      // Constituent identifiers
  int pdgB1 = 0, pdgB2 = 0, pdgB3 = 0;        // Constituent PDG codes

  // Constituent space-time coordinates and momentum components
  float xB1 = 0.f, yB1 = 0.f, zB1 = 0.f, tB1 = 0.f, pxB1 = 0.f, pyB1 = 0.f, pzB1 = 0.f;
  float xB2 = 0.f, yB2 = 0.f, zB2 = 0.f, tB2 = 0.f, pxB2 = 0.f, pyB2 = 0.f, pzB2 = 0.f;
  float xB3 = 0.f, yB3 = 0.f, zB3 = 0.f, tB3 = 0.f, pxB3 = 0.f, pyB3 = 0.f, pzB3 = 0.f;

  // Initialize output objects and histograms
  void init(InitContext const&)
  {
    // Event selection counter
    registry.add("eventCounter", "event Counter", HistType::kTH1F, {{5, 0, 5, "counter"}});
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(1, "Before z-vertex cut");
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(2, "After z-vertex cut");
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(3, "After non-empty constituent lists");
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(4, "After non-empty candidate lists");

    // Output tree for bound-state candidates
    treeBoundState.setObject(new TTree("BoundStateTree", "Tree for coalescence"));

    // Event information
    treeBoundState->Branch("eventID", &eventID, "eventID/L");
    treeBoundState->Branch("nMatterCandidatesPerEvent", &nMatterCandidatesPerEvent, "nMatterCandidatesPerEvent/L");
    treeBoundState->Branch("nAntimatterCandidatesPerEvent", &nAntimatterCandidatesPerEvent, "nAntimatterCandidatesPerEvent/L");

    // Constituent 1: ID, PDG code, space-time coordinates, and momentum components
    treeBoundState->Branch("idB1", &idB1, "idB1/L");
    treeBoundState->Branch("pdgB1", &pdgB1, "pdgB1/I");
    treeBoundState->Branch("xB1", &xB1, "xB1/F");
    treeBoundState->Branch("yB1", &yB1, "yB1/F");
    treeBoundState->Branch("zB1", &zB1, "zB1/F");
    treeBoundState->Branch("tB1", &tB1, "tB1/F");
    treeBoundState->Branch("pxB1", &pxB1, "pxB1/F");
    treeBoundState->Branch("pyB1", &pyB1, "pyB1/F");
    treeBoundState->Branch("pzB1", &pzB1, "pzB1/F");

    // Constituent 2: ID, PDG code, space-time coordinates, and momentum components
    treeBoundState->Branch("idB2", &idB2, "idB2/L");
    treeBoundState->Branch("pdgB2", &pdgB2, "pdgB2/I");
    treeBoundState->Branch("xB2", &xB2, "xB2/F");
    treeBoundState->Branch("yB2", &yB2, "yB2/F");
    treeBoundState->Branch("zB2", &zB2, "zB2/F");
    treeBoundState->Branch("tB2", &tB2, "tB2/F");
    treeBoundState->Branch("pxB2", &pxB2, "pxB2/F");
    treeBoundState->Branch("pyB2", &pyB2, "pyB2/F");
    treeBoundState->Branch("pzB2", &pzB2, "pzB2/F");

    // Constituent 3: ID, PDG code, space-time coordinates, and momentum components (three-body bound states only)
    if (doprocessHypertriton) {
      treeBoundState->Branch("idB3", &idB3, "idB3/L");
      treeBoundState->Branch("pdgB3", &pdgB3, "pdgB3/I");
      treeBoundState->Branch("xB3", &xB3, "xB3/F");
      treeBoundState->Branch("yB3", &yB3, "yB3/F");
      treeBoundState->Branch("zB3", &zB3, "zB3/F");
      treeBoundState->Branch("tB3", &tB3, "tB3/F");
      treeBoundState->Branch("pxB3", &pxB3, "pxB3/F");
      treeBoundState->Branch("pyB3", &pyB3, "pyB3/F");
      treeBoundState->Branch("pzB3", &pzB3, "pzB3/F");
    }
  }

  // Assign mass based on PDG code
  static double massFromPdg(int pdg)
  {
    switch (std::abs(pdg)) {
      case PDG_t::kProton:
        return o2::constants::physics::MassProton;
      case PDG_t::kNeutron:
        return o2::constants::physics::MassNeutron;
      case PDG_t::kLambda0:
        return o2::constants::physics::MassLambda0;
      default:
        return 0.0;
    }
  }

  // Apply a momentum-space skim using Jacobi momenta in the candidate rest frame
  template <typename ReducedPart>
  bool passThreeBodySkim(const ReducedPart& b1, const ReducedPart& b2, const ReducedPart& b3)
  {
    // Constituent masses from PDG codes
    double m1 = massFromPdg(b1.pdgPart);
    double m2 = massFromPdg(b2.pdgPart);
    double m3 = massFromPdg(b3.pdgPart);
    if (m1 <= 0.0 || m2 <= 0.0 || m3 <= 0.0) {
      return false;
    }

    // Constituent four-momenta in the laboratory frame
    auto p4B1 = ROOT::Math::PxPyPzMVector{b1.pxPart, b1.pyPart, b1.pzPart, m1};
    auto p4B2 = ROOT::Math::PxPyPzMVector{b2.pxPart, b2.pyPart, b2.pzPart, m2};
    auto p4B3 = ROOT::Math::PxPyPzMVector{b3.pxPart, b3.pyPart, b3.pzPart, m3};

    // Candidate four-momentum
    auto candidate = p4B1 + p4B2 + p4B3;

    // Boost to the candidate rest frame
    auto betaCandidate = candidate.BoostToCM();
    ROOT::Math::Boost boostToRest(betaCandidate);

    // Constituent momenta in the candidate rest frame
    auto p4B1Star = boostToRest(p4B1);
    auto p4B2Star = boostToRest(p4B2);
    auto p4B3Star = boostToRest(p4B3);
    const ROOT::Math::XYZVector pB1Star = p4B1Star.Vect();
    const ROOT::Math::XYZVector pB2Star = p4B2Star.Vect();
    const ROOT::Math::XYZVector pB3Star = p4B3Star.Vect();

    // Reduced masses entering the Jacobi coordinates
    const double m12 = m1 + m2;
    const double mTot = m1 + m2 + m3;

    // Jacobi momenta in the candidate rest frame
    ROOT::Math::XYZVector pRho = (m2 * pB1Star - m1 * pB2Star) / m12;
    ROOT::Math::XYZVector pLambda = (m3 * (pB1Star + pB2Star) - m12 * pB3Star) / mTot;

    // Reject candidates outside the momentum-space acceptance
    if (pRho.R() > pRhoMax || pLambda.R() > pLambdaMax) {
      return false;
    }

    return true;
  }

  // Group MC particles by their associated MC collision
  Preslice<aod::McParticles> mcParticlesPerMcCollision = o2::aod::mcparticle::mcCollisionId;

  // Process Hypertriton
  void processHypertriton(aod::McCollisions const& collisions, aod::McParticles const& mcParticles)
  {
    // Loop over MC collisions
    for (const auto& collision : collisions) {

      // Event counter before selections
      registry.fill(HIST("eventCounter"), 0.5);

      // Apply z-vertex selection
      if (std::fabs(collision.posZ()) > zVtxMax)
        continue;

      // Event counter after z-vertex cut
      registry.fill(HIST("eventCounter"), 1.5);

      // To be implemented: maybe add INEL>0 selection

      // Get particles in this MC collision
      const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

      // Event ID
      eventID = collision.globalIndex();

      // Containers for candidate constituents
      std::vector<ReducedParticle> protons;
      std::vector<ReducedParticle> neutrons;
      std::vector<ReducedParticle> lambdas;

      // Loop over MC particles
      for (const auto& particle : mcParticlesThisMcColl) {

        // Select physical primary particles only
        if (!particle.isPhysicalPrimary())
          continue;

        // Select particles with |eta|< eta_max
        if (std::fabs(particle.eta()) > etaMax)
          continue;

        // Store protons
        if (std::abs(particle.pdgCode()) == PDG_t::kProton) {
          protons.push_back({particle.globalIndex(), particle.pdgCode(), particle.vx(), particle.vy(), particle.vz(), particle.vt(), particle.px(), particle.py(), particle.pz()});
        }

        // Store neutrons
        if (std::abs(particle.pdgCode()) == PDG_t::kNeutron) {
          neutrons.push_back({particle.globalIndex(), particle.pdgCode(), particle.vx(), particle.vy(), particle.vz(), particle.vt(), particle.px(), particle.py(), particle.pz()});
        }

        // Store lambdas
        if (std::abs(particle.pdgCode()) == PDG_t::kLambda0) {
          lambdas.push_back({particle.globalIndex(), particle.pdgCode(), particle.vx(), particle.vy(), particle.vz(), particle.vt(), particle.px(), particle.py(), particle.pz()});
        }
      } // end of loop over MC particles

      // Reject events that do not contain at least one proton, one neutron, and one lambda
      if (protons.empty() || neutrons.empty() || lambdas.empty())
        continue;

      // Event counter: events containing all three constituent species
      registry.fill(HIST("eventCounter"), 2.5);

      // Count matter and antimatter candidates in the event
      Long64_t nCandidatesMatter(0);
      Long64_t nCandidatesAntiMatter(0);
      for (const auto& proton : protons) {
        for (const auto& neutron : neutrons) {
          for (const auto& lambda : lambdas) {

            const bool isMatter = proton.pdgPart > 0 && neutron.pdgPart > 0 && lambda.pdgPart > 0;
            const bool isAntimatter = proton.pdgPart < 0 && neutron.pdgPart < 0 && lambda.pdgPart < 0;

            // Skip mixed-sign combinations
            if (!isMatter && !isAntimatter)
              continue;
            const bool passSkim = passThreeBodySkim(proton, neutron, lambda);

            if (isMatter && passSkim)
              nCandidatesMatter++;
            if (isAntimatter && passSkim)
              nCandidatesAntiMatter++;
          }
        }
      }

      // Reject events with no accepted candidates
      if (nCandidatesMatter == 0 && nCandidatesAntiMatter == 0)
        continue;

      // Event counter: number of events with at least one candidate
      registry.fill(HIST("eventCounter"), 3.5);

      // Store number of candidates per event
      nMatterCandidatesPerEvent = nCandidatesMatter;
      nAntimatterCandidatesPerEvent = nCandidatesAntiMatter;

      // Loop over accepted triplets and fill the tree
      for (const auto& proton : protons) {
        for (const auto& neutron : neutrons) {
          for (const auto& lambda : lambdas) {

            const bool isMatter = proton.pdgPart > 0 && neutron.pdgPart > 0 && lambda.pdgPart > 0;
            const bool isAntimatter = proton.pdgPart < 0 && neutron.pdgPart < 0 && lambda.pdgPart < 0;

            // Skip mixed-sign combinations
            if (!isMatter && !isAntimatter)
              continue;

            // Apply momentum-space selection
            const bool passSkim = passThreeBodySkim(proton, neutron, lambda);
            if (!passSkim)
              continue;

            // Fill tree
            idB1 = proton.idPart;
            pdgB1 = proton.pdgPart;
            xB1 = proton.xPart;
            yB1 = proton.yPart;
            zB1 = proton.zPart;
            tB1 = proton.tPart;
            pxB1 = proton.pxPart;
            pyB1 = proton.pyPart;
            pzB1 = proton.pzPart;

            idB2 = neutron.idPart;
            pdgB2 = neutron.pdgPart;
            xB2 = neutron.xPart;
            yB2 = neutron.yPart;
            zB2 = neutron.zPart;
            tB2 = neutron.tPart;
            pxB2 = neutron.pxPart;
            pyB2 = neutron.pyPart;
            pzB2 = neutron.pzPart;

            idB3 = lambda.idPart;
            pdgB3 = lambda.pdgPart;
            xB3 = lambda.xPart;
            yB3 = lambda.yPart;
            zB3 = lambda.zPart;
            tB3 = lambda.tPart;
            pxB3 = lambda.pxPart;
            pyB3 = lambda.pyPart;
            pzB3 = lambda.pzPart;

            treeBoundState->Fill();
          }
        }
      }
    } // end of loop over mc collisions
  }
  PROCESS_SWITCH(CoalescenceTreeProducer, processHypertriton, "process hypertriton", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CoalescenceTreeProducer>(cfgc)};
}
