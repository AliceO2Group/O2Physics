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
/// \brief Tree producer for deuteron, triton, helium-3, and hypertriton
/// coalescence studies.
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

using namespace o2;
using namespace o2::framework;

struct CoalescenceTreeProducer {

  // Supported bound-state species
  enum BoundStateSpecies {
    kDeuteron = 0,
    kTriton = 1,
    kHelium3 = 2,
    kHypertriton = 3
  };

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Analysis parameters
  Configurable<float> zVtx{"zVtx", 10.f, "Maximum |z vertex| in cm"};
  Configurable<float> etaMax{"etaMax", 1.5f, "Maximum |eta| for generated particles"};
  Configurable<float> pRhoMax{"pRhoMax", 1.0f, "Maximum Jacobi p_rho in GeV/c"};
  Configurable<float> pLambdaMax{"pLambdaMax", 1.0f, "Maximum Jacobi p_lambda in GeV/c"};
  Configurable<int> boundStateSpecies{"boundStateSpecies", kHypertriton, "Species selector: 0 = deuteron, 1 = triton, 2 = helium3, 3 = hypertriton"};

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  OutputObj<TTree> treeBoundState{"treeBoundState"};

  int64_t eventID = 0;                  // Event ID
  int64_t idB1 = 0, idB2 = 0, idB3 = 0; // MC particle IDs

  int pdgB1 = 0, pdgB2 = 0, pdgB3 = 0;
  int chargeB1 = 0, chargeB2 = 0, chargeB3 = 0;

  // Space-time coordinates and momentum components
  float xB1 = 0.f, yB1 = 0.f, zB1 = 0.f, tB1 = 0.f, pxB1 = 0.f, pyB1 = 0.f, pzB1 = 0.f;
  float xB2 = 0.f, yB2 = 0.f, zB2 = 0.f, tB2 = 0.f, pxB2 = 0.f, pyB2 = 0.f, pzB2 = 0.f;
  float xB3 = 0.f, yB3 = 0.f, zB3 = 0.f, tB3 = 0.f, pxB3 = 0.f, pyB3 = 0.f, pzB3 = 0.f;

  struct Particle {
    int64_t id = 0;
    int pdg = 0;
    int charge = 0;

    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    float t = 0.f;

    float px = 0.f;
    float py = 0.f;
    float pz = 0.f;

    float mass = 0.f;

    Particle() = default;

    template <typename T>
    explicit Particle(T const& p)
      : id(p.globalIndex()),
        pdg(p.pdgCode()),
        charge(chargeFromPdg(pdg)),
        x(p.vx()),
        y(p.vy()),
        z(p.vz()),
        t(p.vt()),
        px(p.px()),
        py(p.py()),
        pz(p.pz()),
        mass(static_cast<float>(massFromPdg(pdg)))
    {
    }

    static int chargeFromPdg(int pdg)
    {
      if (pdg == PDG_t::kProton) {
        return 1;
      }
      if (pdg == PDG_t::kProtonBar) {
        return -1;
      }
      return 0;
    }

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
          return -1.;
      }
    }
  };

  void init(InitContext&)
  {
    registry.add("eventCounter", "event Counter", HistType::kTH1F, {{5, 0, 5, "counter"}});
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(1, "Before z-vertex cut");
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(2, "After z-vertex cut");
    registry.get<TH1>(HIST("eventCounter"))->GetXaxis()->SetBinLabel(3, "After non-empty lists");

    // Tree for pairs/triplets.
    // For deuterons only the first two baryons are stored.
    // For three-body states, the third-baryon branches are also created.
    treeBoundState->Branch("eventID", &eventID);

    treeBoundState->Branch("idB1", &idB1);
    treeBoundState->Branch("idB2", &idB2);

    treeBoundState->Branch("pdgB1", &pdgB1);
    treeBoundState->Branch("pdgB2", &pdgB2);

    treeBoundState->Branch("chargeB1", &chargeB1);
    treeBoundState->Branch("chargeB2", &chargeB2);

    treeBoundState->Branch("xB1", &xB1);
    treeBoundState->Branch("yB1", &yB1);
    treeBoundState->Branch("zB1", &zB1);
    treeBoundState->Branch("tB1", &tB1);
    treeBoundState->Branch("pxB1", &pxB1);
    treeBoundState->Branch("pyB1", &pyB1);
    treeBoundState->Branch("pzB1", &pzB1);

    treeBoundState->Branch("xB2", &xB2);
    treeBoundState->Branch("yB2", &yB2);
    treeBoundState->Branch("zB2", &zB2);
    treeBoundState->Branch("tB2", &tB2);
    treeBoundState->Branch("pxB2", &pxB2);
    treeBoundState->Branch("pyB2", &pyB2);
    treeBoundState->Branch("pzB2", &pzB2);

    if (static_cast<int>(boundStateSpecies) != kDeuteron) {
      treeBoundState->Branch("idB3", &idB3);
      treeBoundState->Branch("pdgB3", &pdgB3);
      treeBoundState->Branch("chargeB3", &chargeB3);

      treeBoundState->Branch("xB3", &xB3);
      treeBoundState->Branch("yB3", &yB3);
      treeBoundState->Branch("zB3", &zB3);
      treeBoundState->Branch("tB3", &tB3);
      treeBoundState->Branch("pxB3", &pxB3);
      treeBoundState->Branch("pyB3", &pyB3);
      treeBoundState->Branch("pzB3", &pzB3);
    }
  }

  ROOT::Math::PxPyPzMVector makeFourVector(Particle const& p)
  {
    return ROOT::Math::PxPyPzMVector{p.px, p.py, p.pz, p.mass};
  }

  void storePair(Particle const& b1, Particle const& b2)
  {
    idB1 = b1.id;
    idB2 = b2.id;

    pdgB1 = b1.pdg;
    pdgB2 = b2.pdg;

    chargeB1 = b1.charge;
    chargeB2 = b2.charge;

    xB1 = b1.x;
    yB1 = b1.y;
    zB1 = b1.z;
    tB1 = b1.t;
    pxB1 = b1.px;
    pyB1 = b1.py;
    pzB1 = b1.pz;

    xB2 = b2.x;
    yB2 = b2.y;
    zB2 = b2.z;
    tB2 = b2.t;
    pxB2 = b2.px;
    pyB2 = b2.py;
    pzB2 = b2.pz;

    treeBoundState->Fill();
  }

  void storeTriplet(Particle const& b1, Particle const& b2, Particle const& b3)
  {
    idB1 = b1.id;
    idB2 = b2.id;
    idB3 = b3.id;

    pdgB1 = b1.pdg;
    pdgB2 = b2.pdg;
    pdgB3 = b3.pdg;

    chargeB1 = b1.charge;
    chargeB2 = b2.charge;
    chargeB3 = b3.charge;

    xB1 = b1.x;
    yB1 = b1.y;
    zB1 = b1.z;
    tB1 = b1.t;
    pxB1 = b1.px;
    pyB1 = b1.py;
    pzB1 = b1.pz;

    xB2 = b2.x;
    yB2 = b2.y;
    zB2 = b2.z;
    tB2 = b2.t;
    pxB2 = b2.px;
    pyB2 = b2.py;
    pzB2 = b2.pz;

    xB3 = b3.x;
    yB3 = b3.y;
    zB3 = b3.z;
    tB3 = b3.t;
    pxB3 = b3.px;
    pyB3 = b3.py;
    pzB3 = b3.pz;

    treeBoundState->Fill();
  }

  bool passTwoBodySkim(Particle const& b1, Particle const& b2)
  {
    auto p4B1 = makeFourVector(b1);
    auto p4B2 = makeFourVector(b2);

    auto candidate = p4B1 + p4B2;

    // BoostToCM() returns the boost vector needed to go from the lab frame
    // to the candidate center-of-mass frame.
    auto betaCandidate = candidate.BoostToCM();
    ROOT::Math::Boost boostToRest(betaCandidate);

    auto p4B1Star = boostToRest(p4B1);
    auto p4B2Star = boostToRest(p4B2);

    /*
    // Space-time propagation block.
    // Not needed for the present skim because the selection is based only on
    // the relative momentum pRho in the candidate rest frame.
    //
    // This block can be re-enabled later if one wants to apply coordinate-space
    // or full phase-space coalescence selections, e.g. cuts on relative distance
    // after propagation to a common time.

    ROOT::Math::XYZTVector x4B1{b1.x, b1.y, b1.z, b1.t};
    ROOT::Math::XYZTVector x4B2{b2.x, b2.y, b2.z, b2.t};

    auto x4B1Star = boostToRest(x4B1);
    auto x4B2Star = boostToRest(x4B2);

    const double tMax = std::max(x4B1Star.T(), x4B2Star.T());

    ROOT::Math::XYZVector rB1{x4B1Star.X(), x4B1Star.Y(), x4B1Star.Z()};
    ROOT::Math::XYZVector rB2{x4B2Star.X(), x4B2Star.Y(), x4B2Star.Z()};

    ROOT::Math::XYZVector vB1 = p4B1Star.Vect() / p4B1Star.E();
    ROOT::Math::XYZVector vB2 = p4B2Star.Vect() / p4B2Star.E();

    rB1 += vB1 * (tMax - x4B1Star.T());
    rB2 += vB2 * (tMax - x4B2Star.T());
    */

    const ROOT::Math::XYZVector pB1Star = p4B1Star.Vect();
    const ROOT::Math::XYZVector pB2Star = p4B2Star.Vect();

    const double m12 = b1.mass + b2.mass;

    ROOT::Math::XYZVector pRho = (b2.mass * pB1Star - b1.mass * pB2Star) / m12;

    if (pRho.R() > pRhoMax) {
      return false;
    }
    return true;
  }

  bool passThreeBodySkim(Particle const& b1, Particle const& b2, Particle const& b3)
  {
    auto p4B1 = makeFourVector(b1);
    auto p4B2 = makeFourVector(b2);
    auto p4B3 = makeFourVector(b3);

    auto candidate = p4B1 + p4B2 + p4B3;

    // BoostToCM() returns the boost vector needed to go from the lab frame
    // to the candidate center-of-mass frame.
    auto betaCandidate = candidate.BoostToCM();
    ROOT::Math::Boost boostToRest(betaCandidate);

    auto p4B1Star = boostToRest(p4B1);
    auto p4B2Star = boostToRest(p4B2);
    auto p4B3Star = boostToRest(p4B3);

    /*
    // Space-time propagation block.
    // Not needed for the present skim because the selection is based only on
    // the Jacobi relative momenta pRho and pLambda in the candidate rest frame.
    //
    // This block can be re-enabled later if one wants to apply coordinate-space
    // or full phase-space coalescence selections, e.g. cuts on relative distances
    // after propagation to a common time.

    ROOT::Math::XYZTVector x4B1{b1.x, b1.y, b1.z, b1.t};
    ROOT::Math::XYZTVector x4B2{b2.x, b2.y, b2.z, b2.t};
    ROOT::Math::XYZTVector x4B3{b3.x, b3.y, b3.z, b3.t};

    auto x4B1Star = boostToRest(x4B1);
    auto x4B2Star = boostToRest(x4B2);
    auto x4B3Star = boostToRest(x4B3);

    const double tMax = std::max({x4B1Star.T(), x4B2Star.T(), x4B3Star.T()});

    ROOT::Math::XYZVector rB1{x4B1Star.X(), x4B1Star.Y(), x4B1Star.Z()};
    ROOT::Math::XYZVector rB2{x4B2Star.X(), x4B2Star.Y(), x4B2Star.Z()};
    ROOT::Math::XYZVector rB3{x4B3Star.X(), x4B3Star.Y(), x4B3Star.Z()};

    ROOT::Math::XYZVector vB1 = p4B1Star.Vect() / p4B1Star.E();
    ROOT::Math::XYZVector vB2 = p4B2Star.Vect() / p4B2Star.E();
    ROOT::Math::XYZVector vB3 = p4B3Star.Vect() / p4B3Star.E();

    rB1 += vB1 * (tMax - x4B1Star.T());
    rB2 += vB2 * (tMax - x4B2Star.T());
    rB3 += vB3 * (tMax - x4B3Star.T());
    */

    const ROOT::Math::XYZVector pB1Star = p4B1Star.Vect();
    const ROOT::Math::XYZVector pB2Star = p4B2Star.Vect();
    const ROOT::Math::XYZVector pB3Star = p4B3Star.Vect();

    const double m12 = b1.mass + b2.mass;
    const double mTot = b1.mass + b2.mass + b3.mass;

    ROOT::Math::XYZVector pRho = (b2.mass * pB1Star - b1.mass * pB2Star) / m12;
    ROOT::Math::XYZVector pLambda = (b3.mass * (pB1Star + pB2Star) - m12 * pB3Star) / mTot;

    if (pRho.R() > pRhoMax) {
      return false;
    }

    if (pLambda.R() > pLambdaMax) {
      return false;
    }

    return true;
  }

  void buildPairs(std::vector<Particle> const& b1Candidates,
                  std::vector<Particle> const& b2Candidates)
  {
    if (b1Candidates.empty() || b2Candidates.empty()) {
      return;
    }

    for (const auto& b1 : b1Candidates) {
      for (const auto& b2 : b2Candidates) {
        if (!passTwoBodySkim(b1, b2)) {
          continue;
        }

        storePair(b1, b2);
      }
    }
  }

  void buildTriplets(std::vector<Particle> const& b1Candidates,
                     std::vector<Particle> const& b2Candidates,
                     std::vector<Particle> const& b3Candidates,
                     bool identicalB1B2 = false,
                     bool identicalB2B3 = false)
  {
    if (b1Candidates.empty() || b2Candidates.empty() || b3Candidates.empty()) {
      return;
    }

    if (identicalB1B2) {
      for (size_t iB1 = 0; iB1 < b1Candidates.size(); ++iB1) {
        for (size_t iB2 = iB1 + 1; iB2 < b2Candidates.size(); ++iB2) {
          for (const auto& b3 : b3Candidates) {
            const auto& b1 = b1Candidates[iB1];
            const auto& b2 = b2Candidates[iB2];

            if (!passThreeBodySkim(b1, b2, b3)) {
              continue;
            }

            storeTriplet(b1, b2, b3);
          }
        }
      }
      return;
    }

    if (identicalB2B3) {
      for (const auto& b1 : b1Candidates) {
        for (size_t iB2 = 0; iB2 < b2Candidates.size(); ++iB2) {
          for (size_t iB3 = iB2 + 1; iB3 < b3Candidates.size(); ++iB3) {
            const auto& b2 = b2Candidates[iB2];
            const auto& b3 = b3Candidates[iB3];

            if (!passThreeBodySkim(b1, b2, b3)) {
              continue;
            }

            storeTriplet(b1, b2, b3);
          }
        }
      }
      return;
    }

    for (const auto& b1 : b1Candidates) {
      for (const auto& b2 : b2Candidates) {
        for (const auto& b3 : b3Candidates) {
          if (!passThreeBodySkim(b1, b2, b3)) {
            continue;
          }

          storeTriplet(b1, b2, b3);
        }
      }
    }
  }

  template <typename McParticles>
  void selectGeneratedParticles(McParticles const& mcParticlesThisMcColl,
                                std::vector<Particle>& protons,
                                std::vector<Particle>& antiProtons,
                                std::vector<Particle>& neutrons,
                                std::vector<Particle>& antiNeutrons,
                                std::vector<Particle>& lambdas,
                                std::vector<Particle>& antiLambdas)
  {
    for (const auto& particle : mcParticlesThisMcColl) {
      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      if (std::fabs(particle.eta()) > etaMax) {
        continue;
      }

      const int pdg = particle.pdgCode();

      if (pdg == PDG_t::kProton) {
        protons.emplace_back(particle);
      } else if (pdg == PDG_t::kProtonBar) {
        antiProtons.emplace_back(particle);
      } else if (pdg == PDG_t::kNeutron) {
        neutrons.emplace_back(particle);
      } else if (pdg == PDG_t::kNeutronBar) {
        antiNeutrons.emplace_back(particle);
      } else if (pdg == PDG_t::kLambda0) {
        lambdas.emplace_back(particle);
      } else if (pdg == PDG_t::kLambda0Bar) {
        antiLambdas.emplace_back(particle);
      }
    }
  }

  void fillDeuteron(std::vector<Particle> const& protons,
                    std::vector<Particle> const& antiProtons,
                    std::vector<Particle> const& neutrons,
                    std::vector<Particle> const& antiNeutrons)
  {
    // Deuteron: p + n
    buildPairs(protons, neutrons);

    // anti-deuteron: anti-p + anti-n
    buildPairs(antiProtons, antiNeutrons);
  }

  void fillHelium3(std::vector<Particle> const& protons,
                   std::vector<Particle> const& antiProtons,
                   std::vector<Particle> const& neutrons,
                   std::vector<Particle> const& antiNeutrons)
  {
    // 3He: p + p + n
    buildTriplets(protons, protons, neutrons, true, false);

    // anti-3He: anti-p + anti-p + anti-n
    buildTriplets(antiProtons, antiProtons, antiNeutrons, true, false);
  }

  void fillTriton(std::vector<Particle> const& protons,
                  std::vector<Particle> const& antiProtons,
                  std::vector<Particle> const& neutrons,
                  std::vector<Particle> const& antiNeutrons)
  {
    // Triton: p + n + n
    buildTriplets(protons, neutrons, neutrons, false, true);

    // anti-triton: anti-p + anti-n + anti-n
    buildTriplets(antiProtons, antiNeutrons, antiNeutrons, false, true);
  }

  void fillHypertriton(std::vector<Particle> const& protons,
                       std::vector<Particle> const& antiProtons,
                       std::vector<Particle> const& neutrons,
                       std::vector<Particle> const& antiNeutrons,
                       std::vector<Particle> const& lambdas,
                       std::vector<Particle> const& antiLambdas)
  {
    // Hypertriton: p + n + Lambda
    buildTriplets(protons, neutrons, lambdas);

    // anti-hypertriton: anti-p + anti-n + anti-Lambda
    buildTriplets(antiProtons, antiNeutrons, antiLambdas);
  }

  bool hasDeuteronCandidates(std::vector<Particle> const& protons,
                             std::vector<Particle> const& antiProtons,
                             std::vector<Particle> const& neutrons,
                             std::vector<Particle> const& antiNeutrons) const
  {
    return (!protons.empty() && !neutrons.empty()) ||
           (!antiProtons.empty() && !antiNeutrons.empty());
  }

  bool hasHelium3Candidates(std::vector<Particle> const& protons,
                            std::vector<Particle> const& antiProtons,
                            std::vector<Particle> const& neutrons,
                            std::vector<Particle> const& antiNeutrons) const
  {
    constexpr std::size_t minimumSizeContainer = 2;
    return (protons.size() >= minimumSizeContainer && !neutrons.empty()) ||
           (antiProtons.size() >= minimumSizeContainer && !antiNeutrons.empty());
  }

  bool hasTritonCandidates(std::vector<Particle> const& protons,
                           std::vector<Particle> const& antiProtons,
                           std::vector<Particle> const& neutrons,
                           std::vector<Particle> const& antiNeutrons) const
  {
    constexpr std::size_t minimumSizeContainer = 2;
    return (!protons.empty() && neutrons.size() >= minimumSizeContainer) ||
           (!antiProtons.empty() && antiNeutrons.size() >= minimumSizeContainer);
  }

  bool hasHypertritonCandidates(std::vector<Particle> const& protons,
                                std::vector<Particle> const& antiProtons,
                                std::vector<Particle> const& neutrons,
                                std::vector<Particle> const& antiNeutrons,
                                std::vector<Particle> const& lambdas,
                                std::vector<Particle> const& antiLambdas) const
  {
    return (!protons.empty() && !neutrons.empty() && !lambdas.empty()) ||
           (!antiProtons.empty() && !antiNeutrons.empty() && !antiLambdas.empty());
  }

  void processCoalescence(aod::McCollision const& collision, aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("eventCounter"), 0.5);

    if (std::fabs(collision.posZ()) > zVtx) {
      return;
    }

    registry.fill(HIST("eventCounter"), 1.5);

    eventID = collision.globalIndex();

    std::vector<Particle> protons;
    std::vector<Particle> antiProtons;
    std::vector<Particle> neutrons;
    std::vector<Particle> antiNeutrons;
    std::vector<Particle> lambdas;
    std::vector<Particle> antiLambdas;

    const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

    selectGeneratedParticles(mcParticlesThisMcColl,
                             protons,
                             antiProtons,
                             neutrons,
                             antiNeutrons,
                             lambdas,
                             antiLambdas);

    switch (boundStateSpecies) {
      case kDeuteron:
        if (!hasDeuteronCandidates(protons, antiProtons, neutrons, antiNeutrons)) {
          return;
        }
        registry.fill(HIST("eventCounter"), 2.5);
        fillDeuteron(protons, antiProtons, neutrons, antiNeutrons);
        break;
      case kHelium3:
        if (!hasHelium3Candidates(protons, antiProtons, neutrons, antiNeutrons)) {
          return;
        }
        registry.fill(HIST("eventCounter"), 2.5);
        fillHelium3(protons, antiProtons, neutrons, antiNeutrons);
        break;
      case kTriton:
        if (!hasTritonCandidates(protons, antiProtons, neutrons, antiNeutrons)) {
          return;
        }
        registry.fill(HIST("eventCounter"), 2.5);
        fillTriton(protons, antiProtons, neutrons, antiNeutrons);
        break;
      case kHypertriton:
        if (!hasHypertritonCandidates(protons, antiProtons, neutrons, antiNeutrons, lambdas, antiLambdas)) {
          return;
        }
        registry.fill(HIST("eventCounter"), 2.5);
        fillHypertriton(protons, antiProtons, neutrons, antiNeutrons, lambdas, antiLambdas);
        break;
      default:
        LOGF(fatal, "Unknown boundStateSpecies=%d. Use 0 = deuteron, 1 = triton, 2 = helium3, 3 = hypertriton", static_cast<int>(boundStateSpecies));
    }
  }
  PROCESS_SWITCH(CoalescenceTreeProducer, processCoalescence, "process coalescence", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CoalescenceTreeProducer>(cfgc)};
}
