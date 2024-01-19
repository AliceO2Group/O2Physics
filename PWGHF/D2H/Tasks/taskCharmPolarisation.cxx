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

/// \file taskCharmPolarisation.cxx
/// \brief Analysis task for non-scalar charm hadron polarisation
///
/// \author F. Grosa (CERN) fabrizio.grosa@cern.ch
/// \author S. Kundu (CERN) sourav.kundu@cern.ch

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

// #include "Common/Core/EventPlaneHelper.h"
// #include "Common/DataModel/Qvectors.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace charm_polarisation
{
enum DecayChannel : uint8_t {
  DstarToDzeroPi = 0,
  LcToPKPi,
  LcToPK0S
};
}
} // namespace o2::aod

struct TaskPolarisationCharmHadrons {
  using CandDstarWSelFlag = soa::Join<aod::HfCandDstar, aod::HfSelDstarToD0Pi>;

  float massPi{0.f};
  float massProton{0.f};
  float massDstar{0.f};

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 Pi"};

  ConfigurableAxis thnAxisInvMass{"thnAxisInvMass", {100, 1.78, 2.05}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis thnAxisPt{"thnAxisPt", {100, 0., 100.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis thnAxisPz{"thnAxisPz", {100, 0., 100.}, "|#it{p}_{z}| (GeV/#it{c})"};
  ConfigurableAxis thnAxisY{"thnAxisY", {20, -1., 1.}, "#it{y}"};
  ConfigurableAxis thnAxisCosThetaStarHelicity{"thnAxisCosThetaStarHelicity", {20, -1., 1.}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis thnAxisCosThetaStarProduction{"thnAxisCosThetaStarProduction", {20, -1., 1.}, "cos(#vartheta_{production})"};
  ConfigurableAxis thnAxisCosThetaStarRandom{"thnAxisCosThetaStarRandom", {20, -1., 1.}, "cos(#vartheta_{random})"};
  ConfigurableAxis thnAxisCosThetaStarBeam{"thnAxisCosThetaStarBeam", {20, -1., 1.}, "cos(#vartheta_{beam})"};
  ConfigurableAxis thnConfigAxisMlBkg{"thnConfigAxisMlBkg", {100, 0., 1.}, "ML bkg"};
  ConfigurableAxis thnConfigAxisMlPrompt{"thnConfigAxisMlPrompt", {100, 0., 1.}, "ML prompt"};
  ConfigurableAxis thnConfigAxisMlNonPrompt{"thnConfigAxisMlNonPrompt", {100, 0., 1.}, "ML non-prompt"};
  // ConfigurableAxis thnAxisCent{"thnAxisCent", {102, -1., 101.}, "centrality (%)"};

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;

  HfHelper hfHelper;
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (doprocessDstar && doprocessDstarWithMl) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    }

    massPi = o2::constants::physics::MassPiPlus;
    massProton = o2::constants::physics::MassProton;
    massDstar = o2::constants::physics::MassDStar;

    if (doprocessDstarWithMl) {
      registry.add("hSparseCharmPolarisation", "THn for polarisation studies", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisCosThetaStarProduction, thnAxisCosThetaStarBeam, thnAxisCosThetaStarRandom, thnConfigAxisMlBkg, thnConfigAxisMlPrompt, thnConfigAxisMlNonPrompt});
    } else {
      registry.add("hSparseCharmPolarisation", "THn for polarisation studies", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisCosThetaStarProduction, thnAxisCosThetaStarBeam, thnAxisCosThetaStarRandom});
    }
  }; // end init

  /// \param candidates are the selected candidates
  /// \param tracks are the tracks
  template <uint8_t channel, bool withMl, typename Cand>
  void runPolarisationAnalysis(Cand const& candidate, Tracks const& tracks)
  {
    float pxDau{-1000.}, pyDau{-1000.}, pzDau{-1000.};
    float pxCharmHad{-1000.}, pyCharmHad{-1000.}, pzCharmHad{-1000.};
    float massDau{0.}, invMassCharmHad{0.}, invMassCharmHadForSparse{0.};
    float rapidity{-999.};
    std::array<float, 3> outputMl{};
    if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
      pxDau = candidate.pxSoftPi();
      pyDau = candidate.pySoftPi();
      pzDau = candidate.pzSoftPi();
      pxCharmHad = candidate.pxDstar();
      pyCharmHad = candidate.pyDstar();
      pzCharmHad = candidate.pzDstar();
      massDau = massPi;
      auto prongSoftPi = candidate.template prongPi_as<aod::Tracks>();
      invMassCharmHad = (prongSoftPi.sign() > 0) ? candidate.invMassDstar() : candidate.invMassAntiDstar();
      invMassCharmHadForSparse = (prongSoftPi.sign() > 0) ? (invMassCharmHad - candidate.invMassD0()) : (invMassCharmHad - candidate.invMassD0Bar()); // different for D*
      rapidity = candidate.y(massDstar);
      if constexpr (withMl) {
        outputMl[0] = -1.; // not yet implemented in the selector
        outputMl[1] = -1.; // not yet implemented in the selector
        outputMl[2] = -1.; // not yet implemented in the selector
      }
    }

    float phiRandom = gRandom->Uniform(0., constants::math::TwoPI);
    float thetaRandom = gRandom->Uniform(0., constants::math::PI);
    ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
    ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, invMassCharmHad);
    ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
    ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
    ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

    ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
    ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0., 0., 1.);
    ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
    ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.);

    float cosThetaStarHelicity = std::abs(helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2()));
    float cosThetaStarProduction = std::abs(normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2()));
    float cosThetaStarBeam = std::abs(beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()));
    float cosThetaStarRandom = std::abs(randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()));

    if constexpr (withMl) {
      registry.fill(HIST("hSparseCharmPolarisation"), invMassCharmHadForSparse, candidate.pt(), pzCharmHad, rapidity, cosThetaStarHelicity, cosThetaStarProduction, cosThetaStarBeam, cosThetaStarRandom, outputMl[0], outputMl[1], outputMl[2]);
    } else {
      registry.fill(HIST("hSparseCharmPolarisation"), invMassCharmHadForSparse, candidate.pt(), pzCharmHad, rapidity, cosThetaStarHelicity, cosThetaStarProduction, cosThetaStarBeam, cosThetaStarRandom);
    }
  }

  // Dstar with rectangular cuts
  void processDstar(soa::Filtered<CandDstarWSelFlag>::iterator const& dstarCandidate, Tracks const& tracks)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false>(dstarCandidate, tracks);
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with rectangular cuts
  void processDstarWithMl(soa::Filtered<CandDstarWSelFlag>::iterator const&, Tracks const&)
  {
    // DUMMY
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarWithMl, "Process Dstar candidates with ML (DUMMY)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPolarisationCharmHadrons>(cfgc)};
}
