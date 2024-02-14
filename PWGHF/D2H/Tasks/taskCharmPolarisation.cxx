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
/// \author M. Faggin (CERN) mattia.faggin@cern.ch

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
enum MassHyposLcToPKPi : uint8_t {
  PKPi = 0,
  PiKP,
  NMassHypoLcToPKPi
};
} // namespace charm_polarisation
} // namespace o2::aod

struct TaskPolarisationCharmHadrons {
  using CandDstarWSelFlag = soa::Join<aod::HfCandDstar, aod::HfSelDstarToD0Pi>;
  using CandLcToPKPiWSelFlag = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;

  float massPi{0.f};
  float massProton{0.f};
  float massDstar{0.f};
  float massLc{0.f};

  uint8_t nMassHypos{0u};

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 Pi"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc decay to P K Pi"};

  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {200, 0.139, 0.179}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0., 100.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPz{"configThnAxisPz", {100, -50., 50.}, "#it{p}_{z} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisY{"configThnAxisY", {20, -1., 1.}, "#it{y}"};
  ConfigurableAxis configThnAxisCosThetaStarHelicity{"configThnAxisCosThetaStarHelicity", {20, -1., 1.}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis configThnAxisCosThetaStarProduction{"configThnAxisCosThetaStarProduction", {20, -1., 1.}, "cos(#vartheta_{production})"};
  ConfigurableAxis configThnAxisCosThetaStarRandom{"configThnAxisCosThetaStarRandom", {20, -1., 1.}, "cos(#vartheta_{random})"};
  ConfigurableAxis configThnAxisCosThetaStarBeam{"configThnAxisCosThetaStarBeam", {20, -1., 1.}, "cos(#vartheta_{beam})"};
  ConfigurableAxis configThnAxisMlBkg{"configThnAxisMlBkg", {100, 0., 1.}, "ML bkg"};
  // ConfigurableAxis configThnAxisMlPrompt{"configThnAxisMlPrompt", {100, 0., 1.}, "ML prompt"};
  ConfigurableAxis configThnAxisMlNonPrompt{"configThnAxisMlNonPrompt", {100, 0., 1.}, "ML non-prompt"};
  // ConfigurableAxis configThnAxisCent{"configThnAxisCent", {102, -1., 101.}, "centrality (%)"};

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Filter filterSelectLcToPKPiCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi) || (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPKPi);

  HfHelper hfHelper;
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<int, 4> processes = {doprocessDstar, doprocessDstarWithMl, doprocessLcToPKPi, doprocessLcToPKPiWithMl};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    massPi = o2::constants::physics::MassPiPlus;
    massProton = o2::constants::physics::MassProton;
    massDstar = o2::constants::physics::MassDStar;
    massLc = o2::constants::physics::MassLambdaCPlus;

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPz{configThnAxisPz, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec thnAxisY{configThnAxisY, "#it{y}"};
    const AxisSpec thnAxisCosThetaStarHelicity{configThnAxisCosThetaStarHelicity, "cos(#vartheta_{helicity})"};
    const AxisSpec thnAxisCosThetaStarProduction{configThnAxisCosThetaStarProduction, "cos(#vartheta_{production})"};
    const AxisSpec thnAxisCosThetaStarRandom{configThnAxisCosThetaStarRandom, "cos(#vartheta_{random})"};
    const AxisSpec thnAxisCosThetaStarBeam{configThnAxisCosThetaStarBeam, "cos(#vartheta_{beam})"};
    const AxisSpec thnAxisMlBkg{configThnAxisMlBkg, "ML bkg"};
    const AxisSpec thnAxisMlNonPrompt{configThnAxisMlNonPrompt, "ML non-prompt"};

    if (doprocessDstarWithMl || doprocessLcToPKPiWithMl) {
      registry.add("hSparseCharmPolarisation", "THn for polarisation studies", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisCosThetaStarProduction, thnAxisCosThetaStarBeam, thnAxisCosThetaStarRandom, thnAxisMlBkg, /*thnConfigAxisMlPrompt,*/ thnAxisMlNonPrompt});
    } else {
      registry.add("hSparseCharmPolarisation", "THn for polarisation studies", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisCosThetaStarProduction, thnAxisCosThetaStarBeam, thnAxisCosThetaStarRandom});
    }

    // inv. mass hypothesis to loop over
    // e.g.: Lc->pKpi has the ambiguity pKpi vs. piKp
    if (doprocessLcToPKPi || doprocessLcToPKPiWithMl) {
      nMassHypos = charm_polarisation::MassHyposLcToPKPi::NMassHypoLcToPKPi;
    } else {
      // D*, Lc->pK0s
      nMassHypos = 1;
    }

  }; // end init

  /// \param candidates are the selected candidates
  template <uint8_t channel, bool withMl, typename Cand>
  void runPolarisationAnalysis(Cand const& candidate)
  {

    // loop over mass hypotheses
    for (uint8_t iMass = 0u; iMass < nMassHypos; iMass++) {

      // variable definition
      float pxDau{-1000.}, pyDau{-1000.}, pzDau{-1000.};
      float pxCharmHad{-1000.}, pyCharmHad{-1000.}, pzCharmHad{-1000.};
      float massDau{0.}, invMassCharmHad{0.}, invMassCharmHadForSparse{0.};
      float rapidity{-999.};
      std::array<float, 3> outputMl{-1., -1., -1.};

      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        // Dstar analysis
        // polarization measured from the soft-pion daughter (*)

        pxDau = candidate.pxSoftPi();
        pyDau = candidate.pySoftPi();
        pzDau = candidate.pzSoftPi();
        pxCharmHad = candidate.pxDstar();
        pyCharmHad = candidate.pyDstar();
        pzCharmHad = candidate.pzDstar();
        massDau = massPi; // (*)
        invMassCharmHad = (candidate.signSoftPi() > 0) ? candidate.invMassDstar() : candidate.invMassAntiDstar();
        invMassCharmHadForSparse = (candidate.signSoftPi() > 0) ? (invMassCharmHad - candidate.invMassD0()) : (invMassCharmHad - candidate.invMassD0Bar()); // different for D*
        rapidity = candidate.y(massDstar);
        if constexpr (withMl) {
          outputMl[0] = -1.; // not yet implemented in the selector
          outputMl[1] = -1.; // not yet implemented in the selector
          outputMl[2] = -1.; // not yet implemented in the selector
        }
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        // Lc->pKpi analysis
        // polarization measured from the proton daughter (*)

        if (iMass == charm_polarisation::MassHyposLcToPKPi::PKPi && candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
          // reconstructed as pKpi
          pxDau = candidate.pxProng0();
          pyDau = candidate.pyProng0();
          pzDau = candidate.pzProng0();
          invMassCharmHadForSparse = hfHelper.invMassLcToPKPi(candidate);
          if constexpr (withMl) {
            if (candidate.mlProbLcToPKPi().size() == 3) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              outputMl[0] = candidate.mlProbLcToPKPi()[0];
              outputMl[1] = candidate.mlProbLcToPKPi()[1];
              outputMl[2] = candidate.mlProbLcToPKPi()[2];
            }
          }
        } else if (iMass == charm_polarisation::MassHyposLcToPKPi::PiKP && candidate.isSelLcToPiKP() >= selectionFlagLcToPKPi) {
          // reconstructed as piKp
          pxDau = candidate.pxProng2();
          pyDau = candidate.pyProng2();
          pzDau = candidate.pzProng2();
          invMassCharmHadForSparse = hfHelper.invMassLcToPiKP(candidate);
          if constexpr (withMl) {
            if (candidate.mlProbLcToPiKP().size() == 3) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              outputMl[0] = candidate.mlProbLcToPiKP()[0];
              outputMl[1] = candidate.mlProbLcToPiKP()[1];
              outputMl[2] = candidate.mlProbLcToPiKP()[2];
            }
          }
        }

        // NB: no need to check cases in which candidate.isSelLcToPKPi() and candidate.isSelLcToPiKP() are both false, because they are rejected already by the Filter

        // hypothesis-independent variables
        pxCharmHad = candidate.px();
        pyCharmHad = candidate.py();
        pzCharmHad = candidate.pz();
        massDau = massProton; // (*)
        rapidity = candidate.y(massLc);

      } // Lc->pKpi

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

      float cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
      float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
      float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

      if constexpr (withMl) {
        registry.fill(HIST("hSparseCharmPolarisation"), invMassCharmHadForSparse, candidate.pt(), pzCharmHad, rapidity, cosThetaStarHelicity, cosThetaStarProduction, cosThetaStarBeam, cosThetaStarRandom, outputMl[0], /*outputMl[1],*/ outputMl[2]);
      } else {
        registry.fill(HIST("hSparseCharmPolarisation"), invMassCharmHadForSparse, candidate.pt(), pzCharmHad, rapidity, cosThetaStarHelicity, cosThetaStarProduction, cosThetaStarBeam, cosThetaStarRandom);
      }
    } /// end loop over mass hypotheses
  }

  /////////////////////////
  //   Dstar analysis   ///
  /////////////////////////

  // Dstar with rectangular cuts
  void processDstar(soa::Filtered<CandDstarWSelFlag>::iterator const& dstarCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false>(dstarCandidate);
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with ML cuts
  void processDstarWithMl(soa::Filtered<CandDstarWSelFlag>::iterator const&)
  {
    // DUMMY
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarWithMl, "Process Dstar candidates with ML (DUMMY)", false);

  ////////////////////////////
  //   Lc->pKpi analysis   ///
  ////////////////////////////

  // Lc->pKpi with rectangular cuts
  void processLcToPKPi(soa::Filtered<CandLcToPKPiWSelFlag>::iterator const& lcCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false>(lcCandidate);
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPi, "Process Lc candidates without ML", false);

  // Lc->pKpi with ML cuts
  void processLcToPKPiWithMl(soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfMlLcToPKPi>>::iterator const& lcCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true>(lcCandidate);
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiWithMl, "Process Lc candidates with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPolarisationCharmHadrons>(cfgc)};
}
