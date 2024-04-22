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
enum CosThetaStarType : uint8_t {
  Helicity = 0,
  Production,
  Beam,
  Random,
  NTypes
};
enum DecayChannel : uint8_t {
  DstarToDzeroPi = 0,
  LcToPKPi,
  LcToPK0S,
  NChannels
};
enum MassHyposLcToPKPi : uint8_t {
  PKPi = 0,
  PiKP,
  NMassHypoLcToPKPi
};
} // namespace charm_polarisation
} // namespace o2::aod

struct TaskPolarisationCharmHadrons {
  using CandDstarWSelFlag = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  using CandLcToPKPiWSelFlag = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;

  float massPi{0.f};
  float massProton{0.f};
  float massKaon{0.f};
  float massDstar{0.f};
  float massLc{0.f};
  float bkgRotationAngleStep{0.f};

  uint8_t nMassHypos{0u};

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 Pi"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc decay to P K Pi"};

  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {200, 0.139f, 0.179f}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.f, 100.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPz{"configThnAxisPz", {100, -50.f, 50.f}, "#it{p}_{z} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisY{"configThnAxisY", {20, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis configThnAxisCosThetaStarHelicity{"configThnAxisCosThetaStarHelicity", {20, -1.f, 1.f}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis configThnAxisCosThetaStarProduction{"configThnAxisCosThetaStarProduction", {20, -1.f, 1.f}, "cos(#vartheta_{production})"};
  ConfigurableAxis configThnAxisCosThetaStarRandom{"configThnAxisCosThetaStarRandom", {20, -1.f, 1.f}, "cos(#vartheta_{random})"};
  ConfigurableAxis configThnAxisCosThetaStarBeam{"configThnAxisCosThetaStarBeam", {20, -1.f, 1.f}, "cos(#vartheta_{beam})"};
  ConfigurableAxis configThnAxisMlBkg{"configThnAxisMlBkg", {100, 0.f, 1.f}, "ML bkg"};
  ConfigurableAxis configThnAxisInvMassD0{"configThnAxisInvMassD0", {250, 1.65f, 2.15f}, "#it{M}(D^{0}) (GeV/#it{c}^{2})"}; // only for D*+
  // ConfigurableAxis configThnAxisMlPrompt{"configThnAxisMlPrompt", {100, 0.f, 1.f}, "ML prompt"};
  ConfigurableAxis configThnAxisMlNonPrompt{"configThnAxisMlNonPrompt", {100, 0.f, 1.f}, "ML non-prompt"};
  // ConfigurableAxis configThnAxisCent{"configThnAxisCent", {102, -1.f, 101.f}, "centrality (%)"};

  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 0, "Number of rotated copies (background) per each original candidate"};

  /// output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", true, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", true, "Activate the THnSparse with cosThStar w.r.t. beam axis"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", true, "Activate the THnSparse with cosThStar w.r.t. random axis"};

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Filter filterSelectLcToPKPiCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi) || (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPKPi);

  HfHelper hfHelper;
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    /// check process functions
    std::array<int, 8> processes = {doprocessDstar, doprocessDstarWithMl, doprocessLcToPKPi, doprocessLcToPKPiWithMl, doprocessDstarMc, doprocessDstarMcWithMl, doprocessLcToPKPiMc, doprocessLcToPKPiMcWithMl};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    /// check output THnSparses
    std::array<int, 4> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarProduction, activateTHnSparseCosThStarBeam, activateTHnSparseCosThStarRandom};
    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (activateTHnSparseCosThStarHelicity) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (activateTHnSparseCosThStarProduction) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (activateTHnSparseCosThStarBeam) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active.");
      }
      if (activateTHnSparseCosThStarRandom) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }

    // check bkg rotation for MC (not supported currently)
    if (nBkgRotations > 0 && (doprocessDstarMc || doprocessDstarMcWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl)) {
      LOGP(fatal, "No background rotation supported for MC.");
    }

    massPi = o2::constants::physics::MassPiPlus;
    massProton = o2::constants::physics::MassProton;
    massKaon = o2::constants::physics::MassKaonCharged;
    massDstar = o2::constants::physics::MassDStar;
    massLc = o2::constants::physics::MassLambdaCPlus;
    bkgRotationAngleStep = constants::math::TwoPI / (nBkgRotations + 1); // nBkgRotations==0: 2π (no rotation); nBkgRotations==1: π; nBkgRotations==2: 2π/3, 4π/3; ...

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassD0{configThnAxisInvMassD0, "#it{M}(D^{0}) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisPz{configThnAxisPz, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec thnAxisY{configThnAxisY, "#it{y}"};
    const AxisSpec thnAxisCosThetaStarHelicity{configThnAxisCosThetaStarHelicity, "cos(#vartheta_{helicity})"};
    const AxisSpec thnAxisCosThetaStarProduction{configThnAxisCosThetaStarProduction, "cos(#vartheta_{production})"};
    const AxisSpec thnAxisCosThetaStarRandom{configThnAxisCosThetaStarRandom, "cos(#vartheta_{random})"};
    const AxisSpec thnAxisCosThetaStarBeam{configThnAxisCosThetaStarBeam, "cos(#vartheta_{beam})"};
    const AxisSpec thnAxisMlBkg{configThnAxisMlBkg, "ML bkg"};
    const AxisSpec thnAxisMlNonPrompt{configThnAxisMlNonPrompt, "ML non-prompt"};
    const AxisSpec thnAxisIsRotatedCandidate{2, -0.5f, 1.5f, "rotated bkg"};

    if (doprocessDstarWithMl || doprocessDstarMcWithMl) {
      /// analysis for D*+ meson with ML, w/o rot. background axis
      if (doprocessDstarWithMl) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
      }
    } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl) {
      /// analysis for Lc+ baryon with ML, w/ rot. background axis (for data only)
      if (doprocessLcToPKPiWithMl) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt});
        }
      }
    } else if (doprocessDstar || doprocessDstarMc) {
      /// analysis for D*+ meson, w/o rot. background axis
      if (doprocessDstar) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom});
        }
      }
    } else if (doprocessLcToPKPi || doprocessLcToPKPiMc) {
      /// analysis for Lc+ baryon, rot. background axis (for data only)
      if (doprocessLcToPKPi) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom});
        }
      }
    }

    // MC Gen histos
    if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl) {
      if (activateTHnSparseCosThStarHelicity) {
        registry.add("hGenPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity});
        registry.add("hGenNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarHelicity});
      }
      if (activateTHnSparseCosThStarProduction) {
        registry.add("hGenPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction});
        registry.add("hGenNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarProduction});
      }
      if (activateTHnSparseCosThStarBeam) {
        registry.add("hGenPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam});
        registry.add("hGenNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarBeam});
      }
      if (activateTHnSparseCosThStarRandom) {
        registry.add("hGenPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom});
        registry.add("hGenNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisPz, thnAxisY, thnAxisCosThetaStarRandom});
      }
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

  /// \param invMassCharmHad is the invariant-mass of the candidate
  /// \param ptCharmHad is the pt of the candidate
  /// \param pzCharmHad is the pz of the candidate
  /// \param rapCharmHad is the rapidity of the candidate
  /// \param invMassD0 is the invariant-mass of the D0 daugher (only for D*+)
  /// \param cosThetaStar is the cosThetaStar of the candidate
  /// \param outputMl is the array with ML output scores
  /// \param isRotatedCandidate is a flag that keeps the info of the rotation of the candidate for bkg studies
  /// \param origin is the MC origin
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, charm_polarisation::CosThetaStarType cosThetaStarType>
  void fillRecoHistos(float invMassCharmHad, float ptCharmHad, float pzCharmHad, float rapCharmHad, float invMassD0, float cosThetaStar, std::array<float, 3> outputMl, int isRotatedCandidate, int8_t origin)
  {
    if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if constexpr (!doMc) {                                                            // data
        if constexpr (withMl) {                                                         // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {  // D*+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if constexpr (!doMc) {                                                                     // data
        if constexpr (withMl) {                                                                  // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {           // D*+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if constexpr (!doMc) {                                                               // data
        if constexpr (withMl) {                                                            // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {     // D*+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if constexpr (!doMc) {                                                                 // data
        if constexpr (withMl) {                                                              // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {       // D*+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, invMassD0, cosThetaStar);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
            }
          }
        }
      }
    }
  }

  /// \param ptCharmHad is the pt of the particle
  /// \param pzCharmHad is the pz of the particle
  /// \param rapCharmHad is the rapidity of the particle
  /// \param cosThetaStar is the cosThetaStar of the particle
  /// \param origin is the MC origin
  template <charm_polarisation::CosThetaStarType cosThetaStarType>
  void fillGenHistos(float ptCharmHad, float pzCharmHad, float rapCharmHad, float cosThetaStar, int8_t origin)
  {
    if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if (origin == RecoDecay::OriginType::Prompt) {                                    // prompt
        registry.fill(HIST("hGenPromptHelicity"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptHelicity"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if (origin == RecoDecay::OriginType::Prompt) {                                             // prompt
        registry.fill(HIST("hGenPromptProduction"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptProduction"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if (origin == RecoDecay::OriginType::Prompt) {                                       // prompt
        registry.fill(HIST("hGenPromptBeam"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptBeam"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if (origin == RecoDecay::OriginType::Prompt) {                                         // prompt
        registry.fill(HIST("hGenPromptRandom"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptRandom"), ptCharmHad, pzCharmHad, rapCharmHad, cosThetaStar);
      }
    }
  }

  /// \param candidates are the selected candidates
  /// \param bkgRotationId is the id for the background rotation
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, typename Cand>
  void runPolarisationAnalysis(Cand const& candidate, int bkgRotationId = 0)
  {
    int8_t origin{RecoDecay::OriginType::None};
    int8_t massHypoMcTruth{-1};
    if constexpr (doMc) {
      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        if (!TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // this candidate is not signal, skip
          return;
        }
        origin = candidate.originMcRec();
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if (!TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_3prong::DecayType::LcToPKPi)) { // this candidate is not signal, skip
          return;
        }
        origin = candidate.originMcRec();
        if (candidate.isCandidateSwapped()) {
          massHypoMcTruth = charm_polarisation::MassHyposLcToPKPi::PiKP;
        } else {
          massHypoMcTruth = charm_polarisation::MassHyposLcToPKPi::PKPi;
        }
      }
    }

    // loop over mass hypotheses
    for (uint8_t iMass = 0u; iMass < nMassHypos; iMass++) {

      // variable definition
      float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
      float pxCharmHad{-1000.f}, pyCharmHad{-1000.f}, pzCharmHad{-1000.f};
      float massDau{0.f}, invMassCharmHad{0.f}, invMassCharmHadForSparse{0.f}, invMassD0{0.f};
      float rapidity{-999.f};
      std::array<float, 3> outputMl{-1.f, -1.f, -1.f};
      int isRotatedCandidate = 0; // currently meaningful only for Lc->pKpi

      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        // Dstar analysis
        // polarization measured from the soft-pion daughter (*)

        massDau = massPi; // (*)
        const float bkgRotAngle = bkgRotationAngleStep * bkgRotationId;
        std::array<float, 3> threeVecSoftPi{candidate.pxSoftPi() * std::cos(bkgRotAngle) - candidate.pySoftPi() * std::sin(bkgRotAngle), candidate.pxSoftPi() * std::sin(bkgRotAngle) + candidate.pySoftPi() * std::cos(bkgRotAngle), candidate.pzSoftPi()}; // we rotate the soft pion
        std::array<float, 3> threeVecD0Prong0{candidate.pxProng0(), candidate.pxProng0(), candidate.pzProng0()};
        std::array<float, 3> threeVecD0Prong1{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()};
        if (bkgRotationId > 0) {
          isRotatedCandidate = 1;
          pxDau = threeVecSoftPi[0];
          pyDau = threeVecSoftPi[1];
          pzDau = threeVecSoftPi[2];
          pxCharmHad = threeVecSoftPi[0] + threeVecD0Prong0[0] + threeVecD0Prong1[0];
          pyCharmHad = threeVecSoftPi[1] + threeVecD0Prong0[1] + threeVecD0Prong1[1];
          pzCharmHad = threeVecSoftPi[2] + threeVecD0Prong0[2] + threeVecD0Prong1[2];
          if (candidate.signSoftPi() > 0) {
            invMassCharmHad = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1, threeVecSoftPi}, std::array{massPi, massKaon, massPi});
            invMassD0 = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1}, std::array{massPi, massKaon});
          } else {
            invMassCharmHad = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1, threeVecSoftPi}, std::array{massKaon, massPi, massPi});
            invMassD0 = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1}, std::array{massKaon, massPi});
          }
          rapidity = RecoDecay::y(std::array{pxCharmHad, pyCharmHad, pzCharmHad}, massDstar);
        } else {
          isRotatedCandidate = 0;
          pxDau = candidate.pxSoftPi();
          pyDau = candidate.pySoftPi();
          pzDau = candidate.pzSoftPi();
          pxCharmHad = candidate.pxDstar();
          pyCharmHad = candidate.pyDstar();
          pzCharmHad = candidate.pzDstar();
          if (candidate.signSoftPi() > 0) {
            invMassCharmHad = candidate.invMassDstar();
            invMassD0 = candidate.invMassD0();
          } else {
            invMassCharmHad = candidate.invMassAntiDstar();
            invMassD0 = candidate.invMassD0Bar();
          }
          rapidity = candidate.y(massDstar);
        }
        invMassCharmHadForSparse = invMassCharmHad - invMassD0;

        if constexpr (withMl) {
          outputMl[0] = candidate.mlProbDstarToD0Pi()[0];
          outputMl[1] = candidate.mlProbDstarToD0Pi()[1];
          outputMl[2] = candidate.mlProbDstarToD0Pi()[2];
        }
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        // Lc->pKpi analysis
        // polarization measured from the proton daughter (*)

        if constexpr (doMc) { // we keep only the good hypo in the MC
          if ((iMass == charm_polarisation::MassHyposLcToPKPi::PiKP && massHypoMcTruth == charm_polarisation::MassHyposLcToPKPi::PKPi) || (iMass == charm_polarisation::MassHyposLcToPKPi::PKPi && massHypoMcTruth == charm_polarisation::MassHyposLcToPKPi::PiKP)) {
            continue;
          }
        }

        /// mass-hypothesis-independent variables
        /// daughters momenta
        const float bkgRotAngle = bkgRotationAngleStep * bkgRotationId;
        std::array<float, 3> threeVecLcProng0{candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()};
        std::array<float, 3> threeVecLcRotatedProng1{candidate.pxProng1() * std::cos(bkgRotAngle) - candidate.pyProng1() * std::sin(bkgRotAngle), candidate.pxProng1() * std::sin(bkgRotAngle) + candidate.pyProng1() * std::cos(bkgRotAngle), candidate.pzProng1()};
        std::array<float, 3> threeVecLcProng2{candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()};
        if (bkgRotationId > 0) {
          /// rotational background - pt of the kaon track rotated
          /// update candidate momentum
          isRotatedCandidate = 1;
          pxCharmHad = threeVecLcProng0[0] + threeVecLcRotatedProng1[0] + threeVecLcProng2[0];
          pyCharmHad = threeVecLcProng0[1] + threeVecLcRotatedProng1[1] + threeVecLcProng2[1];
          pzCharmHad = threeVecLcProng0[2] + threeVecLcRotatedProng1[2] + threeVecLcProng2[2];
        } else {
          /// original candidate (kaon track not rotated)
          isRotatedCandidate = 0;
          pxCharmHad = candidate.px();
          pyCharmHad = candidate.py();
          pzCharmHad = candidate.pz();
        }
        massDau = massProton; // (*)
        rapidity = RecoDecay::y(candidate.pVector(), massLc);

        /// mass-hypothesis-dependent variables
        if (iMass == charm_polarisation::MassHyposLcToPKPi::PKPi && candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
          // reconstructed as pKpi
          pxDau = candidate.pxProng0();
          pyDau = candidate.pyProng0();
          pzDau = candidate.pzProng0();
          if (bkgRotationId) {
            /// rotational background - pt of the kaon track rotated
            invMassCharmHad = RecoDecay::m(std::array{threeVecLcProng0, threeVecLcRotatedProng1, threeVecLcProng2}, std::array{massProton, massKaon, massPi});
            invMassCharmHadForSparse = invMassCharmHad;
          } else {
            /// original candidate (kaon track not rotated)
            invMassCharmHad = hfHelper.invMassLcToPKPi(candidate);
            invMassCharmHadForSparse = hfHelper.invMassLcToPKPi(candidate);
          }
          if constexpr (withMl) {
            if (candidate.mlProbLcToPKPi().size() == 3) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              // !!! NB: each rotated candidates inherits the BDT scores of the original candidate, even if the candidate pt changed after the rotation of the kaon-track pt !!!
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
          if (bkgRotationId) {
            /// rotational background - pt of the kaon track rotated
            invMassCharmHad = RecoDecay::m(std::array{threeVecLcProng0, threeVecLcRotatedProng1, threeVecLcProng2}, std::array{massPi, massKaon, massProton});
            invMassCharmHadForSparse = invMassCharmHad;
          } else {
            /// original candidate (kaon track not rotated)
            invMassCharmHad = hfHelper.invMassLcToPiKP(candidate);
            invMassCharmHadForSparse = hfHelper.invMassLcToPiKP(candidate);
          }
          if constexpr (withMl) {
            if (candidate.mlProbLcToPiKP().size() == 3) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              // !!! NB: each rotated candidates inherits the BDT scores of the original candidate, even if the candidate pt changed after the rotation of the kaon-track pt !!!
              outputMl[0] = candidate.mlProbLcToPiKP()[0];
              outputMl[1] = candidate.mlProbLcToPiKP()[1];
              outputMl[2] = candidate.mlProbLcToPiKP()[2];
            }
          }
        } else {
          // NB: no need to check cases in which candidate.isSelLcToPKPi() and candidate.isSelLcToPiKP() are both false, because they are rejected already by the Filter
          // ... but we need to put this protections here!
          // Otherwise, a candidate selected as pKpi only has invMassCharmHad==0 when iMass == charm_polarisation::MassHyposLcToPKPi::PiKP and viceversa
          continue;
        }

      } // Lc->pKpi

      float phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      float thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
      ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
      ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, invMassCharmHad);
      ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
      ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
      ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

      float ptCharmHad = std::sqrt(pxCharmHad * pxCharmHad + pyCharmHad * pyCharmHad); // this definition is valid for both rotated and original candidates

      if (activateTHnSparseCosThStarHelicity) {
        ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
        float cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Helicity>(invMassCharmHadForSparse, ptCharmHad, pzCharmHad, rapidity, invMassD0, cosThetaStarHelicity, outputMl, isRotatedCandidate, origin);
      }
      if (activateTHnSparseCosThStarProduction) {
        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
        float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Production>(invMassCharmHadForSparse, ptCharmHad, pzCharmHad, rapidity, invMassD0, cosThetaStarProduction, outputMl, isRotatedCandidate, origin);
      }
      if (activateTHnSparseCosThStarBeam) {
        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Beam>(invMassCharmHadForSparse, ptCharmHad, pzCharmHad, rapidity, invMassD0, cosThetaStarBeam, outputMl, isRotatedCandidate, origin);
      }
      if (activateTHnSparseCosThStarRandom) {
        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Random>(invMassCharmHadForSparse, ptCharmHad, pzCharmHad, rapidity, invMassD0, cosThetaStarRandom, outputMl, isRotatedCandidate, origin);
      }
    } /// end loop over mass hypotheses
  }

  /// \param mcParticle is the Mc particle
  /// \param mcParticles is the table of Mc particles
  template <charm_polarisation::DecayChannel channel, typename Part, typename Particles>
  void runMcGenPolarisationAnalysis(Part const& mcParticle, Particles const& mcParticles)
  {
    int8_t origin{RecoDecay::OriginType::None};
    std::vector<int> listDaughters{};
    float massDau{0.f}, massCharmHad{0.f};
    if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
      if (!TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // this particle is not signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      std::array<int, 2> dauPdgs = {kPiPlus, o2::constants::physics::Pdg::kD0};
      RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 1);
      massDau = massPi;
      massCharmHad = massDstar;
    } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
      if (!TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_3prong::DecayType::LcToPKPi)) { // this particle is not signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      std::array<int, 3> dauPdgs = {kProton, -kKPlus, kPiPlus};
      RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 2);
      massDau = massProton;
      massCharmHad = massLc;
    }

    float rapidity = mcParticle.y();
    if (std::abs(rapidity) > 1.f) { // we do not keep particles with |y| > 1
      return;
    }

    float pxCharmHad = mcParticle.px();
    float pyCharmHad = mcParticle.py();
    float pzCharmHad = mcParticle.pz();
    float ptCharmHad = mcParticle.pt();

    float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
    for (const auto& dauIdx : listDaughters) {
      auto dauPart = mcParticles.rawIteratorAt(dauIdx);
      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        if (std::abs(dauPart.pdgCode()) == kPiPlus) {
          pxDau = dauPart.px();
          pyDau = dauPart.py();
          pzDau = dauPart.pz();
          break;
        }
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if (std::abs(dauPart.pdgCode()) == kProton) {
          pxDau = dauPart.px();
          pyDau = dauPart.py();
          pzDau = dauPart.pz();
          break;
        }
      }
    }

    float phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
    float thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
    ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
    ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, massCharmHad);
    ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
    ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
    ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

    if (activateTHnSparseCosThStarHelicity) {
      ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
      float cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Helicity>(ptCharmHad, pzCharmHad, rapidity, cosThetaStarHelicity, origin);
    }
    if (activateTHnSparseCosThStarProduction) {
      ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
      float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Production>(ptCharmHad, pzCharmHad, rapidity, cosThetaStarProduction, origin);
    }
    if (activateTHnSparseCosThStarBeam) {
      ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Beam>(ptCharmHad, pzCharmHad, rapidity, cosThetaStarBeam, origin);
    }
    if (activateTHnSparseCosThStarRandom) {
      ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Random>(ptCharmHad, pzCharmHad, rapidity, cosThetaStarRandom, origin);
    }
  }

  /////////////////////////
  //   Dstar analysis   ///
  /////////////////////////

  // Dstar with rectangular cuts
  void processDstar(soa::Filtered<CandDstarWSelFlag>::iterator const& dstarCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(dstarCandidate);

    for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(dstarCandidate, iRotation);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with ML cuts
  void processDstarWithMl(soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi>>::iterator const& dstarCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate);

    /// rotational background
    for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate, iRotation);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarWithMl, "Process Dstar candidates with ML", false);

  // Dstar in MC with rectangular cuts
  void processDstarMc(soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfCandDstarMcRec>> const& dstarCandidates, soa::Join<aod::McParticles, aod::HfCandDstarMcGen> const& mcParticles)
  {
    for (const auto& dstarCandidate : dstarCandidates) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, true>(dstarCandidate);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi>(mcParticle, mcParticles);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarMc, "Process Dstar candidates in MC without ML", false);

  // Dstar in MC with ML cuts
  void processDstarMcWithMl(soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfCandDstarMcRec, aod::HfMlDstarToD0Pi>> const& dstarCandidates, soa::Join<aod::McParticles, aod::HfCandDstarMcGen> const& mcParticles)
  {
    for (const auto& dstarCandidate : dstarCandidates) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, true>(dstarCandidate);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi>(mcParticle, mcParticles);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarMcWithMl, "Process Dstar candidates in MC with ML", false);

  ////////////////////////////
  //   Lc->pKpi analysis   ///
  ////////////////////////////

  // Lc->pKpi with rectangular cuts
  void processLcToPKPi(soa::Filtered<CandLcToPKPiWSelFlag>::iterator const& lcCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, false>(lcCandidate);

    /// rotational background
    for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, false>(lcCandidate, iRotation);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPi, "Process Lc candidates without ML", false);

  // Lc->pKpi with ML cuts
  void processLcToPKPiWithMl(soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfMlLcToPKPi>>::iterator const& lcCandidate)
  {
    runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate);

    /// rotational background
    for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate, iRotation);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiWithMl, "Process Lc candidates with ML", false);

  // Lc->pKpi in MC with rectangular cuts
  void processLcToPKPiMc(soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfCand3ProngMcRec>> const& lcCandidates, soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles)
  {
    for (const auto& lcCandidate : lcCandidates) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, true>(lcCandidate);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiMc, "Process Lc candidates in MC without ML", false);

  // Lc->pKpi in MC with ML cuts
  void processLcToPKPiMcWithMl(soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfMlLcToPKPi, aod::HfCand3ProngMcRec>> const& lcCandidates, soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles)
  {
    for (const auto& lcCandidate : lcCandidates) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, true>(lcCandidate);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiMcWithMl, "Process Lc candidates in MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPolarisationCharmHadrons>(cfgc)};
}
