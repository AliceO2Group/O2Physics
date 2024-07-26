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

/// columns for table to study the Lc->PKPi background
DECLARE_SOA_COLUMN(MassLc, massLc, float);
DECLARE_SOA_COLUMN(PtLc, ptLc, float);
DECLARE_SOA_COLUMN(RapidityLc, rapidityLc, float);
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);
DECLARE_SOA_COLUMN(PdgMotherProng0, pdgMotherProng0, int);
DECLARE_SOA_COLUMN(PdgMotherProng1, pdgMotherProng1, int);
DECLARE_SOA_COLUMN(PdgMotherProng2, pdgMotherProng2, int);
DECLARE_SOA_COLUMN(MassKPi, massKPi, float);
DECLARE_SOA_COLUMN(MassKProton, massKProton, float);
DECLARE_SOA_COLUMN(MassPiProton, massPiProton, float);
DECLARE_SOA_COLUMN(BdtBkgScore, bdtBkgScore, float);
DECLARE_SOA_COLUMN(IsRealPKPi, isRealPKPi, int8_t);
DECLARE_SOA_COLUMN(IsRealLcPKPi, isRealLcPKPi, int8_t);
DECLARE_SOA_COLUMN(IsReflected, isReflected, int8_t);

} // namespace charm_polarisation

/// table to study the Lc->PKPi background
DECLARE_SOA_TABLE(HfLcPolBkg, "AOD", "HFLCPOLBKG",
                  charm_polarisation::MassLc,
                  charm_polarisation::PtLc,
                  charm_polarisation::RapidityLc,
                  charm_polarisation::CosThetaStar,
                  charm_polarisation::PdgMotherProng0,
                  charm_polarisation::PdgMotherProng1,
                  charm_polarisation::PdgMotherProng2,
                  charm_polarisation::MassKPi,
                  charm_polarisation::MassKProton,
                  charm_polarisation::MassPiProton,
                  charm_polarisation::BdtBkgScore,
                  charm_polarisation::IsRealPKPi,
                  charm_polarisation::IsRealLcPKPi,
                  charm_polarisation::IsReflected);

} // namespace o2::aod

struct TaskPolarisationCharmHadrons {
  Produces<o2::aod::HfLcPolBkg> rowCandLcBkg;

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
  ConfigurableAxis configThnAxisY{"configThnAxisY", {20, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis configThnAxisCosThetaStarHelicity{"configThnAxisCosThetaStarHelicity", {20, -1.f, 1.f}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis configThnAxisCosThetaStarProduction{"configThnAxisCosThetaStarProduction", {20, -1.f, 1.f}, "cos(#vartheta_{production})"};
  ConfigurableAxis configThnAxisCosThetaStarRandom{"configThnAxisCosThetaStarRandom", {20, -1.f, 1.f}, "cos(#vartheta_{random})"};
  ConfigurableAxis configThnAxisCosThetaStarBeam{"configThnAxisCosThetaStarBeam", {20, -1.f, 1.f}, "cos(#vartheta_{beam})"};
  ConfigurableAxis configThnAxisMlBkg{"configThnAxisMlBkg", {100, 0.f, 1.f}, "ML bkg"};
  ConfigurableAxis configThnAxisInvMassD0{"configThnAxisInvMassD0", {250, 1.65f, 2.15f}, "#it{M}(D^{0}) (GeV/#it{c}^{2})"};                           // only for D*+
  ConfigurableAxis configThnAxisInvMassKPiLc{"configThnAxisInvMassKPiLc", {120, 0.65f, 1.25f}, "#it{M}(K#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"}; // only for Lc+->pKpi
  // ConfigurableAxis configThnAxisMlPrompt{"configThnAxisMlPrompt", {100, 0.f, 1.f}, "ML prompt"};
  ConfigurableAxis configThnAxisMlNonPrompt{"configThnAxisMlNonPrompt", {100, 0.f, 1.f}, "ML non-prompt"};
  // ConfigurableAxis configThnAxisCent{"configThnAxisCent", {102, -1.f, 101.f}, "centrality (%)"};
  ConfigurableAxis configThnAxisNumPvContributors{"configThnAxisNumPvContributors", {300, -0.5f, 299.5f}, "num PV contributors"};
  ConfigurableAxis configThnAxisPtB{"configThnAxisPtB", {3000, 0.f, 300.f}, "#it{p}_{T}(B mother) (GeV/#it{c})"};
  ConfigurableAxis configThnAxisAbsEtaTrackMin{"configThnAxisEtaTrackMin", {3, 0.f, 0.3f}, "min |#it{#eta_{track}}|"};
  ConfigurableAxis configThnAxisNumItsClsMin{"configThnAxisNumItsClsMin", {4, 3.5f, 7.5f}, "min #it{N}_{cls ITS}"};
  ConfigurableAxis configThnAxisNumTpcClsMin{"configThnAxisNumTpcClsMin", {3, 79.5f, 140.5f}, "min #it{N}_{cls TPC}"};

  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 0, "Number of rotated copies (background) per each original candidate"};
  Configurable<float> minRotAngleMultByPi{"minRotAngleMultByPi", 5. / 6, "Minimum angle rotation for track rotation, to be multiplied by pi"};
  Configurable<float> maxRotAngleMultByPi{"maxRotAngleMultByPi", 7. / 6, "Maximum angle rotation for track rotation, to be multiplied by pi"};

  // activate study of systematic uncertainties of tracking
  Configurable<bool> activateTrackingSys{"activateTrackingSys", false, "Activate the study of systematic uncertainties of tracking (only w/o ML)"};

  /// output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", true, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", true, "Activate the THnSparse with cosThStar w.r.t. beam axis"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", true, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  float minInvMass{0.f};
  float maxInvMass{1000.f};

  /// table for Lc->pKpi background studies in MC
  Configurable<int> cosThStarAxisLcPKPiBkgMc{"cosThStarAxisLcPKPiBkgMc", 1, "cos(Theta*) axis for background studies (1 = helicity; 2 = production; 3 = beam; 4 = random)"};

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Filter filterSelectLcToPKPiCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi) || (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPKPi);

  using CollisionsWithMcLabels = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels>>;
  using TracksWithMcLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

  using McParticlesDstarMatched = soa::Join<aod::McParticles, aod::HfCandDstarMcGen>;
  using McParticles3ProngMatched = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;

  using CandDstarWSelFlag = soa::Join<aod::HfCandDstars, aod::HfSelDstarToD0Pi>;
  using CandLcToPKPiWSelFlag = soa::Join<aod::HfCand3Prong, aod::HfSelLc>;

  using FilteredCandDstarWSelFlag = soa::Filtered<CandDstarWSelFlag>;
  using FilteredCandDstarWSelFlagAndMl = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi>>;
  using FilteredCandDstarWSelFlagAndMc = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfCandDstarMcRec>>;
  using FilteredCandDstarWSelFlagAndMcAndMl = soa::Filtered<soa::Join<CandDstarWSelFlag, aod::HfMlDstarToD0Pi, aod::HfCandDstarMcRec>>;

  using FilteredCandLcToPKPiWSelFlag = soa::Filtered<CandLcToPKPiWSelFlag>;
  using FilteredCandLcToPKPiWSelFlagAndMl = soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfMlLcToPKPi>>;
  using FilteredCandLcToPKPiWSelFlagAndMc = soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfCand3ProngMcRec>>;
  using FilteredCandLcToPKPiWSelFlagAndMcAndMl = soa::Filtered<soa::Join<CandLcToPKPiWSelFlag, aod::HfCand3ProngMcRec, aod::HfMlLcToPKPi>>;

  SliceCache cache;
  Preslice<FilteredCandDstarWSelFlag> dstarPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMl> dstarWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMc> dstarWithMcPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMcAndMl> dstarWithMcAndMlPerCollision = aod::hf_cand::collisionId;

  Preslice<FilteredCandLcToPKPiWSelFlag> lcToPKPiPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMl> lcToPKPiWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMc> lcToPKPiWithMcPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMcAndMl> lcToPKPiWithMcAndMlPerCollision = aod::hf_cand::collisionId;

  HfHelper hfHelper;
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    /// check process functions
    std::array<int, 9> processes = {doprocessDstar, doprocessDstarWithMl, doprocessLcToPKPi, doprocessLcToPKPiWithMl, doprocessDstarMc, doprocessDstarMcWithMl, doprocessLcToPKPiMc, doprocessLcToPKPiMcWithMl, doprocessLcToPKPiBackgroundMcWithMl};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    if (activateTrackingSys && (doprocessDstarWithMl || doprocessLcToPKPiWithMl || doprocessDstarMcWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl)) {
      LOGP(fatal, "Tracking systematic study only available in case of no ML!");
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
    if (nBkgRotations > 0 && (doprocessDstarMc || doprocessDstarMcWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl)) {
      LOGP(fatal, "No background rotation supported for MC.");
    }

    massPi = o2::constants::physics::MassPiPlus;
    massProton = o2::constants::physics::MassProton;
    massKaon = o2::constants::physics::MassKaonCharged;
    massDstar = o2::constants::physics::MassDStar;
    massLc = o2::constants::physics::MassLambdaCPlus;
    bkgRotationAngleStep = (nBkgRotations > 1) ? (maxRotAngleMultByPi - minRotAngleMultByPi) * constants::math::PI / (nBkgRotations - 1) : 0.;

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassD0{configThnAxisInvMassD0, "#it{M}(D^{0}) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassKPiLc{configThnAxisInvMassKPiLc, "#it{M}(K#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisY{configThnAxisY, "#it{y}"};
    const AxisSpec thnAxisCosThetaStarHelicity{configThnAxisCosThetaStarHelicity, "cos(#vartheta_{helicity})"};
    const AxisSpec thnAxisCosThetaStarProduction{configThnAxisCosThetaStarProduction, "cos(#vartheta_{production})"};
    const AxisSpec thnAxisCosThetaStarRandom{configThnAxisCosThetaStarRandom, "cos(#vartheta_{random})"};
    const AxisSpec thnAxisCosThetaStarBeam{configThnAxisCosThetaStarBeam, "cos(#vartheta_{beam})"};
    const AxisSpec thnAxisMlBkg{configThnAxisMlBkg, "ML bkg"};
    const AxisSpec thnAxisMlNonPrompt{configThnAxisMlNonPrompt, "ML non-prompt"};
    const AxisSpec thnAxisIsRotatedCandidate{2, -0.5f, 1.5f, "rotated bkg"};
    const AxisSpec thnAxisNumPvContributors{configThnAxisNumPvContributors, "num PV contributors"};
    const AxisSpec thnAxisPtB{configThnAxisPtB, "#it{p}_{T}(B mother) (GeV/#it{c})"};
    const AxisSpec thnAxisDausAcc{2, -0.5f, 1.5f, "daughters in acceptance"};
    const AxisSpec thnAxisResoChannelLc{4, -0.5, 3.5, "0: direct  1,2,3: resonant"}; // 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
    const AxisSpec thnAxisAbsEtaTrackMin{configThnAxisAbsEtaTrackMin, "min |#it{#eta_{track}}|"};
    const AxisSpec thnAxisNumItsClsMin{configThnAxisNumItsClsMin, "min #it{N}_{cls ITS}"};
    const AxisSpec thnAxisNumTpcClsMin{configThnAxisNumTpcClsMin, "min #it{N}_{cls TPC}"};

    auto invMassBins = thnAxisInvMass.binEdges;
    minInvMass = invMassBins.front();
    maxInvMass = invMassBins.back();

    registry.add("hNumPvContributorsAll", "Number of PV contributors for all events ;num. PV contributors; counts", HistType::kTH1F, {thnAxisNumPvContributors});
    registry.add("hNumPvContributorsCand", "Number of PV contributors for events with candidates;num. PV contributors; counts", HistType::kTH1F, {thnAxisNumPvContributors});
    registry.add("hNumPvContributorsCandInMass", "Number of PV contributors for events with candidates in the signal region;num. PV contributors; counts", HistType::kTH1F, {thnAxisNumPvContributors});

    if (doprocessDstarWithMl || doprocessDstarMcWithMl) {
      /// analysis for D*+ meson with ML, w/o rot. background axis
      if (doprocessDstarWithMl) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisPtB});
        }
      }
    } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
      /// analysis for Lc+ baryon with ML, w/ rot. background axis (for data only)
      if (doprocessLcToPKPiWithMl) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisMlBkg, thnAxisMlNonPrompt, thnAxisResoChannelLc});
        }
      }
    } else if (doprocessDstar || doprocessDstarMc) {
      /// analysis for D*+ meson, w/o rot. background axis
      if (doprocessDstar) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarHelicity, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarProduction, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarBeam, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisPtB});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStarRandom, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisPtB});
        }
      }
    } else if (doprocessLcToPKPi || doprocessLcToPKPiMc) {
      /// analysis for Lc+ baryon, rot. background axis (for data only)
      if (doprocessLcToPKPi) {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
      } else {
        if (activateTHnSparseCosThStarHelicity) {
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarProduction) {
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarBeam) {
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
        if (activateTHnSparseCosThStarRandom) {
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- reco non-prompt signal", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom, thnAxisResoChannelLc, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate});
        }
      }
    }

    // MC Gen histos
    if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
      if (activateTHnSparseCosThStarHelicity) {
        registry.add("hGenPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisDausAcc, thnAxisResoChannelLc});
        registry.add("hGenNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc});
      }
      if (activateTHnSparseCosThStarProduction) {
        registry.add("hGenPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarProduction, thnAxisDausAcc, thnAxisResoChannelLc});
        registry.add("hGenNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarProduction, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc});
      }
      if (activateTHnSparseCosThStarBeam) {
        registry.add("hGenPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarBeam, thnAxisDausAcc, thnAxisResoChannelLc});
        registry.add("hGenNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarBeam, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc});
      }
      if (activateTHnSparseCosThStarRandom) {
        registry.add("hGenPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- gen prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarRandom, thnAxisDausAcc, thnAxisResoChannelLc});
        registry.add("hGenNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis -- gen non-prompt signal", HistType::kTHnSparseF, {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarRandom, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc});
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
  /// \param numPvContributors is the number of PV contributors
  /// \param rapCharmHad is the rapidity of the candidate
  /// \param invMassD0 is the invariant-mass of the D0 daugher (only for D*+)
  /// \param invMassKPiLc is the invariant-mass of the K-pi pair (only for Lc+)
  /// \param cosThetaStar is the cosThetaStar of the candidate
  /// \param outputMl is the array with ML output scores
  /// \param isRotatedCandidate is a flag that keeps the info of the rotation of the candidate for bkg studies
  /// \param origin is the MC origin
  /// \param ptBhadMother is the pt of the b-hadron mother (only in case of non-prompt)
  /// \param resoChannelLc indicates the Lc decay channel (direct, resonant)
  /// \param absEtaMin is the minimum absolute eta of the daughter tracks
  /// \param numItsClsMin is the minimum number of ITS clusters of the daughter tracks
  /// \param numTpcClsMin is the minimum number of TPC clusters of the daughter tracks
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, charm_polarisation::CosThetaStarType cosThetaStarType>
  void fillRecoHistos(float invMassCharmHad, float ptCharmHad, int numPvContributors, float rapCharmHad, float invMassD0, float invMassKPiLc, float cosThetaStar, std::array<float, 3> outputMl, int isRotatedCandidate, int8_t origin, float ptBhadMother, int8_t resoChannelLc, float absEtaMin, int numItsClsMin, int numTpcClsMin)
  {
    if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if constexpr (!doMc) {                                                            // data
        if constexpr (withMl) {                                                         // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {  // D*+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if constexpr (!doMc) {                                                                     // data
        if constexpr (withMl) {                                                                  // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {           // D*+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if constexpr (!doMc) {                                                               // data
        if constexpr (withMl) {                                                            // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {     // D*+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if constexpr (!doMc) {                                                                 // data
        if constexpr (withMl) {                                                              // with ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {       // D*+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
          }
        } else {                                                                       // without ML
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          } else {                                                                       // non-prompt
            if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother);
            } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin);
            }
          }
        }
      }
    }
  }

  /// \param ptCharmHad is the pt of the particle
  /// \param numPvContributors is the number of PV contributors
  /// \param rapCharmHad is the rapidity of the particle
  /// \param cosThetaStar is the cosThetaStar of the particle
  /// \param origin is the MC origin
  /// \param ptBhadMother is the pt of the b-hadron mother (only in case of non-prompt)
  /// \param areDausInAcc is a flag indicating whether the daughters are in acceptance or not
  /// \param resoChannelLc indicates the Lc decay channel (direct, resonant)
  template <charm_polarisation::CosThetaStarType cosThetaStarType>
  void fillGenHistos(float ptCharmHad, int numPvContributors, float rapCharmHad, float cosThetaStar, int8_t origin, float ptBhadMother, bool areDausInAcc, uint8_t resoChannelLc)
  {
    if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if (origin == RecoDecay::OriginType::Prompt) {                                    // prompt
        registry.fill(HIST("hGenPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if (origin == RecoDecay::OriginType::Prompt) {                                             // prompt
        registry.fill(HIST("hGenPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if (origin == RecoDecay::OriginType::Prompt) {                                       // prompt
        registry.fill(HIST("hGenPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc);
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if (origin == RecoDecay::OriginType::Prompt) {                                         // prompt
        registry.fill(HIST("hGenPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc);
      } else { // non-prompt
        registry.fill(HIST("hGenNonPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc);
      }
    }
  }

  /// \param numPvContributors is the number of PV contributors
  /// \param nCands is the number of candidates associated to a collision
  /// \param nCandsInMass is the number of candidates in the signal mass region associated to a colslision
  void fillMultHistos(int numPvContributors, int nCands, int nCandsInMass)
  {
    registry.fill(HIST("hNumPvContributorsAll"), numPvContributors);
    if (nCands > 0) {
      registry.fill(HIST("hNumPvContributorsCand"), numPvContributors);
    }
    if (nCandsInMass > 0) {
      registry.fill(HIST("hNumPvContributorsCandInMass"), numPvContributors);
    }
  }

  /// \param invMass is the invariant mass
  /// \return true if candidate in signal region
  template <charm_polarisation::DecayChannel channel>
  bool isInSignalRegion(float invMass)
  {
    if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
      if (0.142f < invMass && invMass < 0.15f) {
        return true;
      }
    } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+ (to be tuned!)
      if (2.25f < invMass && invMass < 2.35f) {
        return true;
      }
    }

    return false;
  }

  /// \param daughter is the daughter particle
  /// \param ptMin is the minimum pt
  /// \param etaMax is the maximum eta
  /// \return true if daughter is in acceptance
  template <typename Part>
  bool isDaughterInAcceptance(Part const& daughter, float ptMin, float etaMax)
  {
    if (daughter.pt() < ptMin) {
      return false;
    }
    if (std::abs(daughter.eta()) > etaMax) {
      return false;
    }

    return true;
  }

  /// \param prongTrack is the track we want to find the mother of
  /// \param idMothers is the vector containing the mother IDs
  /// \param particles are the MC particles
  template <typename Trk, typename Part>
  void searchFirstLevelMother(Trk const& prongTrack, std::vector<int>& idMothers, Part const& /*particles*/)
  {
    /// particle associated to the prong track
    if (!prongTrack.has_mcParticle()) {
      return;
    }
    auto prongParticle = prongTrack.template mcParticle_as<Part>();
    /// leave the vector of mother indices empty if the currect paticle has no mothers
    if (!prongParticle.has_mothers()) {
      return;
    }
    // loop over the mother particles of the analysed particle
    for (auto iMother = prongParticle.mothersIds().front(); iMother <= prongParticle.mothersIds().back(); ++iMother) {
      idMothers.push_back(iMother);
      break; // we keep only the first one
    }
  };

  /// prongTracks is the vector of daughter tracks
  /// etaMin is the minimum eta
  /// nItsClsMin is the minumum number of clusters in ITS
  /// nTpcClsMin is the minumum number of clusters in TPC
  template <typename Trk>
  void getTrackingInfos(std::vector<Trk> const& prongTracks, float& etaMin, int& nItsClsMin, int& nTpcClsMin)
  {
    etaMin = -1.f;
    nItsClsMin = -1;
    nTpcClsMin = -1;

    for (const auto& track : prongTracks) {
      if (std::abs(track.eta()) < etaMin) {
        etaMin = std::abs(track.eta());
      }
      if (track.itsNCls() < nItsClsMin) {
        nItsClsMin = track.itsNCls();
      }
      if (track.tpcNClsCrossedRows() < nTpcClsMin) {
        nTpcClsMin = track.tpcNClsCrossedRows();
      }
    }
  }

  /// \param candidates are the selected candidates
  /// \param bkgRotationId is the id for the background rotation
  /// \param numPvContributors is the number of PV contributors
  /// \param particles are the generated particles
  /// \param tracks are the reconstructed tracks
  /// \return true if candidate in signal region
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, bool studyLcPKPiBkgMc = false, typename Cand, typename Part, typename Trk>
  bool runPolarisationAnalysis(Cand const& candidate, int bkgRotationId, int numPvContributors, Part const& particles, Trk const& /*tracks*/)
  {
    bool isCandidateInSignalRegion{false};
    int8_t origin{RecoDecay::OriginType::None};
    int8_t massHypoMcTruth{-1};
    float ptBhadMother{-1.f};
    int8_t resoChannelLc = -1;
    if constexpr (doMc) {
      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        if (!TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // this candidate is not signal, skip
          return isCandidateInSignalRegion;
        }
        origin = candidate.originMcRec();
        ptBhadMother = candidate.ptBhadMotherPart();
        int pdgBhadMother = candidate.pdgBhadMotherPart();
        // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
        if (origin == RecoDecay::OriginType::NonPrompt && (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 553 || pdgBhadMother == 557)) {
          return isCandidateInSignalRegion;
        }
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if constexpr (!studyLcPKPiBkgMc) {                                                                // skip this if studyLcPKPiBkgMc is true, since we are interested in background
          if (!TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_3prong::DecayType::LcToPKPi)) { // this candidate is not signal, skip
            return isCandidateInSignalRegion;
          }
          origin = candidate.originMcRec();
          if (candidate.isCandidateSwapped()) {
            massHypoMcTruth = charm_polarisation::MassHyposLcToPKPi::PiKP;
          } else {
            massHypoMcTruth = charm_polarisation::MassHyposLcToPKPi::PKPi;
          }
          resoChannelLc = candidate.flagMcDecayChanRec(); /// 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
        }
      }
    }

    // loop over mass hypotheses
    for (uint8_t iMass = 0u; iMass < nMassHypos; iMass++) {

      // variable definition
      float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
      float pxCharmHad{-1000.f}, pyCharmHad{-1000.f}, pzCharmHad{-1000.f};
      float massDau{0.f}, invMassCharmHad{0.f}, invMassCharmHadForSparse{0.f}, invMassD0{0.f}, invMassKPiLc{0.f};
      float rapidity{-999.f};
      std::array<float, 3> outputMl{-1.f, -1.f, -1.f};
      int isRotatedCandidate = 0; // currently meaningful only for Lc->pKpi

      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        // Dstar analysis
        // polarization measured from the soft-pion daughter (*)

        massDau = massPi; // (*)
        const float bkgRotAngle = (bkgRotationId > 0) ? minRotAngleMultByPi * constants::math::PI + bkgRotationAngleStep * (bkgRotationId - 1) : 0;

        std::array<float, 3> threeVecSoftPi{candidate.pxSoftPi() * std::cos(bkgRotAngle) - candidate.pySoftPi() * std::sin(bkgRotAngle), candidate.pxSoftPi() * std::sin(bkgRotAngle) + candidate.pySoftPi() * std::cos(bkgRotAngle), candidate.pzSoftPi()}; // we rotate the soft pion
        std::array<float, 3> threeVecD0Prong0{candidate.pVectorProng0()};
        std::array<float, 3> threeVecD0Prong1{candidate.pVectorProng1()};
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
        const float bkgRotAngle = (bkgRotationId > 0) ? minRotAngleMultByPi * constants::math::PI + bkgRotationAngleStep * (bkgRotationId - 1) : 0;

        std::array<float, 3> threeVecLcProng0{candidate.pVectorProng0()};
        std::array<float, 3> threeVecLcRotatedProng1{candidate.pxProng1() * std::cos(bkgRotAngle) - candidate.pyProng1() * std::sin(bkgRotAngle), candidate.pxProng1() * std::sin(bkgRotAngle) + candidate.pyProng1() * std::cos(bkgRotAngle), candidate.pzProng1()};
        std::array<float, 3> threeVecLcProng2{candidate.pVectorProng2()};
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
          // invariant mass of the KPi pair
          invMassKPiLc = hfHelper.invMassKPiPairLcToPKPi(candidate);
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
          // invariant mass of the KPi pair
          invMassKPiLc = hfHelper.invMassKPiPairLcToPiKP(candidate);
        } else {
          // NB: no need to check cases in which candidate.isSelLcToPKPi() and candidate.isSelLcToPiKP() are both false, because they are rejected already by the Filter
          // ... but we need to put this protections here!
          // Otherwise, a candidate selected as pKpi only has invMassCharmHad==0 when iMass == charm_polarisation::MassHyposLcToPKPi::PiKP and viceversa
          continue;
        }

      } // Lc->pKpi

      if (invMassCharmHadForSparse < minInvMass || invMassCharmHadForSparse > maxInvMass) {
        continue;
      }

      float phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      float thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
      ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
      ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, invMassCharmHad);
      ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
      ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
      ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();

      float ptCharmHad = std::sqrt(pxCharmHad * pxCharmHad + pyCharmHad * pyCharmHad); // this definition is valid for both rotated and original candidates

      if (!isCandidateInSignalRegion) { // it could be that only one mass hypothesis is in signal region
        isCandidateInSignalRegion = isInSignalRegion<channel>(invMassCharmHadForSparse);
      }

      float absEtaTrackMin{-1.f};
      int numItsClsMin{-1}, numTpcClsMin{-1};

      if constexpr (!withMl) {
        if (activateTrackingSys) {
          auto trackProng0 = candidate.template prong0_as<Trk>();
          auto trackProng1 = candidate.template prong1_as<Trk>();
          if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
            auto trackProng2 = candidate.template prongPi_as<Trk>();
            getTrackingInfos(std::vector{trackProng0, trackProng1, trackProng2}, absEtaTrackMin, numItsClsMin, numTpcClsMin);
          } else if (channel == charm_polarisation::DecayChannel::LcToPKPi) {
            auto trackProng2 = candidate.template prong2_as<Trk>();
            getTrackingInfos(std::vector{trackProng0, trackProng1, trackProng2}, absEtaTrackMin, numItsClsMin, numTpcClsMin);
          } 
        }
      }

      float cosThetaStarHelicity = -10.f;
      float cosThetaStarProduction = -10.f;
      float cosThetaStarBeam = -10.f;
      float cosThetaStarRandom = -10.f;
      if (activateTHnSparseCosThStarHelicity) {
        ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
        cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Helicity>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarHelicity, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin);
      }
      if (activateTHnSparseCosThStarProduction) {
        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
        cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Production>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarProduction, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin);
      }
      if (activateTHnSparseCosThStarBeam) {
        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Beam>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarBeam, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin);
      }
      if (activateTHnSparseCosThStarRandom) {
        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Random>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarRandom, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin);
      }

      /// Table for Lc->pKpi background studies
      /// Defined only in MC simulations, to study resonances and reflected signal
      if constexpr (doMc && channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if constexpr (studyLcPKPiBkgMc) {
          /******************************************************************************************
          The code below can work only without grouping on "mcCollision".
          In fact, grouping by "mcCollision" introduces the following inconsistencies:

            1) the particle getters "track.template mcParticle_as<Part>()" retrieve the daughter particles quering the full particle table in the dataframe.
               In other words, even if the 3-prong candidate is reconstructed in a completely wrong reco. collision due to the track-to-collision associator,
               therefore this collision points to a "mcCollision" different from the current one and the daughter particles are associated to this different "mcCollision",
               the getter "mcParticle_as<Part>" works anyway, because it works with unbound tables ignoring the fact that "particles" is grouped;

            2) when we look for the mother index from the daughter particles of the previous point, but the daughter particles belong to a "mcCollision" different from the current one,
               then also the mother particle belongs to this different "mcCollision". This means that the mother index ( - "particles.offset()") is outside the "particles.size()",
               because the table "particles" is grouped w.r.t. the current "mcCollision".
          *******************************************************************************************/

          /// check if the tracks are associated to a pion + a kaon + a proton
          int8_t isRealPKPi = 0; /// true only if the triplet is formed by a MC pion + a MC kaon + a MC proton
          bool isGenPKPi = false;
          bool isGenPiKP = false;
          auto trackProng0 = candidate.template prong0_as<Trk>();
          auto trackProng1 = candidate.template prong1_as<Trk>();
          auto trackProng2 = candidate.template prong2_as<Trk>();
          int pdgProng0 = 0;
          int pdgProng1 = 0;
          int pdgProng2 = 0;
          if (trackProng0.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng0 = trackProng0.template mcParticle_as<Part>();
            pdgProng0 = particleProng0.pdgCode();
          }
          if (trackProng1.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng1 = trackProng1.template mcParticle_as<Part>();
            pdgProng1 = particleProng1.pdgCode();
          }
          if (trackProng2.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng2 = trackProng2.template mcParticle_as<Part>();
            pdgProng2 = particleProng2.pdgCode();
          }
          isGenPKPi = std::abs(pdgProng0) == kProton && std::abs(pdgProng1) == kKPlus && std::abs(pdgProng2) == kPiPlus;
          isGenPiKP = std::abs(pdgProng0) == kPiPlus && std::abs(pdgProng1) == kKPlus && std::abs(pdgProng2) == kProton;
          if (isGenPKPi || isGenPiKP) {
            isRealPKPi = 1;
          }

          /// check if the triplet is reflected or not
          /// i.e. generated as pKpi but reconstructed as piKp, or viceversa
          int8_t isReflected = 0;
          if (isRealPKPi && ((iMass == charm_polarisation::MassHyposLcToPKPi::PKPi && candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi && isGenPiKP) || (iMass == charm_polarisation::MassHyposLcToPKPi::PiKP && candidate.isSelLcToPiKP() >= selectionFlagLcToPKPi && isGenPKPi))) {
            isReflected = 1;
          }

          /// check if the pKpi triplet is a Lc->pKpi
          int8_t isRealLcPKPi = 0;
          if (isRealPKPi && TESTBIT(std::abs(candidate.flagMcMatchRec()), aod::hf_cand_3prong::DecayType::LcToPKPi)) {
            isRealLcPKPi = 1;
          }

          /// look for daughters' mothers (1st level only)
          std::vector<int> idMothersProng0 = {};
          std::vector<int> idMothersProng1 = {};
          std::vector<int> idMothersProng2 = {};
          searchFirstLevelMother(trackProng0, idMothersProng0, particles);
          searchFirstLevelMother(trackProng1, idMothersProng1, particles);
          searchFirstLevelMother(trackProng2, idMothersProng2, particles);

          /// check if daughter pairs have the same mother
          /// it should be enough to check the 1st one only (only particles from partonic events, or interactions with material, can have more than 1 mother)
          int pdgMotherProng0 = -1;
          int pdgMotherProng1 = -1;
          int pdgMotherProng2 = -1;
          if (idMothersProng0.size() > 0 && idMothersProng1.size() > 0 && idMothersProng0.at(0) == idMothersProng1.at(0)) {
            /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
            /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
            auto mother = particles.rawIteratorAt(idMothersProng0.at(0) - particles.offset());
            pdgMotherProng0 = std::abs(mother.pdgCode()); // PDG code of the mother
            pdgMotherProng1 = pdgMotherProng0;
          }
          if (idMothersProng1.size() > 0 && idMothersProng2.size() > 0 && idMothersProng1.at(0) == idMothersProng2.at(0)) {
            /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
            /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
            auto mother = particles.rawIteratorAt(idMothersProng1.at(0) - particles.offset());
            pdgMotherProng1 = std::abs(mother.pdgCode()); // PDG code of the mother
            pdgMotherProng2 = pdgMotherProng1;
          }
          if (idMothersProng0.size() > 0 && idMothersProng2.size() > 0 && idMothersProng0.at(0) == idMothersProng2.at(0)) {
            /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
            /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
            auto mother = particles.rawIteratorAt(idMothersProng0.at(0) - particles.offset());
            pdgMotherProng0 = std::abs(mother.pdgCode()); // PDG code of the mother
            pdgMotherProng2 = pdgMotherProng0;
          }

          /// calculate inv. masses for pairs, depending on mass hypothesis
          std::array<float, 3> pVecPion = {};
          std::array<float, 3> pVecKaon = candidate.pVectorProng1();
          std::array<float, 3> pVecProton = {};
          if (iMass == charm_polarisation::MassHyposLcToPKPi::PKPi && candidate.isSelLcToPKPi() >= selectionFlagLcToPKPi) {
            pVecProton = candidate.pVectorProng0();
            pVecPion = candidate.pVectorProng2();
          } else if (iMass == charm_polarisation::MassHyposLcToPKPi::PiKP && candidate.isSelLcToPiKP() >= selectionFlagLcToPKPi) {
            pVecProton = candidate.pVectorProng2();
            pVecPion = candidate.pVectorProng0();
          }
          const float massKPi = RecoDecay::m(std::array{pVecKaon, pVecPion}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
          const float massKProton = RecoDecay::m(std::array{pVecKaon, pVecProton}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassProton});
          const float massPiProton = RecoDecay::m(std::array{pVecPion, pVecProton}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassProton});

          /// Fill the table for selected candidates
          /// No need to check explicitly if candidates are selected, since the Filter is applied
          float cosThetaStarForTable = -10.f;
          switch (cosThStarAxisLcPKPiBkgMc) {
            case 1:
              cosThetaStarForTable = cosThetaStarHelicity;
              break;
            case 2:
              cosThetaStarForTable = cosThetaStarProduction;
              break;
            case 3:
              cosThetaStarForTable = cosThetaStarBeam;
              break;
            case 4:
              cosThetaStarForTable = cosThetaStarRandom;
              break;
            default:
              LOG(fatal) << "cosThStarAxisLcPKPiBkgMc must be between 1 and 4 (1: helicity; 2: production; 3: beam; 4: random), but cosThStarAxisLcPKPiBkgMc = " << cosThStarAxisLcPKPiBkgMc << ". Fix it!";
              break;
          }
          rowCandLcBkg(invMassCharmHadForSparse, ptCharmHad, rapidity,
                       cosThetaStarForTable,
                       pdgMotherProng0, pdgMotherProng1, pdgMotherProng2,
                       massKPi, massKProton, massPiProton,
                       outputMl.at(0),
                       isRealPKPi, isRealLcPKPi, isReflected);
        } // end studyLcPKPiBkgMc
      }   // end table for Lc->pKpi background studies

    } /// end loop over mass hypotheses

    return isCandidateInSignalRegion;
  }

  /// \param mcParticle is the Mc particle
  /// \param mcParticles is the table of Mc particles
  /// \param numPvContributors is the number of PV contributors in the associated reco collision
  template <charm_polarisation::DecayChannel channel, typename Part, typename Particles>
  void runMcGenPolarisationAnalysis(Part const& mcParticle, Particles const& mcParticles, int numPvContributors)
  {
    int8_t origin{RecoDecay::OriginType::None};
    std::vector<int> listDaughters{};
    float massDau{0.f}, massCharmHad{0.f};
    float ptBhadMother{-1.f};
    bool areDauInAcc{true};
    int8_t resoChannelLc = -1;
    if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
      if (!TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_dstar::DecayType::DstarToD0Pi)) { // this particle is not signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(mcParticle.idxBhadMotherPart() - mcParticles.offset());
        int pdgBhadMother = std::abs(bHadMother.pdgCode());
        // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
        if (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 553 || pdgBhadMother == 557) {
          return;
        }
        ptBhadMother = bHadMother.pt();
      }

      std::array<int, 2> dauPdgs = {kPiPlus, o2::constants::physics::Pdg::kD0};
      RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 1);
      massDau = massPi;
      massCharmHad = massDstar;
    } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
      if (!TESTBIT(std::abs(mcParticle.flagMcMatchGen()), aod::hf_cand_3prong::DecayType::LcToPKPi)) { // this particle is not signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      resoChannelLc = mcParticle.flagMcDecayChanGen();
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
      auto dauPart = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        if (std::abs(dauPart.pdgCode()) == kPiPlus) {
          pxDau = dauPart.px();
          pyDau = dauPart.py();
          pzDau = dauPart.pz();
          if (areDauInAcc) {
            areDauInAcc = isDaughterInAcceptance(dauPart, 0.1, 0.8);
          }
        } else if (areDauInAcc) { // check also D0 daughters
          std::vector<int> listDaughtersD0{};
          std::array<int, 2> dauPdgsD0 = {kPiPlus, -kKPlus};
          RecoDecay::getDaughters(mcParticle, &listDaughtersD0, dauPdgsD0, 1);
          for (const auto& dauIdxD0 : listDaughtersD0) {
            auto dauPartD0 = mcParticles.rawIteratorAt(dauIdxD0 - mcParticles.offset());
            if (areDauInAcc) {
              areDauInAcc = isDaughterInAcceptance(dauPartD0, 0.3, 0.8);
            }
          }
        }
      } else if constexpr (channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if (std::abs(dauPart.pdgCode()) == kProton) {
          pxDau = dauPart.px();
          pyDau = dauPart.py();
          pzDau = dauPart.pz();
        }
        if (areDauInAcc) {
          areDauInAcc = isDaughterInAcceptance(dauPart, 0.3, 0.8);
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
      fillGenHistos<charm_polarisation::CosThetaStarType::Helicity>(ptCharmHad, numPvContributors, rapidity, cosThetaStarHelicity, origin, ptBhadMother, areDauInAcc, resoChannelLc);
    }
    if (activateTHnSparseCosThStarProduction) {
      ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
      float cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Production>(ptCharmHad, numPvContributors, rapidity, cosThetaStarProduction, origin, ptBhadMother, areDauInAcc, resoChannelLc);
    }
    if (activateTHnSparseCosThStarBeam) {
      ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      float cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Beam>(ptCharmHad, numPvContributors, rapidity, cosThetaStarBeam, origin, ptBhadMother, areDauInAcc, resoChannelLc);
    }
    if (activateTHnSparseCosThStarRandom) {
      ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      float cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Random>(ptCharmHad, numPvContributors, rapidity, cosThetaStarRandom, origin, ptBhadMother, areDauInAcc, resoChannelLc);
    }
  }

  /////////////////////////
  //   Dstar analysis   ///
  /////////////////////////

  // Dstar with rectangular cuts
  void processDstar(aod::Collisions const& collisions,
                    FilteredCandDstarWSelFlag const& dstarCandidates,
                    TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }

        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(dstarCandidate, iRotation, numPvContributors, -1 /*MC particles*/, tracks);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with ML cuts
  void processDstarWithMl(aod::Collisions const& collisions,
                          FilteredCandDstarWSelFlagAndMl const& dstarCandidates)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/)) {
          nCandsInSignalRegion++;
        }

        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate, iRotation, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarWithMl, "Process Dstar candidates with ML", false);

  // Dstar in MC with rectangular cuts
  void processDstarMc(aod::McCollisions::iterator const&,
                      McParticlesDstarMatched const& mcParticles,
                      CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                      FilteredCandDstarWSelFlagAndMc const& dstarCandidates,
                      TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMcPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarMc, "Process Dstar candidates in MC without ML", false);

  // Dstar in MC with ML cuts
  void processDstarMcWithMl(aod::McCollisions::iterator const&,
                            McParticlesDstarMatched const& mcParticles,
                            CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                            FilteredCandDstarWSelFlagAndMcAndMl const& dstarCandidates)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processDstarMcWithMl, "Process Dstar candidates in MC with ML", false);

  ////////////////////////////
  //   Lc->pKpi analysis   ///
  ////////////////////////////

  // Lc->pKpi with rectangular cuts
  void processLcToPKPi(aod::Collisions const& collisions,
                       FilteredCandLcToPKPiWSelFlag const& lcCandidates,
                       TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, false>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }

        /// rotational background
        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, false>(lcCandidate, iRotation, numPvContributors, -1 /*MC particles*/, tracks);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPi, "Process Lc candidates without ML", false);

  // Lc->pKpi with ML cuts
  void processLcToPKPiWithMl(aod::Collisions const& collisions,
                             FilteredCandLcToPKPiWSelFlagAndMl const& lcCandidates)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiWithMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/)) {
          nCandsInSignalRegion++;
        }

        /// rotational background
        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate, iRotation, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiWithMl, "Process Lc candidates with ML", false);

  // Lc->pKpi in MC with rectangular cuts
  void processLcToPKPiMc(aod::McCollisions::iterator const&,
                         McParticles3ProngMatched const& mcParticles,
                         CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                         FilteredCandLcToPKPiWSelFlagAndMc const& lcCandidates,
                         TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, false, true>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiMc, "Process Lc candidates in MC without ML", false);

  // Lc->pKpi in MC with ML cuts
  void processLcToPKPiMcWithMl(aod::McCollisions::iterator const&,
                               McParticles3ProngMatched const& mcParticles,
                               CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroups
                               FilteredCandLcToPKPiWSelFlagAndMcAndMl const& lcCandidates)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, true>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, -1 /*tracks*/)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiMcWithMl, "Process Lc candidates in MC with ML", false);

  // Lc->pKpi in MC with ML cuts w/o mcCollision grouping (to study Lc background)
  void processLcToPKPiBackgroundMcWithMl(McParticles3ProngMatched const& mcParticles,
                                         FilteredCandLcToPKPiWSelFlagAndMcAndMl const& lcCandidates,
                                         TracksWithMcLabels const& tracks)
  {
    for (const auto& lcCandidate : lcCandidates) {
      runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, true, true>(lcCandidate, 0, /*numPvContributors*/ -1, mcParticles, tracks);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles, /*numPvContributorsGen*/ -1);
    }
  }
  PROCESS_SWITCH(TaskPolarisationCharmHadrons, processLcToPKPiBackgroundMcWithMl, "Process Lc candidates in MC with ML w/o mcCollision grouping", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskPolarisationCharmHadrons>(cfgc)};
}
