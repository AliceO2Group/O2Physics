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
/// \author M. Li (CCNU) mingze.li@cern.ch

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/Utils/utilsFlow.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include <CommonConstants/MathConstants.h>
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
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h> // IWYU pragma: keep (do not replace with Math/Vector3Dfwd.h)
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::analysis::hf_flow_utils;

namespace o2::aod
{
namespace charm_polarisation
{
enum CosThetaStarType : uint8_t {
  Helicity = 0,
  Production,
  Beam,
  Random,
  EP,
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
enum QvecEstimator : uint8_t {
  FV0A = 0,
  FT0M,
  FT0C,
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
DECLARE_SOA_COLUMN(BdtNonPromptScore, bdtNonPromptScore, float);
DECLARE_SOA_COLUMN(IsRealPKPi, isRealPKPi, int8_t);
DECLARE_SOA_COLUMN(IsRealLcPKPi, isRealLcPKPi, int8_t);
DECLARE_SOA_COLUMN(IsReflected, isReflected, int8_t);
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(Origin, origin, int8_t);

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
                  charm_polarisation::BdtNonPromptScore,
                  charm_polarisation::IsRealPKPi,
                  charm_polarisation::IsRealLcPKPi,
                  charm_polarisation::IsReflected,
                  charm_polarisation::Charge,
                  charm_polarisation::Origin);

} // namespace o2::aod

struct HfTaskCharmPolarisation {
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

  // Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qVecDetector{"qVecDetector", 2, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0C: 2)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimator ((None: 0, FT0C: 2, FT0M: 3))"};
  Configurable<int> centralityMin{"centralityMin", 30, "Minimum centrality (0-100) to be considered in the analysis"};
  Configurable<int> centralityMax{"centralityMax", 50, "Maximum centrality (0-100) to be considered in the analysis"};

  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 0, "Number of rotated copies (background) per each original candidate"};
  Configurable<float> minRotAngleMultByPi{"minRotAngleMultByPi", 5. / 6, "Minimum angle rotation for track rotation, to be multiplied by pi"};
  Configurable<float> maxRotAngleMultByPi{"maxRotAngleMultByPi", 7. / 6, "Maximum angle rotation for track rotation, to be multiplied by pi"};

  // activate study of systematic uncertainties of tracking
  Configurable<bool> activateTrackingSys{"activateTrackingSys", false, "Activate the study of systematic uncertainties of tracking"};

  /// output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", true, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", true, "Activate the THnSparse with cosThStar w.r.t. beam axis"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", true, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  Configurable<bool> activateTHnSparseCosThStarEP{"activateTHnSparseCosThStarEP", false, "Activate the THnSparse with cosThStar w.r.t. reaction plane axis"};
  Configurable<bool> activatePartRecoDstar{"activatePartRecoDstar", false, "Activate the study of partly reconstructed D*+ -> D0 (-> KPiPi0) Pi decays"};
  float minInvMass{0.f};
  float maxInvMass{1000.f};

  /// table for Lc->pKpi background studies in MC
  Configurable<int> cosThStarAxisLcPKPiBkgMc{"cosThStarAxisLcPKPiBkgMc", 1, "cos(Theta*) axis for background studies (1 = helicity; 2 = production; 3 = beam; 4 = random)"};

  /// veto conditions for Lc->pKpi analysis
  struct : ConfigurableGroup {
    Configurable<bool> applyLcBkgVeto{"applyLcBkgVeto", false, "Flag to enable the veto on D+ and Ds+ background for Lc->pKpi analysis"};
    /// background from D+->K-pi+pi+
    Configurable<bool> enableLcBkgVetoDplusKPiPi{"enableLcBkgVetoDplusKPiPi", false, "Flag to enable the veto on D+->K-pi+pi+ for Lc->pKpi analysis"};
    Configurable<float> massDplusKPiPiMinVeto{"massDplusKPiPiMinVeto", 1.85, "Min. value for D+->K-pi+pi+ veto"};
    Configurable<float> massDplusKPiPiMaxVeto{"massDplusKPiPiMaxVeto", 1.90, "Max. value for D+->K-pi+pi+ veto"};
    /// background from D+->K+K-pi+
    Configurable<bool> enableLcBkgVetoDplusKKPi{"enableLcBkgVetoDplusKKPi", false, "Flag to enable the veto on D+->K+K-pi+ for Lc->pKpi analysis"};
    Configurable<float> massDplusKKPiMinVeto{"massDplusKKPiMinVeto", 1.85, "Min. value for D+->K+K-pi+ veto"}; // one can use also massDplusKPiPiMinVeto, but this allows more flexibility in analysis
    Configurable<float> massDplusKKPiMaxVeto{"massDplusKKPiMaxVeto", 1.90, "Max. value for D+->K+K-pi+ veto"}; // one can use also massDplusKPiPiMaxVeto, but this allows more flexibility in analysis
    /// background from Ds+->K+K-pi+
    Configurable<bool> enableLcBkgVetoDsKKPi{"enableLcBkgVetoDsKKPi", false, "Flag to enable the veto on Ds+->K+K-pi+ for Lc->pKpi analysis"};
    Configurable<float> massDsKKPiMinVeto{"massDsKKPiMinVeto", 1.94, "Min. value for Ds+->K+K-pi+ veto"};
    Configurable<float> massDsKKPiMaxVeto{"massDsKKPiMaxVeto", 2.00, "Max. value for Ds+->K+K-pi+ veto"};
  } lcBkgVeto;
  struct : ConfigurableGroup {
    /// monitoring histograms (Dalitz plot)
    Configurable<bool> activateTHnLcChannelMonitor{"activateTHnLcChannelMonitor", false, "Flag to switch on the monitoring THnSparse of M2(Kpi), M2(pK), M2(ppi), pt correlation for Lc -> pKpi"};
    ConfigurableAxis configThnAxisInvMass2KPiLcMonitoring{"configThnAxisInvMass2KPiLcMonitoring", {200, 0.3f, 2.3f}, "#it{M}^{2}(K#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisInvMass2PKLcMonitoring{"configThnAxisInvMass2PKLcMonitoring", {320, 2.f, 5.2f}, "#it{M}^{2}(pK) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisInvMass2PPiLcMonitoring{"configThnAxisInvMass2PPiLcMonitoring", {400, 1.f, 5.f}, "#it{M}^{2}(p#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};

    /// veto conditions on Lc->pKpi signals
    Configurable<bool> applyLcSignalVeto{"applyLcSignalVeto", false, "Flag to enable the veto on Lc->pKpi resonant channels"};
    Configurable<float> mass2PPiLcMinVeto{"mass2PPiLcMinVeto", 1.f, "Min. value for Delta++(<-Lc) mass veto"};
    Configurable<float> mass2PPiLcMaxVeto{"mass2PPiLcMaxVeto", 1.6f, "Max. value for Delta++(<-Lc) mass veto"};

  } lcPKPiChannels;

  /// Monitoring of phi Euler angle
  Configurable<bool> activateTHnEulerPhiMonitor{"activateTHnEulerPhiMonitor", false, "Flag to switch on the monitoring THnSparse vs. Euler angle phi (Lc -> pKpi)"};

  /// Application of rapidity cut for reconstructed candidates
  Configurable<float> rapidityCut{"rapidityCut", 999.f, "Max. value of reconstructed candidate rapidity (abs. value)"};

  SliceCache cache;
  EventPlaneHelper epHelper; // event plane helper
  HfEventSelection hfEvSel;  // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  using CollisionsWithMcLabels = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels>>;
  using CollisionsWithMcLabelsAndCent = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Ms, aod::CentFT0Cs>>;
  using CollsWithQVecs = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::QvectorBTots, aod::CentFT0Ms, aod::CentFT0Cs>;
  using TracksWithMcLabels = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>;
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

  Filter filterSelectDstarCandidates = aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi;
  Filter filterSelectLcToPKPiCandidates = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi) || (aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPKPi);

  Preslice<FilteredCandDstarWSelFlag> dstarPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMl> dstarWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMc> dstarWithMcPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandDstarWSelFlagAndMcAndMl> dstarWithMcAndMlPerCollision = aod::hf_cand::collisionId;

  Preslice<FilteredCandLcToPKPiWSelFlag> lcToPKPiPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMl> lcToPKPiWithMlPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMc> lcToPKPiWithMcPerCollision = aod::hf_cand::collisionId;
  Preslice<FilteredCandLcToPKPiWSelFlagAndMcAndMl> lcToPKPiWithMcAndMlPerCollision = aod::hf_cand::collisionId;

  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

  ConfigurableAxis configTHnAxisEulerPhi{"configTHnAxisEulerPhi", {24, -o2::constants::math::PI, o2::constants::math::PI}, "Euler polar angle #phi"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {200, 0.139f, 0.179f}, "#it{M} (GeV/#it{c}^{2})"}; // o2-linter: disable=pdg/explicit-mass (false positive)
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.f, 100.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisY{"configThnAxisY", {20, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis configThnAxisCosThetaStarHelicity{"configThnAxisCosThetaStarHelicity", {20, -1.f, 1.f}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis configThnAxisCosThetaStarProduction{"configThnAxisCosThetaStarProduction", {20, -1.f, 1.f}, "cos(#vartheta_{production})"};
  ConfigurableAxis configThnAxisCosThetaStarRandom{"configThnAxisCosThetaStarRandom", {20, -1.f, 1.f}, "cos(#vartheta_{random})"};
  ConfigurableAxis configThnAxisCosThetaStarBeam{"configThnAxisCosThetaStarBeam", {20, -1.f, 1.f}, "cos(#vartheta_{beam})"};
  ConfigurableAxis configThnAxisCosThetaStarEP{"configThnAxisCosThetaStarEP", {20, -1.f, 1.f}, "cos(#vartheta_{EP})"};
  ConfigurableAxis configThnAxisMlBkg{"configThnAxisMlBkg", {100, 0.f, 1.f}, "ML bkg"};
  ConfigurableAxis configThnAxisInvMassD0{"configThnAxisInvMassD0", {250, 1.65f, 2.15f}, "#it{M}(D^{0}) (GeV/#it{c}^{2})"};                           // only for D*+
  ConfigurableAxis configThnAxisInvMassKPiLc{"configThnAxisInvMassKPiLc", {120, 0.65f, 1.25f}, "#it{M}(K#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"}; // only for Lc+->pKpi
  // ConfigurableAxis configThnAxisMlPrompt{"configThnAxisMlPrompt", {100, 0.f, 1.f}, "ML prompt"};
  ConfigurableAxis configThnAxisMlNonPrompt{"configThnAxisMlNonPrompt", {100, 0.f, 1.f}, "ML non-prompt"};
  ConfigurableAxis configThnAxisCent{"configThnAxisCent", {102, -1.f, 101.f}, "centrality (%)"};
  ConfigurableAxis configThnAxisNumPvContributors{"configThnAxisNumPvContributors", {300, -0.5f, 299.5f}, "num PV contributors"};
  ConfigurableAxis configThnAxisPtB{"configThnAxisPtB", {3000, 0.f, 300.f}, "#it{p}_{T}(B mother) (GeV/#it{c})"};
  ConfigurableAxis configThnAxisAbsEtaTrackMin{"configThnAxisAbsEtaTrackMin", {3, 0.f, 0.3f}, "min |#it{#eta_{track}}|"};
  ConfigurableAxis configThnAxisNumItsClsMin{"configThnAxisNumItsClsMin", {4, 3.5f, 7.5f}, "min #it{N}_{cls ITS}"};
  ConfigurableAxis configThnAxisNumTpcClsMin{"configThnAxisNumTpcClsMin", {3, 79.5f, 140.5f}, "min #it{N}_{cls TPC}"};
  ConfigurableAxis configThnAxisCharge{"configThnAxisCharge", {2, -2.f, 2.f}, "electric charge"};
  ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {100, 0.f, 100.f}, "centrality (%)"};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    /// check process functions
    std::array<int, 13> processes = {doprocessDstar, doprocessDstarWithMl, doprocessLcToPKPi, doprocessLcToPKPiWithMl, doprocessDstarMc, doprocessDstarMcWithMl, doprocessLcToPKPiMc, doprocessLcToPKPiMcWithMl, doprocessLcToPKPiBackgroundMcWithMl, doprocessDstarInPbPb, doprocessDstarWithMlInPbPb, doprocessDstarMcInPbPb, doprocessDstarMcWithMlInPbPb};
    const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    } else if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    /// check output THnSparses
    std::array<int, 5> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarProduction, activateTHnSparseCosThStarBeam, activateTHnSparseCosThStarRandom, activateTHnSparseCosThStarEP};
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
      if (activateTHnSparseCosThStarEP) {
        LOGP(info, "THnSparse with cosThStar w.r.t. event plane axis active.");
      }
    }

    if (activatePartRecoDstar && !(doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb)) {
      LOGP(fatal, "Check on partly reconstructed D* mesons only possible for processDstarMc and processDstarMcWithMl");
    }

    // check bkg rotation for MC (not supported currently)
    if (nBkgRotations > 0 && (doprocessDstarMc || doprocessDstarMcWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb)) {
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
    const AxisSpec thnAxisCosThetaStarEP{configThnAxisCosThetaStarEP, "cos(#vartheta_{EP})"}; // reaction plane
    const AxisSpec thnAxisMlBkg{configThnAxisMlBkg, "ML bkg"};
    const AxisSpec thnAxisMlNonPrompt{configThnAxisMlNonPrompt, "ML non-prompt"};
    const AxisSpec thnAxisIsRotatedCandidate{2, -0.5f, 1.5f, "rotated bkg"};
    const AxisSpec thnAxisNumPvContributors{configThnAxisNumPvContributors, "num PV contributors"};
    const AxisSpec thnAxisPtB{configThnAxisPtB, "#it{p}_{T}(B mother) (GeV/#it{c})"};
    const AxisSpec thnAxisDausAcc{2, -0.5f, 1.5f, "daughters in acceptance"};
    const AxisSpec thnAxisDauToMuons{4, -0.5f, 3.5f, "daughters decayed to muons"};
    const AxisSpec thnAxisResoChannelLc{4, -0.5, 3.5, "0: direct  1,2,3: resonant"}; // 0: direct; 1: Λc± → p± K*; 2: Λc± → Δ(1232)±± K∓; 3: Λc± → Λ(1520) π±
    const AxisSpec thnAxisAbsEtaTrackMin{configThnAxisAbsEtaTrackMin, "min |#it{#eta_{track}}|"};
    const AxisSpec thnAxisNumItsClsMin{configThnAxisNumItsClsMin, "min #it{N}_{cls ITS}"};
    const AxisSpec thnAxisNumTpcClsMin{configThnAxisNumTpcClsMin, "min #it{N}_{cls TPC}"};
    const AxisSpec thnAxisCharge{configThnAxisCharge, "charge"};
    const AxisSpec thnAxisInvMass2KPiLcMonitoring{lcPKPiChannels.configThnAxisInvMass2KPiLcMonitoring, "#it{M}^{2}(K#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2}"};
    const AxisSpec thnAxisInvMass2PKLcMonitoring{lcPKPiChannels.configThnAxisInvMass2PKLcMonitoring, "#it{M}^{2}(pK) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMass2PPiLcMonitoring{lcPKPiChannels.configThnAxisInvMass2PPiLcMonitoring, "#it{M}^{2}(p#pi) from #Lambda_{c}^{+} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisTHnAxisEulerPhi{configTHnAxisEulerPhi, "Euler polar angle #phi"};
    const AxisSpec thnAxisCentrality{configThnAxisCentrality, "centrality (%)"};

    auto invMassBins = thnAxisInvMass.binEdges;
    minInvMass = invMassBins.front();
    maxInvMass = invMassBins.back();

    registry.add("hNumPvContributorsAll", "Number of PV contributors for all events ;num. PV contributors; counts", HistType::kTH1D, {thnAxisNumPvContributors});
    registry.add("hNumPvContributorsCand", "Number of PV contributors for events with candidates;num. PV contributors; counts", HistType::kTH1D, {thnAxisNumPvContributors});
    registry.add("hNumPvContributorsCandInMass", "Number of PV contributors for events with candidates in the signal region;num. PV contributors; counts", HistType::kTH1D, {thnAxisNumPvContributors});
    if (doprocessDstarInPbPb || doprocessDstarMcInPbPb || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
      registry.add("hCentrality", "Centrality distribution for D*+ candidates;centrality (%); counts", HistType::kTH1D, {thnAxisCentrality});
    }

    if (activateTHnSparseCosThStarHelicity) {
      std::vector<AxisSpec> hHelicityaxes = {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY};
      if (doprocessDstar || doprocessDstarMc || doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarInPbPb || doprocessDstarMcInPbPb || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
        hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisInvMassD0, thnAxisCosThetaStarHelicity});
        if (doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
          hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTrackingSys) {
          hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
        }
        if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb) {
          std::vector<AxisSpec> hRecoPromptHelicityAxes(hHelicityaxes);
          hRecoPromptHelicityAxes.insert(hRecoPromptHelicityAxes.end(), {thnAxisDauToMuons});
          std::vector<AxisSpec> hRecoNonPromptHelicityAxes(hHelicityaxes);
          hRecoNonPromptHelicityAxes.insert(hRecoNonPromptHelicityAxes.end(), {thnAxisDauToMuons, thnAxisPtB});
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptHelicityAxes);
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptHelicityAxes);
          if (activatePartRecoDstar) {
            registry.add("hPartRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for partially reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptHelicityAxes);
            registry.add("hPartRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for partially reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptHelicityAxes);
          }
        } else {
          if (nBkgRotations > 0) {
            hHelicityaxes.push_back(thnAxisIsRotatedCandidate);
          }
          registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, hHelicityaxes);
        }
      } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPi || doprocessLcToPKPiMc) { // Lc->pKpi
        hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisInvMassKPiLc, thnAxisCosThetaStarHelicity});
        if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) { // Lc->pKpi with ML
          hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) { // Lc->pKpi MC
          std::vector<AxisSpec> hRecoHelicityAxes(hHelicityaxes);
          if (doprocessLcToPKPiMc) { // Lc->pKpi MC without ML, have one more axis for rotated candidates
            hRecoHelicityAxes.insert(hRecoHelicityAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
          } else {
            hRecoHelicityAxes.insert(hRecoHelicityAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisCharge});
          }
          registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for reconstructed prompt Lc+ candidates", HistType::kTHnSparseF, hRecoHelicityAxes);
          registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for reconstructed non-prompt Lc+ candidates", HistType::kTHnSparseF, hRecoHelicityAxes);
        }
        hHelicityaxes.insert(hHelicityaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
        registry.add("hHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores", HistType::kTHnSparseF, hHelicityaxes);

        if (activateTHnEulerPhiMonitor) {
          std::vector<AxisSpec> hEulerPhiAxes = {thnAxisInvMass, thnAxisPt, thnAxisInvMassKPiLc, thnAxisTHnAxisEulerPhi};
          if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
          }
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisResoChannelLc});
          }
          hEulerPhiAxes.push_back(thnAxisCharge);
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            registry.add("hRecPromptEulerPhiHelicity", "THn for polarisation studies with Euler phi w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
            registry.add("hRecNonPromptEulerPhiHelicity", "THn for polarisation studies with Euler phi w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
          } else {
            registry.add("hEulerPhiHelicity", "THn for polarisation studies with Euler phi w.r.t. helicity axis", HistType::kTHnSparseF, hEulerPhiAxes);
          }
        }
      }
      if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
        std::vector<AxisSpec> const hgenPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        std::vector<AxisSpec> const hgenNonPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarHelicity, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        registry.add("hGenPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
        registry.add("hGenNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hGenPartRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for partially reconstructed generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
          registry.add("hGenPartRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores for partially reconstructed generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        }
      }
    }

    if (activateTHnSparseCosThStarProduction) {
      std::vector<AxisSpec> hProductionaxes = {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY};
      if (doprocessDstar || doprocessDstarMc || doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarInPbPb || doprocessDstarMcInPbPb || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {

        hProductionaxes.insert(hProductionaxes.end(), {thnAxisInvMassD0, thnAxisCosThetaStarProduction});
        if (doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
          hProductionaxes.insert(hProductionaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTrackingSys) {
          hProductionaxes.insert(hProductionaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
        }
        if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb) {
          std::vector<AxisSpec> hRecoPromptProductionAxes(hProductionaxes);
          hRecoPromptProductionAxes.insert(hRecoPromptProductionAxes.end(), {thnAxisDauToMuons});
          std::vector<AxisSpec> hRecoNonPromptProductionAxes(hProductionaxes);
          hRecoNonPromptProductionAxes.insert(hRecoNonPromptProductionAxes.end(), {thnAxisDauToMuons, thnAxisPtB});
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptProductionAxes);
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptProductionAxes);
          if (activatePartRecoDstar) {
            registry.add("hPartRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for partially reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptProductionAxes);
            registry.add("hPartRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for partially reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptProductionAxes);
          }
        } else {
          if (nBkgRotations > 0) {
            hProductionaxes.push_back(thnAxisIsRotatedCandidate);
          }
          registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, hProductionaxes);
        }
      } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPi || doprocessLcToPKPiMc) {
        hProductionaxes.insert(hProductionaxes.end(), {thnAxisInvMassKPiLc, thnAxisCosThetaStarProduction});
        if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          hProductionaxes.insert(hProductionaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          std::vector<AxisSpec> hRecoProductionAxes(hProductionaxes);
          if (doprocessLcToPKPiMc) {
            hRecoProductionAxes.insert(hRecoProductionAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
          } else {
            hRecoProductionAxes.insert(hRecoProductionAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisCharge});
          }
          registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for reconstructed prompt Lc+ candidates", HistType::kTHnSparseF, hRecoProductionAxes);
          registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for reconstructed non-prompt Lc+ candidates", HistType::kTHnSparseF, hRecoProductionAxes);
        }
        hProductionaxes.insert(hProductionaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
        registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores", HistType::kTHnSparseF, hProductionaxes);
        if (activateTHnEulerPhiMonitor) {
          std::vector<AxisSpec> hEulerPhiAxes = {thnAxisInvMass, thnAxisPt, thnAxisInvMassKPiLc, thnAxisTHnAxisEulerPhi};
          if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
          }
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisResoChannelLc});
          }
          hEulerPhiAxes.push_back(thnAxisCharge);
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            registry.add("hRecPromptEulerPhiProduction", "THn for polarisation studies with Euler phi w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
            registry.add("hRecNonPromptEulerPhiProduction", "THn for polarisation studies with Euler phi w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
          } else {
            registry.add("hEulerPhiProduction", "THn for polarisation studies with Euler phi w.r.t. production axis", HistType::kTHnSparseF, hEulerPhiAxes);
          }
        }
      }
      if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
        std::vector<AxisSpec> const hgenPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarProduction, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        std::vector<AxisSpec> const hgenNonPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarProduction, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        registry.add("hGenPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
        registry.add("hGenNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hGenPartRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for partially reconstructed generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
          registry.add("hGenPartRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores for partially reconstructed generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        }
      }
    }

    if (activateTHnSparseCosThStarBeam) {
      std::vector<AxisSpec> hBeamaxes = {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY};
      if (doprocessDstar || doprocessDstarMc || doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarInPbPb || doprocessDstarMcInPbPb || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
        hBeamaxes.insert(hBeamaxes.end(), {thnAxisInvMassD0, thnAxisCosThetaStarBeam});
        if (doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
          hBeamaxes.insert(hBeamaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTrackingSys) {
          hBeamaxes.insert(hBeamaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
        }
        if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb) {
          std::vector<AxisSpec> hRecoPromptBeamAxes(hBeamaxes);
          hRecoPromptBeamAxes.insert(hRecoPromptBeamAxes.end(), {thnAxisDauToMuons});
          std::vector<AxisSpec> hRecoNonPromptBeamAxes(hBeamaxes);
          hRecoNonPromptBeamAxes.insert(hRecoNonPromptBeamAxes.end(), {thnAxisDauToMuons, thnAxisPtB});
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptBeamAxes);
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptBeamAxes);
          if (activatePartRecoDstar) {
            registry.add("hPartRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for partially reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptBeamAxes);
            registry.add("hPartRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for partially reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptBeamAxes);
          }
        } else {
          if (nBkgRotations > 0) {
            hBeamaxes.push_back(thnAxisIsRotatedCandidate);
          }
          registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, hBeamaxes);
        }
      } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPi || doprocessLcToPKPiMc) {
        hBeamaxes.insert(hBeamaxes.end(), {thnAxisInvMassKPiLc, thnAxisCosThetaStarBeam});
        if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          hBeamaxes.insert(hBeamaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          std::vector<AxisSpec> hRecoBeamAxes(hBeamaxes);
          if (doprocessLcToPKPiMc) {
            hRecoBeamAxes.insert(hRecoBeamAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
          } else {
            hRecoBeamAxes.insert(hRecoBeamAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisCharge});
          }
          registry.add("hRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for reconstructed prompt Lc+ candidates", HistType::kTHnSparseF, hRecoBeamAxes);
          registry.add("hRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for reconstructed non-prompt Lc+ candidates", HistType::kTHnSparseF, hRecoBeamAxes);
        }
        hBeamaxes.insert(hBeamaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
        registry.add("hBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores", HistType::kTHnSparseF, hBeamaxes);
        if (activateTHnEulerPhiMonitor) {
          std::vector<AxisSpec> hEulerPhiAxes = {thnAxisInvMass, thnAxisPt, thnAxisInvMassKPiLc, thnAxisTHnAxisEulerPhi};
          if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
          }
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            hEulerPhiAxes.insert(hEulerPhiAxes.end(), {thnAxisResoChannelLc});
          }
          hEulerPhiAxes.push_back(thnAxisCharge);
          if (doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPiMc) {
            registry.add("hRecPromptEulerPhiBeam", "THn for polarisation studies with Euler phi w.r.t. beam axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
            registry.add("hRecNonPromptEulerPhiBeam", "THn for polarisation studies with Euler phi w.r.t. beam axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, hEulerPhiAxes);
          } else {
            registry.add("hEulerPhiBeam", "THn for polarisation studies with Euler phi w.r.t. beam axis", HistType::kTHnSparseF, hEulerPhiAxes);
          }
        }
      }
      if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
        std::vector<AxisSpec> const hgenPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarBeam, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        std::vector<AxisSpec> const hgenNonPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarBeam, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        registry.add("hGenPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
        registry.add("hGenNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hGenPartRecoPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for partially reconstructed generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
          registry.add("hGenPartRecoNonPromptBeam", "THn for polarisation studies with cosThStar w.r.t. beam axis and BDT scores for partially reconstructed generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        }
      }
    }

    if (activateTHnSparseCosThStarRandom) {
      std::vector<AxisSpec> hRandomaxes = {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY};
      if (doprocessDstar || doprocessDstarMc || doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarInPbPb || doprocessDstarMcInPbPb || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {

        hRandomaxes.insert(hRandomaxes.end(), {thnAxisInvMassD0, thnAxisCosThetaStarRandom});
        if (doprocessDstarWithMl || doprocessDstarMcWithMl || doprocessDstarWithMlInPbPb || doprocessDstarMcWithMlInPbPb) {
          hRandomaxes.insert(hRandomaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTrackingSys) {
          hRandomaxes.insert(hRandomaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
        }
        if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb) {
          std::vector<AxisSpec> hRecoPromptRandomAxes(hRandomaxes);
          hRecoPromptRandomAxes.insert(hRecoPromptRandomAxes.end(), {thnAxisDauToMuons});
          std::vector<AxisSpec> hRecoNonPromptRandomAxes(hRandomaxes);
          hRecoNonPromptRandomAxes.insert(hRecoNonPromptRandomAxes.end(), {thnAxisDauToMuons, thnAxisPtB});
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptRandomAxes);
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptRandomAxes);
          if (activatePartRecoDstar) {
            registry.add("hPartRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for partially reconstructed prompt D*+ candidates", HistType::kTHnSparseF, hRecoPromptRandomAxes);
            registry.add("hPartRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for partially reconstructed non-prompt D*+ candidates", HistType::kTHnSparseF, hRecoNonPromptRandomAxes);
          }
        } else {
          if (nBkgRotations > 0) {
            hRandomaxes.push_back(thnAxisIsRotatedCandidate);
          }
          registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, hRandomaxes);
        }
      } else if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl || doprocessLcToPKPi || doprocessLcToPKPiMc) {
        hRandomaxes.insert(hRandomaxes.end(), {thnAxisInvMassKPiLc, thnAxisCosThetaStarRandom});
        if (doprocessLcToPKPiWithMl || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          hRandomaxes.insert(hRandomaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
          std::vector<AxisSpec> hRecoRandomAxes(hRandomaxes);
          if (doprocessLcToPKPiMc) {
            hRecoRandomAxes.insert(hRecoRandomAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
          } else {
            hRecoRandomAxes.insert(hRecoRandomAxes.end(), {thnAxisResoChannelLc, thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisCharge});
          }
          registry.add("hRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for reconstructed prompt Lc+ candidates", HistType::kTHnSparseF, hRecoRandomAxes);
          registry.add("hRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for reconstructed non-prompt Lc+ candidates", HistType::kTHnSparseF, hRecoRandomAxes);
        }
        hRandomaxes.insert(hRandomaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin, thnAxisIsRotatedCandidate, thnAxisCharge});
        registry.add("hRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores", HistType::kTHnSparseF, hRandomaxes);
      }
      if (doprocessDstarMc || doprocessDstarMcWithMl || doprocessDstarMcInPbPb || doprocessDstarMcWithMlInPbPb || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl || doprocessLcToPKPiBackgroundMcWithMl) {
        std::vector<AxisSpec> const hgenPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarRandom, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        std::vector<AxisSpec> const hgenNonPromptAxes = {thnAxisPt, thnAxisNumPvContributors, thnAxisY, thnAxisCosThetaStarRandom, thnAxisPtB, thnAxisDausAcc, thnAxisResoChannelLc, thnAxisCharge};
        registry.add("hGenPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
        registry.add("hGenNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hGenPartRecoPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for partially reconstructed generated prompt D*+ candidates", HistType::kTHnSparseF, hgenPromptAxes);
          registry.add("hGenPartRecoNonPromptRandom", "THn for polarisation studies with cosThStar w.r.t. random axis and BDT scores for partially reconstructed generated non-prompt D*+ candidates", HistType::kTHnSparseF, hgenNonPromptAxes);
        }
      }
    }

    if (activateTHnSparseCosThStarEP && !(doprocessDstarInPbPb || doprocessDstarWithMlInPbPb)) {
      LOGP(fatal, "THnSparse with cosThStar w.r.t. event plane axis is not supported for pp analysis, please check the configuration!");
    } else if (activateTHnSparseCosThStarEP) {
      std::vector<AxisSpec> hEPaxes = {thnAxisInvMass, thnAxisPt, thnAxisNumPvContributors, thnAxisY};
      if (doprocessDstarInPbPb || doprocessDstarWithMlInPbPb) {
        hEPaxes.insert(hEPaxes.end(), {thnAxisInvMassD0, thnAxisCosThetaStarEP});
        if (doprocessDstarWithMlInPbPb) {
          hEPaxes.insert(hEPaxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
        }
        if (activateTrackingSys) {
          hEPaxes.insert(hEPaxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
        }
        if (nBkgRotations > 0) {
          hEPaxes.push_back(thnAxisIsRotatedCandidate);
        }
        registry.add("hEP", "THn for polarisation studies with cosThStar w.r.t. event plane axis and BDT scores", HistType::kTHnSparseF, hEPaxes);
      }
    }

    /// control plots for Lc->pKPi
    if ((doprocessLcToPKPi || doprocessLcToPKPiWithMl || doprocessLcToPKPiMc || doprocessLcToPKPiMcWithMl) && lcPKPiChannels.activateTHnLcChannelMonitor) {
      registry.add("hMass2PairsLcPKPi", "THnSparse to monitor M2(Kpi), M2(pK), M2(ppi), pt correlation for Lc -> pKpi", HistType::kTHnSparseF, {thnAxisInvMass2KPiLcMonitoring, thnAxisInvMass2PKLcMonitoring, thnAxisInvMass2PPiLcMonitoring, thnAxisPt});
    }

    /// Event-plane related histograms
    if (doprocessResolEventPlane) {
      const AxisSpec axisCosDeltaPhi{1000, -1., 1., "cos(2(#Psi_{2}(A) #minus #Psi_{2}(B)))"};
      const AxisSpec axisPsi{180, -o2::constants::math::PIHalf, o2::constants::math::PIHalf, "#Psi_{2}"};

      registry.add("resolEvPlane/hEvPlaneAngleFV0A", ";centrality;#Psi_{2} (FV0A)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleFT0A", ";centrality;#Psi_{2} (FT0A)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleFT0C", ";centrality;#Psi_{2} (FT0C)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleFT0M", ";centrality;#Psi_{2} (FT0M)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleTPCpos", ";centrality;#Psi_{2} (TPC pos)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleTPCneg", ";centrality;#Psi_{2} (TPC neg)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});
      registry.add("resolEvPlane/hEvPlaneAngleTPCtot", ";centrality;#Psi_{2} (TPC)", {HistType::kTH2F, {thnAxisCentrality, axisPsi}});

      registry.add("resolEvPlane/hResolEvPlaneFT0CFT0A", "hResolEvPlaneFT0CFT0A; centrality; cos(2(#Psi_{2}(FT0C) #minus #Psi_{2}(FT0A)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0CFV0A", "hResolEvPlaneFT0CFV0A; centrality; cos(2(#Psi_{2}(FT0C) #minus #Psi_{2}(FV0A)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0CTPCpos", "hResolEvPlaneFT0CTPCpos; centrality; cos(2(#Psi_{2}(FT0C) #minus #Psi_{2}(TPC pos)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0CTPCneg", "hResolEvPlaneFT0CTPCneg; centrality; cos(2(#Psi_{2}(FT0C) #minus #Psi_{2}(TPC neg)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0CTPCtot", "hResolEvPlaneFT0CTPCtot; centrality; cos(2(#Psi_{2}(FT0C) #minus #Psi_{2}(TPC tot)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0AFV0A", "hResolEvPlaneFT0AFV0A; centrality; cos(2(#Psi_{2}(FT0A) #minus #Psi_{2}(FV0A)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0ATPCpos", "hResolEvPlaneFT0ATPCpos; centrality; cos(2(#Psi_{2}(FT0A) #minus #Psi_{2}(TPC pos)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0ATPCneg", "hResolEvPlaneFT0ATPCneg; centrality; cos(2(#Psi_{2}(FT0A) #minus #Psi_{2}(TPC neg)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0ATPCtot", "hResolEvPlaneFT0ATPCtot; centrality; cos(2(#Psi_{2}(FT0A) #minus #Psi_{2}(TPC)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0MFV0A", "hResolEvPlaneFT0MFV0A; centrality; cos(2(#Psi_{2}(FV0A) #minus #Psi_{2}(FV0A)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0MTPCpos", "hResolEvPlaneFT0MTPCpos; centrality; cos(2(#Psi_{2}(FT0M) #minus #Psi_{2}(TPC pos)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0MTPCneg", "hResolEvPlaneFT0MTPCneg; centrality; cos(2(#Psi_{2}(FT0M) #minus #Psi_{2}(TPC neg)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFT0MTPCtot", "hResolEvPlaneFT0MTPCtot; centrality; cos(2(#Psi_{2}(FT0M) #minus #Psi_{2}(TPC tot)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFV0ATPCpos", "hResolEvPlaneFV0ATPCpos; centrality; cos(2(#Psi_{2}(FV0A) #minus #Psi_{2}(TPC pos)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFV0ATPCneg", "hResolEvPlaneFV0ATPCneg; centrality; cos(2(#Psi_{2}(FV0A) #minus #Psi_{2}(TPC neg)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneFV0ATPCtot", "hResolEvPlaneFV0ATPCtot; centrality; cos(2(#Psi_{2}(FV0A) #minus #Psi_{2}(TPC)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
      registry.add("resolEvPlane/hResolEvPlaneTPCposTPCneg", "hResolEvPlaneTPCposTPCneg; centrality; cos(2(#Psi_{2}(TPC pos) #minus #Psi_{2}(TPC neg)))", {HistType::kTH2F, {thnAxisCentrality, axisCosDeltaPhi}});
    }

    // inv. mass hypothesis to loop over
    // e.g.: Lc->pKpi has the ambiguity pKpi vs. piKp
    if (doprocessLcToPKPi || doprocessLcToPKPiWithMl) {
      nMassHypos = charm_polarisation::MassHyposLcToPKPi::NMassHypoLcToPKPi;
    } else {
      // D*, Lc->pK0s
      nMassHypos = 1;
    }

    if (doprocessResolEventPlane) {
      int dummyVariable;
      hfEvSel.init(registry, dummyVariable);
      ccdb->setURL("http://alice-ccdb.cern.ch");
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
    }
  }; // end init

  /// \param invMassCharmHad is the invariant-mass of the candidate
  /// \param ptCharmHad is the pt of the candidate
  /// \param numPvContributors is the number of PV contributors
  /// \param rapCharmHad is the rapidity of the candidate
  /// \param invMassD0 is the invariant-mass of the D0 daugher (only for D*+)
  /// \param invMassKPiLc is the invariant-mass of the K-pi pair (only for Lc+)
  /// \param invMassPKLc is the invariant-mass of the p-K pair (only for Lc+)
  /// \param invMassPPiLc is the invariant-mass of the p-pi pair (only for Lc+)
  /// \param cosThetaStar is the cosThetaStar of the candidate
  /// \param outputMl is the array with ML output scores
  /// \param isRotatedCandidate is a flag that keeps the info of the rotation of the candidate for bkg studies
  /// \param origin is the MC origin
  /// \param ptBhadMother is the pt of the b-hadron mother (only in case of non-prompt)
  /// \param resoChannelLc indicates the Lc decay channel (direct, resonant)
  /// \param absEtaMin is the minimum absolute eta of the daughter tracks
  /// \param numItsClsMin is the minimum number of ITS clusters of the daughter tracks
  /// \param numTpcClsMin is the minimum number of TPC clusters of the daughter tracks
  /// \param charge is the charge of the hadron
  /// \param nMuons is the number of muons from daughter decays
  /// \param isPartRecoDstar is a flag indicating if it is a partly reconstructed Dstar meson (MC only)
  template <charm_polarisation::DecayChannel Channel, bool WithMl, bool DoMc, charm_polarisation::CosThetaStarType CosThetaStarType>
  void fillRecoHistos(float invMassCharmHad, float ptCharmHad, int numPvContributors, float rapCharmHad, float invMassD0, float invMassKPiLc, float cosThetaStar, float phiEuler, std::array<float, 3> outputMl, int isRotatedCandidate, int8_t origin, float ptBhadMother, int8_t resoChannelLc, float absEtaMin, int numItsClsMin, int numTpcClsMin, int8_t charge, int8_t nMuons, bool isPartRecoDstar)
  {
    if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if constexpr (!DoMc) {                                                            // data
        if constexpr (WithMl) {                                                         // with ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
              } else {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], charge);
            }
          }
        } else {                                                                       // without ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, isRotatedCandidate);
              } else {
                registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, charge);
            }
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (WithMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiHelicity"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          }
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if constexpr (!DoMc) {                                                                     // data
        if constexpr (WithMl) {                                                                  // with ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {           // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
              } else {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], charge);
            }
          }
        } else {                                                                       // without ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, isRotatedCandidate);
              } else {
                registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, charge);
            }
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (WithMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiProduction"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          }
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if constexpr (!DoMc) {                                                               // data
        if constexpr (WithMl) {                                                            // with ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {     // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
              } else {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], charge);
            }
          }
        } else {                                                                       // without ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, isRotatedCandidate);
              } else {
                registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hBeam"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
            if (activateTHnEulerPhiMonitor) {
              registry.fill(HIST("hEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, charge);
            }
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (WithMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, charge);
              }
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecPromptEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptBeam"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
              if (activateTHnEulerPhiMonitor) {
                registry.fill(HIST("hRecNonPromptEulerPhiBeam"), invMassCharmHad, ptCharmHad, invMassKPiLc, phiEuler, resoChannelLc, charge);
              }
            }
          }
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if constexpr (!DoMc) {                                                                 // data
        if constexpr (WithMl) {                                                              // with ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {       // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
              } else {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
          }
        } else {                                                                       // without ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, isRotatedCandidate);
              } else {
                registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar);
              }
            }
          } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
            registry.fill(HIST("hRandom"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassKPiLc, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, charge);
          }
        }
      } else {                                                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (WithMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                }
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
            }
          }
        } else {                                                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
              } else {
                registry.fill(HIST("hPartRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
              } else {
                registry.fill(HIST("hPartRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
              }
            } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+
              registry.fill(HIST("hRecoNonPromptRandom"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassKPiLc, cosThetaStar, resoChannelLc, absEtaMin, numItsClsMin, numTpcClsMin, charge);
            }
          }
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::EP) { // EP
      if constexpr (!DoMc) {                                                             // data
        if constexpr (WithMl) {                                                          // with ML
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {   // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], isRotatedCandidate);
              } else {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2]);
              }
            }
          }
        } else {
          if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
            if (activateTrackingSys) {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate);
              } else {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin);
              }
            } else {
              if (nBkgRotations > 0) {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar, isRotatedCandidate);
              } else {
                registry.fill(HIST("hEP"), invMassCharmHad, ptCharmHad, numPvContributors, std::abs(rapCharmHad), invMassD0, cosThetaStar);
              }
            }
          }
        }
      } else {
        if constexpr (WithMl) {                                                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) {                                 // prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                } else {
                  registry.fill(HIST("hPartRecoPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons);
                }
              }
            }
          } else {                                                                       // non-prompt
            if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
              if (activateTrackingSys) {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], absEtaMin, numItsClsMin, numTpcClsMin, nMuons, ptBhadMother);
                }
              } else {
                if (!isPartRecoDstar) {
                  registry.fill(HIST("hRecoNonPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                } else {
                  registry.fill(HIST("hPartRecoNonPromptEP"), invMassCharmHad, ptCharmHad, numPvContributors, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, ptBhadMother);
                }
              }
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
  /// \param isPartRecoDstar is a flag indicating if it is a partly reconstructed Dstar->D0pi->Kpipipi0 meson (MC only)
  template <charm_polarisation::CosThetaStarType CosThetaStarType>
  void fillGenHistos(float ptCharmHad, int numPvContributors, float rapCharmHad, float cosThetaStar, int8_t origin, float ptBhadMother, bool areDausInAcc, uint8_t resoChannelLc, int8_t charge, bool isPartRecoDstar)
  {
    if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if (origin == RecoDecay::OriginType::Prompt) {                                    // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptHelicity"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if (origin == RecoDecay::OriginType::Prompt) {                                             // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptProduction"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Beam) { // Beam
      if (origin == RecoDecay::OriginType::Prompt) {                                       // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptBeam"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::Random) { // Random
      if (origin == RecoDecay::OriginType::Prompt) {                                         // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptRandom"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        }
      }
    } else if constexpr (CosThetaStarType == charm_polarisation::CosThetaStarType::EP) { // EP
      if (origin == RecoDecay::OriginType::Prompt) {                                     // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptEP"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoPromptEP"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, areDausInAcc, resoChannelLc, charge);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptEP"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptEP"), ptCharmHad, numPvContributors, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, resoChannelLc, charge);
        }
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
  template <charm_polarisation::DecayChannel Channel>
  bool isInSignalRegion(float invMass)
  {
    float invMassMin;
    float invMassMax;
    if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
      invMassMin = 0.142f;
      invMassMax = 0.15f;
      if (invMassMin < invMass && invMass < invMassMax) {
        return true;
      }
    } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) { // Lc+ (to be tuned!)
      invMassMin = 2.25f;
      invMassMax = 2.35f;
      if (invMassMin < invMass && invMass < invMassMax) {
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
    etaMin = 10.f;
    nItsClsMin = 10;
    nTpcClsMin = 1000;

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
  template <charm_polarisation::DecayChannel Channel, bool WithMl, bool DoMc, bool StudyLcPkPiBkgMc = false, bool WithEp = false, typename Cand, typename Part, typename Trk, typename QVecs = void>
  bool runPolarisationAnalysis(Cand const& candidate, int bkgRotationId, int numPvContributors, Part const& particles, Trk const& /*tracks*/, QVecs const* qVecs = nullptr)
  {
    if constexpr (WithEp) {
      assert(qVecs && "EP analysis requested but qVecs == nullptr");
    }

    constexpr std::size_t NScores{3u};

    bool isCandidateInSignalRegion{false};
    int8_t origin{RecoDecay::OriginType::None};
    int8_t massHypoMcTruth{-1};
    float ptBhadMother{-1.f};
    int8_t resoChannelLc = -1;
    int8_t charge = -99;
    bool partRecoDstar{false};
    if constexpr (DoMc) {
      if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        partRecoDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPiPi0 && std::abs(candidate.flagMcMatchRecD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0;
        bool const signalDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi && std::abs(candidate.flagMcMatchRecD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
        if (!signalDstar && (!partRecoDstar || !activatePartRecoDstar)) { // this candidate is not signal and not partially reconstructed signal, skip
          return isCandidateInSignalRegion;
        }
        origin = candidate.originMcRec();
        ptBhadMother = candidate.ptBhadMotherPart();
        int const pdgBhadMother = candidate.pdgBhadMotherPart();
        // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
        if (origin == RecoDecay::OriginType::NonPrompt && (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 557)) { // o2-linter: disable=pdg/explicit-code, magic-number (constants not in the PDG header)
          return isCandidateInSignalRegion;
        }
      } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if constexpr (!StudyLcPkPiBkgMc) {                                                                    // skip this if studyLcPKPiBkgMc is true, since we are interested in background
          if (std::abs(candidate.flagMcMatchRec()) != hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) { // this candidate is not signal, skip
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

        /// Lc electric charge from MC truth
        /// This is checked when the reconstructed 3-prong candidate is matched to MC with RecoDecay::getMatchedMCRec
        int8_t const flagMc = candidate.flagMcMatchRec();
        charge = std::abs(flagMc) > 0 ? flagMc / std::abs(flagMc) : 0; /// 0 should never happen, debug protection
      }
    } else {
      /// data
      if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
        /// Calculate the electric charge from reconstructed daughter tracks
        /// Lc charge == first daughter charge
        auto trackProng0 = candidate.template prong0_as<Trk>();
        charge = static_cast<int8_t>(trackProng0.sign());
      }
    }

    // loop over mass hypotheses
    for (uint8_t iMass = 0u; iMass < nMassHypos; iMass++) {

      // variable definition
      float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
      float pxCharmHad{-1000.f}, pyCharmHad{-1000.f}, pzCharmHad{-1000.f};
      double massDau{0.}, invMassCharmHad{0.}, invMassCharmHadForSparse{0.}, invMassD0{0.}, invMassKPiLc{0.}, invMassPKLc{0.}, invMassPPiLc{0.};
      float rapidity{-999.f};
      std::array<float, 3> outputMl{-1.f, -1.f, -1.f};
      int isRotatedCandidate = 0; // currently meaningful only for Lc->pKpi

      if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
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

        if constexpr (WithMl) {
          outputMl[0] = candidate.mlProbDstarToD0Pi()[0];
          outputMl[1] = candidate.mlProbDstarToD0Pi()[1];
          outputMl[2] = candidate.mlProbDstarToD0Pi()[2];
        }
      } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
        // Lc->pKpi analysis
        // polarization measured from the proton daughter (*)

        if constexpr (DoMc) { // we keep only the good hypo in the MC
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
        float invMassPiKPi = 0.f; // bkg. from D+ -> K+pi-pi-
        float invMassKKPi = 0.f;  // bkg. from D+, Ds+ -> K+K-pi+ (1st mass hypothesis)
        float invMassPiKK = 0.f;  // bkg. from D+, Ds+ -> pi+K-K+ (2nd mass hypothesis)
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
            invMassCharmHad = HfHelper::invMassLcToPKPi(candidate);
            invMassCharmHadForSparse = HfHelper::invMassLcToPKPi(candidate);
          }
          if constexpr (WithMl) {
            if (candidate.mlProbLcToPKPi().size() == NScores) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              // !!! NB: each rotated candidates inherits the BDT scores of the original candidate, even if the candidate pt changed after the rotation of the kaon-track pt !!!
              outputMl[0] = candidate.mlProbLcToPKPi()[0];
              outputMl[1] = candidate.mlProbLcToPKPi()[1];
              outputMl[2] = candidate.mlProbLcToPKPi()[2];
            }
          }
          // invariant mass of the KPi pair
          invMassKPiLc = HfHelper::invMassKPiPairLcToPKPi(candidate);
          invMassPKLc = HfHelper::invMassPKPairLcToPKPi(candidate);
          invMassPPiLc = HfHelper::invMassPPiPairLcToPKPi(candidate);

          // D+ and Ds+ invariant mass values, to put a veto on background sources
          invMassPiKPi = HfHelper::invMassDplusToPiKPi(candidate); // bkg. from D+ -> K+pi-pi-
          invMassKKPi = HfHelper::invMassDsToKKPi(candidate);      // bkg. from D+, Ds+ -> K+K-pi+ (1st mass hypothesis)
          invMassPiKK = HfHelper::invMassDsToPiKK(candidate);      // bkg. from D+, Ds+ -> pi+K-K+ (2nd mass hypothesis)

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
            invMassCharmHad = HfHelper::invMassLcToPiKP(candidate);
            invMassCharmHadForSparse = HfHelper::invMassLcToPiKP(candidate);
          }
          if constexpr (WithMl) {
            if (candidate.mlProbLcToPiKP().size() == NScores) {
              // protect from empty vectors
              // the BDT output score might be empty if no preselections were enabled (selectionFlag null)
              // !!! NB: each rotated candidates inherits the BDT scores of the original candidate, even if the candidate pt changed after the rotation of the kaon-track pt !!!
              outputMl[0] = candidate.mlProbLcToPiKP()[0];
              outputMl[1] = candidate.mlProbLcToPiKP()[1];
              outputMl[2] = candidate.mlProbLcToPiKP()[2];
            }
          }
          // invariant mass of the KPi pair
          invMassKPiLc = HfHelper::invMassKPiPairLcToPiKP(candidate);
          invMassPKLc = HfHelper::invMassPKPairLcToPiKP(candidate);
          invMassPPiLc = HfHelper::invMassPPiPairLcToPiKP(candidate);

          // D+ and Ds+ invariant mass values, to put a veto on background sources
          invMassPiKPi = HfHelper::invMassDplusToPiKPi(candidate); // bkg. from D+ -> K+pi-pi-
          invMassKKPi = HfHelper::invMassDsToKKPi(candidate);      // bkg. from D+, Ds+ -> K+K-pi+ (1st mass hypothesis)
          invMassPiKK = HfHelper::invMassDsToPiKK(candidate);      // bkg. from D+, Ds+ -> pi+K-K+ (2nd mass hypothesis)

        } else {
          // NB: no need to check cases in which candidate.isSelLcToPKPi() and candidate.isSelLcToPiKP() are both false, because they are rejected already by the Filter
          // ... but we need to put this protections here!
          // Otherwise, a candidate selected as pKpi only has invMassCharmHad==0 when iMass == charm_polarisation::MassHyposLcToPKPi::PiKP and viceversa
          continue;
        }

        /// put veto on D+, Ds+ inv. masses, to reduce the background
        if (lcBkgVeto.applyLcBkgVeto && ((lcBkgVeto.enableLcBkgVetoDplusKPiPi && lcBkgVeto.massDplusKPiPiMinVeto < invMassPiKPi && invMassPiKPi < lcBkgVeto.massDplusKPiPiMaxVeto) /*bkg. from D+ -> K+pi-pi-*/ ||
                                         (lcBkgVeto.enableLcBkgVetoDplusKKPi && lcBkgVeto.massDplusKKPiMinVeto < invMassKKPi && invMassKKPi < lcBkgVeto.massDplusKKPiMaxVeto) /*bkg. from D+ -> K+K-pi+ (1st mass hypothesis)*/ ||
                                         (lcBkgVeto.enableLcBkgVetoDplusKKPi && lcBkgVeto.massDplusKKPiMinVeto < invMassPiKK && invMassPiKK < lcBkgVeto.massDplusKKPiMaxVeto) /*bkg. from D+ -> K+K-pi+ (2nd mass hypothesis)*/ ||
                                         (lcBkgVeto.enableLcBkgVetoDsKKPi && lcBkgVeto.massDsKKPiMinVeto < invMassKKPi && invMassKKPi < lcBkgVeto.massDsKKPiMaxVeto) /*bkg. from Ds+ -> K+K-pi+ (1st mass hypothesis)*/ ||
                                         (lcBkgVeto.enableLcBkgVetoDsKKPi && lcBkgVeto.massDsKKPiMinVeto < invMassPiKK && invMassPiKK < lcBkgVeto.massDsKKPiMaxVeto)) /*bkg. from Ds+ -> K+K-pi+ (2nd mass hypothesis)*/) {
          /// this candidate has D+ and/or Ds+ in the veto range, let's reject it
          continue;
        }

        /// control plots on pair masses
        double const invMass2KPiLc = invMassKPiLc * invMassKPiLc;
        double const invMass2PKLc = invMassPKLc * invMassPKLc;
        double const invMass2PPiLc = invMassPPiLc * invMassPPiLc;
        if (lcPKPiChannels.activateTHnLcChannelMonitor && bkgRotationId == 0) {
          /// fill Dalitz plot only for genuine candidates (i.e. non-rotated)
          registry.fill(HIST("hMass2PairsLcPKPi"), invMass2KPiLc, invMass2PKLc, invMass2PPiLc, candidate.pt());
        }

        /// veto cut on pair masses
        if (lcPKPiChannels.applyLcSignalVeto && lcPKPiChannels.mass2PPiLcMinVeto < invMass2PPiLc && invMass2PPiLc < lcPKPiChannels.mass2PPiLcMaxVeto) {
          /// this candidate has a significant contribution from Lc+ -> Delta++ K-, let's reject it
          continue;
        }

      } // Lc->pKpi

      if (invMassCharmHadForSparse < minInvMass || invMassCharmHadForSparse > maxInvMass) {
        continue;
      }

      /// apply rapidity selection on the reconstructed candidate
      if (std::abs(rapidity) > rapidityCut) {
        continue;
      }

      float const phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
      float const thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
      ROOT::Math::PxPyPzMVector const fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
      ROOT::Math::PxPyPzMVector const fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, invMassCharmHad);
      ROOT::Math::Boost const boost{fourVecMother.BoostToCM()};
      ROOT::Math::PxPyPzMVector const fourVecDauCM = boost(fourVecDau);
      ROOT::Math::XYZVector const threeVecDauCM = fourVecDauCM.Vect();

      float const ptCharmHad = std::sqrt(pxCharmHad * pxCharmHad + pyCharmHad * pyCharmHad); // this definition is valid for both rotated and original candidates

      if (!isCandidateInSignalRegion) { // it could be that only one mass hypothesis is in signal region
        isCandidateInSignalRegion = isInSignalRegion<Channel>(invMassCharmHadForSparse);
      }

      float absEtaTrackMin{-1.f};
      int numItsClsMin{-1}, numTpcClsMin{-1};

      if (activateTrackingSys) {
        auto trackProng0 = candidate.template prong0_as<Trk>();
        auto trackProng1 = candidate.template prong1_as<Trk>();
        if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
          auto trackProng2 = candidate.template prongPi_as<Trk>();
          getTrackingInfos(std::vector{trackProng0, trackProng1, trackProng2}, absEtaTrackMin, numItsClsMin, numTpcClsMin);
        } else if (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
          auto trackProng2 = candidate.template prong2_as<Trk>();
          getTrackingInfos(std::vector{trackProng0, trackProng1, trackProng2}, absEtaTrackMin, numItsClsMin, numTpcClsMin);
        }
      }

      // helicity
      ROOT::Math::XYZVector const helicityVec = fourVecMother.Vect();
      float cosThetaStarHelicity = -10.f;
      float phiHelicity = -10.f;
      // production
      ROOT::Math::XYZVector const normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
      float cosThetaStarProduction = -10.f;
      float phiProduction = -10.f;
      // beam
      ROOT::Math::XYZVector const beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      float cosThetaStarBeam = -10.f;
      float phiBeam = -10.f;
      // random
      float cosThetaStarRandom = -10.f;

      int8_t nMuons{0u};
      if constexpr (DoMc) {
        nMuons = candidate.nTracksDecayed();
      }

      if constexpr (WithEp && !DoMc) {
        /// EP analysis
        float const xQvec = (*qVecs).at(0);
        float const yQvec = (*qVecs).at(1);
        ROOT::Math::XYZVector const qVecNorm = ROOT::Math::XYZVector(yQvec, -xQvec, 0.f);
        float const phiEP = -99.f;

        if (activateTHnSparseCosThStarEP) {
          // EP
          float cosThetaStarEP = qVecNorm.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(qVecNorm.Mag2());
          fillRecoHistos<Channel, WithMl, DoMc, charm_polarisation::CosThetaStarType::EP>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarEP, phiEP, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin, charge, nMuons, partRecoDstar);
        }
      }

      if (activateTHnSparseCosThStarHelicity) {
        // helicity
        cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
        phiHelicity = std::atan2(beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()), normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2())));
        fillRecoHistos<Channel, WithMl, DoMc, charm_polarisation::CosThetaStarType::Helicity>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarHelicity, phiHelicity, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin, charge, nMuons, partRecoDstar);
      }
      if (activateTHnSparseCosThStarProduction) {
        // production
        cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
        phiProduction = std::atan2(normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2())), helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2())));
        fillRecoHistos<Channel, WithMl, DoMc, charm_polarisation::CosThetaStarType::Production>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarProduction, phiProduction, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin, charge, nMuons, partRecoDstar);
      }
      if (activateTHnSparseCosThStarBeam) {
        // beam
        cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        phiBeam = std::atan2(helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2())), beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()));
        fillRecoHistos<Channel, WithMl, DoMc, charm_polarisation::CosThetaStarType::Beam>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarBeam, phiBeam, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin, charge, nMuons, partRecoDstar);
      }
      if (activateTHnSparseCosThStarRandom) {
        // random
        ROOT::Math::XYZVector const randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
        fillRecoHistos<Channel, WithMl, DoMc, charm_polarisation::CosThetaStarType::Random>(invMassCharmHadForSparse, ptCharmHad, numPvContributors, rapidity, invMassD0, invMassKPiLc, cosThetaStarRandom, -99.f, outputMl, isRotatedCandidate, origin, ptBhadMother, resoChannelLc, absEtaTrackMin, numItsClsMin, numTpcClsMin, charge, nMuons, partRecoDstar);
      }

      /// Table for Lc->pKpi background studies
      /// Defined only in MC simulations, to study resonances and reflected signal
      if constexpr (DoMc && Channel == charm_polarisation::DecayChannel::LcToPKPi) {
        if constexpr (StudyLcPkPiBkgMc) {
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
          int8_t originProng0 = -1;
          int8_t originProng1 = -1;
          int8_t originProng2 = -1;
          std::vector<int> idxBhadMothersProng0{};
          std::vector<int> idxBhadMothersProng1{};
          std::vector<int> idxBhadMothersProng2{};
          if (trackProng0.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng0 = trackProng0.template mcParticle_as<Part>();
            pdgProng0 = particleProng0.pdgCode();
            originProng0 = RecoDecay::getCharmHadronOrigin(particles, particleProng0, false, &idxBhadMothersProng0);
          }
          if (trackProng1.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng1 = trackProng1.template mcParticle_as<Part>();
            pdgProng1 = particleProng1.pdgCode();
            originProng1 = RecoDecay::getCharmHadronOrigin(particles, particleProng1, false, &idxBhadMothersProng1);
          }
          if (trackProng2.has_mcParticle()) {
            /// BEWARE: even when grouping by mcCollision, mcParticle_as<> gets the mcParticle even if it belongs to a different mcCollision
            /// because _as<> works with unbound tables. (*)
            auto particleProng2 = trackProng2.template mcParticle_as<Part>();
            pdgProng2 = particleProng2.pdgCode();
            originProng2 = RecoDecay::getCharmHadronOrigin(particles, particleProng2, false, &idxBhadMothersProng2);
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

          /// check the origin (prompt, non-prompt of the triplet)
          /// need to check each prong, since they might come from combinatorial background
          /// convention:
          ///  - all 3 prongs from the same B hadron: non-prompt
          ///  - all 3 prongs claimed to be prompt: prompt  -->  check on same charm mother done offline with prong PDG daughters (more difficult for beauty, due to more intermediate resonances) and checking that the distribution peaks somehow (otherwise: combinatorial background)
          ///  - otherwise: none
          int8_t originTriplet = RecoDecay::OriginType::None;
          if (originProng0 == RecoDecay::OriginType::Prompt && originProng1 == RecoDecay::OriginType::Prompt && originProng2 == RecoDecay::OriginType::Prompt) {
            /// we claim this triplet as prong w/o checking if all triplets have the same mother
            originTriplet = RecoDecay::OriginType::Prompt;
          } else if (originProng0 == RecoDecay::OriginType::NonPrompt && originProng1 == RecoDecay::OriginType::NonPrompt && originProng2 == RecoDecay::OriginType::NonPrompt) {
            /// check if the three particles share the same B-hadron id. If yes: claim the triplet as "non-prompt"
            int const idBMotherProng0 = idxBhadMothersProng0.at(0);
            int const idBMotherProng1 = idxBhadMothersProng1.at(0);
            int const idBMotherProng2 = idxBhadMothersProng2.at(0);
            if (idBMotherProng0 == idBMotherProng1 && idBMotherProng1 == idBMotherProng2) {
              originTriplet = RecoDecay::OriginType::NonPrompt;
            }
          }

          /// check if the pKpi triplet is a Lc->pKpi
          int8_t isRealLcPKPi = 0;
          if (isRealPKPi && (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi)) {
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
          bool const atLeast2ProngsFromSameMother = (!idMothersProng0.empty() && !idMothersProng1.empty() && idMothersProng0.at(0) == idMothersProng1.at(0)) ||
                                                    (!idMothersProng1.empty() && !idMothersProng2.empty() && idMothersProng1.at(0) == idMothersProng2.at(0)) ||
                                                    (!idMothersProng0.empty() && !idMothersProng2.empty() && idMothersProng0.at(0) == idMothersProng2.at(0));
          if (atLeast2ProngsFromSameMother) {
            if (!idMothersProng0.empty()) {
              /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
              /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
              auto mother = particles.rawIteratorAt(idMothersProng0.at(0) - particles.offset());
              pdgMotherProng0 = std::abs(mother.pdgCode()); // PDG code of the mother
            }
            if (!idMothersProng1.empty()) {
              /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
              /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
              auto mother = particles.rawIteratorAt(idMothersProng1.at(0) - particles.offset());
              pdgMotherProng1 = std::abs(mother.pdgCode()); // PDG code of the mother
            }
            if (!idMothersProng2.empty()) {
              /// BEWARE: in case of mcCollision grouping, the idMother can anyway point to a particle in another collision (*)
              /// therefore the rawIteratorAt call might crash the code because one goes above the (grouped) particles table size
              auto mother = particles.rawIteratorAt(idMothersProng2.at(0) - particles.offset());
              pdgMotherProng2 = std::abs(mother.pdgCode()); // PDG code of the mother
            }
          }

          /// calculate inv. masses for pairs, depending on mass hypothesis
          std::array<float, 3> pVecPion = {};
          std::array<float, 3> const pVecKaon = candidate.pVectorProng1();
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
                       outputMl.at(0), outputMl.at(2),
                       isRealPKPi, isRealLcPKPi, isReflected,
                       charge, originTriplet);
        } // end studyLcPKPiBkgMc
      } // end table for Lc->pKpi background studies

    } /// end loop over mass hypotheses

    return isCandidateInSignalRegion;
  }

  /// \param mcParticle is the Mc particle
  /// \param mcParticles is the table of Mc particles
  /// \param numPvContributors is the number of PV contributors in the associated reco collision
  template <charm_polarisation::DecayChannel Channel, bool WithCent = false, typename Part, typename Particles, typename Cent = void>
  void runMcGenPolarisationAnalysis(Part const& mcParticle, Particles const& mcParticles, int numPvContributors, Cent const* centrality = nullptr)
  {
    if constexpr (WithCent) {
      assert(qVecs && "Centrality analysis requested but Cent == nullptr");
    }
    if constexpr (WithCent) {
      if (*centrality < centralityMin || *centrality > centralityMax) {
        return; // skip this collision if outside of the centrality range
      }
    }

    int8_t origin{RecoDecay::OriginType::None};
    std::vector<int> listDaughters{};
    float massDau{0.f}, massCharmHad{0.f};
    float ptBhadMother{-1.f};
    bool areDauInAcc{true};
    int8_t resoChannelLc = -1;
    int8_t charge = -99;
    bool partRecoDstar{false};
    if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
      partRecoDstar = (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPiPi0) && (std::abs(mcParticle.flagMcMatchGenD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0);
      bool const signalDstar = (std::abs(mcParticle.flagMcMatchGen()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi) && (std::abs(mcParticle.flagMcMatchGenD0()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);

      if (!signalDstar && (!activatePartRecoDstar || !partRecoDstar)) { // this particle is not signal and not partially reconstructed signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      if (origin == RecoDecay::OriginType::NonPrompt) {
        auto bHadMother = mcParticles.rawIteratorAt(mcParticle.idxBhadMotherPart() - mcParticles.offset());
        int const pdgBhadMother = std::abs(bHadMother.pdgCode());
        // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
        if (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 557) { // o2-linter: disable=pdg/explicit-code, magic-number (constants not in the PDG header)
          return;
        }
        ptBhadMother = bHadMother.pt();
      }

      std::array<int, 2> const dauPdgs = {kPiPlus, o2::constants::physics::Pdg::kD0};
      RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 1);
      massDau = massPi;
      massCharmHad = massDstar;
    } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
      if (std::abs(mcParticle.flagMcMatchGen()) != hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) { // this particle is not signal, skip
        return;
      }
      origin = mcParticle.originMcGen();
      resoChannelLc = mcParticle.flagMcDecayChanGen();
      std::array<int, 3> const dauPdgs = {kProton, -kKPlus, kPiPlus};
      RecoDecay::getDaughters(mcParticle, &listDaughters, dauPdgs, 2);
      massDau = massProton;
      massCharmHad = massLc;

      /// electric charge from PDG code
      int const pdgCode = mcParticle.pdgCode();
      charge = static_cast<int8_t>(pdgCode / std::abs(pdgCode));
    }

    float const rapidity = mcParticle.y();
    if (std::abs(rapidity) > 1.f) { // we do not keep particles with |y| > 1
      return;
    }

    float const pxCharmHad = mcParticle.px();
    float const pyCharmHad = mcParticle.py();
    float const pzCharmHad = mcParticle.pz();
    float const ptCharmHad = mcParticle.pt();

    float pxDau{-1000.f}, pyDau{-1000.f}, pzDau{-1000.f};
    for (const auto& dauIdx : listDaughters) {
      auto dauPart = mcParticles.rawIteratorAt(dauIdx - mcParticles.offset());
      if constexpr (Channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        if (std::abs(dauPart.pdgCode()) == kPiPlus) {
          pxDau = dauPart.px();
          pyDau = dauPart.py();
          pzDau = dauPart.pz();
          if (areDauInAcc) {
            areDauInAcc = isDaughterInAcceptance(dauPart, 0.1, 0.8);
          }
        } else if (areDauInAcc) { // check also D0 daughters
          std::vector<int> listDaughtersD0{};
          std::array<int, 2> const dauPdgsD0 = {kPiPlus, -kKPlus};
          RecoDecay::getDaughters(mcParticle, &listDaughtersD0, dauPdgsD0, 1);
          for (const auto& dauIdxD0 : listDaughtersD0) {
            auto dauPartD0 = mcParticles.rawIteratorAt(dauIdxD0 - mcParticles.offset());
            areDauInAcc = isDaughterInAcceptance(dauPartD0, 0.3, 0.8);
          }
        }
      } else if constexpr (Channel == charm_polarisation::DecayChannel::LcToPKPi) {
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

    float const phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
    float const thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
    ROOT::Math::PxPyPzMVector const fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
    ROOT::Math::PxPyPzMVector const fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, massCharmHad);
    ROOT::Math::Boost const boost{fourVecMother.BoostToCM()};
    ROOT::Math::PxPyPzMVector const fourVecDauCM = boost(fourVecDau);
    ROOT::Math::XYZVector const threeVecDauCM = fourVecDauCM.Vect();

    if (activateTHnSparseCosThStarHelicity) {
      ROOT::Math::XYZVector const helicityVec = fourVecMother.Vect();
      float const cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Helicity>(ptCharmHad, numPvContributors, rapidity, cosThetaStarHelicity, origin, ptBhadMother, areDauInAcc, resoChannelLc, charge, partRecoDstar);
    }
    if (activateTHnSparseCosThStarProduction) {
      ROOT::Math::XYZVector const normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
      float const cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Production>(ptCharmHad, numPvContributors, rapidity, cosThetaStarProduction, origin, ptBhadMother, areDauInAcc, resoChannelLc, charge, partRecoDstar);
    }
    if (activateTHnSparseCosThStarBeam) {
      ROOT::Math::XYZVector const beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
      float const cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Beam>(ptCharmHad, numPvContributors, rapidity, cosThetaStarBeam, origin, ptBhadMother, areDauInAcc, resoChannelLc, charge, partRecoDstar);
    }
    if (activateTHnSparseCosThStarRandom) {
      ROOT::Math::XYZVector const randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
      float const cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());
      fillGenHistos<charm_polarisation::CosThetaStarType::Random>(ptCharmHad, numPvContributors, rapidity, cosThetaStarRandom, origin, ptBhadMother, areDauInAcc, resoChannelLc, charge, partRecoDstar);
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
      int const numPvContributors = collision.numContrib();
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
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with ML cuts
  void processDstarWithMl(aod::Collisions const& collisions,
                          FilteredCandDstarWSelFlagAndMl const& dstarCandidates,
                          TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }

        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(dstarCandidate, iRotation, numPvContributors, -1 /*MC particles*/, tracks);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarWithMl, "Process Dstar candidates with ML", false);

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
      int const numPvContributors = collision.numContrib();
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
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarMc, "Process Dstar candidates in MC without ML", false);

  // Dstar in MC with ML cuts
  void processDstarMcWithMl(aod::McCollisions::iterator const&,
                            McParticlesDstarMatched const& mcParticles,
                            CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                            FilteredCandDstarWSelFlagAndMcAndMl const& dstarCandidates,
                            TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarMcWithMl, "Process Dstar candidates in MC with ML", false);

  void processDstarInPbPb(CollsWithQVecs const& collisions,
                          FilteredCandDstarWSelFlag const& dstarCandidates,
                          TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      const auto centrality = getCentralityColl(collision, centEstimator);
      if (centrality < centralityMin || centrality > centralityMax) {
        continue; // skip this collision if outside of the centrality range
      }
      registry.fill(HIST("hCentrality"), centrality);

      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      std::array<float, 3> const qVecs = getQvec(collision, qVecDetector.value);

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false, false, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks, &qVecs)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarInPbPb, "Process Dstar candidates in PbPb collisions", false);

  void processDstarWithMlInPbPb(CollsWithQVecs const& collisions,
                                FilteredCandDstarWSelFlagAndMl const& dstarCandidates,
                                TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      const auto centrality = getCentralityColl(collision, centEstimator);
      if (centrality < centralityMin || centrality > centralityMax) {
        continue; // skip this collision if outside of the centrality range
      }
      registry.fill(HIST("hCentrality"), centrality);

      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      std::array<float, 3> const qVecs = getQvec(collision, qVecDetector.value);

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false, false, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks, &qVecs)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarWithMlInPbPb, "Process Dstar candidates with ML in PbPb collisions", false);

  void processDstarMcInPbPb(aod::McCollisions::iterator const&,
                            McParticlesDstarMatched const& mcParticles,
                            CollisionsWithMcLabelsAndCent const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                            FilteredCandDstarWSelFlagAndMc const& dstarCandidates,
                            TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};

    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      const auto centrality = getCentralityColl(collision, centEstimator);
      if (centrality < centralityMin || centrality > centralityMax) {
        continue; // skip this collision if outside of the centrality range
      }
      registry.fill(HIST("hCentrality"), centrality);

      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
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
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcParticle.mcCollision().globalIndex());
      const auto cent = getCentralityGenColl(recoCollsPerMcColl, centEstimator);
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true>(mcParticle, mcParticles, numPvContributorsGen, &cent);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarMcInPbPb, "Process Dstar candidates in PbPb MC without ML", false);

  void processDstarMcWithMlInPbPb(aod::McCollisions::iterator const&,
                                  McParticlesDstarMatched const& mcParticles,
                                  CollisionsWithMcLabelsAndCent const& collisions, // this is grouped with SmallGroupsCollisionsWithMcLabels const& collisions,
                                  FilteredCandDstarWSelFlagAndMcAndMl const& dstarCandidates,
                                  TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};

    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      const auto centrality = getCentralityColl(collision, centEstimator);
      if (centrality < centralityMin || centrality > centralityMax) {
        continue; // skip this collision if outside of the centrality range
      }
      registry.fill(HIST("hCentrality"), centrality);

      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedDstarCandidates = dstarCandidates.sliceBy(dstarWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& dstarCandidate : groupedDstarCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, true>(dstarCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
    for (const auto& mcParticle : mcParticles) {
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcParticle.mcCollision().globalIndex());
      const auto cent = getCentralityGenColl(recoCollsPerMcColl, centEstimator);
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true>(mcParticle, mcParticles, numPvContributorsGen, &cent);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processDstarMcWithMlInPbPb, "Process Dstar candidates in PbPb MC with ML", false);

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
      int const numPvContributors = collision.numContrib();
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
  PROCESS_SWITCH(HfTaskCharmPolarisation, processLcToPKPi, "Process Lc candidates without ML", false);

  // Lc->pKpi with ML cuts
  void processLcToPKPiWithMl(aod::Collisions const& collisions,
                             FilteredCandLcToPKPiWSelFlagAndMl const& lcCandidates,
                             TracksWithExtra const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiWithMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }

        /// rotational background
        for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, false>(lcCandidate, iRotation, numPvContributors, -1 /*MC particles*/, tracks);
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processLcToPKPiWithMl, "Process Lc candidates with ML", false);

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
      int const numPvContributors = collision.numContrib();
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
  PROCESS_SWITCH(HfTaskCharmPolarisation, processLcToPKPiMc, "Process Lc candidates in MC without ML", false);

  // Lc->pKpi in MC with ML cuts
  void processLcToPKPiMcWithMl(aod::McCollisions::iterator const&,
                               McParticles3ProngMatched const& mcParticles,
                               CollisionsWithMcLabels const& collisions, // this is grouped with SmallGroups
                               FilteredCandLcToPKPiWSelFlagAndMcAndMl const& lcCandidates,
                               TracksWithExtra const& tracks)
  {
    int numPvContributorsGen{0};
    for (const auto& collision : collisions) { // loop over reco collisions associated to this gen collision
      auto thisCollId = collision.globalIndex();
      int const numPvContributors = collision.numContrib();
      auto groupedLcCandidates = lcCandidates.sliceBy(lcToPKPiWithMcAndMlPerCollision, thisCollId);
      int nCands{0}, nCandsInSignalRegion{0};

      if (numPvContributors > numPvContributorsGen) { // we take the associated reconstructed collision with higher number of PV contributors
        numPvContributorsGen = numPvContributors;
      }

      for (const auto& lcCandidate : groupedLcCandidates) {
        nCands++;
        if (runPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi, true, true>(lcCandidate, 0, numPvContributors, -1 /*MC particles*/, tracks)) {
          nCandsInSignalRegion++;
        }
      }
      fillMultHistos(numPvContributors, nCands, nCandsInSignalRegion);
    }

    for (const auto& mcParticle : mcParticles) {
      runMcGenPolarisationAnalysis<charm_polarisation::DecayChannel::LcToPKPi>(mcParticle, mcParticles, numPvContributorsGen);
    }
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processLcToPKPiMcWithMl, "Process Lc candidates in MC with ML", false);

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
  PROCESS_SWITCH(HfTaskCharmPolarisation, processLcToPKPiBackgroundMcWithMl, "Process Lc candidates in MC with ML w/o mcCollision grouping", false);

  // Event-plane resolution
  void processResolEventPlane(CollsWithQVecs::iterator const& collision,
                              aod::BCsWithTimestamps const&)
  {
    float centrality{-1.f};
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentralityEstimator::FT0C, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    if (rejectionMask != 0) {
      return;
    }
    centrality = getCentralityColl(collision, CentralityEstimator::FT0C);

    float const psiFT0a = epHelper.GetEventPlane(collision.qvecFT0ARe(), collision.qvecFT0AIm(), 2);
    float const psiFT0c = epHelper.GetEventPlane(collision.qvecFT0CRe(), collision.qvecFT0CIm(), 2);
    float const psiFT0m = epHelper.GetEventPlane(collision.qvecFT0MRe(), collision.qvecFT0MIm(), 2);
    float const psiFV0a = epHelper.GetEventPlane(collision.qvecFV0ARe(), collision.qvecFV0AIm(), 2);
    float const psiBPoss = epHelper.GetEventPlane(collision.qvecBPosRe(), collision.qvecBPosIm(), 2);
    float const psiBNegs = epHelper.GetEventPlane(collision.qvecBNegRe(), collision.qvecBNegIm(), 2);
    float const psiBTots = epHelper.GetEventPlane(collision.qvecBTotRe(), collision.qvecBTotIm(), 2);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleFV0A"), centrality, psiFV0a);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleFT0A"), centrality, psiFT0a);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleFT0C"), centrality, psiFT0c);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleFT0M"), centrality, psiFT0m);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleTPCpos"), centrality, psiBPoss);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleTPCneg"), centrality, psiBNegs);
    registry.fill(HIST("resolEvPlane/hEvPlaneAngleTPCtot"), centrality, psiBTots);

    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0CFT0A"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0c, psiFT0a, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0CFV0A"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0c, psiFV0a, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0CTPCpos"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0c, psiBPoss, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0CTPCneg"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0c, psiBNegs, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0CTPCtot"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0c, psiBTots, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0AFV0A"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0a, psiFV0a, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0ATPCpos"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0a, psiBPoss, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0ATPCneg"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0a, psiBNegs, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0ATPCtot"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0a, psiBTots, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0MFV0A"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0m, psiFV0a, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0MTPCpos"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0m, psiBPoss, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0MTPCneg"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0m, psiBNegs, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFT0MTPCtot"), centrality, std::cos(2 * getDeltaPsiInRange(psiFT0m, psiBTots, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFV0ATPCpos"), centrality, std::cos(2 * getDeltaPsiInRange(psiFV0a, psiBPoss, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFV0ATPCneg"), centrality, std::cos(2 * getDeltaPsiInRange(psiFV0a, psiBNegs, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneFV0ATPCtot"), centrality, std::cos(2 * getDeltaPsiInRange(psiFV0a, psiBTots, 2)));
    registry.fill(HIST("resolEvPlane/hResolEvPlaneTPCposTPCneg"), centrality, std::cos(2 * getDeltaPsiInRange(psiBPoss, psiBNegs, 2)));
  }
  PROCESS_SWITCH(HfTaskCharmPolarisation, processResolEventPlane, "Process event-plane resolution", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmPolarisation>(cfgc)};
}
