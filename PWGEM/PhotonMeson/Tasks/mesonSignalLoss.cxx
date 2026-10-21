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
// ========================
//
// Calculate signal and event loss for neutral mesons and photons
// Joshua König, joshua.konig@cern.ch

#include "PWGEM/PhotonMeson/DataModel/EventTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
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
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct MesonSignalLoss {

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  using RecoCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using EventAccBits = o2::aod::pmevsel::EventAcceptanceBits;

  enum MultiplicityEstimator { kCentralBarrel,
                               kFT0A,
                               kFT0C };

  // stage kStageHasReco + k means "passed event cut k of EventAccBits"
  enum Stage { kStageAllMC = 0,
               kStageTVX,
               kStageHasReco };

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  Preslice<aod::McCollisions> mcCollisionsPerBC = aod::mccollision::bcId;
  Preslice<aod::McParticles> particlesPerMCCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<RecoCollisions> recoCollisionsPerMCCollision = aod::mccollisionlabel::mcCollisionId;

  struct : ConfigurableGroup {
    std::string prefix = "axisConfig";
    ConfigurableAxis cfgPtAxis{"cfgPtAxis", {400, 0.f, 20.f}, "pT axis for generated particles"};
    ConfigurableAxis cfgMultiplicityAxis{"cfgMultiplicityAxis", {200, -0.5f, 199.5f}, "axis for true charged-particle multiplicity"};
  } axisConfig;

  struct : ConfigurableGroup {
    std::string prefix = "particleConfig";
    Configurable<float> cfgRapidity{"cfgRapidity", 0.8, "rapidity selection"};
    Configurable<float> cfgSecondaryDecayLength{"cfgSecondaryDecayLength", 1.f, "Minimum ancestor decay length (cm) to classify a particle as secondary"};
    Configurable<bool> cfgEnableMultiplicityHistograms{"cfgEnableMultiplicityHistograms", false, "enable generated-particle and event-stage histograms versus true charged-particle multiplicity"};
    Configurable<std::string> cfgMultiplicityEstimator{"cfgMultiplicityEstimator", "centralBarrel", "eta range for true multiplicity: centralBarrel, FT0A, or FT0C"};
  } particleConfig;

  struct : ConfigurableGroup {
    std::string prefix = "rctConfig";
    Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", false, "require good detector flag in run condtion table"};
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for PbPb"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } rctConfig;

  struct : ConfigurableGroup {
    std::string prefix = "eventCuts";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -1e+10, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 1e+10, "max. Zvtx"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", false, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireTVXinEMC{"cfgRequireTVXinEMC", false, "require kTVXinEMC (only for EMC analyses)"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", false, "require sel8 condition"};
  } eventCuts;

  Service<o2::framework::O2DatabasePDG> pdgDB{};
  o2::aod::rctsel::RCTFlagsChecker rctChecker;
  MultiplicityEstimator multiplicityEstimator{kCentralBarrel};

  /// \brief Creates the histograms used for signal-loss and event-selection analysis.
  void Inithistograms()
  {
    registry.add("hBCCounter", "hBCCounter;;#it{N}_{BC}", kTH1F, {{3, -0.5, 3 - 0.5}});
    auto h = registry.get<TH1>(HIST("hBCCounter"));
    h->GetXaxis()->SetBinLabel(1, "All BCs");
    h->GetXaxis()->SetBinLabel(2, "kTVX fired");
    h->GetXaxis()->SetBinLabel(3, "kTVX fired && has MC collision");

    registry.add("hMCCollStage", "hMCCollStage;;#it{N}_{MC coll}", kTH1F, {{20, -0.5, 20 - 0.5}});
    auto hStage = registry.get<TH1>(HIST("hMCCollStage"));
    hStage->GetXaxis()->SetBinLabel(kStageAllMC + 1, "all MC coll");
    hStage->GetXaxis()->SetBinLabel(kStageTVX + 1, "kTVX fired");
    hStage->GetXaxis()->SetBinLabel(kStageHasReco + 1, "has reco coll");
    for (int cut = 1; cut < EventAccBits::kSize; ++cut) {
      hStage->GetXaxis()->SetBinLabel(kStageHasReco + cut + 1, o2::aod::pmevsel::vecEventAccBitNames[cut].c_str());
    }

    registry.add("QC/hNRecoPerMCColl", "reco collisions per MC collision;N_{reco};N_{MC coll}", kTH1F, {{10, -0.5, 9.5}});
    registry.add("QC/hNSelectedRecoPerMCColl", "selected reco collisions per MC collision;N_{sel};N_{MC coll}", kTH1F, {{10, -0.5, 9.5}});
    registry.add("QC/hNContribVsOrder", "N contributors;order in matched list;N_{contrib}", kTH2F, {{10, -0.5, 9.5}, {200, 0, 200}});

    const AxisSpec ptAxis{axisConfig.cfgPtAxis, "p_{T} (GeV/c)"};
    const AxisSpec multAxis{axisConfig.cfgMultiplicityAxis, "N_{ch}^{true}"};
    registry.add("Generated/Pi0/hGeneratedParticlePt", "Generated Pi0;p_{T} (GeV/c);selection stage", kTH2F, {ptAxis, {20, -0.5, 20 - 0.5}});
    registry.add("Generated/Eta/hGeneratedParticlePt", "Generated Eta;p_{T} (GeV/c);selection stage", kTH2F, {ptAxis, {20, -0.5, 20 - 0.5}});
    registry.add("Generated/Gamma/hGeneratedParticlePt", "Generated Gamma;p_{T} (GeV/c);selection stage", kTH2F, {ptAxis, {20, -0.5, 20 - 0.5}});
    if (particleConfig.cfgEnableMultiplicityHistograms) {
      registry.add("Generated/Multiplicity/hMCCollStage", "MC collision stages vs. true multiplicity;N_{ch}^{true};selection stage", kTH2F, {multAxis, {20, -0.5, 20 - 0.5}});
      registry.add("Generated/Multiplicity/Pi0/hGeneratedParticlePt", "Generated Pi0;p_{T} (GeV/c);N_{ch}^{true};selection stage", kTH3F, {ptAxis, multAxis, {20, -0.5, 20 - 0.5}});
      registry.add("Generated/Multiplicity/Eta/hGeneratedParticlePt", "Generated Eta;p_{T} (GeV/c);N_{ch}^{true};selection stage", kTH3F, {ptAxis, multAxis, {20, -0.5, 20 - 0.5}});
      registry.add("Generated/Multiplicity/Gamma/hGeneratedParticlePt", "Generated Gamma;p_{T} (GeV/c);N_{ch}^{true};selection stage", kTH3F, {ptAxis, multAxis, {20, -0.5, 20 - 0.5}});
    }
  }

  /// \brief Configures the multiplicity estimator and initializes histograms.
  /// \param[in] context Workflow initialization context.
  void init([[maybe_unused]] InitContext& context)
  {
    // Init multiplicity estimator (only relevant for multiplicity dependent analyses)
    if (particleConfig.cfgMultiplicityEstimator.value == "centralBarrel") {
      multiplicityEstimator = kCentralBarrel;
    } else if (particleConfig.cfgMultiplicityEstimator.value == "FT0A") {
      multiplicityEstimator = kFT0A;
    } else if (particleConfig.cfgMultiplicityEstimator.value == "FT0C") {
      multiplicityEstimator = kFT0C;
    } else if (particleConfig.cfgEnableMultiplicityHistograms) {
      LOGF(fatal, "Unknown cfgMultiplicityEstimator '%s'; use centralBarrel, FT0A, or FT0C", particleConfig.cfgMultiplicityEstimator.value.c_str());
    }
    Inithistograms();
  }

  /// \brief Checks whether eta lies in the configured multiplicity-estimator range.
  /// \param[in] eta Pseudorapidity to check.
  /// \return True if eta is inside the configured range, otherwise false.
  bool isInMultiplicityEtaRange(float eta) const
  {
    bool decision = false;
    if (multiplicityEstimator == kCentralBarrel) {
      decision = std::abs(eta) < 0.8f; // o2-linter: disable=magic-number (rapidity range of central barrel)
    } else if (multiplicityEstimator == kFT0A) {
      decision = 3.3f < eta && eta < 5.0f; // o2-linter: disable=magic-number (rapidity range of FT0)
    } else if (multiplicityEstimator == kFT0C) {
      decision = -3.3f < eta && eta < -2.1f; // o2-linter: disable=magic-number (rapidity range of FT0)
    }
    return decision;
  }

  /// \brief Counts charged physical-primary particles in the configured eta range.
  /// \param[in] particles Particles to include in the multiplicity count.
  /// \return Number of charged physical-primary particles in range.
  template <typename TParticles>
  int countTrueMultiplicity(TParticles const& particles) const
  {
    int multiplicity = 0;
    for (const auto& particle : particles) {
      if (!particle.isPhysicalPrimary() || !isInMultiplicityEtaRange(particle.eta())) {
        continue;
      }
      const auto* pdgParticle = pdgDB->GetParticle(particle.pdgCode());
      if (pdgParticle && std::abs(pdgParticle->Charge()) >= 1.) { // o2-linter: disable=magic-number (check that particle has electric charge)
        ++multiplicity;
      }
    }
    return multiplicity;
  }

  /// \brief Determines whether a particle has a secondary-production origin.
  /// \param[in] particle Particle whose origin is checked.
  /// \param[in] mcParticles MC-particle table used to trace the particle ancestry.
  /// \return True if the particle is classified as secondary, otherwise false.
  bool isSecondary(aod::McParticle const& particle, aod::McParticles const& mcParticles) const
  {
    std::vector<int64_t> particlesToCheck{particle.globalIndex()};
    std::unordered_set<int64_t> checkedParticles;

    while (!particlesToCheck.empty()) {
      const auto particleIndex = particlesToCheck.back();
      particlesToCheck.pop_back();
      if (particleIndex < 0 || particleIndex >= mcParticles.size() || !checkedParticles.insert(particleIndex).second) {
        continue;
      }

      const auto currentParticle = mcParticles.iteratorAt(particleIndex);
      const auto process = currentParticle.getProcess();
      const bool fromDecay = currentParticle.producedByGenerator() || process == TMCProcess::kPDecay || process == TMCProcess::kPRadDecay;
      if (!currentParticle.producedByGenerator() && !fromDecay) {
        return true;
      }

      for (const auto motherIndex : currentParticle.mothersIds()) {
        if (motherIndex < 0 || motherIndex >= mcParticles.size()) {
          continue;
        }
        const auto mother = mcParticles.iteratorAt(motherIndex);
        const float decayLength = std::sqrt(std::pow(currentParticle.vx() - mother.vx(), 2) + std::pow(currentParticle.vy() - mother.vy(), 2) + std::pow(currentParticle.vz() - mother.vz(), 2));
        if (fromDecay && decayLength > particleConfig.cfgSecondaryDecayLength.value) {
          return true;
        }
        particlesToCheck.push_back(motherIndex);
      }
    }
    return false;
  }

  /// \brief Fills generated-particle histograms for the requested selection stages.
  /// \param[in] particlesInMCCollision Particles associated with the MC collision.
  /// \param[in] mcParticles MC-particle table used to inspect particle ancestry.
  /// \param[in] selectionStage First selection stage to fill.
  /// \param[in] selectionStageMax Last selection stage to fill.
  /// \param[in] trueMultiplicity True charged-particle multiplicity of the collision.
  template <typename TParticles>
  void fillGeneratedParticles(TParticles const& particlesInMCCollision, aod::McParticles const& mcParticles, int selectionStage, int selectionStageMax, int trueMultiplicity)
  {
    if (selectionStageMax < selectionStage) {
      selectionStageMax = selectionStage;
    }
    for (const auto& particle : particlesInMCCollision) {
      if (std::abs(particle.y()) > particleConfig.cfgRapidity) {
        continue;
      }
      if (isSecondary(particle, mcParticles)) {
        continue;
      }
      for (int iSel = selectionStage; iSel <= selectionStageMax; ++iSel) {
        switch (particle.pdgCode()) {
          case PDG_t::kPi0:
            registry.fill(HIST("Generated/Pi0/hGeneratedParticlePt"), particle.pt(), iSel); // also add rapdidity??
            if (particleConfig.cfgEnableMultiplicityHistograms) {
              registry.fill(HIST("Generated/Multiplicity/Pi0/hGeneratedParticlePt"), particle.pt(), trueMultiplicity, iSel);
            }
            break;
          case o2::constants::physics::Pdg::kEta:
            registry.fill(HIST("Generated/Eta/hGeneratedParticlePt"), particle.pt(), iSel); // also add rapdidity??
            if (particleConfig.cfgEnableMultiplicityHistograms) {
              registry.fill(HIST("Generated/Multiplicity/Eta/hGeneratedParticlePt"), particle.pt(), trueMultiplicity, iSel);
            }
            break;
          case PDG_t::kGamma:
            registry.fill(HIST("Generated/Gamma/hGeneratedParticlePt"), particle.pt(), iSel); // also add rapdidity??
            if (particleConfig.cfgEnableMultiplicityHistograms) {
              registry.fill(HIST("Generated/Multiplicity/Gamma/hGeneratedParticlePt"), particle.pt(), trueMultiplicity, iSel);
            }
            break;
          default:
            break;
        }
      }
    }
  }

  /// \brief Returns the event-acceptance bit for the first failed collision selection.
  /// \param[in] collision Reconstructed collision to evaluate.
  /// \return Event-acceptance bit identifying the first failed cut, or kSize if accepted.
  template <typename TCollision>
  int isSelectedEvent(TCollision const& collision)
  {
    if (!collision.has_mcCollision()) {
      return EventAccBits::kHasMCColl;
    }

    if (collision.posZ() < eventCuts.cfgZvtxMin || eventCuts.cfgZvtxMax < collision.posZ()) {
      return EventAccBits::kGoodZVtx;
    }

    if (eventCuts.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return EventAccBits::kIsFT0AND;
    }

    if (eventCuts.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return EventAccBits::kNoTFB;
    }

    if (eventCuts.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return EventAccBits::kITSROFB;
    }

    if (eventCuts.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return EventAccBits::kNoSameBunchPileUp;
    }

    if (eventCuts.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return EventAccBits::kGoodZVtxFTOPV;
    }

    if (eventCuts.cfgRequireNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return EventAccBits::kNoCollInTimeRange;
    }

    if (!(eventCuts.cfgTrackOccupancyMin <= collision.trackOccupancyInTimeRange() && collision.trackOccupancyInTimeRange() < eventCuts.cfgTrackOccupancyMax)) {
      return EventAccBits::kGoodTrackOccupancy;
    }

    if (!(eventCuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventCuts.cfgFT0COccupancyMax)) {
      return EventAccBits::kGoodFT0Occupancy;
    }

    if (eventCuts.cfgRequireTVXinEMC && !collision.alias_bit(triggerAliases::kTVXinEMC)) {
      return EventAccBits::kTVXInEMC;
    }

    // if constexpr (std::is_same_v<std::decay_t<TCollision>, MyCollisions_Cent::iterator>) {
    //   std::array<float, 3> centralities = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
    //   if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
    //     return false;
    //   }
    // }

    if (rctConfig.cfgRequireGoodRCT && !rctChecker.checkTable(collision)) {
      return EventAccBits::kGoodRCT;
    }

    if (eventCuts.cfgRequireSel8 && !collision.sel8()) {
      return EventAccBits::kGoodSel8;
    }

    return EventAccBits::kSize;
  }

  /// \brief Processes MC bunch crossings and fills generated signal-loss histograms.
  /// \param[in] bcs Bunch crossings to process.
  /// \param[in] collisions Reconstructed collisions and their MC labels.
  /// \param[in] mcCollisions MC collisions to process.
  /// \param[in] mcParticles MC particles used for multiplicity and ancestry checks.
  void processMBMC(MyBCs const& bcs, RecoCollisions const& collisions, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    // loop over all bcs
    for (const auto& bc : bcs) {
      registry.fill(HIST("hBCCounter"), 0);

      // Check if TVX trigger fired in that bc
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        registry.fill(HIST("hBCCounter"), 1);
      }

      auto mcCollisionsInBC = mcCollisions.sliceBy(mcCollisionsPerBC, bc.globalIndex());
      // check if has at least one mc collision
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX) && mcCollisionsInBC.size() > 0) {
        registry.fill(HIST("hBCCounter"), 2);
      }

      // loop over all mc collisions in the current bc (should be only 1 mc collision in most cases)
      for (const auto& mcCollision : mcCollisionsInBC) {
        int lastStage = kStageAllMC;
        auto particlesInMCCollision = mcParticles.sliceBy(particlesPerMCCollision, mcCollision.globalIndex());
        const int trueMultiplicity = particleConfig.cfgEnableMultiplicityHistograms ? countTrueMultiplicity(particlesInMCCollision) : 0;

        if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
          lastStage = kStageTVX;

          auto matchedRecoCollisions = collisions.sliceBy(recoCollisionsPerMCCollision, mcCollision.globalIndex());
          registry.fill(HIST("QC/hNRecoPerMCColl"), matchedRecoCollisions.size());

          // if at least one reco collision is found for the current mc collision, loop over all reco collisions and check if the collision is selected
          if (matchedRecoCollisions.size() > 0) {

            int bestReject = EventAccBits::kHasMCColl;
            int nSelected = 0;
            int order = 0;
            for (const auto& recColl : matchedRecoCollisions) {
              // this function returns if the collision is accepted or not. If not, the reason is encoded in the return value
              const int reject = isSelectedEvent(recColl);
              // furthest-passing matched collision decides, so split vertices are not double counted
              bestReject = std::max(bestReject, reject);
              if (reject == EventAccBits::kSize) {
                ++nSelected;
              }
              // Fill how many contributers to the primary vertex are in the current collision (QA info)
              registry.fill(HIST("QC/hNContribVsOrder"), order++, recColl.numContrib());
            }
            // Fill number of collisions that are assigned to the current mc collision
            registry.fill(HIST("QC/hNSelectedRecoPerMCColl"), nSelected);
            // this variable is used to keep track at which stage the collision was rejected or if it was accepted. It is used later for filling the generated meson distributions
            lastStage = kStageHasReco + (bestReject - 1);
          }
        }
        // Fill the particles of interested as a function of the reject reason (lastStage)
        fillGeneratedParticles(particlesInMCCollision, mcParticles, kStageAllMC, lastStage, trueMultiplicity);
        for (int stage = kStageAllMC; stage <= lastStage; ++stage) {
          registry.fill(HIST("hMCCollStage"), stage);
          if (particleConfig.cfgEnableMultiplicityHistograms) {
            registry.fill(HIST("Generated/Multiplicity/hMCCollStage"), trueMultiplicity, stage);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(MesonSignalLoss, processMBMC, "Process MC bunch crossings for signal and event loss", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<MesonSignalLoss>(context)};
}
