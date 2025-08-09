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

/// \file taskDstarPolarisationInJet.cxx
/// \brief Analysis task for Dstar spin alignment in jet (rho00 vs. z||)
///
/// \author F. Grosa (CERN) fabrizio.grosa@cern.ch
/// \author Mingze Li (CCNU) Mingze.li@cern.ch

#include "PWGHF/Core/DecayChannels.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetHFUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <string>
#include <vector>

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
  JetAxis,
  Production,
  Beam,
  Random,
  NTypes
};
enum DecayChannel : uint8_t {
  DstarToDzeroPi = 0,
  NChannels
};
} // namespace charm_polarisation
} // namespace o2::aod

struct HfTaskDstarPolarisationInJet {

  float massPi{0.f};
  float massKaon{0.f};
  float massProton{0.f};
  float massDstar{0.f};

  float bkgRotationAngleStep{0.f};

  uint8_t nMassHypos{0u};

  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 Pi"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  /// activate rotational background
  Configurable<int> nBkgRotations{"nBkgRotations", 0, "Number of rotated copies (background) per each original candidate"};
  Configurable<float> minRotAngleMultByPi{"minRotAngleMultByPi", 5. / 6, "Minimum angle rotation for track rotation, to be multiplied by pi"};
  Configurable<float> maxRotAngleMultByPi{"maxRotAngleMultByPi", 7. / 6, "Maximum angle rotation for track rotation, to be multiplied by pi"};

  // activate study of systematic uncertainties of tracking
  Configurable<bool> activateTrackingSys{"activateTrackingSys", false, "Activate the study of systematic uncertainties of tracking"};

  /// output THnSparses
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarJetAxis{"activateTHnSparseCosThStarJetAxis", true, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", true, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activatePartRecoDstar{"activatePartRecoDstar", false, "Activate the study of partly reconstructed D*+ -> D0 (-> KPiPi0) Pi decays"};
  float invMassMin{0.f};
  float invMassMax{1000.f};

  /// Application of rapidity cut for reconstructed candidates
  Configurable<float> maxAbsRapidityCut{"maxAbsRapidityCut", 999.f, "Max. value of reconstructed candidate rapidity (abs. value)"};

  // Tables for MC jet matching
  using DstarJets = soa::Join<aod::DstarChargedJets, aod::DstarChargedJetConstituents>;
  using JetMCDTable = soa::Join<aod::DstarChargedMCDetectorLevelJets, aod::DstarChargedMCDetectorLevelJetConstituents, aod::DstarChargedMCDetectorLevelJetsMatchedToDstarChargedMCParticleLevelJets>;
  using TracksWithExtra = soa::Join<aod::JTracks, aod::TracksExtra>;

  // slices for accessing proper HF mcdjets collision associated to mccollisions
  Preslice<JetMCDTable> dstarMCDJetsPerCollisionPreslice = aod::jet::collisionId;
  Preslice<DstarJets> dstarJetsPerCollision = aod::jet::collisionId;
  PresliceUnsorted<aod::JetCollisionsMCD> collisionsPerMCCollisionPreslice = aod::jmccollisionlb::mcCollisionId;

  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {200, 0.139f, 0.179f}, "#it{M} (GeV/#it{c}^{2})"}; // o2-linter: disable=pdg/explicit-mass (false positive)
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.f, 100.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisY{"configThnAxisY", {20, -1.f, 1.f}, "#it{y}"};
  ConfigurableAxis configThnAxisCosThetaStar{"configThnAxisCosThetaStar", {20, -1.f, 1.f}, "cos(#vartheta_{helicity})"};
  ConfigurableAxis configThnAxisMlBkg{"configThnAxisMlBkg", {100, 0.f, 1.f}, "ML bkg"};
  ConfigurableAxis configThnAxisInvMassD0{"configThnAxisInvMassD0", {250, 1.65f, 2.15f}, "#it{M}(D^{0}) (GeV/#it{c}^{2})"}; // only for D*+
  // ConfigurableAxis configThnAxisMlPrompt{"configThnAxisMlPrompt", {100, 0.f, 1.f}, "ML prompt"};
  ConfigurableAxis configThnAxisMlNonPrompt{"configThnAxisMlNonPrompt", {100, 0.f, 1.f}, "ML non-prompt"};
  ConfigurableAxis configThnAxisNumCandidates{"configThnAxisNumCandidates", {1, -0.5f, 0.5f}, "num candidates"};
  ConfigurableAxis configThnAxisPtB{"configThnAxisPtB", {3000, 0.f, 300.f}, "#it{p}_{T}(B mother) (GeV/#it{c})"};
  ConfigurableAxis configThnAxisAbsEtaTrackMin{"configThnAxisAbsEtaTrackMin", {3, 0.f, 0.3f}, "min |#it{#eta_{track}}|"};
  ConfigurableAxis configThnAxisNumItsClsMin{"configThnAxisNumItsClsMin", {4, 3.5f, 7.5f}, "min #it{N}_{cls ITS}"};
  ConfigurableAxis configThnAxisNumTpcClsMin{"configThnAxisNumTpcClsMin", {3, 79.5f, 140.5f}, "min #it{N}_{cls TPC}"};
  ConfigurableAxis configThnAxisCharge{"configThnAxisCharge", {2, -2.f, 2.f}, "electric charge"};
  ConfigurableAxis configThnAxisProjection{"configThnAxisProjection", {300, 0.f, 10.f}, "z^{D^{*+},jet}_{||}"};
  ConfigurableAxis configThnAxisJetPt{"configThnAxisJetPt", {100, 0.f, 100.f}, "#it{p}_{T} (GeV/#it{c})"};

  HistogramRegistry registry{"registry", {}};

  std::vector<int> eventSelectionBits;

  void init(InitContext&)
  {
    // initialise event selection:
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(eventSelections.value);

    /// check process functions
    const int nProcesses =
      static_cast<int>(doprocessDstar) +
      static_cast<int>(doprocessDstarWithMl) +
      static_cast<int>(doprocessDstarMc) +
      static_cast<int>(doprocessDstarMcWithMl);
    // std::array<int, 4> processes = {doprocessDstar, doprocessDstarWithMl, doprocessDstarMc, doprocessDstarMcWithMl};
    // const int nProcesses = std::accumulate(processes.begin(), processes.end(), 0);
    if (nProcesses > 1) {
      LOGP(fatal, "Only one process function should be enabled at a time, please check your configuration");
    }
    if (nProcesses == 0) {
      LOGP(fatal, "No process function enabled");
    }

    /// check output THnSparses
    std::array<int, 3> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarJetAxis, activateTHnSparseCosThStarProduction};
    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    }
    if (activateTHnSparseCosThStarHelicity) {
      LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
    }
    if (activateTHnSparseCosThStarJetAxis) {
      LOGP(info, "THnSparse with cosThStar w.r.t. jet axis active.");
    }
    if (activateTHnSparseCosThStarProduction) {
      LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
    }

    if (activatePartRecoDstar && !(doprocessDstarMc || doprocessDstarMcWithMl)) {
      LOGP(fatal, "Check on partly reconstructed D* mesons only possible for processDstarMc and processDstarMcWithMl");
    }

    // check bkg rotation for MC (not supported currently)
    if (nBkgRotations > 0 && (doprocessDstarMc || doprocessDstarMcWithMl)) {
      LOGP(fatal, "No background rotation supported for MC.");
    }

    massPi = o2::constants::physics::MassPiPlus;
    massProton = o2::constants::physics::MassProton;
    massKaon = o2::constants::physics::MassKaonCharged;
    massDstar = o2::constants::physics::MassDStar;
    bkgRotationAngleStep = (nBkgRotations > 1) ? (maxRotAngleMultByPi - minRotAngleMultByPi) * constants::math::PI / (nBkgRotations - 1) : 0.;

    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassD0{configThnAxisInvMassD0, "#it{M}(D^{0}) (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisY{configThnAxisY, "#it{y}"};
    const AxisSpec thnAxisCosThetaStar{configThnAxisCosThetaStar, "cos(#vartheta)"};
    const AxisSpec thnAxisMlBkg{configThnAxisMlBkg, "ML bkg"};
    const AxisSpec thnAxisMlNonPrompt{configThnAxisMlNonPrompt, "ML non-prompt"};
    const AxisSpec thnAxisIsRotatedCandidate{2, -0.5f, 1.5f, "rotated bkg"};
    const AxisSpec thnAxisNumcandidates{configThnAxisNumCandidates, "num candidates"};
    const AxisSpec thnAxisPtB{configThnAxisPtB, "#it{p}_{T}(B mother) (GeV/#it{c})"};
    const AxisSpec thnAxisDausAcc{2, -0.5f, 1.5f, "daughters in acceptance"};
    const AxisSpec thnAxisDauToMuons{4, -0.5f, 3.5f, "daughters decayed to muons"};
    const AxisSpec thnAxisAbsEtaTrackMin{configThnAxisAbsEtaTrackMin, "min |#it{#eta_{track}}|"};
    const AxisSpec thnAxisNumItsClsMin{configThnAxisNumItsClsMin, "min #it{N}_{cls ITS}"};
    const AxisSpec thnAxisNumTpcClsMin{configThnAxisNumTpcClsMin, "min #it{N}_{cls TPC}"};
    const AxisSpec thnAxisCharge{configThnAxisCharge, "charge"};
    const AxisSpec thnAxisProjection{configThnAxisProjection, "z_{||}"};
    const AxisSpec thnAxisJetPt{configThnAxisJetPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisSelFlag{2, -0.5f, 1.5f, "Sel flag"};

    auto invMassBins = thnAxisInvMass.binEdges;
    invMassMin = invMassBins.front();
    invMassMax = invMassBins.back();

    std::vector<AxisSpec> thnRecDataAxes = {thnAxisInvMass, thnAxisPt, thnAxisY, thnAxisInvMassD0, thnAxisCosThetaStar};
    if (activateTrackingSys) {
      thnRecDataAxes.insert(thnRecDataAxes.end(), {thnAxisAbsEtaTrackMin, thnAxisNumItsClsMin, thnAxisNumTpcClsMin});
    }
    if (nBkgRotations > 0) {
      thnRecDataAxes.push_back(thnAxisIsRotatedCandidate);
    }
    if (doprocessDstarWithMl || doprocessDstarMcWithMl) {
      thnRecDataAxes.insert(thnRecDataAxes.end(), {thnAxisMlBkg, thnAxisMlNonPrompt});
    }
    std::vector<AxisSpec> thnRecPromptAxes = thnRecDataAxes;
    std::vector<AxisSpec> thnRecNonPromptAxes = thnRecDataAxes;

    if (doprocessDstar || doprocessDstarWithMl) {
      thnRecDataAxes.insert(thnRecDataAxes.end(), {thnAxisProjection, thnAxisJetPt});
      // Data histos
      if (activateTHnSparseCosThStarHelicity) {
        registry.add("hHelicity", "THn for polarisation studies with cosThStar ", HistType::kTHnSparseF, thnRecDataAxes);
      }
      if (activateTHnSparseCosThStarProduction) {
        registry.add("hProduction", "THn for polarisation studies with cosThStar w.r.t. production axis", HistType::kTHnSparseF, thnRecDataAxes);
      }
      if (activateTHnSparseCosThStarJetAxis) {
        registry.add("hJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis", HistType::kTHnSparseF, thnRecDataAxes);
      }
    } else if (doprocessDstarMc || doprocessDstarMcWithMl) {
      // MC Reco histos
      thnRecPromptAxes.insert(thnRecPromptAxes.end(), {thnAxisDauToMuons, thnAxisProjection, thnAxisJetPt});
      thnRecNonPromptAxes.insert(thnRecNonPromptAxes.end(), {thnAxisPtB, thnAxisDauToMuons, thnAxisProjection, thnAxisJetPt});
      if (activateTHnSparseCosThStarHelicity) {
        registry.add("hRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
        registry.add("hRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hPartRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- partially reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
          registry.add("hPartRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis and BDT scores -- partially reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        }
      }
      if (activateTHnSparseCosThStarProduction) {
        registry.add("hRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
        registry.add("hRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hPartRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- partially reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
          registry.add("hPartRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis and BDT scores -- partially reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        }
      }
      if (activateTHnSparseCosThStarJetAxis) {
        registry.add("hRecoPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis and BDT scores -- reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
        registry.add("hRecoNonPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis and BDT scores -- reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        if (activatePartRecoDstar) {
          registry.add("hPartRecoPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis and BDT scores -- partially reco prompt signal", HistType::kTHnSparseF, thnRecPromptAxes);
          registry.add("hPartRecoNonPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis and BDT scores -- partially reco non-prompt signal", HistType::kTHnSparseF, thnRecNonPromptAxes);
        }
      }

      // MC Gen histos
      std::vector<AxisSpec> thnGenPromptAxes = {thnAxisPt, thnAxisY, thnAxisCosThetaStar, thnAxisDausAcc, thnAxisCharge, thnAxisProjection, thnAxisJetPt};
      std::vector<AxisSpec> thnGenNonPromptAxes = {thnAxisPt, thnAxisY, thnAxisCosThetaStar, thnAxisPtB, thnAxisDausAcc, thnAxisCharge, thnAxisProjection, thnAxisJetPt};

      if (activateTHnSparseCosThStarHelicity) {
        registry.add("hGenPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen prompt signal", HistType::kTHnSparseF, thnGenPromptAxes);
        registry.add("hGenNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen non-prompt signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        if (activatePartRecoDstar && (doprocessDstarMc || doprocessDstarMcWithMl)) {
          registry.add("hGenPartRecoPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen prompt partly reco signal", HistType::kTHnSparseF, thnGenPromptAxes);
          registry.add("hGenPartRecoNonPromptHelicity", "THn for polarisation studies with cosThStar w.r.t. helicity axis -- gen non-prompt partly reco signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        }
      }
      if (activateTHnSparseCosThStarProduction) {
        registry.add("hGenPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen prompt signal", HistType::kTHnSparseF, thnGenPromptAxes);
        registry.add("hGenNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen non-prompt signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        if (activatePartRecoDstar && (doprocessDstarMc || doprocessDstarMcWithMl)) {
          registry.add("hGenPartRecoPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen prompt partly reco signal", HistType::kTHnSparseF, thnGenPromptAxes);
          registry.add("hGenPartRecoNonPromptProduction", "THn for polarisation studies with cosThStar w.r.t. production axis -- gen non-prompt partly reco signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        }
      }
      if (activateTHnSparseCosThStarJetAxis) {
        registry.add("hGenPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis -- gen prompt signal", HistType::kTHnSparseF, thnGenPromptAxes);
        registry.add("hGenNonPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis -- gen non-prompt signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        if (activatePartRecoDstar && (doprocessDstarMc || doprocessDstarMcWithMl)) {
          registry.add("hGenPartRecoPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis -- gen prompt partly reco signal", HistType::kTHnSparseF, thnGenPromptAxes);
          registry.add("hGenPartRecoNonPromptJetAxis", "THn for polarisation studies with cosThStar w.r.t. jet axis -- gen non-prompt partly reco signal", HistType::kTHnSparseF, thnGenNonPromptAxes);
        }
      }
    }
    // inv. mass hypothesis to loop over
    nMassHypos = 1;
  }; // end init

  /// \param invMassCharmHad is the invariant-mass of the candidate
  /// \param ptCharmHad is the pt of the candidate
  /// \param rapCharmHad is the rapidity of the candidate
  /// \param invMassD0 is the invariant-mass of the D0 daugher (only for D*+)
  /// \param cosThetaStar is the cosThetaStar of the candidate
  /// \param outputMl is the array with ML output scores
  /// \param isRotatedCandidate is a flag that keeps the info of the rotation of the candidate for bkg studies
  /// \param origin is the MC origin
  /// \param ptBhadMother is the pt of the b-hadron mother (only in case of non-prompt)
  /// \param absEtaMin is the minimum absolute eta of the daughter tracks
  /// \param numItsClsMin is the minimum number of ITS clusters of the daughter tracks
  /// \param numTpcClsMin is the minimum number of TPC clusters of the daughter tracks
  /// \param nMuons is the number of muons from daughter decays
  /// \param isPartRecoDstar is a flag indicating if it is a partly reconstructed Dstar meson (MC only)
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, charm_polarisation::CosThetaStarType cosThetaStarType>
  void fillRecoHistos(float invMassCharmHad, float ptCharmHad, float rapCharmHad, float invMassD0, float cosThetaStar, std::array<float, 3> outputMl, int isRotatedCandidate, int8_t origin, float ptBhadMother, float absEtaMin, int numItsClsMin, int numTpcClsMin, int8_t nMuons, bool isPartRecoDstar, float zParallel, float jetPt)
  {

    if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if constexpr (!doMc) {                                                            // data
        if constexpr (withMl) {                                                         // with ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          }
        } else { // without ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, zParallel, jetPt);
            }
          }
        }
      } else {                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            }
          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        } else {                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              }
            }

          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptHelicity"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if constexpr (!doMc) {                                                                     // data
        if constexpr (withMl) {                                                                  // with ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          }
        } else { // without ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, zParallel, jetPt);
            }
          }
        }
      } else {                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            }
          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        } else {                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              }
            }

          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptProduction"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        }
      }
    } else if constexpr (cosThetaStarType == charm_polarisation::CosThetaStarType::JetAxis) { // JetAxis
      if constexpr (!doMc) {                                                                  // data
        if constexpr (withMl) {                                                               // with ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            } else {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], zParallel, jetPt);
            }
          }
        } else { // without ML
          if (activateTrackingSys) {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, zParallel, jetPt);
            }
          } else {
            if (nBkgRotations > 0) {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, isRotatedCandidate, zParallel, jetPt);
            } else {
              registry.fill(HIST("hJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, zParallel, jetPt);
            }
          }
        }
      } else {                                           // MC --> no distinction among channels, since rotational bkg not supported
        if constexpr (withMl) {                          // with ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], nMuons, zParallel, jetPt);
              }
            }
          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, outputMl[0], /*outputMl[1],*/ outputMl[2], ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        } else {                                         // without ML
          if (origin == RecoDecay::OriginType::Prompt) { // prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, nMuons, zParallel, jetPt);
              }
            }

          } else { // non-prompt
            if (activateTrackingSys) {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, absEtaMin, numItsClsMin, numTpcClsMin, ptBhadMother, nMuons, zParallel, jetPt);
              }
            } else {
              if (!isPartRecoDstar) {
                registry.fill(HIST("hRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              } else {
                registry.fill(HIST("hPartRecoNonPromptJetAxis"), invMassCharmHad, ptCharmHad, rapCharmHad, invMassD0, cosThetaStar, ptBhadMother, nMuons, zParallel, jetPt);
              }
            }
          }
        }
      }
    }
  }

  /// \param ptCharmHad is the pt of the particle
  /// \param rapCharmHad is the rapidity of the particle
  /// \param cosThetaStar is the cosThetaStar of the particle
  /// \param origin is the MC origin
  /// \param ptBhadMother is the pt of the b-hadron mother (only in case of non-prompt)
  /// \param areDausInAcc is a flag indicating whether the daughters are in acceptance or not
  /// \param isPartRecoDstar is a flag indicating if it is a partly reconstructed Dstar->D0pi->Kpipipi0 meson (MC only)
  void fillGenHistos(charm_polarisation::CosThetaStarType cosThetaStarType, float ptCharmHad, float rapCharmHad, float cosThetaStar, int8_t origin, float ptBhadMother, bool areDausInAcc, int8_t charge, bool isPartRecoDstar, float zParallel, float jetPt)
  {
    if (cosThetaStarType == charm_polarisation::CosThetaStarType::Helicity) { // Helicity
      if (origin == RecoDecay::OriginType::Prompt) {                          // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptHelicity"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoPromptHelicity"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptHelicity"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptHelicity"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        }
      }
    } else if (cosThetaStarType == charm_polarisation::CosThetaStarType::Production) { // Production
      if (origin == RecoDecay::OriginType::Prompt) {                                   // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptProduction"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoPromptProduction"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptProduction"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptProduction"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        }
      }
    } else if (cosThetaStarType == charm_polarisation::CosThetaStarType::JetAxis) { // JetAxis
      if (origin == RecoDecay::OriginType::Prompt) {                                // prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenPromptJetAxis"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoPromptJetAxis"), ptCharmHad, rapCharmHad, cosThetaStar, areDausInAcc, charge, zParallel, jetPt);
        }
      } else { // non-prompt
        if (!isPartRecoDstar) {
          registry.fill(HIST("hGenNonPromptJetAxis"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        } else {
          registry.fill(HIST("hGenPartRecoNonPromptJetAxis"), ptCharmHad, rapCharmHad, cosThetaStar, ptBhadMother, areDausInAcc, charge, zParallel, jetPt);
        }
      }
    }
  }

  /// \param invMass is the invariant mass
  /// \return true if candidate in signal region
  template <charm_polarisation::DecayChannel channel>
  bool isInSignalRegion(float invMass)
  {
    if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) { // D*+
      float invMassSigMin = 0.142f;
      float invMassSigMax = 0.15f;
      if (invMassSigMin < invMass && invMass < invMassSigMax) {
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

  /// \param candidates are the selected candidates
  /// \param bkgRotationId is the id for the background rotation
  /// \param particles are the generated particles
  /// \param tracks are the reconstructed tracks
  /// \return true if candidate in signal region
  template <charm_polarisation::DecayChannel channel, bool withMl, bool doMc, typename Jet, typename Cand>
  bool runPolarisationAnalysis(Jet const& jet, Cand const& candidate, int bkgRotationId)
  {
    bool isCandidateInSignalRegion{false};
    int8_t origin{RecoDecay::OriginType::None};
    float ptBhadMother{-1.f};
    bool partRecoDstar{false};
    if constexpr (doMc) {
      if constexpr (channel == charm_polarisation::DecayChannel::DstarToDzeroPi) {
        partRecoDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPiPi0 && std::abs(candidate.flagMcMatchRecCharm()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiKPi0;
        bool signalDstar = std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_dstar::DecayChannelMain::DstarToPiKPi && std::abs(candidate.flagMcMatchRecCharm()) == hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK;
        if (!signalDstar && (!partRecoDstar || !activatePartRecoDstar)) { // this candidate is not signal and not partially reconstructed signal, skip
          return isCandidateInSignalRegion;
        }
        origin = candidate.originMcRec();
        ptBhadMother = candidate.ptBhadMotherPart();
        int pdgBhadMother = candidate.pdgBhadMotherPart();
        // For unknown reasons there are charm hadrons coming directly from beauty diquarks without an intermediate B-hadron which have an unreasonable correlation between the pT of the charm hadron and the beauty mother. We also remove charm hadrons from quarkonia.
        if (origin == RecoDecay::OriginType::NonPrompt && (pdgBhadMother == 5101 || pdgBhadMother == 5103 || pdgBhadMother == 5201 || pdgBhadMother == 5203 || pdgBhadMother == 5301 || pdgBhadMother == 5303 || pdgBhadMother == 5401 || pdgBhadMother == 5403 || pdgBhadMother == 5503 || pdgBhadMother == 553 || pdgBhadMother == 555 || pdgBhadMother == 553 || pdgBhadMother == 557)) { // o2-linter: disable=pdg/explicit-code, magic-number (constants not in the PDG header)
          return isCandidateInSignalRegion;
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
        const float bkgRotAngle = (bkgRotationId > 0) ? minRotAngleMultByPi * constants::math::PI + bkgRotationAngleStep * (bkgRotationId - 1) : 0;

        std::array<float, 3> threeVecSoftPi{candidate.pxProng1() * std::cos(bkgRotAngle) - candidate.pyProng1() * std::sin(bkgRotAngle), candidate.pxProng1() * std::sin(bkgRotAngle) + candidate.pyProng1() * std::cos(bkgRotAngle), candidate.pzProng1()}; // we rotate the soft pion
        std::array<float, 3> threeVecD0Prong0{candidate.pxProng0Charm(), candidate.pyProng0Charm(), candidate.pzProng0Charm()};
        std::array<float, 3> threeVecD0Prong1{candidate.pxProng1Charm(), candidate.pyProng1Charm(), candidate.pzProng1Charm()};
        if (bkgRotationId > 0) {
          isRotatedCandidate = 1;
          pxDau = threeVecSoftPi[0];
          pyDau = threeVecSoftPi[1];
          pzDau = threeVecSoftPi[2];
          std::array<float, 3> threeVecCand = RecoDecay::pVec(threeVecSoftPi, threeVecD0Prong0, threeVecD0Prong1);
          pxCharmHad = threeVecCand[0];
          pyCharmHad = threeVecCand[1];
          pzCharmHad = threeVecCand[2];
          if (candidate.signProng1() > 0) {
            invMassCharmHad = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1, threeVecSoftPi}, std::array{massPi, massKaon, massPi});
            invMassD0 = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1}, std::array{massPi, massKaon});
          } else {
            invMassCharmHad = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1, threeVecSoftPi}, std::array{massKaon, massPi, massPi});
            invMassD0 = RecoDecay::m(std::array{threeVecD0Prong0, threeVecD0Prong1}, std::array{massKaon, massPi});
          }
          rapidity = RecoDecay::y(threeVecCand, massDstar);
        } else {
          isRotatedCandidate = 0;
          pxDau = candidate.pxProng1();
          pyDau = candidate.pyProng1();
          pzDau = candidate.pzProng1();
          std::array<float, 3> threeVecCand = RecoDecay::pVec(std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()},
                                                              std::array{candidate.pyProng0Charm(), candidate.pxProng0Charm(), candidate.pzProng0Charm()},
                                                              std::array{candidate.pxProng1Charm(), candidate.pyProng1Charm(), candidate.pzProng1Charm()});
          pxCharmHad = threeVecCand[0];
          pyCharmHad = threeVecCand[1];
          pzCharmHad = threeVecCand[2];
          invMassCharmHad = candidate.m();
          invMassD0 = candidate.invMassCharm();
          rapidity = candidate.y();
        }
        invMassCharmHadForSparse = invMassCharmHad - invMassD0;

        if constexpr (withMl) {
          std::copy_n(candidate.mlScores().begin(), outputMl.size(), outputMl.begin());
        }
      }
      if (invMassCharmHadForSparse < invMassMin || invMassCharmHadForSparse > invMassMax) {
        continue;
      }

      /// apply rapidity selection on the reconstructed candidate
      if (std::abs(rapidity) > maxAbsRapidityCut) {
        continue;
      }

      ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(pxDau, pyDau, pzDau, massDau);
      ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(pxCharmHad, pyCharmHad, pzCharmHad, invMassCharmHad);
      ROOT::Math::Boost boost{fourVecMother.BoostToCM()};
      ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);
      ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();
      ROOT::Math::PxPyPzMVector fourVecJet = ROOT::Math::PxPyPzMVector(jet.px(), jet.py(), jet.pz(), jet.mass());
      ROOT::Math::PxPyPzMVector fourVecJetCM = boost(fourVecJet);

      float ptCharmHad = std::sqrt(pxCharmHad * pxCharmHad + pyCharmHad * pyCharmHad); // this definition is valid for both rotated and original candidates

      if (!isCandidateInSignalRegion) { // it could be that only one mass hypothesis is in signal region
        isCandidateInSignalRegion = isInSignalRegion<channel>(invMassCharmHadForSparse);
      }
      float absEtaTrackMin = candidate.absEtaTrackMin();
      int numItsClsMin = candidate.numItsClsMin();
      int numTpcClsMin = candidate.numTpcClsMin();

      // helicity
      ROOT::Math::XYZVector helicityVec = fourVecMother.Vect();
      float cosThetaStarHelicity = -10.f;
      // production
      ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(pyCharmHad, -pxCharmHad, 0.f);
      float cosThetaStarProduction = -10.f;
      // jet axis
      ROOT::Math::XYZVector jetaxisVec = fourVecJetCM.Vect();
      float cosThetaStarJet = -10.f;

      float zParallel = helicityVec.Dot(jetaxisVec) / std::sqrt(jetaxisVec.Mag2());
      float jetPt = jet.pt();

      int8_t nMuons{0u};
      if constexpr (doMc) {
        nMuons = candidate.nTracksDecayed();
      }
      if (activateTHnSparseCosThStarHelicity) {
        // helicity
        cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(helicityVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Helicity>(invMassCharmHadForSparse, ptCharmHad, rapidity, invMassD0, cosThetaStarHelicity, outputMl, isRotatedCandidate, origin, ptBhadMother, absEtaTrackMin, numItsClsMin, numTpcClsMin, nMuons, partRecoDstar, zParallel, jetPt);
      }
      if (activateTHnSparseCosThStarProduction) {
        // production
        cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(normalVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::Production>(invMassCharmHadForSparse, ptCharmHad, rapidity, invMassD0, cosThetaStarProduction, outputMl, isRotatedCandidate, origin, ptBhadMother, absEtaTrackMin, numItsClsMin, numTpcClsMin, nMuons, partRecoDstar, zParallel, jetPt);
      }
      if (activateTHnSparseCosThStarJetAxis) {
        // jet axis
        cosThetaStarJet = jetaxisVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2()) / std::sqrt(jetaxisVec.Mag2());
        fillRecoHistos<channel, withMl, doMc, charm_polarisation::CosThetaStarType::JetAxis>(invMassCharmHadForSparse, ptCharmHad, rapidity, invMassD0, cosThetaStarJet, outputMl, isRotatedCandidate, origin, ptBhadMother, absEtaTrackMin, numItsClsMin, numTpcClsMin, nMuons, partRecoDstar, zParallel, jetPt);
      }
    } /// end loop over mass hypotheses

    return isCandidateInSignalRegion;
  }

  /////////////////////////
  //   Dstar analysis   ///
  /////////////////////////

  // Dstar with rectangular cuts
  void processDstar(aod::JetCollisions const& collisions,
                    DstarJets const& jets,
                    aod::CandidatesDstarData const&)
  {
    for (const auto& collision : collisions) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedDstarjets = jets.sliceBy(dstarJetsPerCollision, thisCollId);
      for (const auto& jet : groupedDstarjets) {
        for (const auto& dstarCandidate : jet.candidates_as<aod::CandidatesDstarData>()) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(jet, dstarCandidate, 0);
          for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
            runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, false>(jet, dstarCandidate, iRotation);
          }
          break; // hf jet should have only one Dstar candidate but for safety
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDstarPolarisationInJet, processDstar, "Process Dstar candidates without ML", true);

  // Dstar with ML cuts
  void processDstarWithMl(aod::JetCollisions const& collisions,
                          DstarJets const& jets,
                          aod::CandidatesDstarData const&)
  {
    for (const auto& collision : collisions) {
      if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
        continue;
      }
      auto thisCollId = collision.globalIndex();
      auto groupedDstarJets = jets.sliceBy(dstarJetsPerCollision, thisCollId);
      for (const auto& jet : groupedDstarJets) {
        for (const auto& dstarCandidate : jet.candidates_as<aod::CandidatesDstarData>()) {
          runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(jet, dstarCandidate, 0);
          for (int iRotation{1}; iRotation <= nBkgRotations; ++iRotation) {
            runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, false>(jet, dstarCandidate, iRotation);
          }
          break; // hf jet should have only one Dstar candidate but for safety
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDstarPolarisationInJet, processDstarWithMl, "Process Dstar candidates with ML", false);

  // Dstar in MC with rectangular cuts
  void processDstarMc(aod::JetMcCollisions const& mcCollisions,
                      aod::JetCollisions const& collisions,
                      JetMCDTable const& mcdJets,
                      aod::CandidatesDstarMCD const&)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mcCollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
          continue;
        }
        const auto dstarMcDJetsPerCollision = mcdJets.sliceBy(dstarMCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdJet : dstarMcDJetsPerCollision) {
          for (const auto& mcdDstarCand : mcdJet.candidates_as<aod::CandidatesDstarMCD>()) {
            runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, false, true>(mcdJet, mcdDstarCand, 0);
            break;
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDstarPolarisationInJet, processDstarMc, "Process Dstar candidates in MC without ML", false);

  // Dstar in MC with ML cuts
  void processDstarMcWithMl(aod::JetMcCollisions const& mcCollisions,
                            aod::JetCollisions const& collisions,
                            JetMCDTable const& mcdJets,
                            aod::CandidatesDstarMCD const&)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto collisionsPerMCCollision = collisions.sliceBy(collisionsPerMCCollisionPreslice, mcCollision.globalIndex());
      for (const auto& collision : collisionsPerMCCollision) {
        if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
          continue;
        }
        const auto dstarMcdJetsPerCollision = mcdJets.sliceBy(dstarMCDJetsPerCollisionPreslice, collision.globalIndex());
        for (const auto& mcdJet : dstarMcdJetsPerCollision) {
          for (const auto& mcdDstarCand : mcdJet.candidates_as<aod::CandidatesDstarMCD>()) {
            runPolarisationAnalysis<charm_polarisation::DecayChannel::DstarToDzeroPi, true, true>(mcdJet, mcdDstarCand, 0);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDstarPolarisationInJet, processDstarMcWithMl, "Process Dstar candidates in MC with ML", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskDstarPolarisationInJet>(cfgc)};
}
