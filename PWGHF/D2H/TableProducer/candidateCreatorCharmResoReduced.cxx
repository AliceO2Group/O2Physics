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

/// \file candidateCreatorCharmResoReduced.cxx
/// \brief Reconstruction of Resonance candidates
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, Università degli Studi di Torino
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, INFN Bari

#include "PWGHF/D2H/Core/SelectorCutsRedDataFormat.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsAnalysis.h"
#include "PWGHF/Utils/utilsMcMatching.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

enum Selections : uint8_t {
  NoSel = 0,
  DSel,
  BachSel,
  NSelSteps
};

enum D0Sel : uint8_t {
  SelectedD0 = 0,
  SelectedD0Bar
};

enum DMesonType : uint8_t {
  Dplus = 1,
  Dstar,
  D0
};

enum BachelorType : uint8_t {
  V0 = 1,
  Track
};

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};

enum TrackType : uint8_t {
  Pion = 0,
  Kaon,
  Proton
};

struct HfCandidateCreatorCharmResoReduced {
  // Produces: Tables with resonance info
  Produces<aod::HfCandCharmReso> rowCandidateReso;
  // Tables with bachelors indices for MC matching and Task
  Produces<aod::Hf3PrV0Ids> rowCandidateResoIndices3PrV0s;
  Produces<aod::HfDstarV0Ids> rowCandidateResoIndicesDstarV0s;
  Produces<aod::Hf2PrV0Ids> rowCandidateResoIndices2PrV0s;
  Produces<aod::Hf3PrTrkIds> rowCandidateResoIndices3PrTrks;
  Produces<aod::HfDstarTrkIds> rowCandidateResoIndicesDstarTrks;
  Produces<aod::Hf2PrTrkIds> rowCandidateResoIndices2PrTrks;
  // D Configurables
  struct : ConfigurableGroup {
    std::string prefix = "dmesonsCuts";
    Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{hf_cuts_d_daughter::vecBinsPt}, "pT bin limits for D daughter cuts"};
    Configurable<LabeledArray<double>> cutsD{"cutsD", {hf_cuts_d_daughter::Cuts[0], hf_cuts_d_daughter::NBinsPt, hf_cuts_d_daughter::NCutVars, hf_cuts_d_daughter::labelsPt, hf_cuts_d_daughter::labelsCutVar}, "D daughter selections"};
    Configurable<bool> keepSideBands{"keepSideBands", false, "flag to keep events from D meson sidebands for backgorund estimation"};
  } cfgDmesCuts;
  // V0 cuts configurables
  struct : ConfigurableGroup {
    std::string prefix = "v0Cuts";
    Configurable<LabeledArray<double>> cutsV0{"cutsV0", {hf_cuts_v0_daughter::Cuts[0], hf_cuts_v0_daughter::NBinsPt, hf_cuts_v0_daughter::NCutVars, hf_cuts_v0_daughter::labelsPt, hf_cuts_v0_daughter::labelsCutVar}, "V0 daughter selections"};
    Configurable<std::vector<double>> binsPtV0{"binsPtV0", std::vector<double>{hf_cuts_v0_daughter::vecBinsPt}, "pT bin limits for V0 daughter cuts"};
    Configurable<int> v0Type{"v0Type", 0, "V0 type to be selected (0: K0s, 1: Lambda"};
  } cfgV0Cuts;
  // Track cuts configurables
  struct : ConfigurableGroup {
    std::string prefix = "trackCuts";
    Configurable<LabeledArray<double>> cutsTrk{"cutsTrk", {hf_cuts_track_daughter::Cuts[0], hf_cuts_track_daughter::NCutVars, hf_cuts_track_daughter::labelsCutVar}, "Track daughter selections, set to -1 to disable cuts"};
    Configurable<int> massHypo{"massHypo", 1, "Mass Hypothesis for the track daughters (0: pion, 1: kaon, 2: proton)"};
  } cfgTrackCuts;
  // Mixed Event configurables
  struct : ConfigurableGroup {
    std::string prefix = "mixedEvent";
    Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
    Configurable<int> numberEventsToSkip{"numberEventsToSkip", -1, "Number of events to Skip in ME process"};
    ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 45., 60., 75., 95, 250}, "event multiplicity pools (PV contributors for now)"};
    ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -4, -1, 1, 4, 10.0}, "z vertex position pools"};
  } cfgMixedEvent;
  struct : ConfigurableGroup {
    std::string prefix = "trackRotation";
    Configurable<bool> enable{"enable", false, "enable rotation of track/V0 for background estimation"};
    Configurable<int> numRotations{"numRotations", 12, "number of track/V0 rotations"};
    Configurable<float> minRotAngleMultByPi{"minRotAngleMultByPi", 5. / 6, "Minimum angle rotation for track rotation, to be multiplied by pi"};
    Configurable<float> maxRotAngleMultByPi{"maxRotAngleMultByPi", 7. / 6, "Maximum angle rotation for track rotation, to be multiplied by pi"};
  } cfgTrackRotation;
  // Histogram axes configurables
  struct : ConfigurableGroup {
    std::string prefix = "histAxes";
    ConfigurableAxis axisPtD{"axisPtD", {100, 0., 50}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisPtV0{"axisPtV0", {100, 0., 50}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisPtReso{"axisPtReso", {100, 0., 50}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisMassD{"axisMassD", {100, 1.7f, 2.1f}, "inv. mass (D) (GeV/#it{c}^{2})"};
    ConfigurableAxis axisMassV0{"axisMassV0", {100, 0.45f, 0.55f}, "inv. mass (V_{0}) (GeV/#it{c}^{2})"};
    ConfigurableAxis axisMassDsj{"axisMassDsj", {400, 0.49f, 0.89f}, "inv. mass (DV_{0}) (GeV/#it{c}^{2})"};
  } cfgHistAxes;
  // Other Configurables
  Configurable<bool> rejectPairsWithCommonDaughter{"rejectPairsWithCommonDaughter", true, "flag to reject the pairs that share a daughter track if not done in the derived data creation"};
  Configurable<bool> useDeltaMass{"useDeltaMass", true, "Use Delta Mass for resonance invariant Mass calculation"};

  SliceCache cache;
  Preslice<aod::HfRedVzeros> candsV0PerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRedTrkNoParams> candsTrackPerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRed3PrNoTrks> cands3PrPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRedDstarNoTrks> candsDstarPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRed2PrNoTrks> cands2PrPerCollision = hf_track_index_reduced::hfRedCollisionId;

  HistogramRegistry registry{"registry"};

  float bkgRotationAngleStep{0.f};
  void init(InitContext const&)
  {
    // histograms
    registry.add("hMassDmesDauVsPt", "D daughter candidates inv. mass", {HistType::kTH2F, {cfgHistAxes.axisMassD, cfgHistAxes.axisPtD}});
    registry.add("hMassV0DauVsPt", "V0 daughter candidates inv. mass", {HistType::kTH2F, {cfgHistAxes.axisMassV0, cfgHistAxes.axisPtV0}});
    registry.add("hMassResoVsPt", "Resonance candidates inv. mass", {HistType::kTH2F, {cfgHistAxes.axisMassDsj, cfgHistAxes.axisPtReso}});
    registry.add("hNPvContCorr", "Collision number of PV contributors ; N contrib ; N contrib", {HistType::kTH2F, {{100, 0, 250}, {100, 0, 250}}});
    registry.add("hZvertCorr", "Collision Z Vtx ; z PV [cm] ; z PV [cm]", {HistType::kTH2F, {{120, -12., 12.}, {120, -12., 12.}}});
    constexpr int kNBinsSelections = Selections::NSelSteps;
    std::string labels[kNBinsSelections];
    labels[Selections::NoSel] = "No selection";
    labels[Selections::DSel] = "D Candidates Selection";
    labels[Selections::BachSel] = "D & other bach. Selection";
    static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
    registry.add("hSelections", "Selections", {HistType::kTH1F, {axisSelections}});
    for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
      registry.get<TH1>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    bkgRotationAngleStep = (cfgTrackRotation.numRotations > 1) ? (cfgTrackRotation.maxRotAngleMultByPi - cfgTrackRotation.minRotAngleMultByPi) * constants::math::PI / (cfgTrackRotation.numRotations - 1) : 0.;
  }

  bool isInMassInterval(float invMass, int ptBin)
  {
    if (!cfgDmesCuts.keepSideBands) {
      return (invMass >= cfgDmesCuts.cutsD->get(ptBin, "invMassSignalLow") && invMass <= cfgDmesCuts.cutsD->get(ptBin, "invMassSignalHigh"));
    }
    return ((invMass >= cfgDmesCuts.cutsD->get(ptBin, "invMassLeftSBLow") && invMass <= cfgDmesCuts.cutsD->get(ptBin, "invMassLeftSBHigh")) ||
            (invMass >= cfgDmesCuts.cutsD->get(ptBin, "invMassRightSBLow") && invMass <= cfgDmesCuts.cutsD->get(ptBin, "invMassRightSBHigh")) ||
            (invMass >= cfgDmesCuts.cutsD->get(ptBin, "invMassSignalLow") && invMass <= cfgDmesCuts.cutsD->get(ptBin, "invMassSignalHigh")));
  }

  /// Basic selection of D candidates
  /// \param candD is the reduced D meson candidate
  /// \return true if selections are passed
  template <DMesonType DType, typename DRedTable>
  uint8_t selctionFlagBachD(DRedTable const& candD)
  {
    uint8_t selection = {BIT(D0Sel::SelectedD0) | BIT(D0Sel::SelectedD0Bar)};
    float invMassD{0.};
    float const ptD = candD.pt();
    int const ptBin = findBin(cfgDmesCuts.binsPtD, ptD);
    if (ptBin == -1) {
      return 0;
    }
    if constexpr (DType == DMesonType::Dplus) {
      invMassD = candD.invMassDplus();
      if (!isInMassInterval(invMassD, ptBin)) {
        return 0;
      }
    } else if constexpr (DType == DMesonType::Dstar) {
      if (candD.sign() > 0) {
        invMassD = candD.invMassDstar() - candD.invMassD0();
      } else {
        invMassD = candD.invMassAntiDstar() - candD.invMassD0Bar();
      }
      if (!isInMassInterval(invMassD, ptBin)) {
        return 0;
      }
    } else if constexpr (DType == DMesonType::D0) {
      if (TESTBIT(candD.selFlagD0(), D0Sel::SelectedD0)) {
        invMassD = candD.invMassD0();
        if (!isInMassInterval(invMassD, ptBin)) {
          CLRBIT(selection, D0Sel::SelectedD0);
        }
      } else {
        CLRBIT(selection, D0Sel::SelectedD0);
      }
      if (TESTBIT(candD.selFlagD0(), D0Sel::SelectedD0Bar)) {
        invMassD = candD.invMassD0Bar();
        if (!isInMassInterval(invMassD, ptBin)) {
          CLRBIT(selection, D0Sel::SelectedD0Bar);
        }
      } else {
        CLRBIT(selection, D0Sel::SelectedD0Bar);
      }
    }
    return selection;
  }

  /// Basic selection of V0 and track candidates
  /// \param candV0 is the reduced V0 candidate
  /// \return true if selections are passed
  template <typename V0RedTable>
  bool isV0Selected(V0RedTable const& candV0)
  {
    int const ptBin = findBin(cfgV0Cuts.binsPtV0, candV0.pt());
    const float invMassLow = cfgV0Cuts.cutsV0->get(ptBin, "invMassLow");
    const float invMassHigh = cfgV0Cuts.cutsV0->get(ptBin, "invMassHigh");
    if (ptBin == -1) {
      return false;
    }
    if (cfgV0Cuts.v0Type == V0Type::K0s) { // K0s
      if (!TESTBIT(candV0.v0Type(), V0Type::K0s)) {
        return false;
      }
      if (candV0.invMassK0s() < invMassLow || candV0.invMassK0s() > invMassHigh) {
        return false;
      }
    } else if (cfgV0Cuts.v0Type == V0Type::Lambda) { // Lambda
      if (!TESTBIT(candV0.v0Type(), V0Type::Lambda) && !TESTBIT(candV0.v0Type(), V0Type::AntiLambda)) {
        return false;
      }

      if (TESTBIT(candV0.v0Type(), V0Type::Lambda)) {
        if (candV0.invMassLambda() < invMassLow || candV0.invMassLambda() > invMassHigh) {
          return false;
        }
      }

      if (TESTBIT(candV0.v0Type(), V0Type::AntiLambda)) {
        if (candV0.invMassAntiLambda() < invMassLow || candV0.invMassAntiLambda() > invMassHigh) {
          return false;
        }
      }
    } else {
      LOG(error) << "Unsupported V0 type for selection: " << cfgV0Cuts.v0Type;
      return false;
    }
    // selection on kinematics and topology
    if (candV0.dca() > cfgV0Cuts.cutsV0->get(ptBin, "dcaMax") || candV0.cpa() < cfgV0Cuts.cutsV0->get(ptBin, "cpaMin") || candV0.v0Radius() < cfgV0Cuts.cutsV0->get(ptBin, "radiusMin")) {
      return false;
    }
    return true;
  }

  // Basic selection of track candidates
  /// \param candTr is the reduced track candidate
  /// \return true if selections are passed
  template <typename TrkRedTable>
  bool isTrackSelected(TrkRedTable const& candTr)
  {
    // pT selection
    if (cfgTrackCuts.cutsTrk->get("ptMin") > 0 && candTr.pt() < cfgTrackCuts.cutsTrk->get("ptMin")) {
      return false;
    }
    // ITS quality selection
    if (cfgTrackCuts.cutsTrk->get("itsNClsMin") > 0 && candTr.itsNCls() < cfgTrackCuts.cutsTrk->get("itsNClsMin")) {
      return false;
    }
    // TPC quality selection
    if (cfgTrackCuts.cutsTrk->get("tpcNCrossedRowsMin") > 0 && candTr.tpcNClsCrossedRows() < cfgTrackCuts.cutsTrk->get("tpcNCrossedRowsMin")) {
      return false;
    }
    if (cfgTrackCuts.cutsTrk->get("tpcChi2Max") > 0 && candTr.tpcChi2NCl() > cfgTrackCuts.cutsTrk->get("tpcChi2Max")) {
      return false;
    }
    // PID selection
    if (cfgTrackCuts.massHypo == TrackType::Pion) { // Pion
      if (cfgTrackCuts.cutsTrk->get("nSigmaTpc") > 0 && candTr.tpcNSigmaPi() > cfgTrackCuts.cutsTrk->get("nSigmaTpc")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaTof") > 0 && candTr.tofNSigmaPi() > cfgTrackCuts.cutsTrk->get("nSigmaTof")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaComb") > 0 && candTr.tpcTofNSigmaPi() > cfgTrackCuts.cutsTrk->get("nSigmaComb")) {
        return false;
      }
    } else if (cfgTrackCuts.massHypo == TrackType::Kaon) { // Kaon
      if (cfgTrackCuts.cutsTrk->get("nSigmaTpc") > 0 && candTr.tpcNSigmaKa() > cfgTrackCuts.cutsTrk->get("nSigmaTpc")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaTof") > 0 && candTr.tofNSigmaKa() > cfgTrackCuts.cutsTrk->get("nSigmaTof")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaComb") > 0 && candTr.tpcTofNSigmaKa() > cfgTrackCuts.cutsTrk->get("nSigmaComb")) {
        return false;
      }
    } else if (cfgTrackCuts.massHypo == TrackType::Proton) { // Proton
      if (cfgTrackCuts.cutsTrk->get("nSigmaTpc") > 0 && candTr.tpcNSigmaPr() > cfgTrackCuts.cutsTrk->get("nSigmaTpc")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaTof") > 0 && candTr.tofNSigmaPr() > cfgTrackCuts.cutsTrk->get("nSigmaTof")) {
        return false;
      }
      if (cfgTrackCuts.cutsTrk->get("nSigmaComb") > 0 && candTr.tpcTofNSigmaPr() > cfgTrackCuts.cutsTrk->get("nSigmaComb")) {
        return false;
      }
    } else {
      LOG(error) << "Unsupported mass hypothesis for track selection: " << cfgTrackCuts.massHypo;
      return false;
    }
    return true;
  }

  /// Fill the output tables with the resonance candidates
  /// \param collision is the collision information
  /// \param candD is the reduced D meson candidate
  /// \param candV0Tr is the reduced V0 or track candidate
  /// \tparam dType is the type of D meson (Dplus, Dstar, D0)
  /// \tparam bachType is the type of bachelor (V0 or Track)
  template <DMesonType DType, BachelorType BachType, typename Coll, typename DRedTable, typename V0TrRedTable>
  void fillOutputTables(Coll const& collision,
                        DRedTable const& candD,
                        V0TrRedTable const& candV0Tr,
                        int selectionFlag)
  {
    std::vector<std::array<float, 3>> pVectorCharmProngs = {candD.pVectorProng0(), candD.pVectorProng1()};
    std::array<float, 3> pVecD = candD.pVector();

    int const numFills = (cfgTrackRotation.enable) ? cfgTrackRotation.numRotations : 1; // number of times we fil the tables: default 1, but more in case of track rotation

    for (int iFill{0}; iFill < numFills; ++iFill) {

      std::array<float, 3> pVecV0Tr = candV0Tr.pVector();
      if (cfgTrackRotation.enable) { // let's rotate
        float const bkgRotAngle = cfgTrackRotation.minRotAngleMultByPi * constants::math::PI + bkgRotationAngleStep * iFill;
        pVecV0Tr = std::array<float, 3>{candV0Tr.px() * std::cos(bkgRotAngle) - candV0Tr.py() * std::sin(bkgRotAngle), candV0Tr.px() * std::sin(bkgRotAngle) + candV0Tr.py() * std::cos(bkgRotAngle), candV0Tr.pz()};
      }

      float invMassReso{-1}, invMassV0Tr{-1}, invMassD{-1};
      int8_t signReso{0}, isWrongSign{0};
      double const ptReso = RecoDecay::pt(RecoDecay::sumOfVec(pVecV0Tr, pVecD));

      if constexpr (DType == DMesonType::Dplus) {
        invMassD = candD.invMassDplus();
        pVectorCharmProngs.push_back(candD.pVectorProng2());
        if constexpr (BachType == BachelorType::V0) {
          if (cfgV0Cuts.v0Type == V0Type::K0s) { // K0s
            invMassV0Tr = candV0Tr.invMassK0s();
            signReso = candD.sign();
            if (useDeltaMass) {
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassK0Short});
            }
          } else if (cfgV0Cuts.v0Type == V0Type::Lambda) { // Lambda
            if (TESTBIT(candV0Tr.v0Type(), V0Type::Lambda)) {
              invMassV0Tr = candV0Tr.invMassLambda();
              signReso = candD.sign();
              isWrongSign = candD.sign() < 0 ? 1 : 0;
            } else if (TESTBIT(candV0Tr.v0Type(), V0Type::AntiLambda)) {
              invMassV0Tr = candV0Tr.invMassAntiLambda();
              signReso = candD.sign();
              isWrongSign = candD.sign() > 0 ? 1 : 0;
            }
            if (useDeltaMass) {
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda}) - invMassD;
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassLambda});
            }
          }
          rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD,
                           invMassV0Tr,
                           signReso,
                           isWrongSign);
          rowCandidateResoIndices3PrV0s(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
          registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
          registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
          registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
        } else if constexpr (BachType == BachelorType::Track) {
          signReso = candD.sign() + candV0Tr.sign();
          isWrongSign = candD.sign() * candV0Tr.sign() > 0 ? 1 : 0;
          if (cfgTrackCuts.massHypo == TrackType::Pion) { // Pion
            invMassV0Tr = MassPiPlus;
            if (useDeltaMass) {
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassPiPlus}) - invMassD;
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassPiPlus});
            }
          } else if (cfgTrackCuts.massHypo == TrackType::Kaon) { // Kaon
            invMassV0Tr = MassKPlus;
            if (useDeltaMass) {
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassKPlus}) - invMassD;
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassKPlus});
            }
          } else if (cfgTrackCuts.massHypo == TrackType::Proton) { // Proton
            invMassV0Tr = MassProton;
            if (useDeltaMass) {
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton}) - invMassD;
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassProton});
            }
          }
          rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD,
                           invMassV0Tr,
                           signReso,
                           isWrongSign);
          rowCandidateResoIndices3PrTrks(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
          registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
          registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
          registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
        }
      } else if constexpr (DType == DMesonType::Dstar) {
        float invMassD0;
        if (candD.sign() > 0) {
          invMassD = candD.invMassDstar();
          invMassD0 = candD.invMassD0();
        } else {
          invMassD = candD.invMassAntiDstar();
          invMassD0 = candD.invMassD0Bar();
        }
        pVectorCharmProngs.push_back(candD.pVectorProng2());
        if constexpr (BachType == BachelorType::V0) {
          signReso = candD.sign();
          if (cfgV0Cuts.v0Type == V0Type::K0s) { // K0s
            invMassV0Tr = candV0Tr.invMassK0s();
            if (useDeltaMass) {
              if (candD.sign() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              }
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassK0Short});
            }
          } else if (cfgV0Cuts.v0Type == V0Type::Lambda) { // Lambda
            if (TESTBIT(candV0Tr.v0Type(), V0Type::Lambda)) {
              invMassV0Tr = candV0Tr.invMassLambda();
            } else if (TESTBIT(candV0Tr.v0Type(), V0Type::AntiLambda)) {
              invMassV0Tr = candV0Tr.invMassAntiLambda();
              isWrongSign = 1;
            }
            if (useDeltaMass) {
              if (candD.sign() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda}) - invMassD;
              }
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassLambda});
            }
          }
          rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD - invMassD0,
                           invMassV0Tr,
                           signReso,
                           isWrongSign);
          rowCandidateResoIndicesDstarV0s(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
          registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
          registry.fill(HIST("hMassDmesDauVsPt"), invMassD - invMassD0, candD.pt());
          registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
        } else if constexpr (BachType == BachelorType::Track) {
          signReso = candD.sign() + candV0Tr.sign();
          isWrongSign = candD.sign() * candV0Tr.sign() > 0 ? 1 : 0;
          if (cfgTrackCuts.massHypo == TrackType::Pion) { // Pion
            invMassV0Tr = MassPiPlus;
            if (useDeltaMass) {
              if (candD.sign() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassPiPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassPiPlus}) - invMassD;
              }
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassPiPlus});
            }
          } else if (cfgTrackCuts.massHypo == TrackType::Kaon) { // Kaon
            invMassV0Tr = MassKPlus;
            if (useDeltaMass) {
              if (candD.sign() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassKPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassKPlus}) - invMassD;
              }
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassKPlus});
            }
          } else if (cfgTrackCuts.massHypo == TrackType::Proton) { // Proton
            invMassV0Tr = MassProton;
            if (useDeltaMass) {
              if (candD.sign() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton}) - invMassD;
              }
            } else {
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassProton});
            }
          }
          rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD - invMassD0,
                           invMassV0Tr,
                           signReso,
                           isWrongSign);
          rowCandidateResoIndicesDstarTrks(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
          registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
          registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
          registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
        }
      } else if constexpr (DType == DMesonType::D0) {
        // D0
        if (TESTBIT(selectionFlag, D0Sel::SelectedD0)) {
          invMassD = candD.invMassD0();
          if constexpr (BachType == BachelorType::V0) {
            signReso = 0;
            if (cfgV0Cuts.v0Type == V0Type::K0s) { // K0s
              invMassV0Tr = candV0Tr.invMassK0s();
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassK0Short}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0, MassK0Short});
              }
            } else if (cfgV0Cuts.v0Type == V0Type::Lambda) { // Lambda
              if (TESTBIT(candV0Tr.v0Type(), V0Type::Lambda)) {
                invMassV0Tr = candV0Tr.invMassLambda();
              } else if (TESTBIT(candV0Tr.v0Type(), V0Type::AntiLambda)) {
                invMassV0Tr = candV0Tr.invMassAntiLambda();
                isWrongSign = 1;
              }
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassLambda}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0, MassLambda});
              }
            }
            rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                             pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                             invMassReso,
                             invMassD,
                             invMassV0Tr,
                             signReso,
                             isWrongSign);
            rowCandidateResoIndices2PrV0s(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
            registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
            registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
            registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
          } else if constexpr (BachType == BachelorType::Track) {
            signReso = candV0Tr.sign();
            isWrongSign = candV0Tr.sign() > 0 ? 0 : 1;
            if (cfgTrackCuts.massHypo == TrackType::Pion) { // Pion
              invMassV0Tr = MassPiPlus;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0, MassPiPlus});
              }
            } else if (cfgTrackCuts.massHypo == TrackType::Kaon) { // Kaon
              invMassV0Tr = MassKPlus;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassKPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0, MassKPlus});
              }
            } else if (cfgTrackCuts.massHypo == TrackType::Proton) { // Proton
              invMassV0Tr = MassProton;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassProton}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0, MassProton});
              }
            }
            rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                             pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                             invMassReso,
                             invMassD,
                             invMassV0Tr,
                             signReso,
                             isWrongSign);
            rowCandidateResoIndices2PrTrks(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
            registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
            registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
            registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
          }
        }
        // D0bar
        if (TESTBIT(selectionFlag, D0Sel::SelectedD0Bar)) {
          invMassD = candD.invMassD0Bar();
          if constexpr (BachType == BachelorType::V0) {
            signReso = 0;
            if (cfgV0Cuts.v0Type == V0Type::K0s) { // K0s
              invMassV0Tr = candV0Tr.invMassK0s();
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0Bar, MassK0Short});
              }
            } else if (cfgV0Cuts.v0Type == V0Type::Lambda) { // Lambda
              if (TESTBIT(candV0Tr.v0Type(), V0Type::Lambda)) {
                invMassV0Tr = candV0Tr.invMassLambda();
                isWrongSign = 1;
              } else if (TESTBIT(candV0Tr.v0Type(), V0Type::AntiLambda)) {
                invMassV0Tr = candV0Tr.invMassAntiLambda();
              }
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassKPlus, MassPiPlus, MassLambda}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0Bar, MassLambda});
              }
            }
            rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                             pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                             invMassReso,
                             invMassD,
                             invMassV0Tr,
                             signReso,
                             isWrongSign);
            rowCandidateResoIndices2PrV0s(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
            registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
            registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
            registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
          } else if constexpr (BachType == BachelorType::Track) {
            signReso = candV0Tr.sign();
            isWrongSign = candV0Tr.sign() > 0 ? 1 : 0;
            if (cfgTrackCuts.massHypo == TrackType::Pion) { // Pion
              invMassV0Tr = MassPiPlus;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassKPlus, MassPiPlus, MassPiPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0Bar, MassPiPlus});
              }
            } else if (cfgTrackCuts.massHypo == TrackType::Kaon) { // Kaon
              invMassV0Tr = MassKPlus;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassKPlus, MassPiPlus, MassKPlus}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0Bar, MassKPlus});
              }
            } else if (cfgTrackCuts.massHypo == TrackType::Proton) { // Proton
              invMassV0Tr = MassProton;
              if (useDeltaMass) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassKPlus, MassPiPlus, MassProton}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassD0Bar, MassProton});
              }
            }
            rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                             pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                             invMassReso,
                             invMassD,
                             invMassV0Tr,
                             signReso,
                             isWrongSign);
            rowCandidateResoIndices2PrTrks(collision.globalIndex(), candD.globalIndex(), candV0Tr.globalIndex());
            registry.fill(HIST("hMassResoVsPt"), invMassReso, ptReso);
            registry.fill(HIST("hMassDmesDauVsPt"), invMassD, candD.pt());
            registry.fill(HIST("hMassV0DauVsPt"), invMassV0Tr, candV0Tr.pt());
          }
        }
      }
    }
  }

  template <DMesonType DType, BachelorType BachType, typename Coll, typename DRedTable, typename V0TrRedTable>
  void runCandidateCreation(Coll const& collision,
                            DRedTable const& candsD,
                            V0TrRedTable const& candsV0Tr)
  {
    // loop on D candidates
    // LOG(info) << "Number of D candidates: " << candsD.size() << ", Number of V0/Track candidates: " << candsV0Tr.size();
    for (const auto& candD : candsD) {
      // selection of D candidates
      registry.fill(HIST("hSelections"), 1);
      uint8_t const selFlagD = selctionFlagBachD<DType>(candD);
      if (selFlagD == 0) {
        continue;
      }
      registry.fill(HIST("hSelections"), 1 + Selections::DSel);
      std::vector<int> dDaughtersIDs = {candD.prong0Id(), candD.prong1Id()};
      if constexpr (DType == DMesonType::Dstar || DType == DMesonType::Dplus) {
        dDaughtersIDs.push_back(candD.prong2Id());
      }
      // loop on V0 or track candidates
      bool alreadyCounted{false};
      for (const auto& candV0Tr : candsV0Tr) {
        if constexpr (BachType == BachelorType::V0) { // Case: V0
          if (rejectPairsWithCommonDaughter && (std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.prong0Id()) != dDaughtersIDs.end() || std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.prong1Id()) != dDaughtersIDs.end())) {
            continue;
          }
          if (!isV0Selected(candV0Tr)) {
            continue;
          }
          if (!alreadyCounted) {
            registry.fill(HIST("hSelections"), 1 + Selections::BachSel);
            alreadyCounted = true;
          }
        } else if constexpr (BachType == BachelorType::Track) { // Case: Track
          if (rejectPairsWithCommonDaughter && std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.trackId()) != dDaughtersIDs.end()) {
            continue;
          }
          if (!isTrackSelected(candV0Tr)) {
            continue;
          }
          if (!alreadyCounted) {
            registry.fill(HIST("hSelections"), 1 + Selections::BachSel);
            alreadyCounted = true;
          }
        }
        // Filling of tables and histograms
        fillOutputTables<DType, BachType>(collision, candD, candV0Tr, selFlagD);
      } // end of loop on V0/Track candidates
    } // end of loop on D candidates
  } // end of function

  // Process data with Mixed Event
  /// \tparam fillMl is a flag to Fill ML scores if present
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param DRedTable is the D bachelors table
  /// \param V0TrRedTable is the V0/Track bachelors table
  template <DMesonType DType, BachelorType BachType, typename Coll, typename DRedTable, typename V0TrRedTable>
  void runCandidateCreationMixedEvent(Coll const& collisions,
                                      DRedTable const& candsD,
                                      V0TrRedTable const& candsV0Tr)
  {
    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
    BinningType const corrBinning{{cfgMixedEvent.zPoolBins, cfgMixedEvent.multPoolBins}, true};
    auto bachTuple = std::make_tuple(candsD, candsV0Tr);
    Pair<Coll, DRedTable, V0TrRedTable, BinningType> const pairs{corrBinning, cfgMixedEvent.numberEventsMixed, cfgMixedEvent.numberEventsToSkip, collisions, bachTuple, &cache};
    for (const auto& [collision1, bachDs, collision2, bachV0Trs] : pairs) {
      registry.fill(HIST("hNPvContCorr"), collision1.numContrib(), collision2.numContrib());
      registry.fill(HIST("hZvertCorr"), collision1.posZ(), collision2.posZ());
      for (const auto& [bachD, bachV0Tr] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(bachDs, bachV0Trs))) {
        // Apply analysis selections on D and V0 bachelors
        uint8_t const selFlagD = selctionFlagBachD<DType>(bachD);
        if (selFlagD == 0) {
          continue;
        }
        if constexpr (BachType == BachelorType::V0) {
          if (!isV0Selected(bachV0Tr)) {
            continue;
          }
        } else if constexpr (BachType == BachelorType::Track) {
          if (!isTrackSelected(bachV0Tr)) {
            continue;
          }
        }
        fillOutputTables<DType, BachType>(collision1, bachD, bachV0Tr, selFlagD);
      }
    }
  } // runCandidateCreationMixedEvent

  // List of Process Functions
  void process3ProngV0s(aod::HfRedCollisions const& collisions,
                        aod::HfRed3PrNoTrks const& candsD,
                        aod::HfRedVzeros const& candsV0)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(cands3PrPerCollision, thisCollId);
      auto v0sThisColl = candsV0.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<DMesonType::Dplus, BachelorType::V0>(collision, candsDThisColl, v0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process3ProngV0s, "Process resonances decaying in a 3 prong D meson and a V0", true);

  void processDstarV0s(aod::HfRedCollisions const& collisions,
                       aod::HfRedDstarNoTrks const& candsD,
                       aod::HfRedVzeros const& candsV0)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDstarPerCollision, thisCollId);
      auto v0sThisColl = candsV0.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<DMesonType::Dstar, BachelorType::V0>(collision, candsDThisColl, v0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarV0s, "Process resonances decaying in a Dstar meson and a V0", false);

  void process2PrV0s(aod::HfRedCollisions const& collisions,
                     aod::HfRed2PrNoTrks const& candsD,
                     aod::HfRedVzeros const& candsV0)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(cands2PrPerCollision, thisCollId);
      auto v0sThisColl = candsV0.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<DMesonType::D0, BachelorType::V0>(collision, candsDThisColl, v0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process2PrV0s, "Process resonances decaying in a 2 prong D meson and a V0", false);

  void process3ProngTracks(aod::HfRedCollisions const& collisions,
                           aod::HfRed3PrNoTrks const& candsD,
                           HfRedTrkNoParams const& candsTr)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(cands3PrPerCollision, thisCollId);
      auto trksThisColl = candsTr.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<DMesonType::Dplus, BachelorType::Track>(collision, candsDThisColl, trksThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process3ProngTracks, "Process resonances decaying in a 3 prong D meson and a Track", false);

  void processDstarTracks(aod::HfRedCollisions const& collisions,
                          aod::HfRedDstarNoTrks const& candsD,
                          HfRedTrkNoParams const& candsTr)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDstarPerCollision, thisCollId);
      auto trksThisColl = candsTr.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<DMesonType::Dstar, BachelorType::Track>(collision, candsDThisColl, trksThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarTracks, "Process resonances decaying in a Dstar meson and a Track", false);

  void process2PrTracks(aod::HfRedCollisions const& collisions,
                        aod::HfRed2PrNoTrks const& candsD,
                        HfRedTrkNoParams const& candsTr)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(cands2PrPerCollision, thisCollId);
      auto trksThisColl = candsTr.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<DMesonType::D0, BachelorType::Track>(collision, candsDThisColl, trksThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process2PrTracks, "Process resonances decaying in a 2 prong D meson and a Track", false);

  // Mixed Event Process Functions
  void process3ProngV0sMixedEvent(aod::HfRedCollisions const& collisions,
                                  aod::HfRed3PrNoTrks const& candsD,
                                  aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<DMesonType::Dplus, BachelorType::V0>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process3ProngV0sMixedEvent, "Process mixed events for resonances decaying in a 3 prong D meson and a V0", false);

  void processDstarV0sMixedEvent(aod::HfRedCollisions const& collisions,
                                 aod::HfRedDstarNoTrks const& candsD,
                                 aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<DMesonType::Dstar, BachelorType::V0>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarV0sMixedEvent, "Process mixed events for resonances decaying in a Dstar meson and a V0", false);

  void process2PrV0sMixedEvent(aod::HfRedCollisions const& collisions,
                               aod::HfRed2PrNoTrks const& candsD,
                               aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<DMesonType::D0, BachelorType::V0>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process2PrV0sMixedEvent, "Process mixed events for resonances decaying in a 2 prong D meson and a V0", false);

  void process3ProngTracksMixedEvent(aod::HfRedCollisions const& collisions,
                                     aod::HfRed3PrNoTrks const& candsD,
                                     HfRedTrkNoParams const& candsTr)
  {
    runCandidateCreationMixedEvent<DMesonType::Dplus, BachelorType::Track>(collisions, candsD, candsTr);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process3ProngTracksMixedEvent, "Process mixed events for resonances decaying in a 3 prong D meson and a Track", false);

  void processDstarTracksMixedEvent(aod::HfRedCollisions const& collisions,
                                    aod::HfRedDstarNoTrks const& candsD,
                                    HfRedTrkNoParams const& candsTr)
  {
    runCandidateCreationMixedEvent<DMesonType::Dstar, BachelorType::Track>(collisions, candsD, candsTr);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarTracksMixedEvent, "Process mixed events for resonances decaying in a Dstar meson and a Track", false);

  void process2PrTracksMixedEvent(aod::HfRedCollisions const& collisions,
                                  aod::HfRed2PrNoTrks const& candsD,
                                  HfRedTrkNoParams const& candsTr)
  {
    runCandidateCreationMixedEvent<DMesonType::D0, BachelorType::Track>(collisions, candsD, candsTr);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, process2PrTracksMixedEvent, "Process mixed events for resonances decaying in a 2 prong D meson and a Track", false);

}; // struct HfCandidateCreatorCharmResoReduced

struct HfCandidateCreatorCharmResoReducedExpressions {

  Produces<aod::HfMcRecRedResos> rowResoMcRec;

  using CandResoWithIndices2PrV0s = soa::Join<aod::HfCandCharmReso, aod::Hf2PrV0Ids>;
  using CandResoWithIndices2PrTrks = soa::Join<aod::HfCandCharmReso, aod::Hf2PrTrkIds>;
  using CandResoWithIndices3PrV0s = soa::Join<aod::HfCandCharmReso, aod::Hf3PrV0Ids>;
  using CandResoWithIndices3PrTrks = soa::Join<aod::HfCandCharmReso, aod::Hf3PrTrkIds>;
  using CandResoWithIndicesDstarV0s = soa::Join<aod::HfCandCharmReso, aod::HfDstarV0Ids>;
  using CandResoWithIndicesDstarTrks = soa::Join<aod::HfCandCharmReso, aod::HfDstarTrkIds>;

  // Configurable axis
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisInvMassReso{"axisInvMassReso", {400, 0.49f, 0.89f}, "inv. mass (DV_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng0{"axisInvMassProng0", {200, 0.14, 0.17}, "inv. mass (D) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng1{"axisInvMassProng1", {200, 0.47, 0.53}, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisDebug{"axisDebug", {16, -0.5, 15.5}, "MC debug flag"};
  ConfigurableAxis axisOrigin{"axisOrigin", {3, -0.5, 2.5}, "MC origin flag"};
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hMassMcMatched", "Reso MC candidates Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcMatchedIncomplete", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcUnmatched", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcNoEntry", "Reso MC candidates w.o. entry in MC Reco table;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
  }

  /// Fill candidate information at MC reconstruction level
  /// \param rowsMcRec MC reco information on DPi pairs
  /// \param candsReso prong global indices of B0 candidates
  template <typename McRec, typename CandResoWithIndices>
  void fillResoMcRec(McRec const& rowsMcRec, CandResoWithIndices const& candsReso)
  {
    for (const auto& candReso : candsReso) {
      bool filledMcInfo{false};
      for (const auto& rowMcRec : rowsMcRec) {
        if ((rowMcRec.prong0Id() != candReso.prong0Id()) || (rowMcRec.prong1Id() != candReso.prong1Id())) {
          continue;
        }
        rowResoMcRec(rowMcRec.flagMcMatchRec(),
                     rowMcRec.flagMcMatchRecD(),
                     rowMcRec.flagMcMatchChanD(),
                     rowMcRec.debugMcRec(),
                     rowMcRec.origin(),
                     rowMcRec.ptGen(),
                     rowMcRec.invMassGen(),
                     rowMcRec.nTracksDecayed());
        filledMcInfo = true;
        if (std::abs(rowMcRec.flagMcMatchRec()) > 0 &&
            !TESTBIT(rowMcRec.debugMcRec(), hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched)) {
          registry.fill(HIST("hMassMcMatched"), candReso.invMass(), candReso.pt());
        } else if (std::abs(rowMcRec.flagMcMatchRec()) > 0 &&
                   TESTBIT(rowMcRec.debugMcRec(), hf_decay::hf_cand_reso::PartialMatchMc::ResoPartlyMatched)) {
          registry.fill(HIST("hMassMcMatchedIncomplete"), candReso.invMass(), candReso.pt());
        } else {
          registry.fill(HIST("hMassMcUnmatched"), candReso.invMass(), candReso.pt());
        }
        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the D-Pi creator
        rowResoMcRec(0, 0, 0, 0, 0, -1.f, -1.f, 0);
        registry.fill(HIST("hMassMcNoEntry"), candReso.invMass(), candReso.pt());
      }
    }
  }

  void processDstarV0Mc(aod::HfDstarV0McRec const& rowsMcRec, CandResoWithIndicesDstarV0s const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processDstarV0Mc, "Process resonances to Dstar V0 MC", false);

  void processDstarTrackMc(aod::HfDstarTrkMcRec const& rowsMcRec, CandResoWithIndicesDstarTrks const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processDstarTrackMc, "Process resonances to Dstar track MC", false);

  void process2PrV0Mc(aod::Hf2PrV0McRec const& rowsMcRec, CandResoWithIndices2PrV0s const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, process2PrV0Mc, "Process resonances to D0 V0 MC", false);

  void process2PrTrackMc(aod::Hf2PrTrkMcRec const& rowsMcRec, CandResoWithIndices2PrTrks const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, process2PrTrackMc, "Process resonances to D0 track MC", false);

  void process3PrV0Mc(aod::Hf3PrV0McRec const& rowsMcRec, CandResoWithIndices3PrV0s const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, process3PrV0Mc, "Process resonances to Dplus V0 MC", false);

  void process3PrTrackMc(aod::Hf3PrTrkMcRec const& rowsMcRec, CandResoWithIndices3PrTrks const& candsReso)
  {
    fillResoMcRec(rowsMcRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, process3PrTrackMc, "Process resonances to Dplus track MC", false);

  void processDummy(aod::HfCandCharmReso const&) {}
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processDummy, "Process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorCharmResoReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorCharmResoReducedExpressions>(cfgc)};
}
