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

#include "Common/Core/RecoDecay.h"

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
#include <cstdint>
#include <numeric>
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
  V0Sel,
  TrackSel,
  NSelSteps
};
enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  Ds2StarToDplusK0s,
  XcToDplusLambda,
  LambdaDminus,
  DstarTrack,
  D0Track
};

enum DType : uint8_t {
  Dplus = 1,
  Dstar,
  D0
};

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};

enum D0SelectionType : uint8_t {
  SelectedD0 = 0,
  SelectedD0Bar
};

enum DecayTypeMc : uint8_t {
  Ds1ToDStarK0ToD0PiK0s = 1,
  Ds2StarToDplusK0sToPiKaPiPiPi,
  Ds1ToDStarK0ToDPlusPi0K0s,
  Ds1ToDStarK0ToD0PiK0sPart,
  Ds1ToDStarK0ToD0NoPiK0sPart,
  Ds1ToDStarK0ToD0PiK0sOneMu,
  Ds2StarToDplusK0sOneMu
};

const int nBinsPt = 7;
constexpr double BinsPt[nBinsPt + 1] = {
  1.,
  2.,
  4.,
  6.,
  8.,
  12.,
  24.,
  1000.};
auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + nBinsPt + 1};

struct HfCandidateCreatorCharmResoReduced {
  // Produces: Tables with resonance info
  Produces<aod::HfCandCharmReso> rowCandidateReso;
  Produces<aod::HfCandChaResTr> rowCandidateResoTrack;
  // Optional daughter ML scores table
  Produces<aod::HfCharmResoMLs> mlScores;
  // Table with candidate indices for MC matching
  Produces<aod::HfResoIndices> rowCandidateResoIndices;

  // Configurables
  Configurable<bool> rejectDV0PairsWithCommonDaughter{"rejectDV0PairsWithCommonDaughter", true, "flag to reject the pairs that share a daughter track if not done in the derived data creation"};
  Configurable<bool> keepSideBands{"keepSideBands", false, "flag to keep events from D meson sidebands for backgorund estimation"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{vecBinsPt}, "Histogram pT bin limits"};
  // Daughters selection cuts
  Configurable<LabeledArray<double>> cutsDDaughter{"cutsDDaughter", {hf_cuts_d_daughter::Cuts[0], hf_cuts_d_daughter::NBinsPt, hf_cuts_d_daughter::NCutVars, hf_cuts_d_daughter::labelsPt, hf_cuts_d_daughter::labelsCutVar}, "D daughter selections"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{hf_cuts_d_daughter::vecBinsPt}, "pT bin limits for D daughter cuts"};
  Configurable<LabeledArray<double>> cutsV0Daughter{"cutsV0Daughter", {hf_cuts_v0_daughter::Cuts[0], hf_cuts_v0_daughter::NBinsPt, hf_cuts_v0_daughter::NCutVars, hf_cuts_v0_daughter::labelsPt, hf_cuts_v0_daughter::labelsCutVar}, "V0 daughter selections"};
  Configurable<std::vector<double>> binsPtV0{"binsPtV0", std::vector<double>{hf_cuts_v0_daughter::vecBinsPt}, "pT bin limits for V0 daughter cuts"};

  // Configurables for ME
  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<int> numberEventsToSkip{"numberEventsToSkip", -1, "Number of events to Skip in ME process"};

  SliceCache cache;

  using HfRed3PrNoTrksWithMl = soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl>;
  using HfRed2PrNoTrksWithMl = soa::Join<aod::HfRed2PrNoTrks, aod::HfRed2ProngsMl>;

  Preslice<aod::HfRedVzeros> candsV0PerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRedTrkNoParams> candsTrackPerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRed3PrNoTrks> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRed2PrNoTrks> candsD0PerCollision = hf_track_index_reduced::hfRedCollisionId;
  Preslice<HfRed3PrNoTrksWithMl> candsDPerCollisionWithMl = hf_track_index_reduced::hfRedCollisionId;
  Preslice<HfRed2PrNoTrksWithMl> candsD0PerCollisionWithMl = hf_track_index_reduced::hfRedCollisionId;

  // Partition of V0 candidates based on v0Type
  Partition<aod::HfRedVzeros> candidatesK0s = aod::hf_reso_v0::v0Type == (uint8_t)1 || aod::hf_reso_v0::v0Type == (uint8_t)3 || aod::hf_reso_v0::v0Type == (uint8_t)5;
  Partition<aod::HfRedVzeros> candidatesLambda = aod::hf_reso_v0::v0Type == (uint8_t)2 || aod::hf_reso_v0::v0Type == (uint8_t)4;

  ConfigurableAxis multPoolBins{"multPoolBins", {VARIABLE_WIDTH, 0., 45., 60., 75., 95, 250}, "event multiplicity pools (PV contributors for now)"};
  ConfigurableAxis zPoolBins{"zPoolBins", {VARIABLE_WIDTH, -10.0, -4, -1, 1, 4, 10.0}, "z vertex position pools"};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // check that only one process function is enabled
    std::array<bool, 16> doprocess{doprocessDs2StarToDplusK0s, doprocessDs2StarToDplusK0sWithMl, doprocessDs1ToDstarK0s, doprocessDs1ToDstarK0sWithMl, doprocessDs1ToDstarK0sMixedEvent, doprocessDs1ToDstarK0sMixedEventWithMl, doprocessDs2StarToDplusK0sMixedEventWithMl,
                                   doprocessXcToDplusLambda, doprocessXcToDplusLambdaWithMl, doprocessLambdaDminus, doprocessLambdaDminusWithMl, doprocessDstarTrack, doprocessDstarTrackWithMl, doprocessD0Track, doprocessD0TrackWithMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }
    // histograms
    const AxisSpec axisPt{(std::vector<double>)vecBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisMassDsj{400, 0.49f, 0.89f, ""};
    registry.add("hMassDs1", "Ds1 candidates;m_{Ds1} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisMassDsj, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDs2Star", "Ds^{*}2 candidates; m_Ds^{*}2 (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {axisMassDsj, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassXcRes", "XcRes candidates; m_XcRes (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{300, 1.1, 1.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassLambdaDminus", "LambdaDminus candidates; m_LambdaDminus (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{300, 1.1, 1.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDstarTrack", "DstarTrack candidates; m_DstarTrack (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 0.9, 1.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0Track", "D0Track candidates; m_D0Track (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 0.8, 1.3}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0BarTrack", "D0Track candidates; m_D0Track (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 0.8, 1.3}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    if (doprocessDs1ToDstarK0sMixedEvent) {
      registry.add("hNPvContCorr", "Collision number of PV contributors ; N contrib ; N contrib", {HistType::kTH2F, {{100, 0, 250}, {100, 0, 250}}});
      registry.add("hZvertCorr", "Collision Z Vtx ; z PV [cm] ; z PV [cm]", {HistType::kTH2F, {{120, -12., 12.}, {120, -12., 12.}}});
    }

    if (activateQA) {
      constexpr int kNBinsSelections = Selections::NSelSteps;
      std::string labels[kNBinsSelections];
      labels[Selections::NoSel] = "No selection";
      labels[Selections::DSel] = "D Candidates Selection";
      labels[Selections::V0Sel] = "D & V0 candidate Selection";
      labels[Selections::TrackSel] = "D & Track candidate Selection";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections", {HistType::kTH1F, {axisSelections}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH1>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }
  }

  /// Basic selection of D candidates
  /// \param candD is the reduced D meson candidate
  /// \return true if selections are passed
  template <DecayChannel channel, typename DRedTable>
  bool isDSelected(DRedTable const& candD)
  {
    float invMassD{0.};
    float ptD = candD.pt();
    int ptBin = findBin(binsPtD, ptD);
    if (ptBin == -1) {
      return false;
    }
    if constexpr (channel == DecayChannel::Ds2StarToDplusK0s || channel == DecayChannel::XcToDplusLambda || channel == DecayChannel::LambdaDminus) {
      invMassD = candD.invMassDplus();
    } else if constexpr (channel == DecayChannel::Ds1ToDstarK0s || channel == DecayChannel::DstarTrack) {
      if (candD.dType() > 0)
        invMassD = candD.invMassDstar() - candD.invMassD0();
      else
        invMassD = candD.invMassAntiDstar() - candD.invMassD0Bar();
    }
    // invariant mass selection
    if (!keepSideBands) {
      if constexpr (channel == DecayChannel::D0Track) {
        if ((candD.invMassD0() < cutsDDaughter->get(ptBin, "invMassSignalLow") || candD.invMassD0() > cutsDDaughter->get(ptBin, "invMassSignalHigh")) &&
            (candD.invMassD0Bar() < cutsDDaughter->get(ptBin, "invMassSignalLow") || candD.invMassD0Bar() > cutsDDaughter->get(ptBin, "invMassSignalHigh"))) {
          return false;
        }
      } else {
        if (invMassD < cutsDDaughter->get(ptBin, "invMassSignalLow") || invMassD > cutsDDaughter->get(ptBin, "invMassSignalHigh")) {
          return false;
        }
      }
    } else {
      if ((invMassD < cutsDDaughter->get(ptBin, "invMassLeftSBLow")) ||
          (invMassD > cutsDDaughter->get(ptBin, "invMassLeftSBHigh") && invMassD < cutsDDaughter->get(ptBin, "invMassSignalLow")) ||
          (invMassD > cutsDDaughter->get(ptBin, "invMassSignalHigh") && invMassD < cutsDDaughter->get(ptBin, "invMassRightSBLow")) ||
          (invMassD > cutsDDaughter->get(ptBin, "invMassRightSBHigh"))) {
        return false;
      }
    }
    return true;
  }

  /// Basic selection of V0 and track candidates
  /// \param candV0 is the reduced V0 candidate
  /// \param candD is the reduced D meson candidate
  /// \return true if selections are passed
  template <DecayChannel channel, typename DRedTable, typename V0RedTable>
  bool isV0Selected(V0RedTable const& candV0, DRedTable const& candD)
  {
    float massV0{0.};
    float invMassV0{0.};
    float ptV0 = candV0.pt();
    int ptBin = findBin(binsPtV0, ptV0);
    if (ptBin == -1) {
      return false;
    }
    if (channel == DecayChannel::Ds2StarToDplusK0s || channel == DecayChannel::Ds1ToDstarK0s) {
      massV0 = MassK0Short;
      invMassV0 = candV0.invMassK0s();
    } else if (channel == DecayChannel::XcToDplusLambda || channel == DecayChannel::LambdaDminus) {
      massV0 = MassLambda;
      int wsFact{1};
      if (channel == DecayChannel::LambdaDminus)
        wsFact = -1;
      uint8_t targetV0Type{0};
      if (wsFact * candD.dType() > 0) {
        invMassV0 = candV0.invMassLambda();
        targetV0Type = V0Type::Lambda;
      } else {
        invMassV0 = candV0.invMassAntiLambda();
        targetV0Type = V0Type::AntiLambda;
      }
      // check skimming cuts
      if (!TESTBIT(candV0.v0Type(), targetV0Type)) {
        return false;
      }
    }
    // selection on V0 candidate mass
    if ((invMassV0 - massV0) > cutsV0Daughter->get(ptBin, "invMassLow") && (massV0 - invMassV0) < cutsV0Daughter->get(ptBin, "invMassLow")) {
      return false;
    }
    // selection on kinematics and topology
    if (candV0.dca() > cutsV0Daughter->get(ptBin, "dcaMax") || candV0.cpa() < cutsV0Daughter->get(ptBin, "cpaMin") || candV0.v0Radius() < cutsV0Daughter->get(ptBin, "radiusMin")) {
      return false;
    }
    return true;
  }

  template <bool fillMl, DecayChannel channel, typename Coll, typename DRedTable, typename V0TrRedTable>
  void runCandidateCreation(Coll const& collision,
                            DRedTable const& candsD,
                            V0TrRedTable const& candsV0Tr)
  {
    // loop on D candidates
    for (const auto& candD : candsD) {
      // selection of D candidates
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1);
      }
      if (!isDSelected<channel>(candD)) {
        continue;
      }
      if (activateQA) {
        registry.fill(HIST("hSelections"), 1 + Selections::DSel);
      }
      float invMassD{0.};
      float invMassD0{0.};
      std::array<std::array<float, 3>, 3> pVectorCharmProngs;
      if constexpr (channel != DecayChannel::D0Track) {
        if (std::abs(candD.dType()) == DType::Dplus)
          invMassD = candD.invMassDplus();
        if (candD.dType() == DType::Dstar) {
          invMassD = candD.invMassDstar();
          invMassD0 = candD.invMassD0();
        }
        if (candD.dType() == (-1) * DType::Dstar) {
          invMassD = candD.invMassAntiDstar();
          invMassD0 = candD.invMassD0Bar();
        }
        pVectorCharmProngs = {candD.pVectorProng0(), candD.pVectorProng1(), candD.pVectorProng2()};
      } else {
        pVectorCharmProngs = {candD.pVectorProng0(), candD.pVectorProng1(), {0.}};
      }
      std::array<float, 3> pVecD = {candD.px(), candD.py(), candD.pz()};

      // loop on V0 or track candidates
      bool alreadyCounted{false};
      for (const auto& candV0Tr : candsV0Tr) {
        if (rejectDV0PairsWithCommonDaughter) {
          if constexpr (channel == DecayChannel::D0Track) {
            const std::array<int, 2> dDaughtersIDs = {candD.prong0Id(), candD.prong1Id()};
            if (std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.globalIndex()) != dDaughtersIDs.end()) {
              continue;
            }
          } else {
            const std::array<int, 3> dDaughtersIDs = {candD.prong0Id(), candD.prong1Id(), candD.prong2Id()};
            if constexpr (channel == DecayChannel::DstarTrack) {
              if (std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.globalIndex()) != dDaughtersIDs.end()) {
                continue;
              }
            } else {
              if (std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.prong0Id()) != dDaughtersIDs.end() || std::find(dDaughtersIDs.begin(), dDaughtersIDs.end(), candV0Tr.prong1Id()) != dDaughtersIDs.end()) {
                continue;
              }
            }
          }
        }
        if constexpr ((channel != DecayChannel::DstarTrack) && (channel != DecayChannel::D0Track)) {
          if (!isV0Selected<channel>(candV0Tr, candD)) {
            continue;
          }
          if (activateQA && !alreadyCounted) {
            registry.fill(HIST("hSelections"), 1 + Selections::V0Sel);
            alreadyCounted = true;
          }
        }

        float invMassReso{0.};
        float invMassV0{0.};
        std::array<float, 3> pVecV0Tr = {candV0Tr.px(), candV0Tr.py(), candV0Tr.pz()};
        float ptReso = RecoDecay::pt(RecoDecay::sumOfVec(pVecV0Tr, pVecD));

        if constexpr (channel == DecayChannel::DstarTrack) {
          if (candD.dType() > 0) {
            invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton});
          } else {
            invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassProton});
          }
          registry.fill(HIST("hMassDstarTrack"), invMassReso - invMassD, ptReso);
        } else if constexpr (channel == DecayChannel::D0Track) {
          if (TESTBIT(candD.selFlagD0(), D0SelectionType::SelectedD0)) {
            invMassD = candD.invMassD0();
            invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassProton});
            registry.fill(HIST("hMassD0Track"), invMassReso - invMassD, ptReso);
          }
          if (TESTBIT(candD.selFlagD0(), D0SelectionType::SelectedD0Bar)) {
            invMassD = candD.invMassD0Bar();
            invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassProton});
            registry.fill(HIST("hMassD0BarTrack"), invMassReso - invMassD, ptReso);
          }
        } else {
          switch (channel) {
            case DecayChannel::Ds1ToDstarK0s:
              invMassV0 = candV0Tr.invMassK0s();
              if (candD.dType() > 0) {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              } else {
                invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[1], pVectorCharmProngs[0], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              }
              registry.fill(HIST("hMassDs1"), invMassReso, ptReso);
              break;
            case DecayChannel::Ds2StarToDplusK0s:
              invMassV0 = candV0Tr.invMassK0s();
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0Short}) - invMassD;
              registry.fill(HIST("hMassDs2Star"), invMassReso, ptReso);
              break;
            case DecayChannel::XcToDplusLambda:
              if (candD.dType() > 0) {
                invMassV0 = candV0Tr.invMassLambda();
              } else {
                invMassV0 = candV0Tr.invMassAntiLambda();
              }
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda}) - invMassD;
              registry.fill(HIST("hMassXcRes"), invMassReso, ptReso);
              break;
            case DecayChannel::LambdaDminus:
              if (candD.dType() < 0) {
                invMassV0 = candV0Tr.invMassLambda();
              } else {
                invMassV0 = candV0Tr.invMassAntiLambda();
              }
              invMassReso = RecoDecay::m(std::array{pVectorCharmProngs[0], pVectorCharmProngs[1], pVectorCharmProngs[2], pVecV0Tr}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassLambda}) - invMassD;
              registry.fill(HIST("hMassLambdaDminus"), invMassReso, ptReso);
              break;
            default:
              break;
          }
        }
        // Filling Output table
        if constexpr (channel == DecayChannel::DstarTrack) {
          rowCandidateResoTrack(pVecD[0], pVecD[1], pVecD[2],
                                candV0Tr.px(), candV0Tr.py(), candV0Tr.pz(),
                                invMassReso,
                                invMassD - invMassD0);
        } else if constexpr (channel == DecayChannel::D0Track) {
          rowCandidateResoTrack(pVecD[0], pVecD[1], pVecD[2],
                                candV0Tr.px(), candV0Tr.py(), candV0Tr.pz(),
                                invMassReso,
                                0);
        } else {
          rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD,
                           invMassV0,
                           candV0Tr.cpa(),
                           candV0Tr.dca(),
                           candV0Tr.v0Radius(),
                           invMassD0);
          rowCandidateResoIndices(collision.globalIndex(),
                                  candD.globalIndex(),
                                  candV0Tr.globalIndex());
        }
        if constexpr (fillMl) {
          mlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
        }
      }
    }
  } // main function
  // Process data with Mixed Event
  /// \tparam fillMl is a flag to Fill ML scores if present
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param DRedTable is the D bachelors table
  /// \param V0TrRedTable is the V0/Track bachelors table
  template <bool fillMl, DecayChannel channel, typename Coll, typename DRedTable, typename V0TrRedTable>
  void runCandidateCreationMixedEvent(Coll const& collisions,
                                      DRedTable const& candsD,
                                      V0TrRedTable const& candsV0Tr)
  {
    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
    BinningType corrBinning{{zPoolBins, multPoolBins}, true};
    auto bachTuple = std::make_tuple(candsD, candsV0Tr);
    Pair<Coll, DRedTable, V0TrRedTable, BinningType> pairs{corrBinning, numberEventsMixed, numberEventsToSkip, collisions, bachTuple, &cache};
    for (const auto& [collision1, bachDs, collision2, bachV0Trs] : pairs) {
      registry.fill(HIST("hNPvContCorr"), collision1.numContrib(), collision2.numContrib());
      registry.fill(HIST("hZvertCorr"), collision1.posZ(), collision2.posZ());
      for (const auto& [bachD, bachV0Tr] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(bachDs, bachV0Trs))) {
        // Apply analysis selections on D and V0 bachelors
        if (!isDSelected<channel>(bachD) || !isV0Selected<channel>(bachV0Tr, bachD)) {
          continue;
        }
        // Retrieve D and V0 informations
        float invMassD{0.};
        float invMassD0{0.};
        if (std::abs(bachD.dType()) == DType::Dplus) {
          invMassD = bachD.invMassDplus();
        }
        if (bachD.dType() == DType::Dstar) {
          invMassD = bachD.invMassDstar();
          invMassD0 = bachD.invMassD0();
        }
        if (bachD.dType() == (-1) * DType::Dstar) {
          invMassD = bachD.invMassAntiDstar();
          invMassD0 = bachD.invMassD0Bar();
        }
        std::array<float, 3> pVecD = {bachD.px(), bachD.py(), bachD.pz()};
        float invMassReso{0.};
        float invMassV0{0.};
        std::array<float, 3> pVecV0Tr = {bachV0Tr.px(), bachV0Tr.py(), bachV0Tr.pz()};
        float ptReso = RecoDecay::pt(RecoDecay::sumOfVec(pVecV0Tr, pVecD));
        switch (channel) {
          case DecayChannel::Ds1ToDstarK0s:
            invMassV0 = bachV0Tr.invMassK0s();
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDStar, MassK0Short});
            registry.fill(HIST("hMassDs1"), invMassReso, ptReso);
            break;
          case DecayChannel::Ds2StarToDplusK0s:
            invMassV0 = bachV0Tr.invMassK0s();
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassK0Short});
            registry.fill(HIST("hMassDs2Star"), invMassReso, ptReso);
            break;
          case DecayChannel::XcToDplusLambda:
            if (bachD.dType() > 0) {
              invMassV0 = bachV0Tr.invMassLambda();
            } else {
              invMassV0 = bachV0Tr.invMassAntiLambda();
            }
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassLambda});
            registry.fill(HIST("hMassXcRes"), invMassReso, ptReso);
            break;
          case DecayChannel::LambdaDminus:
            if (bachD.dType() < 0) {
              invMassV0 = bachV0Tr.invMassLambda();
            } else {
              invMassV0 = bachV0Tr.invMassAntiLambda();
            }
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{MassDPlus, MassLambda});
            registry.fill(HIST("hMassLambdaDminus"), invMassReso, ptReso);
            break;
          default:
            break;
        }
        // Fill output table
        rowCandidateReso(pVecD[0], pVecD[1], pVecD[2],
                         pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                         invMassReso,
                         invMassD,
                         invMassV0,
                         bachV0Tr.cpa(),
                         bachV0Tr.dca(),
                         bachV0Tr.v0Radius(),
                         invMassD0);
        rowCandidateResoIndices(collision1.globalIndex(),
                                bachD.globalIndex(),
                                bachV0Tr.globalIndex());
        if constexpr (fillMl) {
          mlScores(bachD.mlScoreBkgMassHypo0(), bachD.mlScorePromptMassHypo0(), bachD.mlScoreNonpromptMassHypo0());
        }
      }
    }
  } // runCandidateCreationMixedEvent

  // List of Process Functions
  void processDs2StarToDplusK0s(aod::HfRedCollisions const& collisions,
                                aod::HfRed3PrNoTrks const& candsD,
                                aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::Ds2StarToDplusK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0s, "Process Ds2* candidates without ML info", true);

  void processDs2StarToDplusK0sWithMl(aod::HfRedCollisions const& collisions,
                                      soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& candsD,
                                      aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollisionWithMl, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::Ds2StarToDplusK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0sWithMl, "Process Ds2* candidates with Ml info", false);

  void processDs2StarToDplusK0sMixedEvent(aod::HfRedCollisions const& collisions,
                                          aod::HfRed3PrNoTrks const& candsD,
                                          aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<false, DecayChannel::Ds2StarToDplusK0s>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0sMixedEvent, "Process Ds2Star mixed Event without ML", false);

  void processDs2StarToDplusK0sMixedEventWithMl(aod::HfRedCollisions const& collisions,
                                                HfRed3PrNoTrksWithMl const& candsD,
                                                aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<true, DecayChannel::Ds2StarToDplusK0s>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0sMixedEventWithMl, "Process Ds2Star mixed Event with ML", false);

  void processDs1ToDstarK0s(aod::HfRedCollisions const& collisions,
                            aod::HfRed3PrNoTrks const& candsD,
                            aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::Ds1ToDstarK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0s, "Process Ds1 candidates without Ml info", false);

  void processDs1ToDstarK0sWithMl(aod::HfRedCollisions const& collisions,
                                  HfRed3PrNoTrksWithMl const& candsD,
                                  aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollisionWithMl, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::Ds1ToDstarK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0sWithMl, "Process Ds1 candidates with Ml info", false);

  void processDs1ToDstarK0sMixedEvent(aod::HfRedCollisions const& collisions,
                                      aod::HfRed3PrNoTrks const& candsD,
                                      aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<false, DecayChannel::Ds1ToDstarK0s>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0sMixedEvent, "Process Ds1 mixed Event without ML", false);

  void processDs1ToDstarK0sMixedEventWithMl(aod::HfRedCollisions const& collisions,
                                            HfRed3PrNoTrksWithMl const& candsD,
                                            aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreationMixedEvent<true, DecayChannel::Ds1ToDstarK0s>(collisions, candsD, candsV0);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0sMixedEventWithMl, "Process Ds1 mixed Event with ML", false);

  void processXcToDplusLambda(aod::HfRedCollisions const& collisions,
                              aod::HfRed3PrNoTrks const& candsD,
                              aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto lambdaThisColl = candidatesLambda.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::XcToDplusLambda>(collision, candsDThisColl, lambdaThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processXcToDplusLambda, "Process Xc candidates without Ml info", false);

  void processXcToDplusLambdaWithMl(aod::HfRedCollisions const& collisions,
                                    HfRed3PrNoTrksWithMl const& candsD,
                                    aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollisionWithMl, thisCollId);
      auto lambdaThisColl = candidatesLambda.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::XcToDplusLambda>(collision, candsDThisColl, lambdaThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processXcToDplusLambdaWithMl, "Process Xc candidates with Ml info", false);

  void processLambdaDminus(aod::HfRedCollisions const& collisions,
                           aod::HfRed3PrNoTrks const& candsD,
                           aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto lambdaThisColl = candidatesLambda.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::LambdaDminus>(collision, candsDThisColl, lambdaThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processLambdaDminus, "Process LambdaDminus candidates without Ml info", false);

  void processLambdaDminusWithMl(aod::HfRedCollisions const& collisions,
                                 HfRed3PrNoTrksWithMl const& candsD,
                                 aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollisionWithMl, thisCollId);
      auto lambdaThisColl = candidatesLambda.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::LambdaDminus>(collision, candsDThisColl, lambdaThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processLambdaDminusWithMl, "Process LambdaDminus candidates with Ml info", false);
  void processDstarTrack(aod::HfRedCollisions const& collisions,
                         aod::HfRed3PrNoTrks const& candsD,
                         aod::HfRedTrkNoParams const& candidatesTrack)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto trackThisColl = candidatesTrack.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::DstarTrack>(collision, candsDThisColl, trackThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarTrack, "Process DStar candidates without Ml info", false);

  void processDstarTrackWithMl(aod::HfRedCollisions const& collisions,
                               HfRed3PrNoTrksWithMl const& candsD,
                               aod::HfRedTrkNoParams const& candidatesTrack)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollisionWithMl, thisCollId);
      auto trackThisColl = candidatesTrack.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::DstarTrack>(collision, candsDThisColl, trackThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarTrackWithMl, "Process DStar candidates with Ml info", false);

  void processD0Track(aod::HfRedCollisions const& collisions,
                      aod::HfRed2PrNoTrks const& candsD,
                      aod::HfRedTrkNoParams const& candidatesTrack)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsD0PerCollision, thisCollId);
      auto trackThisColl = candidatesTrack.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<false, DecayChannel::D0Track>(collision, candsDThisColl, trackThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processD0Track, "Process D0 candidates without Ml info", false);

  void processD0TrackWithMl(aod::HfRedCollisions const& collisions,
                            HfRed2PrNoTrksWithMl const& candsD,
                            aod::HfRedTrkNoParams const& candidatesTrack)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsD0PerCollisionWithMl, thisCollId);
      auto trackThisColl = candidatesTrack.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::D0Track>(collision, candsDThisColl, trackThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processD0TrackWithMl, "Process D0 candidates with Ml info", false);

}; // struct HfCandidateCreatorCharmResoReduced

struct HfCandidateCreatorCharmResoReducedExpressions {

  Produces<aod::HfMcRecRedResos> rowResoMcRec;

  using CandResoWithIndices = soa::Join<aod::HfCandCharmReso, aod::HfResoIndices>;

  // Configurable axis
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisInvMassReso{"axisInvMassReso", {400, 0.49f, 0.89f}, "inv. mass (DV_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng0{"axisInvMassProng0", {200, 0.14, 0.17}, "inv. mass (D) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng1{"axisInvMassProng1", {200, 0.47, 0.53}, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassD0{"axisInvMassD0", {200, 1.65, 2.05}, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisDebug{"axisDebug", {16, -0.5, 15.5}, "MC debug flag"};
  ConfigurableAxis axisOrigin{"axisOrigin", {3, -0.5, 2.5}, "MC origin flag"};
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hMassMcMatched", "Reso MC candidates Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcMatchedIncomplete", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcUnmatched", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcNoEntry", "Reso MC candidates w.o. entry in MC Reco table;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisInvMassReso, axisPt}});
    registry.add("hMassMcMatchedVsBach0Mass", "Reso MC candidates Matched with generate particle;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassProng0}});
    registry.add("hMassMcUnmatchedVsBach0Mass", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassProng0}});
    registry.add("hMassMcMatchedVsBach1Mass", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassProng1}});
    registry.add("hMassMcUnmatchedVsBach1Mass", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassProng1}});
    registry.add("hMassMcMatchedVsD0Mass", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassD0}});
    registry.add("hMassMcUnmatchedVsD0Mass", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {axisInvMassReso, axisInvMassD0}});
    registry.add("hMassMcUnmatchedVsDebug", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2});debug flag", {HistType::kTH2F, {axisInvMassReso, axisDebug}});
    registry.add("hSparseUnmatchedDebug", "THn for debug of MC matching and Correlated BKG study", HistType::kTHnSparseF, {axisInvMassReso, axisPt, axisInvMassProng0, axisInvMassProng1, axisInvMassD0, axisDebug, axisOrigin});
  }

  /// Fill candidate information at MC reconstruction level
  /// \param rowsDV0McRec MC reco information on DPi pairs
  /// \param candsReso prong global indices of B0 candidates
  template <typename McRec>
  void fillResoMcRec(McRec const& rowsDV0McRec, CandResoWithIndices const& candsReso)
  {
    for (const auto& candReso : candsReso) {
      bool filledMcInfo{false};
      for (const auto& rowDV0McRec : rowsDV0McRec) {
        if ((rowDV0McRec.prong0Id() != candReso.prong0Id()) || (rowDV0McRec.prong1Id() != candReso.prong1Id())) {
          continue;
        }
        rowResoMcRec(rowDV0McRec.flagMcMatchRec(), rowDV0McRec.debugMcRec(), rowDV0McRec.origin(), rowDV0McRec.ptMother());
        filledMcInfo = true;
        if (std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0s || std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi ||
            std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0sOneMu || std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds2StarToDplusK0sOneMu) {
          registry.fill(HIST("hMassMcMatched"), candReso.invMass(), candReso.pt());
          registry.fill(HIST("hMassMcMatchedVsBach0Mass"), candReso.invMass(), candReso.invMassProng0() - candReso.invMassD0());
          registry.fill(HIST("hMassMcMatchedVsBach1Mass"), candReso.invMass(), candReso.invMassProng1());
          registry.fill(HIST("hMassMcMatchedVsD0Mass"), candReso.invMass(), candReso.invMassD0());

        } else if (std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds1ToDStarK0ToD0NoPiK0sPart || std::abs(rowDV0McRec.flagMcMatchRec()) == DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s) {
          registry.fill(HIST("hMassMcMatchedIncomplete"), candReso.invMass(), candReso.pt());
        } else {
          registry.fill(HIST("hMassMcUnmatched"), candReso.invMass(), candReso.pt());
          registry.fill(HIST("hMassMcUnmatchedVsBach0Mass"), candReso.invMass(), candReso.invMassProng0() - candReso.invMassD0());
          registry.fill(HIST("hMassMcUnmatchedVsBach1Mass"), candReso.invMass(), candReso.invMassProng1());
          registry.fill(HIST("hMassMcUnmatchedVsD0Mass"), candReso.invMass(), candReso.invMassD0());
          registry.fill(HIST("hMassMcUnmatchedVsDebug"), candReso.invMass(), rowDV0McRec.debugMcRec());
          registry.fill(HIST("hSparseUnmatchedDebug"), candReso.invMass(), candReso.pt(), candReso.invMassProng0() - candReso.invMassD0(), candReso.invMassProng1(), candReso.invMassD0(), rowDV0McRec.debugMcRec(), rowDV0McRec.origin());
        }

        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the D-Pi creator
        rowResoMcRec(0, -1, -1, -1.f);
        registry.fill(HIST("hMassMcNoEntry"), candReso.invMass(), candReso.pt());
      }
    }
  }

  void processMc(aod::HfMcRecRedDV0s const& rowsDV0McRec, CandResoWithIndices const& candsReso)
  {
    fillResoMcRec(rowsDV0McRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processMc, "Process MC", false);

  void processDummy(CandResoWithIndices const&) {}
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processDummy, "Process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorCharmResoReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorCharmResoReducedExpressions>(cfgc)};
}
