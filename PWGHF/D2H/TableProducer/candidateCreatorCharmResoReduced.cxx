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

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "EventFiltering/PWGHF/HFFilterHelpers.h"

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Core/SelectorCutsRedDataFormat.h"
#include "PWGHF/Utils/utilsAnalysis.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
  DstarTrack
};

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};

enum DecayTypeMc : uint8_t {
  Ds1ToDStarK0ToD0PiK0s = 0,
  Ds2StarToDplusK0sToPiKaPiPiPi,
  Ds1ToDStarK0ToDPlusPi0K0s,
  Ds1ToDStarK0ToDPlusGammaK0s
};

const int nBinsPt = 7;
constexpr double binsPt[nBinsPt + 1] = {
  1.,
  2.,
  4.,
  6.,
  8.,
  12.,
  24.,
  1000.};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

struct HfCandidateCreatorCharmResoReduced {
  // Produces: Tables with resonance info
  Produces<aod::HfCandCharmReso> rowCandidateReso;
  Produces<aod::HfCandChaResTr> rowCandidateResoTrack;
  // Optional daughter ML scores table
  Produces<aod::HfCharmResoMLs> mlScores;

  // Configurables
  Configurable<bool> rejectDV0PairsWithCommonDaughter{"rejectDV0PairsWithCommonDaughter", true, "flag to reject the pairs that share a daughter track if not done in the derived data creation"};
  Configurable<bool> keepSideBands{"keepSideBands", false, "flag to keep events from D meson sidebands for backgorund estimation"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{vecBinsPt}, "Histogram pT bin limits"};
  // Daughters selection cuts
  Configurable<LabeledArray<double>> cutsD{"cutsDdaughter", {hf_cuts_d_daughter::cuts[0], hf_cuts_d_daughter::nBinsPt, hf_cuts_d_daughter::nCutVars, hf_cuts_d_daughter::labelsPt, hf_cuts_d_daughter::labelsCutVar}, "D daughter selections"};
  Configurable<std::vector<double>> binsPtD{"binsPtD", std::vector<double>{hf_cuts_d_daughter::vecBinsPt}, "pT bin limits for D daughter cuts"};
  Configurable<LabeledArray<double>> cutsV0{"cutsV0daughter", {hf_cuts_v0_daughter::cuts[0], hf_cuts_v0_daughter::nBinsPt, hf_cuts_v0_daughter::nCutVars, hf_cuts_v0_daughter::labelsPt, hf_cuts_v0_daughter::labelsCutVar}, "V0 daughter selections"};

  Configurable<std::vector<double>> binsPtV0{"binsPtV0", std::vector<double>{hf_cuts_v0_daughter::vecBinsPt}, "pT bin limits for V0 daughter cuts"};

  using reducedDWithMl = soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl>;

  // Partition of V0 candidates based on v0Type
  Partition<aod::HfRedVzeros> candidatesK0s = aod::hf_reso_v0::v0Type == (uint8_t)1 || aod::hf_reso_v0::v0Type == (uint8_t)3 || aod::hf_reso_v0::v0Type == (uint8_t)5;
  Partition<aod::HfRedVzeros> candidatesLambda = aod::hf_reso_v0::v0Type == (uint8_t)2 || aod::hf_reso_v0::v0Type == (uint8_t)4;

  Preslice<aod::HfRedVzeros> candsV0PerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRedTracks> candsTrackPerCollision = aod::hf_track_index_reduced::hfRedCollisionId;
  Preslice<aod::HfRed3PrNoTrks> candsDPerCollision = hf_track_index_reduced::hfRedCollisionId;

  // Useful constants
  double massK0{0.};
  double massLambda{0.};
  double massProton{0.};
  double massDplus{0.};
  double massDstar{0.};
  double massD0{0.};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // check that only one process function is enabled
    std::array<bool, 10> doprocess{doprocessDs2StarToDplusK0s, doprocessDs2StarToDplusK0sWithMl, doprocessDs1ToDstarK0s, doprocessDs1ToDstarK0sWithMl, doprocessXcToDplusLambda, doprocessXcToDplusLambdaWithMl, doprocessLambdaDminus, doprocessLambdaDminusWithMl, doprocessDstarTrack, doprocessDstarTrackWithMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }
    // histograms
    const AxisSpec axisPt{(std::vector<double>)vecBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    registry.add("hMassDs1", "Ds1 candidates;m_{Ds1} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2.5, 2.7}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDs2Star", "Ds^{*}2 candidates; m_Ds^{*}2 (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 2.4, 2.7}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassXcRes", "XcRes candidates; m_XcRes (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 2.9, 3.3}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassLambdaDminus", "LambdaDminus candidates; m_LambdaDminus (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 2.9, 3.3}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDstarTrack", "DstarTrack candidates; m_DstarTrack (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 0.9, 1.4}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

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
    // mass constants
    massK0 = o2::constants::physics::MassK0Short;
    massLambda = o2::constants::physics::MassLambda;
    massProton = o2::constants::physics::MassProton;
    massDplus = o2::constants::physics::MassDPlus;
    massDstar = o2::constants::physics::MassDStar;
    massD0 = o2::constants::physics::MassD0;
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
    // slection on D candidate mass
    if (channel == DecayChannel::Ds2StarToDplusK0s || channel == DecayChannel::XcToDplusLambda || channel == DecayChannel::LambdaDminus) {
      invMassD = candD.invMassDplus();
    } else if (channel == DecayChannel::Ds1ToDstarK0s || channel == DecayChannel::DstarTrack) {
      if (candD.dType() > 0)
        invMassD = candD.invMassDstar();
      else
        invMassD = candD.invMassAntiDstar();
    }
    // invariant mass selection
    if (!keepSideBands) {
      if (invMassD < cutsD->get(ptBin, "invMassSignalLow") || invMassD > cutsD->get(ptBin, "invMassSignalHigh")) {
        return false;
      }
    } else {
      if ((invMassD < cutsD->get(ptBin, "invMassLeftSBLow")) ||
          (invMassD > cutsD->get(ptBin, "invMassLeftSBHigh") && invMassD < cutsD->get(ptBin, "invMassSignalLow")) ||
          (invMassD > cutsD->get(ptBin, "invMassSignalHigh") && invMassD < cutsD->get(ptBin, "invMassRightSBLow")) ||
          (invMassD > cutsD->get(ptBin, "invMassRightSBHigh"))) {
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
      massV0 = massK0;
      invMassV0 = candV0.invMassK0s();
    } else if (channel == DecayChannel::XcToDplusLambda || channel == DecayChannel::LambdaDminus) {
      massV0 = massLambda;
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
    if ((invMassV0 - massV0) > cutsV0->get(ptBin, "invMassLow") && (massV0 - invMassV0) < cutsV0->get(ptBin, "invMassLow")) {
      return false;
    }
    // selection on kinematics and topology
    if (candV0.dca() > cutsV0->get(ptBin, "dcaMax") || candV0.cpa() < cutsV0->get(ptBin, "cpaMin") || candV0.v0Radius() < cutsV0->get(ptBin, "radiusMin")) {
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
      if (std::abs(candD.dType()) == 1)
        invMassD = candD.invMassDplus();
      if (candD.dType() == 2)
        invMassD = candD.invMassDstar();
      if (candD.dType() == -2)
        invMassD = candD.invMassAntiDstar();
      std::array<float, 3> pVecD = {candD.px(), candD.py(), candD.pz()};
      std::array<int, 3> dDaughtersIds = {candD.prong0Id(), candD.prong1Id(), candD.prong2Id()};

      // loop on V0 or track candidates
      bool alreadyCounted{false};
      for (const auto& candV0Tr : candsV0Tr) {
        if (rejectDV0PairsWithCommonDaughter) {
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
        if constexpr (channel != DecayChannel::DstarTrack) {
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
          invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{massDstar, massProton});
          registry.fill(HIST("hMassDstarTrack"), invMassReso, ptReso);
        } else {
          switch (channel) {
            case DecayChannel::Ds1ToDstarK0s:
              invMassV0 = candV0Tr.invMassK0s();
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{massDstar, massK0});
              registry.fill(HIST("hMassDs1"), invMassReso, ptReso);
              break;
            case DecayChannel::Ds2StarToDplusK0s:
              invMassV0 = candV0Tr.invMassK0s();
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{massDplus, massK0});
              registry.fill(HIST("hMassDs2Star"), invMassReso, ptReso);
              break;
            case DecayChannel::XcToDplusLambda:
              if (candD.dType() > 0) {
                invMassV0 = candV0Tr.invMassLambda();
              } else {
                invMassV0 = candV0Tr.invMassAntiLambda();
              }
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{massDplus, massLambda});
              registry.fill(HIST("hMassXcRes"), invMassReso, ptReso);
              break;
            case DecayChannel::LambdaDminus:
              if (candD.dType() < 0) {
                invMassV0 = candV0Tr.invMassLambda();
              } else {
                invMassV0 = candV0Tr.invMassAntiLambda();
              }
              invMassReso = RecoDecay::m(std::array{pVecD, pVecV0Tr}, std::array{massDplus, massLambda});
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
                                invMassD);
        } else {
          rowCandidateReso(collision.globalIndex(),
                           candD.globalIndex(),
                           candV0Tr.globalIndex(),
                           pVecD[0], pVecD[1], pVecD[2],
                           pVecV0Tr[0], pVecV0Tr[1], pVecV0Tr[2],
                           invMassReso,
                           invMassD,
                           invMassV0,
                           candV0Tr.cpa(),
                           candV0Tr.dca(),
                           candV0Tr.v0Radius());
        }
        if constexpr (fillMl) {
          mlScores(candD.mlScoreBkgMassHypo0(), candD.mlScorePromptMassHypo0(), candD.mlScoreNonpromptMassHypo0());
        }
      }
    }
  } // main function

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
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::Ds2StarToDplusK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0sWithMl, "Process Ds2* candidates with Ml info", false);

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
                                  reducedDWithMl const& candsD,
                                  aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto k0sThisColl = candidatesK0s.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::Ds1ToDstarK0s>(collision, candsDThisColl, k0sThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0sWithMl, "Process Ds1 candidates with Ml info", false);

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
                                    soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& candsD,
                                    aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
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
                                 soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& candsD,
                                 aod::HfRedVzeros const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto lambdaThisColl = candidatesLambda.sliceBy(candsV0PerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::LambdaDminus>(collision, candsDThisColl, lambdaThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processLambdaDminusWithMl, "Process LambdaDminus candidates with Ml info", false);
  void processDstarTrack(aod::HfRedCollisions const& collisions,
                         aod::HfRed3PrNoTrks const& candsD,
                         soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& candidatesTrack)
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
                               soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& candsD,
                               soa::Join<aod::HfRedTrackBases, aod::HfRedTracksCov> const& candidatesTrack)
  //                                    aod::HfRedTracks const&)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsD.sliceBy(candsDPerCollision, thisCollId);
      auto trackThisColl = candidatesTrack.sliceBy(candsTrackPerCollision, thisCollId);
      runCandidateCreation<true, DecayChannel::DstarTrack>(collision, candsDThisColl, trackThisColl);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDstarTrackWithMl, "Process DStar candidates with Ml info", false);

}; // struct HfCandidateCreatorCharmResoReduced

struct HfCandidateCreatorCharmResoReducedExpressions {

  Produces<aod::HfMcRecRedResos> rowResoMcRec;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    const AxisSpec axisPt{(std::vector<double>)vecBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    registry.add("hMassMcMatched", "Reso MC candidates Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2.5, 2.7}, {axisPt}}});
    registry.add("hMassMcMatchedIncomplete", "Reso MC candidates Matched with generate particle w. Invcomplete decay;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2.5, 2.7}, {axisPt}}});
    registry.add("hMassMcUnmatched", "Reso MC candidates NOT Matched with generate particle;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2.5, 2.7}, {axisPt}}});
    registry.add("hMassMcNoEntry", "Reso MC candidates w.o. entry in MC Reco table;m (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 2.5, 2.7}, {axisPt}}});
  }

  /// Fill candidate information at MC reconstruction level
  /// \param rowsDV0McRec MC reco information on DPi pairs
  /// \param candsReso prong global indices of B0 candidates
  template <typename McRec>
  void fillResoMcRec(McRec const& rowsDV0McRec, aod::HfCandCharmReso const& candsReso)
  {
    for (const auto& candReso : candsReso) {
      bool filledMcInfo{false};
      for (const auto& rowDV0McRec : rowsDV0McRec) {
        if ((rowDV0McRec.prong0Id() != candReso.prong0Id()) || (rowDV0McRec.prong1Id() != candReso.prong1Id())) {
          continue;
        }
        rowResoMcRec(rowDV0McRec.flagMcMatchRec(), rowDV0McRec.debugMcRec(), rowDV0McRec.origin(), rowDV0McRec.ptMother());
        filledMcInfo = true;
        if (TESTBIT(rowDV0McRec.flagMcMatchRec(), DecayTypeMc::Ds1ToDStarK0ToD0PiK0s) || TESTBIT(rowDV0McRec.flagMcMatchRec(), DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi)) {
          registry.fill(HIST("hMassMcMatched"), candReso.invMass(), candReso.pt());
        } else if (TESTBIT(rowDV0McRec.flagMcMatchRec(), DecayTypeMc::Ds1ToDStarK0ToDPlusGammaK0s) || TESTBIT(rowDV0McRec.flagMcMatchRec(), DecayTypeMc::Ds1ToDStarK0ToDPlusPi0K0s)) {
          registry.fill(HIST("hMassMcMatchedIncomplete"), candReso.invMass(), candReso.pt());
        } else {
          registry.fill(HIST("hMassMcUnmatched"), candReso.invMass(), candReso.pt());
        }

        break;
      }
      if (!filledMcInfo) { // protection to get same size tables in case something went wrong: we created a candidate that was not preselected in the D-Pi creator
        rowResoMcRec(0, -1, -1, -1.f);
        registry.fill(HIST("hMassMcNoEntry"), candReso.invMass(), candReso.pt());
      }
    }
  }

  void processMc(aod::HfMcRecRedDV0s const& rowsDV0McRec, aod::HfCandCharmReso const& candsReso)
  {
    fillResoMcRec(rowsDV0McRec, candsReso);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReducedExpressions, processMc, "Process MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorCharmResoReduced>(cfgc),
                      adaptAnalysisTask<HfCandidateCreatorCharmResoReducedExpressions>(cfgc)};
}
