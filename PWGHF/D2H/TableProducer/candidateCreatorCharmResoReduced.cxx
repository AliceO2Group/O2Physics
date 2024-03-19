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

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum Selections : uint8_t {
  NoSel = 0,
  DSel,
  V0Sel,
  NSelSteps
};
enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  Ds2StarToDplusK0s,
  XcToDplusLambda
};
enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};
const int nBins = 7;
constexpr double binsPt[nBins + 1] = {
  1.,
  2.,
  4.,
  6.,
  8.,
  12.,
  24.,
  50.};
auto vecBins = std::vector<double>{binsPt, binsPt + nBins + 1};

struct HfCandidateCreatorCharmResoReduced {
  // Produces: Tables with resonance info
  Produces<aod::HfCandCharmReso> rowCandidateReso;

  // Configurables
  Configurable<double> invMassWindowD{"invMassWindowD", 0.5, "invariant-mass window for D candidates (GeV/c2)"};
  Configurable<double> invMassWindowV0{"invMassWindowV0", 0.5, "invariant-mass window for V0 candidates (GeV/c2)"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};
  // Hist Axis
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{vecBins}, "pT bin limits"};

  // Partition of V0 candidates based on v0Type
  Partition<aod::HfRedVzeros> candidatesK0s = aod::hf_reso_cand_reduced::v0Type == (uint8_t)1 || aod::hf_reso_cand_reduced::v0Type == (uint8_t)3 || aod::hf_reso_cand_reduced::v0Type == (uint8_t)5;
  Partition<aod::HfRedVzeros> candidatesLambda = aod::hf_reso_cand_reduced::v0Type == (uint8_t)2 || aod::hf_reso_cand_reduced::v0Type == (uint8_t)4;

  // Useful constants
  double massK0{0.};
  double massLambda{0.};
  double massDplus{0.};
  double massDstar{0.};

  // Histogram registry: if task make it with a THNsparse with all variables you want to save
  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    for (const auto& value : vecBins) {
      LOGF(info, "bin limit %f", value);
    }
    const AxisSpec axisPt{(std::vector<double>)vecBins, "#it{p}_{T} (GeV/#it{c})"};
    // histograms
    registry.add("hMassDs1", "Ds1 candidates;m_{Ds1} - m_{D^{*}} (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{100, 2.4, 2.7}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassDs2Star", "Ds^{*}2 candidates; Ds^{*}2 - m_{D^{#plus}} (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 2.4, 2.7}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassXcRes", "XcRes candidates; XcRes - m_{D^{#plus}} (GeV/#it{c}^{2}) ;entries", {HistType::kTH2F, {{100, 2.9, 3.3}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    if (activateQA) {
      constexpr int kNBinsSelections = Selections::NSelSteps;
      std::string labels[kNBinsSelections];
      labels[Selections::NoSel] = "No selection";
      labels[Selections::DSel] = "D Candidates Selection";
      labels[Selections::V0Sel] = "D & V0 candidate Selection";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections", {HistType::kTH1F, {axisSelections}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH1>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }

    massK0 = o2::constants::physics::MassK0Short;
    massLambda = o2::constants::physics::MassLambda;
    massDplus = o2::constants::physics::MassDPlus;
    massDstar = o2::constants::physics::MassDStar;
  }

  template <DecayChannel channel, typename DRedTable>
  bool isDSelected(DRedTable const& candD)
  {
    float massD{0.};
    // slection on D candidate mass
    if (channel == DecayChannel::Ds2StarToDplusK0s || channel == DecayChannel::XcToDplusLambda) {
      massD = massDplus;
    } else if (channel == DecayChannel::Ds1ToDstarK0s) {
      massD = massDstar;
    }
    if (std::fabs(candD.invMass() - massD) > invMassWindowD) {
      return false;
    }
    return true;
  }

  template <DecayChannel channel, typename DRedTable, typename V0RedTable>
  bool isV0Selected(V0RedTable const& candV0, DRedTable const& candD)
  {
    float massV0{0.};
    float invMassV0{0.};
    // slection on V0 candidate mass
    if (channel == DecayChannel::Ds2StarToDplusK0s || channel == DecayChannel::Ds1ToDstarK0s) {
      massV0 = massK0;
      invMassV0 = candV0.invMassK0s();
    } else if (channel == DecayChannel::XcToDplusLambda) {
      massV0 = massLambda;
      uint8_t targetV0Type{0};
      if (candD.dType() > 0) {
        invMassV0 = candV0.invMassLambda();
        targetV0Type = V0Type::Lambda;
      } else {
        invMassV0 = candV0.invMassAntiLambda();
        targetV0Type = V0Type::AntiLambda;
      }
      if (!TESTBIT(candV0.v0Type(), targetV0Type)) {
        return false;
      }
    }
    if (std::fabs(invMassV0 - massV0) > invMassWindowV0) {
      return false;
    }
    return true;
  }

  template <DecayChannel channel, typename Coll, typename DRedTable, typename V0RedTable>
  void runCandidateCreation(Coll const& collisions,
                            DRedTable const& candsD,
                            V0RedTable const& candsV0)
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
      float invMassD = candD.invMass();
      std::array<float, 3> pVecD = {candD.px(), candD.py(), candD.pz()};
      float ptD = RecoDecay::pt(pVecD);
      ;
      // loop on V0 candidates
      bool alreadyCounted{false};
      for (const auto& candV0 : candsV0) {
        if (!isV0Selected<channel>(candV0, candD)) {
          continue;
        }
        if (activateQA && !alreadyCounted) {
          registry.fill(HIST("hSelections"), 1 + Selections::V0Sel);
          alreadyCounted = true;
        }
        float invMassReso{0.};
        float invMassV0{0.};
        std::array<float, 3> pVecV0 = {candV0.px(), candV0.py(), candV0.pz()};
        float ptV0 = RecoDecay::pt(pVecV0);
        float ptReso = RecoDecay::pt(RecoDecay::sumOfVec(pVecV0, pVecD));
        switch (channel) {
          case DecayChannel::Ds1ToDstarK0s:
            invMassV0 = candV0.invMassK0s();
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0}, std::array{massDstar, massK0});
            registry.fill(HIST("hMassDs1"), invMassReso, ptReso);
            break;
          case DecayChannel::Ds2StarToDplusK0s:
            invMassV0 = candV0.invMassK0s();
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0}, std::array{massDplus, massK0});
            registry.fill(HIST("hMassDs2Star"), invMassReso, ptReso);
            break;
          case DecayChannel::XcToDplusLambda:
            if (candD.dType() > 0) {
              invMassV0 = candV0.invMassLambda();
            } else {
              invMassV0 = candV0.invMassAntiLambda();
            }
            invMassReso = RecoDecay::m(std::array{pVecD, pVecV0}, std::array{massDplus, massLambda});
            registry.fill(HIST("hMassXcRes"), invMassReso, ptReso);
            break;
          default:
            break;
        }
        // Filling Output table
        rowCandidateReso(collisions.globalIndex(),
                         invMassReso,
                         ptReso,
                         invMassD,
                         ptD,
                         invMassV0,
                         ptV0,
                         candV0.cpa(),
                         candV0.dca(),
                         candV0.radius());
      }
    }
  } // main function

  void processDs2StarToDplusK0s(aod::HfRedCollisions::iterator const& collision,
                                aod::HfRed3PrNoTrks const& candsD,
                                aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreation<DecayChannel::Ds2StarToDplusK0s>(collision, candsD, candidatesK0s);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs2StarToDplusK0s, "Process Ds2* candidates without MC info and without ML info", true);

  void processDs1ToDstarK0s(aod::HfRedCollisions::iterator const& collision,
                            aod::HfRed3PrNoTrks const& candsD,
                            aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreation<DecayChannel::Ds1ToDstarK0s>(collision, candsD, candidatesK0s);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processDs1ToDstarK0s, "Process Ds1 candidates without MC info and without ML info", false);

  void processXcToDplusLambda(aod::HfRedCollisions::iterator const& collision,
                              aod::HfRed3PrNoTrks const& candsD,
                              aod::HfRedVzeros const& candsV0)
  {
    runCandidateCreation<DecayChannel::XcToDplusLambda>(collision, candsD, candidatesLambda);
  }
  PROCESS_SWITCH(HfCandidateCreatorCharmResoReduced, processXcToDplusLambda, "Process Xc candidates without MC info and without ML info", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateCreatorCharmResoReduced>(cfgc)};
}
