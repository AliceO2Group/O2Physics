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

/// \file treeCreatorXicToPKPi.cxx
/// \brief Writer of XiC -> pKpi candidates in the form of flat tables to be stored in TTrees.
/// inspired from file treeCreatorLcToPKPi.cxx and to treeCreatorDplusToPiKPi.cxx

/// \author Himanshu Sharma <himanshu.sharma@cern.ch>, INFN Padova
/// \author Cristina Terrevoli <cristina.terrevoli@cern.ch>, INFN Bari

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{

DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(PProng2, pProng2, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(E, e, float);
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);
DECLARE_SOA_COLUMN(NSigTpcPr0, nSigTpcPr0, float);
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);
DECLARE_SOA_COLUMN(NSigTofPr0, nSigTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);
DECLARE_SOA_COLUMN(NSigTpcPr2, nSigTpcPr2, float);
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);
DECLARE_SOA_COLUMN(NSigTofPr2, nSigTofPr2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr0, nSigTpcTofPr0, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr2, nSigTpcTofPr2, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(FlagMc, flagMc, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t);
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int);

// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

DECLARE_SOA_TABLE(HfCandXicLites, "AOD", "HFCANDXICLITE",
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  full::NSigTpcPi0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofPr0,
                  full::NSigTpcKa1,
                  full::NSigTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofPr2,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofPr0,
                  full::NSigTpcTofKa1,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofPr2,
                  hf_sel_candidate_xic::IsSelXicToPKPi,
                  hf_sel_candidate_xic::IsSelXicToPiKP,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::Eta,
                  full::Phi,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped)

DECLARE_SOA_TABLE(HfCandXicFulls, "AOD", "HFCANDXICFULL",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength,
                  hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA,
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::ImpactParameterNormalised0,
                  full::PtProng0,
                  full::PProng0,
                  full::ImpactParameterNormalised1,
                  full::PtProng1,
                  full::PProng1,
                  full::ImpactParameterNormalised2,
                  full::PtProng2,
                  full::PProng2,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::PxProng2,
                  hf_cand::PyProng2,
                  hf_cand::PzProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand::ErrorImpactParameter2,
                  full::NSigTpcPi0,
                  full::NSigTpcPr0,
                  full::NSigTofPi0,
                  full::NSigTofPr0,
                  full::NSigTpcKa1,
                  full::NSigTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcPr2,
                  full::NSigTofPi2,
                  full::NSigTofPr2,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofPr0,
                  full::NSigTpcTofKa1,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofPr2,
                  hf_sel_candidate_xic::IsSelXicToPKPi,
                  hf_sel_candidate_xic::IsSelXicToPiKP,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Cpa,
                  full::CpaXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped);

DECLARE_SOA_TABLE(HfCandXicFullEvs, "AOD", "HFCANDXICFULLEV",
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXicFullPs, "AOD", "HFCANDXICFULLP",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXicToPKPi {
  Produces<o2::aod::HfCandXicFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXicFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXicFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandXicLites> rowCandidateLite;

  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection flag for Xic"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived  tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill  derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of   background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  using CandXicData = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi>>;
  using CandXicMcReco = soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>>;
  using CandXicMcGen = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic;
  Filter filterMcGenMatching = nabs(o2::aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi);

  Partition<CandXicData> selectedXicToPKPiCand = aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic;
  Partition<CandXicData> selectedXicToPiKPCand = aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic;

  Partition<CandXicMcReco> reconstructedCandSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi);
  Partition<CandXicMcReco> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi);

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,
      runNumber);
  }

  /// Fill table accounting for MC information and mass hypothesis
  /// \param doMc true to fill MC information
  /// \param massHypo mass hypothesis considered: 0 = PKPi, 1 = PiKP
  /// \param candidate is candidate
  template <bool DoMc = false, int MassHypo = 0, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    int candSwapped = 0;
    if constexpr (DoMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
      candSwapped = candidate.isCandidateSwapped();
    }

    float invMassXic = 0;
    int selStatusPKPi = candidate.isSelXicToPKPi();
    int selStatusPiKP = candidate.isSelXicToPiKP();

    if constexpr (MassHypo == 0) { // Xic->PKPi
      selStatusPiKP *= -1;
      invMassXic = HfHelper::invMassXicToPKPi(candidate);
    } else if constexpr (MassHypo == 1) { // Xic->PiKP
      selStatusPKPi *= -1;
      invMassXic = HfHelper::invMassXicToPiKP(candidate);
    }
    if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.nSigTpcPi0(),
        candidate.nSigTpcPr0(),
        candidate.nSigTofPi0(),
        candidate.nSigTofPr0(),
        candidate.nSigTpcKa1(),
        candidate.nSigTofKa1(),
        candidate.nSigTpcPi2(),
        candidate.nSigTpcPr2(),
        candidate.nSigTofPi2(),
        candidate.nSigTofPr2(),
        candidate.tpcTofNSigmaPi0(),
        candidate.tpcTofNSigmaPr0(),
        candidate.tpcTofNSigmaKa1(),
        candidate.tpcTofNSigmaPi2(),
        candidate.tpcTofNSigmaPr2(),
        selStatusPKPi,
        selStatusPiKP,
        invMassXic,
        candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.eta(),
        candidate.phi(),
        flagMc,
        originMc,
        candSwapped);

    } else {
      rowCandidateFull(
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.errorDecayLength(),
        candidate.errorDecayLengthXY(),
        candidate.chi2PCA(),
        candidate.rSecondaryVertex(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.impactParameterNormalised0(),
        candidate.ptProng0(),
        RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
        candidate.impactParameterNormalised1(),
        candidate.ptProng1(),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
        candidate.impactParameterNormalised2(),
        candidate.ptProng2(),
        RecoDecay::p(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
        candidate.pxProng0(),
        candidate.pyProng0(),
        candidate.pzProng0(),
        candidate.pxProng1(),
        candidate.pyProng1(),
        candidate.pzProng1(),
        candidate.pxProng2(),
        candidate.pyProng2(),
        candidate.pzProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.errorImpactParameter0(),
        candidate.errorImpactParameter1(),
        candidate.errorImpactParameter2(),
        candidate.nSigTpcPi0(),
        candidate.nSigTpcPr0(),
        candidate.nSigTofPi0(),
        candidate.nSigTofPr0(),
        candidate.nSigTpcKa1(),
        candidate.nSigTofKa1(),
        candidate.nSigTpcPi2(),
        candidate.nSigTpcPr2(),
        candidate.nSigTofPi2(),
        candidate.nSigTofPr2(),
        candidate.tpcTofNSigmaPi0(),
        candidate.tpcTofNSigmaPr0(),
        candidate.tpcTofNSigmaKa1(),
        candidate.tpcTofNSigmaPi2(),
        candidate.tpcTofNSigmaPr2(),
        selStatusPKPi,
        selStatusPiKP,
        invMassXic,
        candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        HfHelper::ctXic(candidate),
        candidate.eta(),
        candidate.phi(),
        HfHelper::yXic(candidate),
        HfHelper::eXic(candidate),
        flagMc,
        originMc,
        candSwapped);
    }
  }

  void processData(aod::Collisions const& collisions,
                   CandXicData const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }
    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(selectedXicToPKPiCand.size() + selectedXicToPiKPCand.size());
    } else {
      rowCandidateFull.reserve(selectedXicToPKPiCand.size() + selectedXicToPiKPCand.size());
    }

    for (const auto& candidate : selectedXicToPKPiCand) {
      if (downSampleBkgFactor < 1.) {
        float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable<false, 0>(candidate);
    }

    for (const auto& candidate : selectedXicToPiKPCand) {
      if (downSampleBkgFactor < 1.) {
        float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable<false, 1>(candidate);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorXicToPKPi, processData, "Process data tree writer", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 CandXicMcReco const&,
                 CandXicMcGen const& mcParticles)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillOnlySignal) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandSig.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandSig.size());
      }

      for (const auto& candidate : reconstructedCandSig) {
        if (candidate.isCandidateSwapped() == 0) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isCandidateSwapped() == 1) {
          fillCandidateTable<true, 1>(candidate);
        }
      }
    } else if (fillOnlyBackground) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandBkg.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandBkg.size());
      }

      for (const auto& candidate : reconstructedCandBkg) {
        if (downSampleBkgFactor < 1.) {
          float const pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
        // Bkg candidates are not matched to MC so rely on selections only
        if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
          fillCandidateTable<true, 1>(candidate);
        }
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandSig.size() + reconstructedCandBkg.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandSig.size() + reconstructedCandBkg.size());
      }

      for (const auto& candidate : reconstructedCandSig) {
        if (candidate.isCandidateSwapped() == 0) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isCandidateSwapped() == 1) {
          fillCandidateTable<true, 1>(candidate);
        }
      }

      for (const auto& candidate : reconstructedCandBkg) {
        // Bkg candidates are not matched to MC so rely on selections only
        if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
          fillCandidateTable<true, 1>(candidate);
        }
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(mcParticles.size());
    for (const auto& particle : mcParticles) {
      rowCandidateFullParticles(
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorXicToPKPi, processMc, "Process MC tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXicToPKPi>(cfgc));
  return workflow;
}
