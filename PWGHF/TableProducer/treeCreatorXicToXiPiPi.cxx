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

/// \file treeCreatorXicToXiPiPi.cxx
/// \brief Writer of Ξc± → Ξ∓ π± π± candidates in the form of flat tables to be stored in TTrees.
///
/// \author Phil Lennart Stahlhut <phil.lennart.stahlhut@cern.ch>, CERN

// include directives (C++, external, ROOT, O2, O2Physics, PWG-non-HF, PWGHF)
// using directives

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace full
{
// track indices
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int);  //! Selection flag of candidate (output of candidateSelector)
DECLARE_SOA_INDEX_COLUMN_FULL(Xi, xi, int, Tracks, "_pi0");
DECLARE_SOA_INDEX_COLUMN_FULL(Pi0, pi0, int, Tracks, "_pi0");
DECLARE_SOA_INDEX_COLUMN_FULL(Pi1, pi1, int, Tracks, "_pi1");
//vertices
DECLARE_SOA_COLUMN(PvX, pvX, float);
DECLARE_SOA_COLUMN(PvY, pvY, float);  
DECLARE_SOA_COLUMN(PvZ, pvZ, float);
DECLARE_SOA_COLUMN(PvErrX, pvErrX, float);
DECLARE_SOA_COLUMN(PvErrY, pvErrY, float);  
DECLARE_SOA_COLUMN(PvErrZ, pvErrZ, float);    
DECLARE_SOA_COLUMN(SvX, svX, float);  
DECLARE_SOA_COLUMN(SvY, svY, float);  
DECLARE_SOA_COLUMN(SvZ, svZ, float);
DECLARE_SOA_COLUMN(Chi2Sv, chi2Sv, float);
DECLARE_SOA_COLUMN(SvErrX, svErrX, float);  
DECLARE_SOA_COLUMN(SvErrY, svErrY, float);  
DECLARE_SOA_COLUMN(SvErrZ, svErrZ, float);
DECLARE_SOA_COLUMN(DecVtxXiX, decVtxXiX, float);  
DECLARE_SOA_COLUMN(DecVtxXiY, decVtxXiY, float);  
DECLARE_SOA_COLUMN(DecVtxXiZ, decVtxXiZ, float);
DECLARE_SOA_COLUMN(Chi2XiVtx, xhi2XiVtx, float);
DECLARE_SOA_COLUMN(DecVtxLamX, decVtxLamX, float);  
DECLARE_SOA_COLUMN(DecVtxLamY, decVtxLamY, float);  
DECLARE_SOA_COLUMN(DecVtxLamZ, decVtxLamZ, float);
DECLARE_SOA_COLUMN(Chi2LamVtx, chi2LamVtx, float);
// properties of XicPlus
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(E, e, float);                                                   //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(M, m, float);                                                   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(P, p, float);                                                   //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                   //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                               //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                               //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(Ct, ct, float);                                                 //! Proper lifetime time ctau of candidate (cm)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! Cosine pointing angle of candidate in transverse plane
// properties of daughter tracks       
DECLARE_SOA_COLUMN(PtXi, ptXi, float);                                                  //! Transverse momentum of Xi (prong0) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterXi, impactParameterXi, float);                        //! Impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedXi, impactParameterNormalisedXi, float);    //! Normalised impact parameter of Xi (prong0)
DECLARE_SOA_COLUMN(PtPi0, ptPi0, float);                                                //! Transverse momentum of Pi0 (prong1) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi0, impactParameterPi0, float);                      //! Impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi0, impactParameterNormalisedPi0, float);  //! Normalised impact parameter of Pi0 (prong1)
DECLARE_SOA_COLUMN(PtPi1, ptPi1, float);                                                //! Transverse momentum of Pi1 (prong2) (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterPi1, impactParameterPi1, float);                      //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(ImpactParameterNormalisedPi1, impactParameterNormalisedPi1, float);  //! Normalised impact parameter of Pi1 (prong2)
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                  //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(CpaXi, cpaXi, float);
DECLARE_SOA_COLUMN(CpaXYXi, cpaXYXi, float);
DECLARE_SOA_COLUMN(CpaLam, cpaLam, float);
DECLARE_SOA_COLUMN(CpaXYLam, cpaXYLam, float);
DECLARE_SOA_COLUMN(DcaPi0Pi1, dcaPi0Pi1, float);
DECLARE_SOA_COLUMN(DcaPi0Xi, dcaPi0Xi, float);
DECLARE_SOA_COLUMN(DcaPi1Xi, dcaPi1Xi, float);
// PID daughters
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);                                 //! TPC Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);                                 //! TOF Nsigma separation for prong1 with pion mass hypothesis
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int); //! Event rejection flag
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);         //! Run number
} // namespace full

DECLARE_SOA_TABLE(HfCandXicLites, "AOD", "HFCANDXICLITE",
                  full::XiId,
                  full::Pi0Id,
                  full::Pi1Id,
                  full::CandidateSelFlag,
                  full::PvX,
                  full::PvY,
                  full::PvZ,
                  full::SvX,
                  full::SvY,
                  full::SvZ,
                  full::Chi2Sv,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterPi1,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  //full::DcaPi0Pi1,
                  //full::DcaPi0Xi,
                  //full::DcaPi1Xi,
                  hf_cand_xictoxipipi::FlagMcMatchRec);

DECLARE_SOA_TABLE(HfCandXicFulls, "AOD", "HFCANDXICFULL",
                  full::XiId,
                  full::Pi0Id,
                  full::Pi1Id,
                  full::CandidateSelFlag,
                  full::PvX,
                  full::PvY,
                  full::PvZ,
                  full::PvErrX,
                  full::PvErrY,
                  full::PvErrZ,
                  full::SvX,
                  full::SvY,
                  full::SvZ,
                  full::Chi2Sv,
                  full::SvErrX,
                  full::SvErrY,
                  full::SvErrZ,
                  full::DecVtxXiX,
                  full::DecVtxXiY,
                  full::DecVtxXiZ,
                  //full::Chi2XiVtx,
                  full::DecVtxLamX,
                  full::DecVtxLamY,
                  full::DecVtxLamZ,
                  //full::Chi2LamVtx,
                  full::Sign,
                  full::E,
                  full::M,
                  full::P,
                  full::Pt,
                  full::Y,
                  full::Eta,
                  full::Phi,
                  full::Ct,
                  full::DecayLength,
                  full::DecayLengthNormalised,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::PtXi,
                  full::PtPi0,
                  full::PtPi1,
                  full::ImpactParameterXi,
                  full::ImpactParameterNormalisedXi,
                  full::ImpactParameterPi0,
                  full::ImpactParameterNormalisedPi0,
                  full::ImpactParameterPi1,
                  full::ImpactParameterNormalisedPi1,
                  full::MaxNormalisedDeltaIP,
                  full::CpaXi,
                  full::CpaXYXi,
                  full::CpaLam,
                  full::CpaXYLam,
                  //full::DcaPi0Pi1,
                  //full::DcaPi0Xi,
                  //full::DcaPi1Xi,
                  hf_cand_xictoxipipi::FlagMcMatchRec);

DECLARE_SOA_TABLE(HfCandXicFullEvs, "AOD", "HFCANDXICFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandXicFullPs, "AOD", "HFCANDXICFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_xictoxipipi::FlagMcMatchGen);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXicToXiPiPi {
  Produces<o2::aod::HfCandXicLites> rowCandidateLite;
  Produces<o2::aod::HfCandXicFulls> rowCandidateFull;
  Produces<o2::aod::HfCandXicFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandXicFullPs> rowCandidateFullParticles;

  Configurable<int> selectionFlagXic{"selectionXic", 1, "Selection Flag for Xic"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  
  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCandXic, aod::HfCandXicMcRec, aod::HfSelXicToXiPiPi>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_xic::isSelXicToXiPiPi >= selectionFlagXic;

  Partition<SelectedCandidatesMc> recSig = nabs(aod::hf_cand_xictoxipipi::flagMcMatchRec) == (int8_t)BIT(aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);
  Partition<SelectedCandidatesMc> recBg = nabs(aod::hf_cand_xictoxipipi::flagMcMatchRec) != (int8_t)BIT(aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi);

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillEvent(const T& collision, int isEventReject, int runNumber)
  {
    rowCandidateFullEvents(
      collision.bcId(),
      collision.numContrib(),
      collision.posX(),
      collision.posY(),
      collision.posZ(),
      isEventReject,  //! filled with 0
      runNumber       //! filled with 1
      );
  }

  template <bool doMc = false, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc = 0;
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
    }
    if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.cascadeId(),
        candidate.pi0Id(),
        candidate.pi1Id(),
        candidate.isSelXicToXiPiPi(),
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.chi2PCA(),
        candidate.sign(),
        candidate.e(o2::constants::physics::MassXiCPlus),
        candidate.invMassXic(),
        candidate.p(),
        candidate.pt(),
        candidate.y(o2::constants::physics::MassXiCPlus),
        candidate.eta(),
        candidate.phi(),
        candidate.ct(o2::constants::physics::MassXiCPlus),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.cosPaXi(),
        candidate.cosPaXYXi(),
        candidate.cosPaLambda(),
        candidate.cosPaXYLambda(),
        //candidate.dcaPi0Pi1(),
        //candidate.dcaPi0Xi(),
        //candidate.dcaPi1Xi(),
        flagMc);
    } else {
      rowCandidateFull(
        candidate.cascadeId(),
        candidate.pi0Id(),
        candidate.pi1Id(),
        candidate.isSelXicToXiPiPi(),
        candidate.posX(),
        candidate.posY(),
        candidate.posZ(),
        candidate.xPvErr(),
        candidate.yPvErr(),
        candidate.zPvErr(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        candidate.chi2PCA(),
        candidate.xSvErr(),
        candidate.ySvErr(),
        candidate.zSvErr(),
        candidate.xDecayVtxXi(),
        candidate.yDecayVtxXi(),
        candidate.zDecayVtxXi(),
        //candidate.kfCascadeChi2(),
        candidate.xDecayVtxLambda(),
        candidate.yDecayVtxLambda(),
        candidate.zDecayVtxLambda(),
        //candidate.kfV0Chi2(),
        candidate.sign(),
        candidate.e(o2::constants::physics::MassXiCPlus),
        candidate.invMassXic(),
        candidate.p(),
        candidate.pt(),
        candidate.y(o2::constants::physics::MassXiCPlus),
        candidate.eta(),
        candidate.phi(),
        candidate.ct(o2::constants::physics::MassXiCPlus),
        candidate.decayLength(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXY(),
        candidate.decayLengthXYNormalised(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameterNormalised0(),
        candidate.impactParameter1(),
        candidate.impactParameterNormalised1(),
        candidate.impactParameter2(),
        candidate.impactParameterNormalised2(),
        candidate.maxNormalisedDeltaIP(),
        candidate.cosPaXi(),
        candidate.cosPaXYXi(),
        candidate.cosPaLambda(),
        candidate.cosPaXYLambda(),
        //candidate.dcaPi0Pi1(),
        //candidate.dcaPi0Xi(),
        //candidate.dcaPi1Xi(),
        flagMc);
    }
  }

  void processData(aod::Collisions const& collisions,
                   soa::Filtered<soa::Join<aod::HfCandXic, aod::HfSelXicToXiPiPi>> const& candidates,
                   TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (fillOnlyBackground && downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (pseudoRndm >= downSampleBkgFactor && candidate.pt() < ptMaxForDownSample) {
          continue;
        }
      }
      fillCandidateTable(candidate);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 SelectedCandidatesMc const& candidates,
                 soa::Join<aod::McParticles, aod::HfCandXicMcGen> const& particles,
                 TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillOnlySignal) {
      rowCandidateFull.reserve(recSig.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recSig.size());
      }
      for (const auto& candidate : recSig) {
        fillCandidateTable<true>(candidate);
      }
    } else if (fillOnlyBackground) {
      rowCandidateFull.reserve(recBg.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(recBg.size());
      }
      for (const auto& candidate : recBg) {
        float pseudoRndm = candidate.ptProng1() * 1000. - (int64_t)(candidate.ptProng1() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
        fillCandidateTable<true>(candidate);
      }
    } else {
      rowCandidateFull.reserve(candidates.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        fillCandidateTable<true>(candidate);
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (TESTBIT(std::abs(particle.flagMcMatchGen()), aod::hf_cand_xictoxipipi::DecayType::XicToXiPiPi)) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassXiCPlus),
          particle.flagMcMatchGen());
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorXicToXiPiPi, processMc, "Process MC", false);
};


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorXicToXiPiPi>(cfgc)};
}
