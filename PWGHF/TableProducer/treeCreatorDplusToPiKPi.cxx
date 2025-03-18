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

/// \file treeCreatorDplusToPiKPi.cxx
/// \brief Writer of D+ → π+ K- π+ candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples or ML training.
///        In this file are defined and filled the output tables
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

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
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);                     //! Radius of secondary vertex (cm)
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                                     //! Transverse momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PProng0, pProng0, float);                                       //! Momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float); //! Normalised impact parameter of prong0
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                                     //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(PProng1, pProng1, float);                                       //! Momentum of prong1 (in GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float); //! Normalised impact parameter of prong1
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);                                     //! Transverse momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(PProng2, pProng2, float);                                       //! Momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float); //! Normalised impact parameter of prong2
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int);                       //! Selection flag of candidate (output of candidateSelector)
DECLARE_SOA_COLUMN(M, m, float);                                                   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                   //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                   //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                               //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                               //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                   //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);                                 //! TPC Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);                                 //! TPC Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);                                 //! TOF Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);                                 //! TOF Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);                                 //! TPC Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);                                 //! TPC Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);                                 //! TOF Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);                                 //! TOF Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);                                 //! TPC Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);                                 //! TPC Nsigma separation for prong2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);                                 //! TOF Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);                                 //! TOF Nsigma separation for prong2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float);                           //! TPC and TOF combined Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa0, nSigTpcTofKa0, float);                           //! TPC and TOF combined Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float);                           //! TPC and TOF combined Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float);                           //! TPC and TOF combined Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float);                           //! TPC and TOF combined Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa2, nSigTpcTofKa2, float);                           //! TPC and TOF combined Nsigma separation for prong2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);             //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(Ct, ct, float);                                                 //! Proper lifetime times c of candidate (cm)
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int); //! Event rejection flag
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);         //! Run number
// ML scores
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float); //! ML score of the first configured index
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float); //! ML score of the second configured index
} // namespace full
DECLARE_SOA_TABLE(HfCandDpMls, "AOD", "HFCANDDPML",
                  full::MlScore0,
                  full::MlScore1)

DECLARE_SOA_TABLE(HfCandDpLites, "AOD", "HFCANDDPLITE",
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
                  hf_cand::ImpactParameterZ0,
                  hf_cand::ImpactParameterZ1,
                  hf_cand::ImpactParameterZ2,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTpcTofPi1,
                  full::NSigTpcTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcKa2,
                  full::NSigTofPi2,
                  full::NSigTofKa2,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofKa2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::FlagMcDecayChanRec)

DECLARE_SOA_TABLE(HfCandDpFulls, "AOD", "HFCANDDPFULL",
                  collision::BCId,
                  collision::NumContrib,
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
                  hf_cand::ImpactParameterZ0,
                  hf_cand::ImpactParameterZ1,
                  hf_cand::ImpactParameterZ2,
                  hf_cand::ErrorImpactParameterZ0,
                  hf_cand::ErrorImpactParameterZ1,
                  hf_cand::ErrorImpactParameterZ2,
                  full::NSigTpcPi0,
                  full::NSigTpcKa0,
                  full::NSigTofPi0,
                  full::NSigTofKa0,
                  full::NSigTpcTofPi0,
                  full::NSigTpcTofKa0,
                  full::NSigTpcPi1,
                  full::NSigTpcKa1,
                  full::NSigTofPi1,
                  full::NSigTofKa1,
                  full::NSigTpcTofPi1,
                  full::NSigTpcTofKa1,
                  full::NSigTpcPi2,
                  full::NSigTpcKa2,
                  full::NSigTofPi2,
                  full::NSigTofKa2,
                  full::NSigTpcTofPi2,
                  full::NSigTpcTofKa2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::FlagMcDecayChanRec);

DECLARE_SOA_TABLE(HfCandDpFullEvs, "AOD", "HFCANDDPFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandDpFullPs, "AOD", "HFCANDDPFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcGen);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorDplusToPiKPi {
  Produces<o2::aod::HfCandDpFulls> rowCandidateFull;
  Produces<o2::aod::HfCandDpFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandDpFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandDpLites> rowCandidateLite;
  Produces<o2::aod::HfCandDpMls> rowCandidateMl;

  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 1, "Selection Flag for Dplus"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlySignalMl{"fillOnlySignalMl", false, "Flag to fill derived tables with MC and ML info"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<std::vector<int>> classMl{"classMlindexes", {0, 2}, "Indexes of ML bkg and non-prompt scores."};

  HfHelper hfHelper;

  using SelectedCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi>>;
  using MatchedGenCandidatesMc = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using SelectedCandidatesMcWithMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus;
  Filter filterMcGenMatching = nabs(o2::aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi));

  Partition<SelectedCandidatesMc> reconstructedCandSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) || nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)); // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±
  Partition<SelectedCandidatesMc> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi));
  Partition<SelectedCandidatesMcWithMl> reconstructedCandSigMl = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DplusToPiKPi)) || nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) || nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DstarToPiKPiBkg)); // DecayType::DsToKKPi is used to flag both Ds± → K± K∓ π± and D± → K± K∓ π±

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
      isEventReject,
      runNumber);
  }

  template <bool doMc = false, bool doMl = false, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    int8_t channelMc = 0;
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
      channelMc = candidate.flagMcDecayChanRec();
    }

    std::vector<float> outputMl = {-999., -999.};
    if constexpr (doMl) {
      for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
        outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
      }
      rowCandidateMl(
        outputMl[0],
        outputMl[1]);
    }

    auto prong0 = candidate.template prong0_as<TracksWPid>();
    auto prong1 = candidate.template prong1_as<TracksWPid>();
    auto prong2 = candidate.template prong2_as<TracksWPid>();

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
        candidate.impactParameterZ0(),
        candidate.impactParameterZ1(),
        candidate.impactParameterZ2(),
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaKa(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaKa(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        prong2.tpcNSigmaPi(),
        prong2.tpcNSigmaKa(),
        prong2.tofNSigmaPi(),
        prong2.tofNSigmaKa(),
        prong2.tpcTofNSigmaPi(),
        prong2.tpcTofNSigmaKa(),
        candidate.isSelDplusToPiKPi(),
        hfHelper.invMassDplusToPiKPi(candidate),
        candidate.pt(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yDplus(candidate),
        flagMc,
        originMc,
        channelMc);
    } else {
      rowCandidateFull(
        candidate.collision().bcId(),
        candidate.collision().numContrib(),
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
        candidate.impactParameterZ0(),
        candidate.impactParameterZ1(),
        candidate.impactParameterZ2(),
        candidate.errorImpactParameterZ0(),
        candidate.errorImpactParameterZ1(),
        candidate.errorImpactParameterZ2(),
        prong0.tpcNSigmaPi(),
        prong0.tpcNSigmaKa(),
        prong0.tofNSigmaPi(),
        prong0.tofNSigmaKa(),
        prong0.tpcTofNSigmaPi(),
        prong0.tpcTofNSigmaKa(),
        prong1.tpcNSigmaPi(),
        prong1.tpcNSigmaKa(),
        prong1.tofNSigmaPi(),
        prong1.tofNSigmaKa(),
        prong1.tpcTofNSigmaPi(),
        prong1.tpcTofNSigmaKa(),
        prong2.tpcNSigmaPi(),
        prong2.tpcNSigmaKa(),
        prong2.tofNSigmaPi(),
        prong2.tofNSigmaKa(),
        prong2.tpcTofNSigmaPi(),
        prong2.tpcTofNSigmaKa(),
        candidate.isSelDplusToPiKPi(),
        hfHelper.invMassDplusToPiKPi(candidate),
        candidate.pt(),
        candidate.p(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        hfHelper.ctDplus(candidate),
        candidate.eta(),
        candidate.phi(),
        hfHelper.yDplus(candidate),
        hfHelper.eDplus(candidate),
        flagMc,
        originMc,
        channelMc);
    }
  }

  void processData(aod::Collisions const& collisions,
                   soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>> const& candidates,
                   TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidates.size());
    } else {
      rowCandidateFull.reserve(candidates.size());
    }
    for (const auto& candidate : candidates) {
      if (downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable(candidate);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDplusToPiKPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 SelectedCandidatesMc const& candidates,
                 MatchedGenCandidatesMc const& particles,
                 SelectedCandidatesMcWithMl const&,
                 TracksWPid const&)
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
        fillCandidateTable<true>(candidate);
      }
    } else if (fillOnlySignalMl) {
      rowCandidateMl.reserve(reconstructedCandSigMl.size());
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandSigMl.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandSigMl.size());
      }
      for (const auto& candidate : reconstructedCandSigMl) {
        if (downSampleBkgFactor < 1.) {
          float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
        fillCandidateTable<true, true>(candidate);
      }
    } else if (fillOnlyBackground) {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(reconstructedCandBkg.size());
      } else {
        rowCandidateFull.reserve(reconstructedCandBkg.size());
      }
      for (const auto& candidate : reconstructedCandBkg) {
        if (downSampleBkgFactor < 1.) {
          float pseudoRndm = candidate.ptProng0() * 1000. - static_cast<int64_t>(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
        fillCandidateTable<true>(candidate);
      }
    } else {
      if (fillCandidateLiteTable) {
        rowCandidateLite.reserve(candidates.size());
      } else {
        rowCandidateFull.reserve(candidates.size());
      }
      for (const auto& candidate : candidates) {
        fillCandidateTable<true>(candidate);
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      rowCandidateFullParticles(
        particle.mcCollision().bcId(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecay::y(particle.pVector(), o2::constants::physics::MassDPlus),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDplusToPiKPi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorDplusToPiKPi>(cfgc)};
}
