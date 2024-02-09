// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file treeCreatorDsToKKPi.cxx
/// \brief Writer of Ds to KKpi candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug, local optimization of analysis on small samples, or ML training.
///        In this file are defined and filled the output tables
///
/// \author Stefano Politanò <stefano.politano@polito.it>, Politecnico & INFN, Torino
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

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
// Decay daughters
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);           //! Transverse momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PProng0, pProng0, float);             //! Momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);           //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(PProng1, pProng1, float);             //! Momentum of prong1 (in GeV/c)
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);           //! Transverse momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(PProng2, pProng2, float);             //! Momentum of prong2 (GeV/c)
DECLARE_SOA_COLUMN(NSigTpcPi0, nSigTpcPi0, float);       //! TPC Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa0, nSigTpcKa0, float);       //! TPC Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi0, nSigTofPi0, float);       //! TOF Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa0, nSigTofKa0, float);       //! TOF Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi0, nSigTpcTofPi0, float); //! TPC and TOF combined Nsigma separation for prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa0, nSigTpcTofKa0, float); //! TPC and TOF combined Nsigma separation for prong0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);       //! TPC Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa1, nSigTpcKa1, float);       //! TPC Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);       //! TOF Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa1, nSigTofKa1, float);       //! TOF Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi1, nSigTpcTofPi1, float); //! TPC and TOF combined Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa1, nSigTpcTofKa1, float); //! TPC and TOF combined Nsigma separation for prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPi2, nSigTpcPi2, float);       //! TPC Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKa2, nSigTpcKa2, float);       //! TPC Nsigma separation for prong2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi2, nSigTofPi2, float);       //! TOF Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKa2, nSigTofKa2, float);       //! TOF Nsigma separation for prong2 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPi2, nSigTpcTofPi2, float); //! TPC and TOF combined Nsigma separation for prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKa2, nSigTpcTofKa2, float); //! TPC and TOF combined Nsigma separation for prong2 with kaon mass hypothesis
// Candidate
DECLARE_SOA_COLUMN(M, m, float);                                             //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                             //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                             //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                         //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                             //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(Ct, ct, float);                                           //! Proper lifetime times c of candidate (cm)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                         //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                     //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);       //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(ImpactParameterXY, impactParameterXY, float);             //! Transverse impact parameter of candidate (cm)
DECLARE_SOA_COLUMN(DeltaMassPhi, deltaMassPhi, float);                       //! Absolute mass difference between kaon-pair and phi-meson invariant mass (Gev/c2)
DECLARE_SOA_COLUMN(AbsCos3PiK, absCos3PiK, float);                           //! Cube of absolute value of the cosine of pion-kaon angle in the phi rest frame
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int); //! Event rejection flag
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);         //! Run number
} // namespace full

DECLARE_SOA_TABLE(HfCandDsLites, "AOD", "HFCANDDSLITE",
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
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
                  hf_sel_candidate_ds::IsSelDsToKKPi,
                  hf_sel_candidate_ds::IsSelDsToPiKK,
                  full::M,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::ImpactParameterXY,
                  full::DeltaMassPhi,
                  full::AbsCos3PiK,
                  hf_cand::Chi2PCA,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::FlagMcDecayChanRec)

DECLARE_SOA_TABLE(HfCandDsFulls, "AOD", "HFCANDDSFULL",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
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
                  hf_sel_candidate_ds::IsSelDsToKKPi,
                  hf_sel_candidate_ds::IsSelDsToPiKK,
                  hf_cand::XSecondaryVertex,
                  hf_cand::YSecondaryVertex,
                  hf_cand::ZSecondaryVertex,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::E,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::Cpa,
                  full::CpaXY,
                  full::MaxNormalisedDeltaIP,
                  full::ImpactParameterXY,
                  full::DeltaMassPhi,
                  full::AbsCos3PiK,
                  hf_cand::Chi2PCA,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_3prong::FlagMcDecayChanRec);

DECLARE_SOA_TABLE(HfCandDsFullEvs, "AOD", "HFCANDDSFULLEV",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandDsFullPs, "AOD", "HFCANDDSFULLP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcGen);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorDsToKKPi {
  Produces<o2::aod::HfCandDsFulls> rowCandidateFull;
  Produces<o2::aod::HfCandDsFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandDsFullPs> rowCandidateFullParticles;
  Produces<o2::aod::HfCandDsLites> rowCandidateLite;

  Configurable<int> decayChannelDs{"decayChannelDs", 1, "Switch between decay channels: 1 for Ds->PhiPi->KKpi, 2 for Ds->K0*K->KKPi"};
  Configurable<int> decayChannelDplusMc{"decayChannelDplusMc", 1, "Switch between decay channels for MC Dplus: 1 for Dplus->PhiPi->KKpi, 2 for Dplus->K0*K->KKPi"};
  Configurable<bool> FillDplusMc{"fillDplusMc", false, "Switch to fill Dplus MC information"};
  Configurable<int> selectionFlagDs{"selectionFlagDs", 1, "Selection flag for Ds"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  // parameters for production of training samples
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillOnlyBackground{"fillOnlyBackground", false, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;

  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDsMcReco = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfCand3ProngMcRec>>;
  using CandDsMcGen = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPiExt, aod::TracksPidKaExt>;

  int offsetDplusDecayChannel = aod::hf_cand_3prong::DecayChannelDToKKPi::DplusToPhiPi - aod::hf_cand_3prong::DecayChannelDToKKPi::DsToPhiPi; // Offset between Dplus and Ds to use the same decay channel. See aod::hf_cand_3prong::DecayChannelDToKKPi

  Filter filterSelectCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;
  Filter filterMcGenMatching = nabs(o2::aod::hf_cand_3prong::flagMcMatchGen) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && (aod::hf_cand_3prong::flagMcDecayChanGen == decayChannelDs || (FillDplusMc && aod::hf_cand_3prong::flagMcDecayChanGen == (decayChannelDplusMc + offsetDplusDecayChannel))); // Do not store Dplus MC if FillDplusMc is false

  Partition<CandDsData> selectedDsToKKPiCand = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlagDs;
  Partition<CandDsData> selectedDsToPiKKCand = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlagDs;

  Partition<CandDsMcReco> reconstructedCandSig = nabs(aod::hf_cand_3prong::flagMcMatchRec) == static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi)) && (aod::hf_cand_3prong::flagMcDecayChanRec == decayChannelDs || (FillDplusMc && aod::hf_cand_3prong::flagMcDecayChanRec == (decayChannelDplusMc + offsetDplusDecayChannel))); // Do not store Dplus MC if FillDplusMc is false
  Partition<CandDsMcReco> reconstructedCandBkg = nabs(aod::hf_cand_3prong::flagMcMatchRec) != static_cast<int8_t>(BIT(aod::hf_cand_3prong::DecayType::DsToKKPi));

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

  /// Fill table accounting for MC information and mass hypothesis
  /// \param doMc true to fill MC information
  /// \param massHypo mass hypothesis considered: 0 = KKPi, 1 = PiKK
  /// \param candidate is candidate
  template <bool doMc = false, int massHypo = 0, typename T>
  void fillCandidateTable(const T& candidate)
  {
    int8_t flagMc = 0;
    int8_t originMc = 0;
    int8_t channelMc = 0;
    float yCand = 0;
    float eCand = 0;
    float ctCand = 0;
    if constexpr (doMc) {
      flagMc = candidate.flagMcMatchRec();
      originMc = candidate.originMcRec();
      channelMc = candidate.flagMcDecayChanRec();
      if (FillDplusMc && candidate.flagMcDecayChanRec() == (decayChannelDplusMc + offsetDplusDecayChannel)) {
        yCand = hfHelper.yDplus(candidate);
        eCand = hfHelper.eDplus(candidate);
        ctCand = hfHelper.ctDplus(candidate);
      } else {
        yCand = hfHelper.yDs(candidate);
        eCand = hfHelper.eDs(candidate);
        ctCand = hfHelper.ctDs(candidate);
      }
    }

    float invMassDs = 0;
    float deltaMassPhiKK = 0;
    float absCos3PiKDs = 0;
    if constexpr (massHypo == 0) {
      invMassDs = hfHelper.invMassDsToKKPi(candidate);
      deltaMassPhiKK = hfHelper.deltaMassPhiDsToKKPi(candidate);
      absCos3PiKDs = std::abs(hfHelper.cos3PiKDsToKKPi(candidate));
    } else if constexpr (massHypo == 1) {
      invMassDs = hfHelper.invMassDsToPiKK(candidate);
      deltaMassPhiKK = hfHelper.deltaMassPhiDsToPiKK(candidate);
      absCos3PiKDs = std::abs(hfHelper.cos3PiKDsToPiKK(candidate));
    }

    auto prong0 = candidate.template prong0_as<TracksWPid>();
    auto prong1 = candidate.template prong1_as<TracksWPid>();
    auto prong2 = candidate.template prong2_as<TracksWPid>();

    if (fillCandidateLiteTable) {
      rowCandidateLite(
        candidate.ptProng0(),
        candidate.ptProng1(),
        candidate.ptProng2(),
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
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
        candidate.isSelDsToKKPi(),
        candidate.isSelDsToPiKK(),
        invMassDs,
        candidate.pt(),
        candidate.eta(),
        candidate.phi(),
        yCand,
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterXY(),
        deltaMassPhiKK,
        absCos3PiKDs,
        candidate.chi2PCA(),
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
        candidate.ptProng0(),
        RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
        candidate.ptProng1(),
        RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
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
        candidate.isSelDsToKKPi(),
        candidate.isSelDsToPiKK(),
        candidate.xSecondaryVertex(),
        candidate.ySecondaryVertex(),
        candidate.zSecondaryVertex(),
        invMassDs,
        candidate.pt(),
        candidate.p(),
        ctCand,
        candidate.eta(),
        candidate.phi(),
        yCand,
        eCand,
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.maxNormalisedDeltaIP(),
        candidate.impactParameterXY(),
        deltaMassPhiKK,
        absCos3PiKDs,
        candidate.chi2PCA(),
        flagMc,
        originMc,
        channelMc);
    }
  }

  void processData(aod::Collisions const& collisions,
                   CandDsData const& candidates,
                   TracksWPid const&)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, 0, 1);
    }

    // Filling candidate properties
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(selectedDsToKKPiCand.size() + selectedDsToPiKKCand.size());
    } else {
      rowCandidateFull.reserve(selectedDsToKKPiCand.size() + selectedDsToPiKKCand.size());
    }

    for (const auto& candidate : selectedDsToKKPiCand) {
      if (downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable<false, 0>(candidate);
    }

    for (const auto& candidate : selectedDsToPiKKCand) {
      if (downSampleBkgFactor < 1.) {
        float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
        if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
          continue;
        }
      }
      fillCandidateTable<false, 1>(candidate);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDsToKKPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 CandDsMcReco const& candidates,
                 CandDsMcGen const& mcParticles,
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
          float pseudoRndm = candidate.ptProng0() * 1000. - (int64_t)(candidate.ptProng0() * 1000);
          if (candidate.pt() < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
            continue;
          }
        }
        // Bkg candidates are not matched to MC so rely on selections only
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
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
        if (candidate.isSelDsToKKPi() >= selectionFlagDs) {
          fillCandidateTable<true, 0>(candidate);
        }
        if (candidate.isSelDsToPiKK() >= selectionFlagDs) {
          fillCandidateTable<true, 1>(candidate);
        }
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(mcParticles.size());
    for (const auto& particle : mcParticles) {
      rowCandidateFullParticles(
        particle.mcCollision().bcId(),
        particle.pt(),
        particle.eta(),
        particle.phi(),
        RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::constants::physics::MassDS),
        particle.flagMcMatchGen(),
        particle.originMcGen());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDsToKKPi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorDsToKKPi>(cfgc)};
}
