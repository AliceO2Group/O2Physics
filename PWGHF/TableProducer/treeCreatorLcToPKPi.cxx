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

/// \file treeCreatorLcToPKPi.cxx
/// \brief Writer of the 3 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
///
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, CERN

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
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
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(MassKPi, massKPi, float); // invariant mass of the candidate Kpi daughters
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
DECLARE_SOA_COLUMN(IsCandidateSwapped, isCandidateSwapped, int8_t);
DECLARE_SOA_INDEX_COLUMN_FULL(Candidate, candidate, int, HfCand3ProngWPidPiKaPr, "_0");
DECLARE_SOA_INDEX_COLUMN(McParticle, mcParticle);
DECLARE_SOA_COLUMN(Channel, channel, int8_t); // direct or resonant
DECLARE_SOA_COLUMN(MlScoreFirstClass, mlScoreFirstClass, float);
DECLARE_SOA_COLUMN(MlScoreSecondClass, mlScoreSecondClass, float);
DECLARE_SOA_COLUMN(MlScoreThirdClass, mlScoreThirdClass, float);
// Events
DECLARE_SOA_INDEX_COLUMN(McCollision, mcCollision);
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(CentFT0A, centFT0A, float);
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, float);
DECLARE_SOA_COLUMN(CentFT0M, centFT0M, float);
DECLARE_SOA_COLUMN(CentFV0A, centFV0A, float);
DECLARE_SOA_COLUMN(CentFDDM, centFDDM, float);
DECLARE_SOA_COLUMN(MultZeqNTracksPV, multZeqNTracksPV, float);
} // namespace full

namespace kf
{
DECLARE_SOA_COLUMN(X, x, float);                                         //! decay vertex X coordinate
DECLARE_SOA_COLUMN(Y, y, float);                                         //! decay vertex Y coordinate
DECLARE_SOA_COLUMN(Z, z, float);                                         //! decay vertex Z coordinate
DECLARE_SOA_COLUMN(ErrX, errX, float);                                   //! decay vertex X coordinate error
DECLARE_SOA_COLUMN(ErrY, errY, float);                                   //! decay vertex Y coordinate error
DECLARE_SOA_COLUMN(ErrZ, errZ, float);                                   //! decay vertex Z coordinate error
DECLARE_SOA_COLUMN(ErrPVX, errPVX, float);                               //! event vertex X coordinate error
DECLARE_SOA_COLUMN(ErrPVY, errPVY, float);                               //! event vertex Y coordinate error
DECLARE_SOA_COLUMN(ErrPVZ, errPVZ, float);                               //! event vertex Z coordinate error
DECLARE_SOA_COLUMN(Chi2PrimProton, chi2PrimProton, float);               //! Chi2 of prong's approach to the PV
DECLARE_SOA_COLUMN(Chi2PrimKaon, chi2PrimKaon, float);                   //! Chi2 of prong's approach to the PV
DECLARE_SOA_COLUMN(Chi2PrimPion, chi2PrimPion, float);                   //! Chi2 of prong's approach to the PV
DECLARE_SOA_COLUMN(DcaProtonKaon, dcaProtonKaon, float);                 //! Distance of closest approach between 2 prongs, cm
DECLARE_SOA_COLUMN(DcaProtonPion, dcaProtonPion, float);                 //! Distance of closest approach between 2 prongs, cm
DECLARE_SOA_COLUMN(DcaPionKaon, dcaPionKaon, float);                     //! Distance of closest approach between 2 prongs, cm
DECLARE_SOA_COLUMN(Chi2GeoProtonKaon, chi2GeoProtonKaon, float);         //! Chi2 of two prongs' approach to each other
DECLARE_SOA_COLUMN(Chi2GeoProtonPion, chi2GeoProtonPion, float);         //! Chi2 of two prongs' approach to each other
DECLARE_SOA_COLUMN(Chi2GeoPionKaon, chi2GeoPionKaon, float);             //! Chi2 of two prongs' approach to each other
DECLARE_SOA_COLUMN(Chi2Geo, chi2Geo, float);                             //! chi2 geo of the full candidate
DECLARE_SOA_COLUMN(Chi2Topo, chi2Topo, float);                           //! chi2 topo of the full candidate (chi2prim of candidate to PV)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                     //! decay length, cm
DECLARE_SOA_COLUMN(DecayLengthError, decayLengthError, float);           //! decay length error
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float); //! decay length over its error
DECLARE_SOA_COLUMN(T, t, float);                                         //! proper lifetime, ps
DECLARE_SOA_COLUMN(ErrT, errT, float);                                   //! lifetime error
DECLARE_SOA_COLUMN(MassInv, massInv, float);                             //! invariant mass
DECLARE_SOA_COLUMN(P, p, float);                                         //! momentum
DECLARE_SOA_COLUMN(Pt, pt, float);                                       //! transverse momentum
DECLARE_SOA_COLUMN(ErrP, errP, float);                                   //! momentum error
DECLARE_SOA_COLUMN(ErrPt, errPt, float);                                 //! transverse momentum error
DECLARE_SOA_COLUMN(IsSelected, isSelected, int);                         //! flag whether candidate was selected in candidateSelectorLc task
DECLARE_SOA_COLUMN(SigBgStatus, sigBgStatus, int);                       //! 0 bg, 1 prompt, 2 non-prompt, 3 wrong order of prongs, -1 default value (impossible, should not be the case), -999 for data
DECLARE_SOA_COLUMN(NSigTpcPi, nSigTpcPi, float);
DECLARE_SOA_COLUMN(NSigTpcKa, nSigTpcKa, float);
DECLARE_SOA_COLUMN(NSigTpcPr, nSigTpcPr, float);
DECLARE_SOA_COLUMN(NSigTofPi, nSigTofPi, float);
DECLARE_SOA_COLUMN(NSigTofKa, nSigTofKa, float);
DECLARE_SOA_COLUMN(NSigTofPr, nSigTofPr, float);
DECLARE_SOA_COLUMN(NSigTpcTofPi, nSigTpcTofPi, float);
DECLARE_SOA_COLUMN(NSigTpcTofKa, nSigTpcTofKa, float);
DECLARE_SOA_COLUMN(NSigTpcTofPr, nSigTpcTofPr, float);
DECLARE_SOA_COLUMN(MultNTracksPV, multNTracksPV, int);
} // namespace kf

namespace kf_collision
{
DECLARE_SOA_COLUMN(PosXErr, posXErr, float); //! PV X coordinate uncertainty
DECLARE_SOA_COLUMN(PosYErr, posYErr, float); //! PV Y coordinate uncertainty
DECLARE_SOA_COLUMN(PosZErr, posZErr, float); //! PV Z coordinate uncertainty
DECLARE_SOA_COLUMN(McPosX, mcPosX, float);   //! PV X coordinate uncertainty
DECLARE_SOA_COLUMN(McPosY, mcPosY, float);   //! PV Y coordinate uncertainty
DECLARE_SOA_COLUMN(McPosZ, mcPosZ, float);   //! PV Z coordinate uncertainty
} // namespace kf_collision

namespace mc_match
{
DECLARE_SOA_COLUMN(P, p, float);           //! Momentum, GeV/c
DECLARE_SOA_COLUMN(Pt, pt, float);         //! Transverse momentum, GeV/c
DECLARE_SOA_COLUMN(XDecay, xDecay, float); //! Secondary (decay) vertex X coordinate, cm
DECLARE_SOA_COLUMN(YDecay, yDecay, float); //! Secondary (decay) vertex Y coordinate, cm
DECLARE_SOA_COLUMN(ZDecay, zDecay, float); //! Secondary (decay) vertex Z coordinate, cm
DECLARE_SOA_COLUMN(LDecay, lDecay, float); //! Decay length, cm (distance between PV and SV, curvature is neglected)
DECLARE_SOA_COLUMN(TDecay, tDecay, float); //! Proper lifetime, ps
DECLARE_SOA_COLUMN(XEvent, xEvent, float); //! Primary (event) vertex X coordinate, cm
DECLARE_SOA_COLUMN(YEvent, yEvent, float); //! Primary (event) vertex Y coordinate, cm
DECLARE_SOA_COLUMN(ZEvent, zEvent, float); //! Primary (event) vertex Z coordinate, cm
} // namespace mc_match

DECLARE_SOA_TABLE(HfCandLcMCs, "AOD", "HFCANDLCMC",
                  mc_match::P, mc_match::Pt,
                  mc_match::XDecay, mc_match::YDecay, mc_match::ZDecay, mc_match::LDecay,
                  mc_match::TDecay,
                  mc_match::XEvent, mc_match::YEvent, mc_match::ZEvent)

DECLARE_SOA_TABLE(HfCandLcKFs, "AOD", "HFCANDLCKF",
                  kf::X, kf::Y, kf::Z, kf::ErrX, kf::ErrY, kf::ErrZ,
                  kf::ErrPVX, kf::ErrPVY, kf::ErrPVZ,
                  kf::Chi2PrimProton, kf::Chi2PrimKaon, kf::Chi2PrimPion,
                  kf::DcaProtonKaon, kf::DcaProtonPion, kf::DcaPionKaon,
                  kf::Chi2GeoProtonKaon, kf::Chi2GeoProtonPion, kf::Chi2GeoPionKaon,
                  kf::Chi2Geo, kf::Chi2Topo, kf::DecayLength, kf::DecayLengthError, kf::DecayLengthNormalised, kf::T, kf::ErrT,
                  kf::MassInv, kf::P, kf::Pt, kf::ErrP, kf::ErrPt,
                  kf::IsSelected, kf::SigBgStatus,
                  kf::MultNTracksPV,
                  kf::NSigTpcPr,
                  kf::NSigTpcKa,
                  kf::NSigTpcPi,
                  kf::NSigTofPr,
                  kf::NSigTofKa,
                  kf::NSigTofPi,
                  kf::NSigTpcTofPr,
                  kf::NSigTpcTofKa,
                  kf::NSigTpcTofPi);

DECLARE_SOA_TABLE(HfCandLcLites, "AOD", "HFCANDLCLITE",
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
                  hf_cand::BitmapProngsContributorsPV,
                  hf_cand::Chi2PCA,
                  full::DecayLength,
                  full::DecayLengthXY,
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
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::Cpa,
                  full::CpaXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcRec,
                  full::IsCandidateSwapped,
                  full::Channel,
                  full::MassKPi,
                  full::MlScoreFirstClass,
                  full::MlScoreSecondClass,
                  full::MlScoreThirdClass);

DECLARE_SOA_TABLE(HfCollIdLCLite, "AOD", "HFCOLLIDLCLITE",
                  full::CollisionId);

DECLARE_SOA_TABLE(HfCandLcFulls, "AOD", "HFCANDLCFULL",
                  full::CollisionId,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  hf_cand::NProngsContributorsPV,
                  hf_cand::BitmapProngsContributorsPV,
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
                  full::CandidateSelFlag,
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
                  full::IsCandidateSwapped,
                  full::CandidateId,
                  full::Channel,
                  full::MassKPi,
                  full::MlScoreFirstClass,
                  full::MlScoreSecondClass,
                  full::MlScoreThirdClass);

DECLARE_SOA_TABLE(HfCandLcFullEvs, "AOD", "HFCANDLCFULLEV",
                  full::CollisionId,
                  full::McCollisionId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  kf_collision::PosXErr,
                  kf_collision::PosYErr,
                  kf_collision::PosZErr,
                  kf_collision::McPosX,
                  kf_collision::McPosY,
                  kf_collision::McPosZ,
                  full::IsEventReject,
                  full::RunNumber,
                  full::CentFT0A,
                  full::CentFT0C,
                  full::CentFT0M,
                  full::CentFV0A,
                  full::CentFDDM,
                  full::MultZeqNTracksPV,
                  kf::MultNTracksPV);

DECLARE_SOA_TABLE(HfCandLcFullPs, "AOD", "HFCANDLCFULLP",
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::FlagMc,
                  full::OriginMcGen,
                  mc_match::P,
                  mc_match::XDecay,
                  mc_match::YDecay,
                  mc_match::ZDecay,
                  mc_match::LDecay,
                  mc_match::TDecay,
                  mc_match::XEvent,
                  mc_match::YEvent,
                  mc_match::ZEvent);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorLcToPKPi {
  Produces<o2::aod::HfCandLcFulls> rowCandidateFull;
  Produces<o2::aod::HfCandLcLites> rowCandidateLite;
  Produces<o2::aod::HfCandLcKFs> rowCandidateKF;
  Produces<o2::aod::HfCandLcMCs> rowCandidateMC;
  Produces<o2::aod::HfCollIdLCLite> rowCollisionId;
  Produces<o2::aod::HfCandLcFullEvs> rowCandidateFullEvents;
  Produces<o2::aod::HfCandLcFullPs> rowCandidateFullParticles;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<bool> fillCandidateLiteTable{"fillCandidateLiteTable", false, "Switch to fill lite table with candidate properties"};
  Configurable<bool> fillCollIdTable{"fillCollIdTable", false, "Fill a single-column table with collision index"};
  Configurable<bool> fillCandidateMcTable{"fillCandidateMcTable", false, "Switch to fill a table with MC particles matched to candidates"};
  Configurable<bool> applyMl{"applyMl", false, "Whether ML was used in candidateSelectorLc"};
  Configurable<bool> keepSignalMc{"keepSignalMc", false, "Fill MC tree with signal candidates"};
  Configurable<bool> keepBkgMc{"keepBkgMc", false, "Fill MC tree with background candidates"};
  Configurable<bool> keepCorrBkgMC{"keepCorrBkgMC", false, "Flag to keep correlated background sources (Λc+ -> p K− π+ π0, p π− π+, p K− K+ and other charm hadrons)"};
  Configurable<double> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of candidates to store in the tree"};
  Configurable<float> downSampleBkgPtMax{"downSampleBkgPtMax", 100.f, "Max. pt for background downsampling"};

  constexpr static float UndefValueFloat = -999.f;
  constexpr static int UndefValueInt = -999;
  constexpr static float NanoToPico = 1000.f;

  using TracksWPid = soa::Join<aod::Tracks, aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidPr, aod::PidTpcTofFullPr>;
  using Cents = soa::Join<aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs>;

  // number showing MC status of the candidate (signal or background, prompt or non-prompt etc.)
  enum SigBgStatus : int {
    Background = 0, // combinatorial background, at least one of the prongs do not originate from the Lc decay
    Prompt,         // signal with Lc produced directly in the event
    NonPrompt,      // signal with Lc produced aftewards the event, e.g. during decay of beauty particle
    WrongOrder,     // all the prongs are from Lc decay, but proton and pion hypothesis are swapped
    Default = -1    // impossible, should not be the case, to catch logical error if any
  };

  /// \brief function which determines if the candidate corresponds to MC-particle or belongs to a combinatorial background
  /// \param candidate candidate to be checked for being signal or background
  /// \param candFlag 0 for PKPi hypothesis and 1 for PiKP hypothesis
  /// \return SigBgStatus enum with value encoding MC status of the candidate
  template <typename CandType>
  SigBgStatus determineSignalBgStatus(const CandType& candidate, int candFlag)
  {
    const int flag = candidate.flagMcMatchRec();
    const int origin = candidate.originMcRec();
    const int swapped = candidate.isCandidateSwapped();

    SigBgStatus status{Default};

    if (std::abs(flag) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
      if (swapped == 0) {
        if (candFlag == 0) {
          if (origin == RecoDecay::OriginType::Prompt) {
            status = Prompt;
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            status = NonPrompt;
          }
        } else {
          status = WrongOrder;
        }
      } else {
        if (candFlag == 1) {
          if (origin == RecoDecay::OriginType::Prompt) {
            status = Prompt;
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            status = NonPrompt;
          }
        } else {
          status = WrongOrder;
        }
      }
    } else {
      status = Background;
    }

    return status;
  }

  void init(InitContext const&)
  {
    std::array<bool, 8> processes = {doprocessDataNoCentralityWithDCAFitterN, doprocessDataWithCentralityWithDCAFitterN, doprocessDataNoCentralityWithKFParticle, doprocessDataWithCentralityWithKFParticle,
                                     doprocessMcNoCentralityWithDCAFitterN, doprocessMcWithCentralityWithDCAFitterN, doprocessMcNoCentralityWithKFParticle, doprocessMcWithCentralityWithKFParticle};
    if (std::accumulate(processes.begin(), processes.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function must be enabled at a time.");
    }
    if ((std::accumulate(processes.begin(), processes.begin() + 4, 0) != 0) && fillCandidateMcTable) {
      LOGP(fatal, "fillCandidateMcTable can be activated only in case of MC processing.");
    }
  }

  /// \brief function to fill event properties
  /// \param collisions Collision table
  template <bool UseCentrality, bool IsMc, typename Colls>
  void fillEventProperties(Colls const& collisions)
  {
    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (const auto& collision : collisions) {

      float centFT0A = -1.f;
      float centFT0C = -1.f;
      float centFT0M = -1.f;
      float centFV0A = -1.f;
      float centFDDM = -1.f;
      if constexpr (UseCentrality) {
        centFT0A = collision.centFT0A();
        centFT0C = collision.centFT0C();
        centFT0M = collision.centFT0M();
        centFV0A = collision.centFV0A();
        centFDDM = collision.centFDDM();
      }

      float mcPosX{UndefValueFloat};
      float mcPosY{UndefValueFloat};
      float mcPosZ{UndefValueFloat};
      int mcCollId{-1};

      if constexpr (IsMc) {
        auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();

        mcPosX = mcCollision.posX();
        mcPosY = mcCollision.posY();
        mcPosZ = mcCollision.posZ();

        mcCollId = collision.mcCollisionId();
      }

      rowCandidateFullEvents(
        collision.globalIndex(),
        mcCollId,
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        std::sqrt(collision.covXX()),
        std::sqrt(collision.covYY()),
        std::sqrt(collision.covZZ()),
        mcPosX,
        mcPosY,
        mcPosZ,
        0,
        collision.bc().runNumber(),
        centFT0A,
        centFT0C,
        centFT0M,
        centFV0A,
        centFDDM,
        collision.multZeqNTracksPV(),
        collision.multNTracksPV());
    }
  }

  /// \brief function to reserve tables size
  /// \param candidatesSize size of the candidates table
  /// \param isMc boolean flag whether MC or data is processed
  template <int ReconstructionType>
  void reserveTables(size_t candidatesSize, bool isMc)
  {
    if (fillCandidateLiteTable) {
      rowCandidateLite.reserve(candidatesSize * 2);
    } else {
      rowCandidateFull.reserve(candidatesSize * 2);
    }
    if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
      rowCandidateKF.reserve(candidatesSize * 2);
    }
    if (fillCollIdTable) {
      /// save also candidate collision indices
      rowCollisionId.reserve(candidatesSize);
    }
    if (isMc && fillCandidateMcTable) {
      rowCandidateMC.reserve(candidatesSize * 2);
    }
  }

  /// \brief function to evaluate invariant mass of the Lc candidate and KPi pair
  /// \param candidate candidate instance
  /// \param candFlag flag indicating if PKPi (0) or PiKP (1) hypothesis is used
  template <typename CandType>
  std::pair<float, float> evaluateInvariantMassesDCAFitter(CandType const& candidate, int candFlag)
  {
    const float invMass = candFlag == 0 ? HfHelper::invMassLcToPKPi(candidate) : HfHelper::invMassLcToPiKP(candidate);
    const float invMassKPi = candFlag == 0 ? HfHelper::invMassKPiPairLcToPKPi(candidate) : HfHelper::invMassKPiPairLcToPiKP(candidate);

    return std::make_pair(invMass, invMassKPi);
  }

  /// \brief function to get ML score values for the current candidate and assign them to input parameters
  /// \param candidate candidate instance
  /// \param candidateMlScore instance of handler of vectors with ML scores associated with the current candidate
  /// \param mlScoreFirstClass ML score for belonging to the first class
  /// \param mlScoreSecondClass ML score for belonging to the second class
  /// \param mlScoreThirdClass ML score for belonging to the third class
  /// \param candFlag flag indicating if PKPi (0) or PiKP (1) hypothesis is used
  void assignMlScores(aod::HfMlLcToPKPi::iterator const& candidateMlScore, float& mlScoreFirstClass, float& mlScoreSecondClass, float& mlScoreThirdClass, int candFlag)
  {
    std::vector<float> mlScores;
    if (candFlag == 0) {
      std::copy(candidateMlScore.mlProbLcToPKPi().begin(), candidateMlScore.mlProbLcToPKPi().end(), std::back_inserter(mlScores));
    } else {
      std::copy(candidateMlScore.mlProbLcToPiKP().begin(), candidateMlScore.mlProbLcToPiKP().end(), std::back_inserter(mlScores));
    }
    constexpr int IndexFirstClass{0};
    constexpr int IndexSecondClass{1};
    constexpr int IndexThirdClass{2};
    if (mlScores.empty()) {
      return; // when candidateSelectorLc rejects a candidate by "usual", non-ML cut, the ml score vector remains empty
    }
    mlScoreFirstClass = mlScores.at(IndexFirstClass);
    mlScoreSecondClass = mlScores.at(IndexSecondClass);
    if (mlScores.size() > IndexThirdClass) {
      mlScoreThirdClass = mlScores.at(IndexThirdClass);
    }
  }

  /// \brief function to fill lite table
  /// \param candidate candidate instance
  /// \param candidateMlScore instance of handler of vectors with ML scores associated with the current candidate
  /// \param candFlag flag indicating if PKPi (0) or PiKP (1) hypothesis is used
  template <bool IsMc, typename CandType>
  void fillLiteTable(CandType const& candidate, aod::HfMlLcToPKPi::iterator const& candidateMlScore, int candFlag)
  {
    auto [functionInvMass, functionInvMassKPi] = evaluateInvariantMassesDCAFitter(candidate, candFlag);
    const float functionCt = HfHelper::ctLc(candidate);
    const float functionY = HfHelper::yLc(candidate);

    int8_t functionFlagMcMatchRec{0};
    int8_t functionOriginMcRec{0};
    int8_t functionIsCandidateSwapped{0};
    int8_t functionFlagMcDecayChanRec{-1};

    if constexpr (IsMc) {
      functionFlagMcMatchRec = candidate.flagMcMatchRec();
      functionOriginMcRec = candidate.originMcRec();
      functionIsCandidateSwapped = candidate.isCandidateSwapped();
      functionFlagMcDecayChanRec = candidate.flagMcDecayChanRec();
    }

    float mlScoreFirstClass{UndefValueFloat};
    float mlScoreSecondClass{UndefValueFloat};
    float mlScoreThirdClass{UndefValueFloat};

    if (applyMl) {
      assignMlScores(candidateMlScore, mlScoreFirstClass, mlScoreSecondClass, mlScoreThirdClass, candFlag);
    }

    rowCandidateLite(
      candidate.posX(),
      candidate.posY(),
      candidate.posZ(),
      candidate.nProngsContributorsPV(),
      candidate.bitmapProngsContributorsPV(),
      candidate.chi2PCA(),
      candidate.decayLength(),
      candidate.decayLengthXY(),
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
      1 << candFlag,
      functionInvMass,
      candidate.pt(),
      candidate.cpa(),
      candidate.cpaXY(),
      functionCt,
      candidate.eta(),
      candidate.phi(),
      functionY,
      functionFlagMcMatchRec,
      functionOriginMcRec,
      functionIsCandidateSwapped,
      functionFlagMcDecayChanRec,
      functionInvMassKPi,
      mlScoreFirstClass,
      mlScoreSecondClass,
      mlScoreThirdClass);

    if (fillCollIdTable) {
      /// save also candidate collision indices
      rowCollisionId(candidate.collisionId());
    }
  }

  /// \brief function to fill lite table
  /// \param candidate candidate instance
  /// \param candidateMlScore instance of handler of vectors with ML scores associated with the current candidate
  /// \param candFlag flag indicating if PKPi (0) or PiKP (1) hypothesis is used
  template <bool IsMc, typename CandType>
  void fillFullTable(CandType const& candidate, aod::HfMlLcToPKPi::iterator const& candidateMlScore, int candFlag)
  {
    auto [functionInvMass, functionInvMassKPi] = evaluateInvariantMassesDCAFitter(candidate, candFlag);
    const float functionCt = HfHelper::ctLc(candidate);
    const float functionY = HfHelper::yLc(candidate);
    const float functionE = HfHelper::eLc(candidate);

    int8_t functionFlagMcMatchRec{0};
    int8_t functionOriginMcRec{0};
    int8_t functionIsCandidateSwapped{0};
    int8_t functionFlagMcDecayChanRec{-1};

    if constexpr (IsMc) {
      functionFlagMcMatchRec = candidate.flagMcMatchRec();
      functionOriginMcRec = candidate.originMcRec();
      functionIsCandidateSwapped = candidate.isCandidateSwapped();
      functionFlagMcDecayChanRec = candidate.flagMcDecayChanRec();
    }

    float mlScoreFirstClass{UndefValueFloat};
    float mlScoreSecondClass{UndefValueFloat};
    float mlScoreThirdClass{UndefValueFloat};

    if (applyMl) {
      assignMlScores(candidateMlScore, mlScoreFirstClass, mlScoreSecondClass, mlScoreThirdClass, candFlag);
    }

    rowCandidateFull(
      candidate.collisionId(),
      candidate.posX(),
      candidate.posY(),
      candidate.posZ(),
      candidate.nProngsContributorsPV(),
      candidate.bitmapProngsContributorsPV(),
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
      1 << candFlag,
      functionInvMass,
      candidate.pt(),
      candidate.p(),
      candidate.cpa(),
      candidate.cpaXY(),
      functionCt,
      candidate.eta(),
      candidate.phi(),
      functionY,
      functionE,
      functionFlagMcMatchRec,
      functionOriginMcRec,
      functionIsCandidateSwapped,
      candidate.globalIndex(),
      functionFlagMcDecayChanRec,
      functionInvMassKPi,
      mlScoreFirstClass,
      mlScoreSecondClass,
      mlScoreThirdClass);
  }

  /// \brief function to fill lite table
  /// \param candidate candidate instance
  /// \param collision collision, to which the candidate belongs
  /// \param candFlag flag indicating if PKPi (0) or PiKP (1) hypothesis is used
  /// \param functionSelection flag indicating if candidate was selected by candidateSelectorLc task
  /// \param sigbgstatus for MC: number indicating if candidate is prompt, non-prompt or background; for data: UndefValueInt
  template <typename CandType, typename CollType>
  void fillKFTable(CandType const& candidate,
                   CollType const& collision,
                   int candFlag,
                   int functionSelection,
                   int sigbgstatus)
  {
    float chi2primProton{};
    float chi2primPion{};
    float dcaProtonKaon{};
    float dcaPionKaon{};
    float chi2GeoProtonKaon{};
    float chi2GeoPionKaon{};
    float mass{};
    float valueTpcNSigmaPr{};
    const float valueTpcNSigmaKa = candidate.nSigTpcKa1();
    float valueTpcNSigmaPi{};
    float valueTofNSigmaPr{};
    const float valueTofNSigmaKa = candidate.nSigTofKa1();
    float valueTofNSigmaPi{};
    float valueTpcTofNSigmaPr{};
    const float valueTpcTofNSigmaKa = candidate.tpcTofNSigmaKa1();
    float valueTpcTofNSigmaPi{};
    if (candFlag == 0) {
      chi2primProton = candidate.kfChi2PrimProng0();
      chi2primPion = candidate.kfChi2PrimProng2();
      dcaProtonKaon = candidate.kfDcaProng0Prong1();
      dcaPionKaon = candidate.kfDcaProng1Prong2();
      chi2GeoProtonKaon = candidate.kfChi2GeoProng0Prong1();
      chi2GeoPionKaon = candidate.kfChi2GeoProng1Prong2();
      mass = candidate.kfMassPKPi();
      valueTpcNSigmaPr = candidate.nSigTpcPr0();
      valueTpcNSigmaPi = candidate.nSigTpcPi2();
      valueTofNSigmaPr = candidate.nSigTofPr0();
      valueTofNSigmaPi = candidate.nSigTofPi2();
      valueTpcTofNSigmaPr = candidate.tpcTofNSigmaPr0();
      valueTpcTofNSigmaPi = candidate.tpcTofNSigmaPi2();
    } else {
      chi2primProton = candidate.kfChi2PrimProng2();
      chi2primPion = candidate.kfChi2PrimProng0();
      dcaProtonKaon = candidate.kfDcaProng1Prong2();
      dcaPionKaon = candidate.kfDcaProng0Prong1();
      chi2GeoProtonKaon = candidate.kfChi2GeoProng1Prong2();
      chi2GeoPionKaon = candidate.kfChi2GeoProng0Prong1();
      mass = candidate.kfMassPiKP();
      valueTpcNSigmaPr = candidate.nSigTpcPr2();
      valueTpcNSigmaPi = candidate.nSigTpcPi0();
      valueTofNSigmaPr = candidate.nSigTofPr2();
      valueTofNSigmaPi = candidate.nSigTofPi0();
      valueTpcTofNSigmaPr = candidate.tpcTofNSigmaPr2();
      valueTpcTofNSigmaPi = candidate.tpcTofNSigmaPi0();
    }
    const float svX = candidate.xSecondaryVertex();
    const float svY = candidate.ySecondaryVertex();
    const float svZ = candidate.zSecondaryVertex();
    const float svErrX = candidate.kfXError();
    const float svErrY = candidate.kfYError();
    const float svErrZ = candidate.kfZError();
    const float pvErrX = candidate.kfXPVError();
    const float pvErrY = candidate.kfYPVError();
    const float pvErrZ = candidate.kfZPVError();
    const float chi2primKaon = candidate.kfChi2PrimProng1();
    const float dcaProtonPion = candidate.kfDcaProng0Prong2();
    const float chi2GeoProtonPion = candidate.kfChi2GeoProng0Prong2();
    const float chi2Geo = candidate.kfChi2Geo();
    const float chi2Topo = candidate.kfChi2Topo();
    const float decayLength = candidate.kfDecayLength();
    const float dl = candidate.kfDecayLengthError();
    const float pt = std::sqrt(candidate.kfPx() * candidate.kfPx() + candidate.kfPy() * candidate.kfPy());
    const float deltaPt = std::sqrt(candidate.kfPx() * candidate.kfPx() * candidate.kfErrorPx() * candidate.kfErrorPx() +
                                    candidate.kfPy() * candidate.kfPy() * candidate.kfErrorPy() * candidate.kfErrorPy()) /
                          pt;
    const float p = std::sqrt(pt * pt + candidate.kfPz() * candidate.kfPz());
    const float deltaP = std::sqrt(pt * pt * deltaPt * deltaPt +
                                   candidate.kfPz() * candidate.kfPz() * candidate.kfErrorPz() * candidate.kfErrorPz()) /
                         p;
    const float lifetime = decayLength * MassLambdaCPlus / LightSpeedCm2PS / p;
    const float deltaT = dl * MassLambdaCPlus / LightSpeedCm2PS / p;
    rowCandidateKF(
      svX, svY, svZ, svErrX, svErrY, svErrZ,
      pvErrX, pvErrY, pvErrZ,
      chi2primProton, chi2primKaon, chi2primPion,
      dcaProtonKaon, dcaProtonPion, dcaPionKaon,
      chi2GeoProtonKaon, chi2GeoProtonPion, chi2GeoPionKaon,
      chi2Geo, chi2Topo, decayLength, dl, decayLength / dl, lifetime, deltaT,
      mass, p, pt, deltaP, deltaPt,
      functionSelection, sigbgstatus,
      collision.multNTracksPV(),
      valueTpcNSigmaPr,
      valueTpcNSigmaKa,
      valueTpcNSigmaPi,
      valueTofNSigmaPr,
      valueTofNSigmaKa,
      valueTofNSigmaPi,
      valueTpcTofNSigmaPr,
      valueTpcTofNSigmaKa,
      valueTpcTofNSigmaPi);
  }

  /// \brief core function to fill tables in MC
  /// \param collisions Collision table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  template <bool UseCentrality, int ReconstructionType, typename Colls, typename CandType>
  void fillTablesMc(Colls const& collisions,
                    aod::McCollisions const&,
                    CandType const& candidates,
                    aod::HfMlLcToPKPi const& candidateMlScores,
                    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                    soa::Join<TracksWPid, o2::aod::McTrackLabels> const&, aod::BCs const&)
  {

    constexpr bool IsMc = true;

    fillEventProperties<UseCentrality, IsMc>(collisions);

    const size_t candidatesSize = candidates.size();
    reserveTables<ReconstructionType>(candidatesSize, IsMc);

    int iCand{0};
    for (const auto& candidate : candidates) {
      auto candidateMlScore = candidateMlScores.rawIteratorAt(iCand);
      ++iCand;
      float ptProng0 = candidate.ptProng0();
      auto collision = candidate.template collision_as<Colls>();
      auto fillTable = [&](int candFlag) {
        double const pseudoRndm = ptProng0 * 1000. - static_cast<int64_t>(ptProng0 * 1000);
        const int functionSelection = candFlag == 0 ? candidate.isSelLcToPKPi() : candidate.isSelLcToPiKP();
        const int sigbgstatus = determineSignalBgStatus(candidate, candFlag);
        const bool isMcCandidateSignal = (sigbgstatus == Prompt) || (sigbgstatus == NonPrompt);
        const bool passSelection = functionSelection >= selectionFlagLc;
        const bool keepAll = !keepSignalMc && !keepBkgMc && !keepCorrBkgMC;
        const int flag = candidate.flagMcMatchRec();
        const bool isCorrBkg = (std::abs(flag) != 0 && std::abs(flag) != o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) || sigbgstatus == WrongOrder;
        const bool notSkippedBkg = isMcCandidateSignal || candidate.pt() > downSampleBkgPtMax || pseudoRndm < downSampleBkgFactor;
        if (passSelection && notSkippedBkg && (keepAll || (keepSignalMc && isMcCandidateSignal) || (keepBkgMc && !isMcCandidateSignal) || (keepCorrBkgMC && isCorrBkg))) {
          if (fillCandidateLiteTable) {
            fillLiteTable<IsMc>(candidate, candidateMlScore, candFlag);
          } else {
            fillFullTable<IsMc>(candidate, candidateMlScore, candFlag);
          }

          if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
            fillKFTable(candidate, collision, candFlag, functionSelection, sigbgstatus);
          }
          if (fillCandidateMcTable) {
            float p{}, pt{}, svX{}, svY{}, svZ{}, pvX{}, pvY{}, pvZ{}, decayLength{}, lifetime{};
            if (!isMcCandidateSignal) {
              p = UndefValueFloat;
              pt = UndefValueFloat;
              svX = UndefValueFloat;
              svY = UndefValueFloat;
              svZ = UndefValueFloat;
              pvX = UndefValueFloat;
              pvY = UndefValueFloat;
              pvZ = UndefValueFloat;
              decayLength = UndefValueFloat;
              lifetime = UndefValueFloat;
            } else {
              auto mcParticleProng0 = candidate.template prong0_as<soa::Join<TracksWPid, o2::aod::McTrackLabels>>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>();
              auto indexMother = RecoDecay::getMother(particles, mcParticleProng0, o2::constants::physics::Pdg::kLambdaCPlus, true);
              auto particleMother = particles.rawIteratorAt(indexMother);
              auto mcCollision = particleMother.template mcCollision_as<aod::McCollisions>();
              p = particleMother.p();
              pt = particleMother.pt();
              const float p2m = p / MassLambdaCPlus;
              const float gamma = std::sqrt(1 + p2m * p2m); // mother's particle Lorentz factor
              pvX = mcCollision.posX();
              pvY = mcCollision.posY();
              pvZ = mcCollision.posZ();
              svX = mcParticleProng0.vx();
              svY = mcParticleProng0.vy();
              svZ = mcParticleProng0.vz();
              decayLength = RecoDecay::distance(std::array<float, 3>{svX, svY, svZ}, std::array<float, 3>{pvX, pvY, pvZ});
              lifetime = mcParticleProng0.vt() * NanoToPico / gamma; // from ns to ps * from lab time to proper time
            }
            rowCandidateMC(
              p, pt,
              svX, svY, svZ, decayLength, lifetime,
              pvX, pvY, pvZ);
          }
        }
      };

      fillTable(0);
      fillTable(1);
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (const auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
        auto mcDaughter0 = particle.template daughters_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>().begin();
        auto mcCollision = particle.template mcCollision_as<aod::McCollisions>();
        auto p = particle.p();
        const float p2m = p / MassLambdaCPlus;
        const float gamma = std::sqrt(1 + p2m * p2m); // mother's particle Lorentz factor
        const float pvX = mcCollision.posX();
        const float pvY = mcCollision.posY();
        const float pvZ = mcCollision.posZ();
        const float svX = mcDaughter0.vx();
        const float svY = mcDaughter0.vy();
        const float svZ = mcDaughter0.vz();
        const float l = RecoDecay::distance(std::array<float, 3>{svX, svY, svZ}, std::array<float, 3>{pvX, pvY, pvZ});
        const float t = mcDaughter0.vt() * NanoToPico / gamma; // from ns to ps * from lab time to proper time
        rowCandidateFullParticles(
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus),
          particle.flagMcMatchGen(),
          particle.originMcGen(),
          p,
          svX, svY, svZ, l, t,
          pvX, pvY, pvZ);
      }
    }
  }

  /// \brief process function for MC w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcNoCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, aod::PVMults> const& collisions,
                                           aod::McCollisions const& mcCollisions,
                                           soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfCand3ProngMcRec, aod::HfSelLc> const& candidates,
                                           aod::HfMlLcToPKPi const& candidateMlScores,
                                           soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                                           soa::Join<TracksWPid, o2::aod::McTrackLabels> const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<false, aod::hf_cand::VertexerType::DCAFitter>(collisions, mcCollisions, candidates, candidateMlScores, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcNoCentralityWithDCAFitterN, "Process MC tree writer w/o centrality with DCAFitterN", false);

  /// \brief process function for MC with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcWithCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, Cents, aod::PVMults> const& collisions,
                                             aod::McCollisions const& mcCollisions,
                                             soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfCand3ProngMcRec, aod::HfSelLc> const& candidates,
                                             aod::HfMlLcToPKPi const& candidateMlScores,
                                             soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                                             soa::Join<TracksWPid, o2::aod::McTrackLabels> const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<true, aod::hf_cand::VertexerType::DCAFitter>(collisions, mcCollisions, candidates, candidateMlScores, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcWithCentralityWithDCAFitterN, "Process MC tree writer with centrality with DCAFitterN", false);

  /// \brief process function for MC w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param particles Generated particle table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcNoCentralityWithKFParticle(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, aod::PVMults> const& collisions,
                                           aod::McCollisions const& mcCollisions,
                                           soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                                           aod::HfMlLcToPKPi const& candidateMlScores,
                                           soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                                           soa::Join<TracksWPid, o2::aod::McTrackLabels> const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<false, aod::hf_cand::VertexerType::KfParticle>(collisions, mcCollisions, candidates, candidateMlScores, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcNoCentralityWithKFParticle, "Process MC tree writer w/o centrality with KFParticle", false);

  /// \brief process function for MC with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param mcCollisions MC collision table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processMcWithCentralityWithKFParticle(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::PVMultZeqs, Cents, aod::PVMults> const& collisions,
                                             aod::McCollisions const& mcCollisions,
                                             soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                                             aod::HfMlLcToPKPi const& candidateMlScores,
                                             soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                                             soa::Join<TracksWPid, o2::aod::McTrackLabels> const& tracks, aod::BCs const& bcs)
  {
    fillTablesMc<true, aod::hf_cand::VertexerType::KfParticle>(collisions, mcCollisions, candidates, candidateMlScores, particles, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processMcWithCentralityWithKFParticle, "Process MC tree writer with centrality with KFParticle", false);

  /// \brief core function to fill tables in data
  /// \param collisions Collision table
  /// \param candidates Lc->pKpi candidate table
  template <bool UseCentrality, int ReconstructionType, typename Colls, typename CandType>
  void fillTablesData(Colls const& collisions,
                      CandType const& candidates,
                      aod::HfMlLcToPKPi const& candidateMlScores,
                      TracksWPid const&, aod::BCs const&)
  {

    constexpr bool IsMc = false;

    fillEventProperties<UseCentrality, IsMc>(collisions);

    const size_t candidatesSize = candidates.size();
    reserveTables<ReconstructionType>(candidatesSize, IsMc);

    // Filling candidate properties

    int iCand{0};
    for (const auto& candidate : candidates) {
      auto candidateMlScore = candidateMlScores.rawIteratorAt(iCand);
      ++iCand;
      float ptProng0 = candidate.ptProng0();
      auto collision = candidate.template collision_as<Colls>();
      auto fillTable = [&](int candFlag) {
        double const pseudoRndm = ptProng0 * 1000. - static_cast<int64_t>(ptProng0 * 1000);
        const int functionSelection = candFlag == 0 ? candidate.isSelLcToPKPi() : candidate.isSelLcToPiKP();
        if (functionSelection >= selectionFlagLc && (candidate.pt() > downSampleBkgPtMax || (pseudoRndm < downSampleBkgFactor && candidate.pt() < downSampleBkgPtMax))) {
          if (fillCandidateLiteTable) {
            fillLiteTable<IsMc>(candidate, candidateMlScore, candFlag);
          } else {
            fillFullTable<IsMc>(candidate, candidateMlScore, candFlag);
          }

          if constexpr (ReconstructionType == aod::hf_cand::VertexerType::KfParticle) {
            fillKFTable(candidate, collision, candFlag, functionSelection, UndefValueInt);
          }
        }
      };

      fillTable(0);
      fillTable(1);
    }
  }

  /// \brief process function for data w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataNoCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::PVMultZeqs, aod::PVMults> const& collisions,
                                             soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc> const& candidates,
                                             aod::HfMlLcToPKPi const& candidateMlScores,
                                             TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<false, aod::hf_cand::VertexerType::DCAFitter>(collisions, candidates, candidateMlScores, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataNoCentralityWithDCAFitterN, "Process data tree writer w/o centrality with DCAFitterN", false);

  /// \brief process function for data with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataWithCentralityWithDCAFitterN(soa::Join<aod::Collisions, aod::PVMultZeqs, Cents, aod::PVMults> const& collisions,
                                               soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc> const& candidates,
                                               aod::HfMlLcToPKPi const& candidateMlScores,
                                               TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<true, aod::hf_cand::VertexerType::DCAFitter>(collisions, candidates, candidateMlScores, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataWithCentralityWithDCAFitterN, "Process data tree writer with centrality with DCAFitterN", true);

  /// \brief process function for data w/o centrality
  /// \param collisions Collision table w/o join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataNoCentralityWithKFParticle(soa::Join<aod::Collisions, aod::PVMultZeqs, aod::PVMults> const& collisions,
                                             soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                                             aod::HfMlLcToPKPi const& candidateMlScores,
                                             TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<false, aod::hf_cand::VertexerType::KfParticle>(collisions, candidates, candidateMlScores, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataNoCentralityWithKFParticle, "Process data tree writer w/o centrality with KFParticle", false);

  /// \brief process function for data with centrality
  /// \param collisions Collision table with join of the centrality table
  /// \param candidates Lc->pKpi candidate table
  /// \param tracks Track table
  /// \param bcs Bunch-crossing table
  void processDataWithCentralityWithKFParticle(soa::Join<aod::Collisions, aod::PVMultZeqs, Cents, aod::PVMults> const& collisions,
                                               soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelLc, aod::HfCand3ProngKF> const& candidates,
                                               aod::HfMlLcToPKPi const& candidateMlScores,
                                               TracksWPid const& tracks, aod::BCs const& bcs)
  {
    fillTablesData<true, aod::hf_cand::VertexerType::KfParticle>(collisions, candidates, candidateMlScores, tracks, bcs);
  }
  PROCESS_SWITCH(HfTreeCreatorLcToPKPi, processDataWithCentralityWithKFParticle, "Process data tree writer with centrality with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLcToPKPi>(cfgc));
  return workflow;
}
