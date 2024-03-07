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

/// \file dataCreatorToXiPi.cxx
/// \brief Writer of the omegac0 or xic0 to Xi Pi candidates in the form of flat tables to be stored self contained derived data
///        In this file are defined and filled the output tables
///
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace full
{
// from creator
DECLARE_SOA_COLUMN(XPv, xPv, float);
DECLARE_SOA_COLUMN(YPv, yPv, float);
DECLARE_SOA_COLUMN(ZPv, zPv, float);
DECLARE_SOA_COLUMN(XDecayVtxCharmBaryon, xDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(YDecayVtxCharmBaryon, yDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(ZDecayVtxCharmBaryon, zDecayVtxCharmBaryon, float);
DECLARE_SOA_COLUMN(XDecayVtxCascade, xDecayVtxCascade, float);
DECLARE_SOA_COLUMN(YDecayVtxCascade, yDecayVtxCascade, float);
DECLARE_SOA_COLUMN(ZDecayVtxCascade, zDecayVtxCascade, float);
DECLARE_SOA_COLUMN(XDecayVtxV0, xDecayVtxV0, float);
DECLARE_SOA_COLUMN(YDecayVtxV0, yDecayVtxV0, float);
DECLARE_SOA_COLUMN(ZDecayVtxV0, zDecayVtxV0, float);
DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t); // sign of pi <- xi
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxPiFromCharmBaryon, pxPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyPiFromCharmBaryon, pyPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzPiFromCharmBaryon, pzPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PxPiFromCasc, pxPiFromCasc, float);
DECLARE_SOA_COLUMN(PyPiFromCasc, pyPiFromCasc, float);
DECLARE_SOA_COLUMN(PzPiFromCasc, pzPiFromCasc, float);
DECLARE_SOA_COLUMN(PxPosV0Dau, pxPosV0Dau, float);
DECLARE_SOA_COLUMN(PyPosV0Dau, pyPosV0Dau, float);
DECLARE_SOA_COLUMN(PzPosV0Dau, pzPosV0Dau, float);
DECLARE_SOA_COLUMN(PxNegV0Dau, pxNegV0Dau, float);
DECLARE_SOA_COLUMN(PyNegV0Dau, pyNegV0Dau, float);
DECLARE_SOA_COLUMN(PzNegV0Dau, pzNegV0Dau, float);
DECLARE_SOA_COLUMN(ImpactParCascXY, impactParCascXY, float);
DECLARE_SOA_COLUMN(ImpactParPiFromCharmBaryonXY, impactParPiFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParPiFromCharmBaryonXY, errImpactParPiFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, double);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, double);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, double);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, double);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, double);
DECLARE_SOA_COLUMN(EtaPiFromCasc, etaPiFromCasc, double);
DECLARE_SOA_COLUMN(EtaPiFromCharmBaryon, etaPiFromCharmBaryon, double);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(IsPionGlbTrkWoDca, isPionGlbTrkWoDca, bool);
DECLARE_SOA_COLUMN(PionItsNCls, pionItsNCls, uint8_t);
// from creator - MC (rec level)
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
// from selector
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCasc, tpcNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCasc, tofNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);

} // namespace full

DECLARE_SOA_TABLE(HfToXiPiEvDatas, "AOD", "HFTOXIPIEVDATA",
                  collision::NumContrib,
                  collision::Chi2,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ);

DECLARE_SOA_TABLE(HfToXiPiDatas, "AOD", "HFTOXIPIDATA",
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxPiFromCasc, full::PyPiFromCasc, full::PzPiFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaPiFromCasc, full::EtaPiFromCharmBaryon,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::ErrorDecayLengthCharmBaryon,
                  full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginRec, full::CollisionMatched);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfDataCreatorToXiPi {

  Produces<o2::aod::HfToXiPiDatas> rowCandidateData;
  Produces<o2::aod::HfToXiPiEvDatas> rowEvData;

  using MyTrackTable = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>;

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillEvent(const T& collision)
  {
    rowEvData(
      collision.numContrib(),
      collision.chi2(),
      collision.posX(),
      collision.posY(),
      collision.posZ());
  }

  template <class TMyTracks, typename T>
  void fillCandidate(const T& candidate, int8_t flagMc, int8_t originMc, bool collisionMatched)
  {
    if(candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {
    
    rowCandidateData(
      candidate.xPv(),
      candidate.yPv(),
      candidate.zPv(),
      candidate.collision().numContrib(),
      candidate.collision().chi2(),
      candidate.xDecayVtxCharmBaryon(),
      candidate.yDecayVtxCharmBaryon(),
      candidate.zDecayVtxCharmBaryon(),
      candidate.xDecayVtxCascade(),
      candidate.yDecayVtxCascade(),
      candidate.zDecayVtxCascade(),
      candidate.xDecayVtxV0(),
      candidate.yDecayVtxV0(),
      candidate.zDecayVtxV0(),
      candidate.signDecay(),
      candidate.pxCharmBaryon(),
      candidate.pyCharmBaryon(),
      candidate.pzCharmBaryon(),
      candidate.pxPiFromCharmBaryon(),
      candidate.pyPiFromCharmBaryon(),
      candidate.pzPiFromCharmBaryon(),
      candidate.pxPiFromCasc(),
      candidate.pyPiFromCasc(),
      candidate.pzPiFromCasc(),
      candidate.pxPosV0Dau(),
      candidate.pyPosV0Dau(),
      candidate.pzPosV0Dau(),
      candidate.pxNegV0Dau(),
      candidate.pyNegV0Dau(),
      candidate.pzNegV0Dau(),
      candidate.impactParCascXY(),
      candidate.impactParPiFromCharmBaryonXY(),
      candidate.errImpactParCascXY(),
      candidate.errImpactParPiFromCharmBaryonXY(),
      candidate.invMassLambda(),
      candidate.invMassCascade(),
      candidate.invMassCharmBaryon(),
      candidate.etaV0PosDau(),
      candidate.etaV0NegDau(),
      candidate.etaPiFromCasc(),
      candidate.etaPiFromCharmBaryon(),
      candidate.dcaXYToPvV0Dau0(),
      candidate.dcaXYToPvV0Dau1(),
      candidate.dcaXYToPvCascDau(),
      candidate.dcaCascDau(),
      candidate.dcaV0Dau(),
      candidate.dcaCharmBaryonDau(),
      candidate.errorDecayLengthCharmBaryon(),
      candidate.template piFromCharmBaryon_as<TMyTracks>().isGlobalTrackWoDCA(),
      candidate.template piFromCharmBaryon_as<TMyTracks>().itsNCls(),
      candidate.pidTpcInfoStored(),
      candidate.pidTofInfoStored(),
      candidate.tpcNSigmaPiFromCharmBaryon(),
      candidate.tpcNSigmaPiFromCasc(),
      candidate.tpcNSigmaPiFromLambda(),
      candidate.tpcNSigmaPrFromLambda(),
      candidate.tofNSigmaPiFromCharmBaryon(),
      candidate.tofNSigmaPiFromCasc(),
      candidate.tofNSigmaPiFromLambda(),
      candidate.tofNSigmaPrFromLambda(),
      flagMc,
      originMc,
      collisionMatched);
    }
  }

  void processData(aod::Collisions const& collisions, MyTrackTable const&,
                   soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEvData.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateData.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<MyTrackTable>(candidate, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfDataCreatorToXiPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions, MyTrackTable const&,
                 soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEvData.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision);
    }

    // Filling candidate properties
    rowCandidateData.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<MyTrackTable>(candidate, candidate.flagMcMatchRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfDataCreatorToXiPi, processMc, "Process MC", false);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfDataCreatorToXiPi>(cfgc)};
}
