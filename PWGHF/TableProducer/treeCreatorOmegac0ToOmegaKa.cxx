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

/// \file treeCreatorOmegac0ToOmegaKa.cxx
/// \brief Writer of the omegac0 to Omega Ka candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

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

namespace o2::aod
{
namespace full
{
// collision info
DECLARE_SOA_COLUMN(IsEventSel8, isEventSel8, bool);
DECLARE_SOA_COLUMN(IsEventSelZ, isEventSelZ, bool);
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
DECLARE_SOA_COLUMN(SignDecay, signDecay, int8_t); // sign of ka <- omega
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxKaFromCharmBaryon, pxKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyKaFromCharmBaryon, pyKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzKaFromCharmBaryon, pzKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PxKaFromCasc, pxKaFromCasc, float);
DECLARE_SOA_COLUMN(PyKaFromCasc, pyKaFromCasc, float);
DECLARE_SOA_COLUMN(PzKaFromCasc, pzKaFromCasc, float);
DECLARE_SOA_COLUMN(PxPosV0Dau, pxPosV0Dau, float);
DECLARE_SOA_COLUMN(PyPosV0Dau, pyPosV0Dau, float);
DECLARE_SOA_COLUMN(PzPosV0Dau, pzPosV0Dau, float);
DECLARE_SOA_COLUMN(PxNegV0Dau, pxNegV0Dau, float);
DECLARE_SOA_COLUMN(PyNegV0Dau, pyNegV0Dau, float);
DECLARE_SOA_COLUMN(PzNegV0Dau, pzNegV0Dau, float);
DECLARE_SOA_COLUMN(ImpactParCascXY, impactParCascXY, float);
DECLARE_SOA_COLUMN(ImpactParKaFromCharmBaryonXY, impactParKaFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParKaFromCharmBaryonXY, errImpactParKaFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaKaFromCasc, etaKaFromCasc, float);
DECLARE_SOA_COLUMN(EtaKaFromCharmBaryon, etaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(NormImpParCascade, normImpParCascade, double);
DECLARE_SOA_COLUMN(NormImpParKaFromCharmBar, normImpParKaFromCharmBar, double);
DECLARE_SOA_COLUMN(IsKaonGlbTrkWoDca, isKaonGlbTrkWoDca, bool);
DECLARE_SOA_COLUMN(KaonItsNCls, kaonItsNCls, uint8_t);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
// from selector
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCharmBaryon, tpcNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCasc, tpcNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCharmBaryon, tofNSigmaKaFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCasc, tofNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);

} // namespace full

DECLARE_SOA_TABLE(HfToOmegaKaEvs, "AOD", "HFTOOMEKAEV",
                  full::IsEventSel8, full::IsEventSelZ);

DECLARE_SOA_TABLE(HfOmegac0ToOmegaKaLites, "AOD", "HFTOOMEKALITE",
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxKaFromCharmBaryon, full::PyKaFromCharmBaryon, full::PzKaFromCharmBaryon,
                  full::PxKaFromCasc, full::PyKaFromCasc, full::PzKaFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParKaFromCharmBaryonXY,
                  full::ErrImpactParCascXY, full::ErrImpactParKaFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaKaFromCasc, full::EtaKaFromCharmBaryon,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::ErrorDecayLengthCharmBaryon, full::NormImpParCascade, full::NormImpParKaFromCharmBar,
                  full::IsKaonGlbTrkWoDca, full::KaonItsNCls,
                  full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaKaFromCharmBaryon, full::TpcNSigmaKaFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaKaFromCharmBaryon, full::TofNSigmaKaFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginMcRec, full::CollisionMatched);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorOmegac0ToOmegaKa {

  Produces<o2::aod::HfOmegac0ToOmegaKaLites> rowCandidateLite;
  Produces<o2::aod::HfToOmegaKaEvs> rowEv;

  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};

  using MyTrackTable = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>;
  using MyEventTable = soa::Join<aod::Collisions, aod::EvSels>;

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillEvent(const T& collision, float cutZPv)
  {
    rowEv(collision.sel8(), std::abs(collision.posZ()) < cutZPv);
  }

  template <typename T>
  void fillCandidateLite(const T& candidate, int8_t flagMc, int8_t originMc, bool collisionMatched)
  {
    if (candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {

      rowCandidateLite(
        candidate.xPv(),
        candidate.yPv(),
        candidate.zPv(),
        candidate.template collision_as<MyEventTable>().numContrib(),
        candidate.template collision_as<MyEventTable>().chi2(),
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
        candidate.pxBachFromCharmBaryon(),
        candidate.pyBachFromCharmBaryon(),
        candidate.pzBachFromCharmBaryon(),
        candidate.pxBachFromCasc(),
        candidate.pyBachFromCasc(),
        candidate.pzBachFromCasc(),
        candidate.pxPosV0Dau(),
        candidate.pyPosV0Dau(),
        candidate.pzPosV0Dau(),
        candidate.pxNegV0Dau(),
        candidate.pyNegV0Dau(),
        candidate.pzNegV0Dau(),
        candidate.impactParCascXY(),
        candidate.impactParBachFromCharmBaryonXY(),
        candidate.errImpactParCascXY(),
        candidate.errImpactParBachFromCharmBaryonXY(),
        candidate.invMassLambda(),
        candidate.invMassCascade(),
        candidate.invMassCharmBaryon(),
        candidate.etaV0PosDau(),
        candidate.etaV0NegDau(),
        candidate.etaBachFromCasc(),
        candidate.etaBachFromCharmBaryon(),
        candidate.dcaXYToPvV0Dau0(),
        candidate.dcaXYToPvV0Dau1(),
        candidate.dcaXYToPvCascDau(),
        candidate.dcaCascDau(),
        candidate.dcaV0Dau(),
        candidate.dcaCharmBaryonDau(),
        candidate.errorDecayLengthCharmBaryon(),
        candidate.impactParCascXY() / candidate.errImpactParCascXY(),
        candidate.impactParBachFromCharmBaryonXY() / candidate.errImpactParBachFromCharmBaryonXY(),
        candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().isGlobalTrackWoDCA(),
        candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().itsNCls(),
        candidate.pidTpcInfoStored(),
        candidate.pidTofInfoStored(),
        candidate.tpcNSigmaKaFromCharmBaryon(),
        candidate.tpcNSigmaKaFromCasc(),
        candidate.tpcNSigmaPiFromLambda(),
        candidate.tpcNSigmaPrFromLambda(),
        candidate.tofNSigmaKaFromCharmBaryon(),
        candidate.tofNSigmaKaFromCasc(),
        candidate.tofNSigmaPiFromLambda(),
        candidate.tofNSigmaPrFromLambda(),
        flagMc,
        originMc,
        collisionMatched);
    }
  }

  void processDataLite(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToOmegaK, aod::HfSelToOmegaKa> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidateLite(candidate, -7, RecoDecay::OriginType::None, true);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaKa, processDataLite, "Process data", true);

  void processMcLite(MyEventTable const& collisions, MyTrackTable const&,
                     soa::Join<aod::HfCandToOmegaK, aod::HfSelToOmegaKa, aod::HfToOmegaKMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidateLite(candidate, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaKa, processMcLite, "Process MC", false);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorOmegac0ToOmegaKa>(cfgc)};
}
