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

/// \file treeCreatorToXiPi.cxx
/// \brief Writer of the omegac0 or xic0 to Xi Pi candidates in the form of flat tables to be stored in TTrees.
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
DECLARE_SOA_COLUMN(Chi2PCACharmBaryon, chi2PCACharmBaryon, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryonXX, covVtxCharmBaryonXX, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryonYY, covVtxCharmBaryonYY, float);
DECLARE_SOA_COLUMN(CovVtxCharmBaryonZZ, covVtxCharmBaryonZZ, float);
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PxPiFromCharmBaryon, pxPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyPiFromCharmBaryon, pyPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzPiFromCharmBaryon, pzPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PxLambda, pxLambda, float);
DECLARE_SOA_COLUMN(PyLambda, pyLambda, float);
DECLARE_SOA_COLUMN(PzLambda, pzLambda, float);
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
DECLARE_SOA_COLUMN(ImpactParCascZ, impactParCascZ, float);
DECLARE_SOA_COLUMN(ImpactParPiFromCharmBaryonZ, impactParPiFromCharmBaryonZ, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParPiFromCharmBaryonXY, errImpactParPiFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, double);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, double);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, double);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, double);
DECLARE_SOA_COLUMN(CosPACharmBaryon, cosPACharmBaryon, double);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, double);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, double);
DECLARE_SOA_COLUMN(CosPAXYCharmBaryon, cosPAXYCharmBaryon, double);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, double);
DECLARE_SOA_COLUMN(CTauOmegac, ctauOmegac, double);
DECLARE_SOA_COLUMN(CTauCascade, ctauCascade, double);
DECLARE_SOA_COLUMN(CTauV0, ctauV0, double);
DECLARE_SOA_COLUMN(CTauXic, ctauXic, double);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, double);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, double);
DECLARE_SOA_COLUMN(EtaPiFromCasc, etaPiFromCasc, double);
DECLARE_SOA_COLUMN(EtaPiFromCharmBaryon, etaPiFromCharmBaryon, double);
DECLARE_SOA_COLUMN(EtaCharmBaryon, etaCharmBaryon, double);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, double);
DECLARE_SOA_COLUMN(EtaV0, etaV0, double);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau0, dcaZToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau1, dcaZToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaZToPvCascDau, dcaZToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(DecLenCharmBaryon, decLenCharmBaryon, double);
DECLARE_SOA_COLUMN(DecLenCascade, decLenCascade, double);
DECLARE_SOA_COLUMN(DecLenV0, decLenV0, double);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYCharmBaryon, errorDecayLengthXYCharmBaryon, float);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
// from selector
DECLARE_SOA_COLUMN(StatusPidLambda, statusPidLambda, bool);
DECLARE_SOA_COLUMN(StatusPidCascade, statusPidCascade, bool);
DECLARE_SOA_COLUMN(StatusPidCharmBaryon, statusPidCharmBaryon, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassCharmBaryon, statusInvMassCharmBaryon, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
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

DECLARE_SOA_TABLE(HfToXiPiFulls, "AOD", "HFTOXIPIFULL",
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay, full::Chi2PCACharmBaryon,
                  full::CovVtxCharmBaryonXX, full::CovVtxCharmBaryonYY, full::CovVtxCharmBaryonZZ,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxCasc, full::PyCasc, full::PzCasc,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxLambda, full::PyLambda, full::PzLambda,
                  full::PxPiFromCasc, full::PyPiFromCasc, full::PzPiFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ImpactParCascZ, full::ImpactParPiFromCharmBaryonZ,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::CosPAV0, full::CosPACharmBaryon, full::CosPACasc, full::CosPAXYV0, full::CosPAXYCharmBaryon, full::CosPAXYCasc,
                  full::CTauOmegac, full::CTauCascade, full::CTauV0, full::CTauXic,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaPiFromCasc, full::EtaPiFromCharmBaryon,
                  full::EtaCharmBaryon, full::EtaCascade, full::EtaV0,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaZToPvV0Dau0, full::DcaZToPvV0Dau1, full::DcaZToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::DecLenCharmBaryon, full::DecLenCascade, full::DecLenV0, full::ErrorDecayLengthCharmBaryon, full::ErrorDecayLengthXYCharmBaryon,
                  full::StatusPidLambda, full::StatusPidCascade, full::StatusPidCharmBaryon,
                  full::StatusInvMassLambda, full::StatusInvMassCascade, full::StatusInvMassCharmBaryon, full::ResultSelections, full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::DebugMcRec, full::OriginRec);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorToXiPi {

  Produces<o2::aod::HfToXiPiFulls> rowCandidateFull;

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillCandidate(const T& candidate, int8_t flagMc, int8_t debugMc, int8_t originMc)
  {
    rowCandidateFull(
      candidate.xPv(),
      candidate.yPv(),
      candidate.zPv(),
      candidate.collision().numContrib(),
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
      candidate.chi2PCACharmBaryon(),
      candidate.covVtxCharmBaryon0(),
      candidate.covVtxCharmBaryon3(),
      candidate.covVtxCharmBaryon5(),
      candidate.pxCharmBaryon(),
      candidate.pyCharmBaryon(),
      candidate.pzCharmBaryon(),
      candidate.pxCasc(),
      candidate.pyCasc(),
      candidate.pzCasc(),
      candidate.pxPiFromCharmBaryon(),
      candidate.pyPiFromCharmBaryon(),
      candidate.pzPiFromCharmBaryon(),
      candidate.pxLambda(),
      candidate.pyLambda(),
      candidate.pzLambda(),
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
      candidate.impactParCascZ(),
      candidate.impactParPiFromCharmBaryonZ(),
      candidate.errImpactParCascXY(),
      candidate.errImpactParPiFromCharmBaryonXY(),
      candidate.invMassLambda(),
      candidate.invMassCascade(),
      candidate.invMassCharmBaryon(),
      candidate.cosPAV0(),
      candidate.cosPACharmBaryon(),
      candidate.cosPACasc(),
      candidate.cosPAXYV0(),
      candidate.cosPAXYCharmBaryon(),
      candidate.cosPAXYCasc(),
      candidate.ctauOmegac(),
      candidate.ctauCascade(),
      candidate.ctauV0(),
      candidate.ctauXic(),
      candidate.etaV0PosDau(),
      candidate.etaV0NegDau(),
      candidate.etaPiFromCasc(),
      candidate.etaPiFromCharmBaryon(),
      candidate.etaCharmBaryon(),
      candidate.etaCascade(),
      candidate.etaV0(),
      candidate.dcaXYToPvV0Dau0(),
      candidate.dcaXYToPvV0Dau1(),
      candidate.dcaXYToPvCascDau(),
      candidate.dcaZToPvV0Dau0(),
      candidate.dcaZToPvV0Dau1(),
      candidate.dcaZToPvCascDau(),
      candidate.dcaCascDau(),
      candidate.dcaV0Dau(),
      candidate.dcaCharmBaryonDau(),
      candidate.decLenCharmBaryon(),
      candidate.decLenCascade(),
      candidate.decLenV0(),
      candidate.errorDecayLengthCharmBaryon(),
      candidate.errorDecayLengthXYCharmBaryon(),
      candidate.statusPidLambda(),
      candidate.statusPidCascade(),
      candidate.statusPidCharmBaryon(),
      candidate.statusInvMassLambda(),
      candidate.statusInvMassCascade(),
      candidate.statusInvMassCharmBaryon(),
      candidate.resultSelections(),
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
      debugMc,
      originMc);
  }

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate(candidate, -7, -7, RecoDecay::OriginType::None);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorToXiPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfToXiPiMCRec> const& candidates)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originRec());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorToXiPi, processMc, "Process MC", false);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorToXiPi>(cfgc)};
}
