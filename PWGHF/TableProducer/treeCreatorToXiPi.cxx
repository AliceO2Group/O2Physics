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
DECLARE_SOA_COLUMN(CovVtxCharmXX, covVtxCharmXX, float);
DECLARE_SOA_COLUMN(CovVtxCharmYY, covVtxCharmYY, float);
DECLARE_SOA_COLUMN(CovVtxCharmZZ, covVtxCharmZZ, float);
DECLARE_SOA_COLUMN(PxCharm, pxCharm, float);
DECLARE_SOA_COLUMN(PyCharm, pyCharm, float);
DECLARE_SOA_COLUMN(PzCharm, pzCharm, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PxPiFromCharm, pxPiFromCharm, float);
DECLARE_SOA_COLUMN(PyPiFromCharm, pyPiFromCharm, float);
DECLARE_SOA_COLUMN(PzPiFromCharm, pzPiFromCharm, float);
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
DECLARE_SOA_COLUMN(DcaXYToPvPiFromCharm, dcaXYToPvPiFromCharm, float);
DECLARE_SOA_COLUMN(ImpactParCascZ, impactParCascZ, float);
DECLARE_SOA_COLUMN(DcaZToPvPiFromCharm, dcaZToPvPiFromCharm, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, double);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, double);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, double);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, double);
DECLARE_SOA_COLUMN(CosPACharm, cosPACharm, double);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, double);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, double);
DECLARE_SOA_COLUMN(CosPAXYCharm, cosPAXYCharm, double);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, double);
DECLARE_SOA_COLUMN(CTauOmegac, ctauOmegac, double);
DECLARE_SOA_COLUMN(CTauCascade, ctauCascade, double);
DECLARE_SOA_COLUMN(CTauV0, ctauV0, double);
DECLARE_SOA_COLUMN(CTauXic, ctauXic, double);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, double);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, double);
DECLARE_SOA_COLUMN(EtaPiFromCasc, etaPiFromCasc, double);
DECLARE_SOA_COLUMN(EtaPiFromCharm, etaPiFromCharm, double);
DECLARE_SOA_COLUMN(EtaCharm, etaCharm, double);
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
DECLARE_SOA_COLUMN(DcaCharmDau, dcaCharmDau, float);
DECLARE_SOA_COLUMN(DecLenCharm, decLenCharm, double);
DECLARE_SOA_COLUMN(DecLenCascade, decLenCascade, double);
DECLARE_SOA_COLUMN(DecLenV0, decLenV0, double);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharm, errorDecayLengthCharm, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYCharm, errorDecayLengthXYCharm, float);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
// from selector
DECLARE_SOA_COLUMN(StatusPidLambda, statusPidLambda, bool);
DECLARE_SOA_COLUMN(StatusPidCascade, statusPidCascade, bool);
DECLARE_SOA_COLUMN(StatusPidCharm, statusPidCharm, bool);
DECLARE_SOA_COLUMN(StatusInvMassLambda, statusInvMassLambda, bool);
DECLARE_SOA_COLUMN(StatusInvMassCascade, statusInvMassCascade, bool);
DECLARE_SOA_COLUMN(StatusInvMassCharmBaryon, statusInvMassCharmBaryon, bool);
DECLARE_SOA_COLUMN(ResultSelections, resultSelections, bool);
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharm, tpcNSigmaPiFromCharm, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCasc, tpcNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharm, tofNSigmaPiFromCharm, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCasc, tofNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);

} // namespace full

DECLARE_SOA_TABLE(HfToXiPiFulls, "AOD", "HFTOXIPIFULL",
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay, full::Chi2PCACharm,
                  full::CovVtxCharmXX, full::CovVtxCharmYY, full::CovVtxCharmZZ,
                  full::PxCharm, full::PyCharm, full::PzCharm,
                  full::PxCasc, full::PyCasc, full::PzCasc,
                  full::PxPiFromCharm, full::PyPiFromCharm, full::PzPiFromCharm,
                  full::PxLambda, full::PyLambda, full::PzLambda,
                  full::PxPiFromCasc, full::PyPiFromCasc, full::PzPiFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::DcaXYToPvPiFromCharm,
                  full::ImpactParCascZ, full::DcaZToPvPiFromCharm,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::CosPAV0, full::CosPACharm, full::CosPACasc, full::CosPAXYV0, full::CosPAXYCharm, full::CosPAXYCasc,
                  full::CTauOmegac, full::CTauCascade, full::CTauV0, full::CTauXic,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaPiFromCasc, full::EtaPiFromCharm,
                  full::EtaCharm, full::EtaCascade, full::EtaV0,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaZToPvV0Dau0, full::DcaZToPvV0Dau1, full::DcaZToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmDau,
                  full::DecLenCharm, full::DecLenCascade, full::DecLenV0, full::ErrorDecayLengthCharm, full::ErrorDecayLengthXYCharm,
                  full::StatusPidLambda, full::StatusPidCascade, full::StatusPidCharm,
                  full::StatusInvMassLambda, full::StatusInvMassCascade, full::StatusInvMassCharmBaryon, full::ResultSelections, full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharm, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharm, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::DebugMcRec);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorToXiPi {

  Produces<o2::aod::HfToXiPiFulls> rowCandidateFull;

  void init(InitContext const&)
  {
  }

  template <typename T>
  void fillCandidate(const T& candidate, int8_t flagMc, int8_t debugMc)
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
      candidate.chi2PCACharm(),
      candidate.covVtxCharm0(),
      candidate.covVtxCharm3(),
      candidate.covVtxCharm5(),
      candidate.pxCharm(),
      candidate.pyCharm(),
      candidate.pzCharm(),
      candidate.pxCasc(),
      candidate.pyCasc(),
      candidate.pzCasc(),
      candidate.pxPiFromCharm(),
      candidate.pyPiFromCharm(),
      candidate.pzPiFromCharm(),
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
      candidate.dcaXYToPvPiFromCharm(),
      candidate.impactParCascZ(),
      candidate.dcaZToPvPiFromCharm(),
      candidate.invMassLambda(),
      candidate.invMassCascade(),
      candidate.invMassCharmBaryon(),
      candidate.cosPAV0(),
      candidate.cosPACharm(),
      candidate.cosPACasc(),
      candidate.cosPAXYV0(),
      candidate.cosPAXYCharm(),
      candidate.cosPAXYCasc(),
      candidate.ctauOmegac(),
      candidate.ctauCascade(),
      candidate.ctauV0(),
      candidate.ctauXic(),
      candidate.etaV0PosDau(),
      candidate.etaV0NegDau(),
      candidate.etaPiFromCasc(),
      candidate.etaPiFromCharm(),
      candidate.etaCharm(),
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
      candidate.dcaCharmDau(),
      candidate.decLenCharm(),
      candidate.decLenCascade(),
      candidate.decLenV0(),
      candidate.errorDecayLengthCharm(),
      candidate.errorDecayLengthXYCharm(),
      candidate.statusPidLambda(),
      candidate.statusPidCascade(),
      candidate.statusPidCharm(),
      candidate.statusInvMassLambda(),
      candidate.statusInvMassCascade(),
      candidate.statusInvMassCharmBaryon(),
      candidate.resultSelections(),
      candidate.pidTpcInfoStored(),
      candidate.pidTofInfoStored(),
      candidate.tpcNSigmaPiFromCharm(),
      candidate.tpcNSigmaPiFromCasc(),
      candidate.tpcNSigmaPiFromLambda(),
      candidate.tpcNSigmaPrFromLambda(),
      candidate.tofNSigmaPiFromCharm(),
      candidate.tofNSigmaPiFromCasc(),
      candidate.tofNSigmaPiFromLambda(),
      candidate.tofNSigmaPrFromLambda(),
      flagMc,
      debugMc);
  }

  void processData(aod::Collisions const& collisions,
                   soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate(candidate, -7, -7);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorToXiPi, processData, "Process data tree writer", true);

  void processMc(aod::Collisions const& collisions,
                 soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfToXiPiMCRec> const& candidates)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorToXiPi, processMc, "Process MC tree writer", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorToXiPi>(cfgc)};
}
