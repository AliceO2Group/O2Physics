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

/// \file treeCreatorToXiPiQa.cxx
/// \brief Writer of the omegac0 or xic0 to Xi Pi candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Jinhyun Park <jinhyun.park@cern.ch>, Pusan National University
/// \author Krista Smith <krista.lizbeth.smith@cern.ch>, Pusan National University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;

// SV Reco method
enum {
  DCAFITTER = 0,
  KFPARTICLE
};

// Table size
enum {
  FULL = 0,
  LITE
};

namespace o2::aod
{
namespace full
{
// collision info
DECLARE_SOA_COLUMN(IsEventSel8, isEventSel8, bool);
DECLARE_SOA_COLUMN(IsEventSelZ, isEventSelZ, bool);
DECLARE_SOA_COLUMN(Centrality, centrality, float);
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
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(CosPACharmBaryon, cosPACharmBaryon, float);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, float);
DECLARE_SOA_COLUMN(CosPAXYV0, cosPAXYV0, float);
DECLARE_SOA_COLUMN(CosPAXYCharmBaryon, cosPAXYCharmBaryon, float);
DECLARE_SOA_COLUMN(CosPAXYCasc, cosPAXYCasc, float);
DECLARE_SOA_COLUMN(CTauOmegac, cTauOmegac, float);
DECLARE_SOA_COLUMN(CTauCascade, cTauCascade, float);
DECLARE_SOA_COLUMN(CTauV0, cTauV0, float);
DECLARE_SOA_COLUMN(CTauXic, cTauXic, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaPiFromCasc, etaPiFromCasc, float);
DECLARE_SOA_COLUMN(EtaPiFromCharmBaryon, etaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCharmBaryon, etaCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, float);
DECLARE_SOA_COLUMN(EtaV0, etaV0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau0, dcaZToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaZToPvV0Dau1, dcaZToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaZToPvCascDau, dcaZToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(DecLenCharmBaryon, decLenCharmBaryon, float);
DECLARE_SOA_COLUMN(DecLenCascade, decLenCascade, float);
DECLARE_SOA_COLUMN(DecLenV0, decLenV0, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthXYCharmBaryon, errorDecayLengthXYCharmBaryon, float);
DECLARE_SOA_COLUMN(NormImpParCascade, normImpParCascade, double);
DECLARE_SOA_COLUMN(NormImpParPiFromCharmBar, normImpParPiFromCharmBar, double);
DECLARE_SOA_COLUMN(NormDecayLenCharmBar, normDecayLenCharmBar, double);
DECLARE_SOA_COLUMN(IsPionGlbTrkWoDca, isPionGlbTrkWoDca, bool);
DECLARE_SOA_COLUMN(PionItsNCls, pionItsNCls, uint8_t);
DECLARE_SOA_COLUMN(NTpcRowsPion, nTpcRowsPion, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsPiFromCasc, nTpcRowsPiFromCasc, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsPosV0Dau, nTpcRowsPosV0Dau, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsNegV0Dau, nTpcRowsNegV0Dau, int16_t);
// from creator KF
DECLARE_SOA_COLUMN(KfDcaXYPiFromXic, kfDcaXYPiFromXic, float);
DECLARE_SOA_COLUMN(KfDcaXYCascToPv, kfDcaXYCascToPv, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCasc, chi2GeoCasc, float);
DECLARE_SOA_COLUMN(Chi2GeoXic, chi2GeoXic, float);
DECLARE_SOA_COLUMN(Chi2MassV0, chi2MassV0, float);
DECLARE_SOA_COLUMN(Chi2MassCasc, chi2MassCasc, float);
DECLARE_SOA_COLUMN(V0ldl, v0ldl, float);
DECLARE_SOA_COLUMN(Cascldl, cascldl, float);
DECLARE_SOA_COLUMN(Xicldl, xicldl, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToPv, chi2TopoCascToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromXicToPv, chi2TopoPiFromXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoXicToPv, chi2TopoXicToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCasc, chi2TopoV0ToCasc, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToXic, chi2TopoCascToXic, float);
DECLARE_SOA_COLUMN(DecayLenXYLambda, decayLenXYLambda, float);
DECLARE_SOA_COLUMN(DecayLenXYCasc, decayLenXYCasc, float);
DECLARE_SOA_COLUMN(DecayLenXYXic, decayLenXYXic, float);
DECLARE_SOA_COLUMN(CosPaV0ToCasc, cosPaV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaV0ToPv, cosPaV0ToPv, float);
DECLARE_SOA_COLUMN(CosPaCascToXic, cosPaCascToXic, float);
DECLARE_SOA_COLUMN(CosPaCascToPv, cosPaCascToPv, float);
DECLARE_SOA_COLUMN(CosPaXicToPv, cosPaXicToPv, float);
DECLARE_SOA_COLUMN(KfRapXic, kfRapXic, float);
DECLARE_SOA_COLUMN(KfptPiFromXic, kfptPiFromXic, float);
DECLARE_SOA_COLUMN(KfptXic, kfptXic, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromXic, cosThetaStarPiFromXic, float);
DECLARE_SOA_COLUMN(CtXic, ctXic, float);
DECLARE_SOA_COLUMN(EtaXic, etaXic, float);
DECLARE_SOA_COLUMN(V0Ndf, v0Ndf, float);
DECLARE_SOA_COLUMN(CascNdf, cascNdf, float);
DECLARE_SOA_COLUMN(XicNdf, xicNdf, float);
DECLARE_SOA_COLUMN(MassV0Ndf, massV0Ndf, float);
DECLARE_SOA_COLUMN(MassCascNdf, massCascNdf, float);
DECLARE_SOA_COLUMN(V0Chi2OverNdf, v0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(CascChi2OverNdf, cascChi2OverNdf, float);
DECLARE_SOA_COLUMN(XicChi2OverNdf, xicChi2OverNdf, float);
DECLARE_SOA_COLUMN(MassV0Chi2OverNdf, massV0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(MassCascChi2OverNdf, massCascChi2OverNdf, float);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
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

DECLARE_SOA_TABLE(HfToXiPiEvs, "AOD", "HFTOXIPIEV",
                  full::IsEventSel8, full::IsEventSelZ);

DECLARE_SOA_TABLE(HfToXiPiFulls, "AOD", "HFTOXIPIFULL",
                  full::XPv, full::YPv, full::ZPv, full::Centrality, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
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
                  full::NormImpParCascade, full::NormImpParPiFromCharmBar, full::NormDecayLenCharmBar, full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::NTpcRowsPion, full::NTpcRowsPiFromCasc, full::NTpcRowsPosV0Dau, full::NTpcRowsNegV0Dau,
                  full::StatusPidLambda, full::StatusPidCascade, full::StatusPidCharmBaryon,
                  full::StatusInvMassLambda, full::StatusInvMassCascade, full::StatusInvMassCharmBaryon, full::ResultSelections,
                  full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::DebugMcRec, full::OriginRec, full::CollisionMatched);

DECLARE_SOA_TABLE(HfToXiPiLites, "AOD", "HFTOXIPILITE",
                  full::XPv, full::YPv, full::ZPv, full::Centrality, collision::NumContrib, collision::Chi2,
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
                  full::ErrorDecayLengthCharmBaryon, full::NormImpParCascade, full::NormImpParPiFromCharmBar,
                  full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::NTpcRowsPion, full::NTpcRowsPiFromCasc, full::NTpcRowsPosV0Dau, full::NTpcRowsNegV0Dau,
                  full::StatusPidLambda, full::StatusPidCascade, full::StatusPidCharmBaryon,
                  full::StatusInvMassLambda, full::StatusInvMassCascade, full::StatusInvMassCharmBaryon, full::ResultSelections,
                  full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginRec, full::CollisionMatched);

DECLARE_SOA_TABLE(HfKfXicFulls, "AOD", "HFKFXICFULL",
                  full::Centrality,
                  // full::StatusPidLambda, full::StatusPidCascade, full::StatusPidCharmBaryon,
                  // full::StatusInvMassLambda, full::StatusInvMassCascade, full::StatusInvMassCharmBaryon,
                  full::ResultSelections,
                  full::TpcNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TofNSigmaPiFromCasc,
                  full::TpcNSigmaPiFromLambda, full::TofNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda, full::TofNSigmaPrFromLambda,
                  full::KfDcaXYPiFromXic, full::DcaCascDau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoXic,
                  full::Chi2MassV0, full::Chi2MassCasc,
                  full::V0ldl, full::Cascldl, // full::Xicldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv, full::Chi2TopoPiFromXicToPv, full::Chi2TopoXicToPv,
                  full::Chi2TopoV0ToCasc, full::Chi2TopoCascToXic,
                  full::DecayLenXYLambda, full::DecayLenXYCasc, full::DecayLenXYXic,
                  full::CosPaV0ToCasc, full::CosPaV0ToPv, full::CosPaCascToXic, full::CosPaCascToPv, // full::CosPaXicToPv,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfRapXic, // full::KfptPiFromXic, full::KfptXic,
                  full::CosThetaStarPiFromXic, full::CtXic, full::EtaXic,
                  full::V0Ndf, full::CascNdf, full::XicNdf,
                  full::MassV0Ndf, full::MassCascNdf,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::XicChi2OverNdf,
                  full::MassV0Chi2OverNdf, full::MassCascChi2OverNdf,
                  full::FlagMcMatchRec, full::DebugMcRec, full::OriginRec, full::CollisionMatched);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorToXiPiQa {

  Produces<o2::aod::HfToXiPiFulls> rowCandidateFull;
  Produces<o2::aod::HfToXiPiLites> rowCandidateLite;
  Produces<o2::aod::HfKfXicFulls> rowKfCandidate;
  Produces<o2::aod::HfToXiPiEvs> rowEv;

  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};

  using MyTrackTable = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>;
  using MyEventTable = soa::Join<aod::Collisions, aod::EvSels>;
  using MyEventTableWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using MyEventTableWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using MyEventTableWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;

  void init(InitContext const&)
  {
    if ((doprocessMcLiteXic0 && doprocessMcLiteOmegac0) || (doprocessMcFullXic0 && doprocessMcFullOmegac0)) {
      LOGF(fatal, "Both Xic0 and Omegac0 MC processes enabled, please choose ONLY one!");
    }
  }

  //////////////////////////////////////////////////////
  //                                                  //
  //       Fill functions to fill in the tables       //
  //                                                  //
  //////////////////////////////////////////////////////

  template <bool useCentrality, typename T>
  void fillEvent(const T& collision, float cutZPv)
  {
    rowEv(collision.sel8(), std::abs(collision.posZ()) < cutZPv);
  }

  template <int svReco, int tableSize, bool useCentrality, typename MyEventTableType, typename T>
  void fillCandidate(const T& candidate, int8_t flagMc, int8_t debugMc, int8_t originMc, bool collisionMatched)
  {
    // save all candidate information
    float centrality = -999.f;
    if constexpr (useCentrality) {
      const auto& collision = candidate.template collision_as<MyEventTableType>();
      centrality = o2::hf_centrality::getCentralityColl(collision);
    }

    if constexpr (svReco == DCAFITTER) {
      if constexpr (tableSize == LITE) {
        rowCandidateLite(candidate.xPv(),
                         candidate.yPv(),
                         candidate.zPv(),
                         centrality,
                         candidate.template collision_as<MyEventTableType>().numContrib(),
                         candidate.template collision_as<MyEventTableType>().chi2(),
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
                         candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template bachelor_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template posTrack_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template negTrack_as<MyTrackTable>().tpcNClsCrossedRows(),
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
                         originMc,
                         collisionMatched);
      } else {
        rowCandidateFull(candidate.xPv(),
                         candidate.yPv(),
                         candidate.zPv(),
                         centrality,
                         candidate.template collision_as<MyEventTableType>().numContrib(),
                         candidate.template collision_as<MyEventTableType>().chi2(),
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
                         candidate.covVtxCharmBaryon0(),
                         candidate.covVtxCharmBaryon3(),
                         candidate.covVtxCharmBaryon5(),
                         candidate.pxCharmBaryon(),
                         candidate.pyCharmBaryon(),
                         candidate.pzCharmBaryon(),
                         candidate.pxCasc(),
                         candidate.pyCasc(),
                         candidate.pzCasc(),
                         candidate.pxBachFromCharmBaryon(),
                         candidate.pyBachFromCharmBaryon(),
                         candidate.pzBachFromCharmBaryon(),
                         candidate.pxLambda(),
                         candidate.pyLambda(),
                         candidate.pzLambda(),
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
                         candidate.impactParCascZ(),
                         candidate.impactParBachFromCharmBaryonZ(),
                         candidate.errImpactParCascXY(),
                         candidate.errImpactParBachFromCharmBaryonXY(),
                         candidate.invMassLambda(),
                         candidate.invMassCascade(),
                         candidate.invMassCharmBaryon(),
                         candidate.cosPAV0(),
                         candidate.cosPACharmBaryon(),
                         candidate.cosPACasc(),
                         candidate.cosPAXYV0(),
                         candidate.cosPAXYCharmBaryon(),
                         candidate.cosPAXYCasc(),
                         candidate.cTauOmegac(),
                         candidate.cTauCascade(),
                         candidate.cTauV0(),
                         candidate.cTauXic(),
                         candidate.etaV0PosDau(),
                         candidate.etaV0NegDau(),
                         candidate.etaBachFromCasc(),
                         candidate.etaBachFromCharmBaryon(),
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
                         candidate.impactParCascXY() / candidate.errImpactParCascXY(),
                         candidate.impactParBachFromCharmBaryonXY() / candidate.errImpactParBachFromCharmBaryonXY(),
                         candidate.decLenCharmBaryon() / candidate.errorDecayLengthCharmBaryon(),
                         candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().isGlobalTrackWoDCA(),
                         candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().itsNCls(),
                         candidate.template bachelorFromCharmBaryon_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template bachelor_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template posTrack_as<MyTrackTable>().tpcNClsCrossedRows(),
                         candidate.template negTrack_as<MyTrackTable>().tpcNClsCrossedRows(),
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
                         originMc,
                         collisionMatched);
      }
    } else {
      if constexpr (tableSize == LITE) {
        // currently, no lite sized table for KFParticle
      } else {
        rowKfCandidate(centrality,
                       candidate.resultSelections(),
                       candidate.tpcNSigmaPiFromCharmBaryon(),
                       candidate.tofNSigmaPiFromCharmBaryon(),
                       candidate.tpcNSigmaPiFromCasc(),
                       candidate.tofNSigmaPiFromCasc(),
                       candidate.tpcNSigmaPiFromLambda(),
                       candidate.tofNSigmaPiFromLambda(),
                       candidate.tpcNSigmaPrFromLambda(),
                       candidate.tofNSigmaPrFromLambda(),
                       candidate.kfDcaXYPiFromXic(),
                       candidate.dcaCascDau(),
                       candidate.dcaCharmBaryonDau(),
                       candidate.kfDcaXYCascToPv(),
                       candidate.dcaXYToPvV0Dau0(),
                       candidate.dcaXYToPvV0Dau1(),
                       candidate.dcaXYToPvCascDau(),
                       candidate.chi2GeoV0(),
                       candidate.chi2GeoCasc(),
                       candidate.chi2GeoXic(),
                       candidate.chi2MassV0(),
                       candidate.chi2MassCasc(),
                       candidate.v0ldl(),
                       candidate.cascldl(),
                       candidate.chi2TopoV0ToPv(),
                       candidate.chi2TopoCascToPv(),
                       candidate.chi2TopoPiFromXicToPv(),
                       candidate.chi2TopoXicToPv(),
                       candidate.chi2TopoV0ToCasc(),
                       candidate.chi2TopoCascToXic(),
                       candidate.decayLenXYLambda(),
                       candidate.decayLenXYCasc(),
                       candidate.decayLenXYXic(),
                       candidate.cosPaV0ToCasc(),
                       candidate.cosPAV0(),
                       candidate.cosPaCascToXic(),
                       candidate.cosPACasc(),
                       candidate.invMassLambda(),
                       candidate.invMassCascade(),
                       candidate.invMassCharmBaryon(),
                       candidate.kfRapXic(),
                       candidate.cosThetaStarPiFromXic(),
                       candidate.cTauXic(),
                       candidate.etaCharmBaryon(),
                       candidate.v0Ndf(),
                       candidate.cascNdf(),
                       candidate.xicNdf(),
                       candidate.massV0Ndf(),
                       candidate.massCascNdf(),
                       candidate.v0Chi2OverNdf(),
                       candidate.cascChi2OverNdf(),
                       candidate.xicChi2OverNdf(),
                       candidate.massV0Chi2OverNdf(),
                       candidate.massCascChi2OverNdf(),
                       flagMc,
                       debugMc,
                       originMc,
                       collisionMatched);
      }
    }
  }

  ////////////////////////////////////
  //                                //
  //       Process functions        //
  //                                //
  ////////////////////////////////////

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~Data with DCAFitter~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processDataFull(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTable>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processDataLite(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTable>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processDataLiteWithFT0M(MyEventTableWithFT0M const& collisions, MyTrackTable const&,
                               soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0M>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processDataLiteWithFT0C(MyEventTableWithFT0C const& collisions, MyTrackTable const&,
                               soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0C>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processDataLiteWithNTracksPV(MyEventTableWithNTracksPV const& collisions, MyTrackTable const&,
                                    soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithNTracksPV>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processDataFull, "Process data with full information w/o centrality", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processDataLite, "Process data and produce lite table version", true);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processDataLiteWithFT0M, "Process data and produce lite table version with FT0M", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processDataLiteWithFT0C, "Process data and produce lite table version with FT0C", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processDataLiteWithNTracksPV, "Process data and produce lite table version with NTracksPV", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~Data with KFParticle~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  void processKfData(MyEventTable const& collisions, MyTrackTable const&,
                     soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, false, MyEventTable>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processKfDataWithFT0M(MyEventTableWithFT0M const& collisions, MyTrackTable const&,
                             soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, false, MyEventTableWithFT0M>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processKfDataWithFT0C(MyEventTableWithFT0C const& collisions, MyTrackTable const&,
                             soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithFT0C>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  void processKfDataWithNTracksPV(MyEventTableWithNTracksPV const& collisions, MyTrackTable const&,
                                  soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithNTracksPV>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfData, "Process KF data, no cent", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithFT0M, "Process KF data, with FT0M", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithFT0C, "Process KF data, with FT0C", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithNTracksPV, "Process KF data, with NTracksPV", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~MC with DCAFitter~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processMcFullXic0(MyEventTable const& collisions, MyTrackTable const&,
                         soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTable>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcFullOmegac0(MyEventTable const& collisions, MyTrackTable const&,
                            soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfOmegacToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTable>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcLiteXic0(MyEventTable const& collisions, MyTrackTable const&,
                         soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTable>(candidate, candidate.flagMcMatchRec(), -7, candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcLiteXic0WithFT0C(MyEventTableWithFT0C const& collisions, MyTrackTable const&,
                                 soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0C>(candidate, candidate.flagMcMatchRec(), -7, candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcLiteXic0WithFT0M(MyEventTableWithFT0M const& collisions, MyTrackTable const&,
                                 soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0M>(candidate, candidate.flagMcMatchRec(), -7, candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcLiteXic0WithNTracksPV(MyEventTableWithNTracksPV const& collisions, MyTrackTable const&,
                                      soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithNTracksPV>(candidate, candidate.flagMcMatchRec(), -7, candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processMcLiteOmegac0(MyEventTable const& collisions, MyTrackTable const&,
                            soa::Join<aod::HfCandToXiPi, aod::HfSelToXiPi, aod::HfOmegacToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTable>(candidate, candidate.flagMcMatchRec(), -7, candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcFullXic0, "Process MC with full information for xic0 w/o centrality", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcFullOmegac0, "Process MC with full information for omegac0", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcLiteXic0, "Process MC and produce lite table version for xic0", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcLiteXic0WithFT0C, "Process MC and produce lite table version for Xic0 with FT0C", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcLiteXic0WithFT0M, "Process MC and produce lite table version for Xic0 with FT0M", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcLiteXic0WithNTracksPV, "Process MC and produce lite table version for Xic0 with NTracksPV", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processMcLiteOmegac0, "Process MC and produce lite table version for omegac0", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~MC with KFParticle~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  void processKfMcXic0(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, false, MyEventTable>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processKfMcXic0WithFT0C(MyEventTableWithFT0C const& collisions, MyTrackTable const&,
                               soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithFT0C>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processKfMcXic0WithFT0M(MyEventTableWithFT0M const& collisions, MyTrackTable const&,
                               soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithFT0M>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  void processKfMcXic0WithNTracksPV(MyEventTableWithNTracksPV const& collisions, MyTrackTable const&,
                                    soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision, zPvCut);
    }

    // Filling candidate table
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithNTracksPV>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originMcRec(), candidate.collisionMatched());
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfMcXic0, "Process MC with information for xic0", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfMcXic0WithFT0C, "Process MC with information for xic0 at FT0C", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfMcXic0WithFT0M, "Process MC with information for xic0 at FT0M", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfMcXic0WithNTracksPV, "Process MC with information for xic0 at Ntrack", false);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorToXiPiQa>(cfgc)};
}
