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
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/CCDB/EventSelectionParams.h"
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

#include <Rtypes.h>

#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
DECLARE_SOA_COLUMN(PtCharmBaryon, ptCharmBaryon, float);
DECLARE_SOA_COLUMN(PtPiFromCharmBaryon, ptPiFromCharmBaryon, float);
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
DECLARE_SOA_COLUMN(PhiCharmBaryon, phiCharmBaryon, float);
DECLARE_SOA_COLUMN(YCharmBaryon, yCharmBaryon, float);
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
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCasc, tpcNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCasc, tofNSigmaPiFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
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
// MC
DECLARE_SOA_COLUMN(ParticlePdg, particlePdg, int);
DECLARE_SOA_COLUMN(PtGenB, ptGenB, float);
DECLARE_SOA_COLUMN(NContribMax, nContribMax, int);
DECLARE_SOA_COLUMN(NRecoColl, nRecoColl, int);
DECLARE_SOA_COLUMN(HasRecoCollTVX, hasRecoCollTVX, bool);
DECLARE_SOA_COLUMN(HasRecoCollSel8, hasRecoCollSel8, bool);
DECLARE_SOA_COLUMN(HasRecoCollSel8Zvtx, hasRecoCollSel8Zvtx, bool);
} // namespace full

DECLARE_SOA_TABLE(HfToXiPiEvs, "AOD", "HFTOXIPIEV",
                  full::IsEventSel8, full::IsEventSelZ);

DECLARE_SOA_TABLE(HfToXiPiFulls, "AOD", "HFTOXIPIFULL",
                  full::XPv, full::YPv, full::ZPv, full::Centrality, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
                  full::PtCharmBaryon, full::PtPiFromCharmBaryon,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxCasc, full::PyCasc, full::PzCasc,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxPiFromCasc, full::PyPiFromCasc, full::PzPiFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ImpactParCascZ, full::ImpactParPiFromCharmBaryonZ,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::CosPAV0, full::CosPACharmBaryon, full::CosPACasc, full::CosPAXYV0, full::CosPAXYCharmBaryon, full::CosPAXYCasc,
                  full::CTauOmegac,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaPiFromCasc, full::EtaPiFromCharmBaryon,
                  full::EtaCharmBaryon, full::EtaCascade, full::EtaV0,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaZToPvV0Dau0, full::DcaZToPvV0Dau1, full::DcaZToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::DecLenCharmBaryon, full::DecLenCascade, full::DecLenV0, full::ErrorDecayLengthCharmBaryon, full::ErrorDecayLengthXYCharmBaryon,
                  full::NormImpParCascade, full::NormImpParPiFromCharmBar, full::NormDecayLenCharmBar, full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::NTpcRowsPion, full::NTpcRowsPiFromCasc, full::NTpcRowsPosV0Dau, full::NTpcRowsNegV0Dau,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginRec);

DECLARE_SOA_TABLE(HfToXiPiLites, "AOD", "HFTOXIPILITE",
                  full::XPv, full::YPv, full::ZPv, full::Centrality, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
                  full::PtCharmBaryon, full::PtPiFromCharmBaryon,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxPiFromCasc, full::PyPiFromCasc, full::PzPiFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::CosPAV0, full::CosPACharmBaryon, full::CosPACasc,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaPiFromCasc, full::EtaPiFromCharmBaryon, full::EtaCharmBaryon,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::ErrorDecayLengthCharmBaryon, full::NormImpParCascade, full::NormImpParPiFromCharmBar,
                  full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::NTpcRowsPion, full::NTpcRowsPiFromCasc, full::NTpcRowsPosV0Dau, full::NTpcRowsNegV0Dau,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginRec);

DECLARE_SOA_TABLE(HfKfXicFulls, "AOD", "HFKFXICFULL",
                  full::Centrality,
                  full::TpcNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TofNSigmaPiFromCasc,
                  full::TpcNSigmaPiFromLambda, full::TofNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda, full::TofNSigmaPrFromLambda,
                  full::KfDcaXYPiFromXic, full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoXic,
                  full::V0ldl, full::Cascldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv, full::Chi2TopoXicToPv,
                  full::Chi2TopoV0ToCasc, full::Chi2TopoCascToXic,
                  full::DecayLenXYLambda, full::DecayLenXYCasc, full::DecayLenXYXic,
                  full::CosPaV0ToCasc, full::CosPaV0ToPv, full::CosPaCascToXic, full::CosPaCascToPv, // full::CosPaXicToPv,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfRapXic, full::KfptPiFromXic, full::KfptXic,
                  full::CosThetaStarPiFromXic, full::CtXic, full::EtaXic,
                  full::V0Ndf, full::CascNdf, full::XicNdf,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::XicChi2OverNdf,
                  full::FlagMcMatchRec, full::OriginRec);

DECLARE_SOA_TABLE(HfCandToXiPiGen, "AOD", "HFCANDTOXIPIGEN",
                  full::InvMassCharmBaryon,
                  full::PtCharmBaryon,
                  full::EtaCharmBaryon,
                  full::PhiCharmBaryon,
                  full::YCharmBaryon,
                  full::FlagMcMatchRec,
                  full::OriginRec,
                  full::ParticlePdg,
                  full::PtGenB,
                  full::NContribMax,
                  full::NRecoColl,
                  full::HasRecoCollTVX,
                  full::HasRecoCollSel8,
                  full::HasRecoCollSel8Zvtx);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorToXiPiQa {

  Produces<o2::aod::HfToXiPiFulls> rowCandidateFull;
  Produces<o2::aod::HfToXiPiLites> rowCandidateLite;
  Produces<o2::aod::HfKfXicFulls> rowKfCandidate;
  Produces<o2::aod::HfCandToXiPiGen> rowCandidateParticles;
  Produces<o2::aod::HfToXiPiEvs> rowEv;

  Configurable<int> selectionFlag{"selectionFlag", 7, "Selection threshold: 1 skims, 3 topology, 7 PID, 15 ML"};
  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};
  Configurable<int8_t> decayChannelSelection{"decayChannelSelection", hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Xic0ToXiPi, "Decay channel used to match particle information"};
  Configurable<bool> fillGenTable{"fillGenTable", true, "Fill generated MC information if requested"};

  using MyTrackTable = soa::Join<aod::TracksWExtra, aod::TrackSelection, aod::TracksPidPi, aod::TracksPidPr>;
  using MyTrackTableLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr>;
  using MyCandTableDca = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfSelToXiPiQa>>;
  using MyCandTableKf = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfCandCascBachKfQa, aod::HfSelToXiPiQa>>;
  using MyCandTableXic0McDca = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfSelToXiPiQa, aod::HfCandCascBachMCRecQa>>;
  using MyCandTableXic0McKf = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfCandCascBachKfQa, aod::HfSelToXiPiQa, aod::HfCandCascBachMCRecQa>>;
  using MyCandTableOmegac0McDca = MyCandTableXic0McDca;
  using MyEventTable = soa::Join<aod::Collisions, aod::EvSels>;
  using MyEventTableWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using MyEventTableWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using MyEventTableWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;
  using MyEventTableWithMcLabels = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MyEventTableWithMcLabelsFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using MyEventTableWithMcLabelsFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using MyEventTableWithMcLabelsNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentNTPVs>;

  using MatchedGenXiPi = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandCascBachMCGenQa>>;

  using McCollisionsCentFT0Cs = soa::Join<aod::McCollisions, aod::McCentFT0Cs>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_to_xi_pi_qa::isSelToXiPi >= selectionFlag;
  Filter filterRecXiPi = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(decayChannelSelection);
  Filter filterGenXiPi = nabs(aod::hf_cand_mc_flag::flagMcMatchGen) == static_cast<int8_t>(decayChannelSelection);

  PresliceUnsorted<MyEventTableWithMcLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

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
  void fillEvent(const T& collision)
  {
    rowEv(collision.sel8(), std::abs(collision.posZ()) < zPvCut);
  }

  template <int svReco, int tableSize, bool useCentrality, typename MyEventTableType, typename T>
  void fillCandidate(const T& candidate, MyTrackTable const& tracks, MyTrackTableLf const& lfTracks, int8_t flagMc, int8_t originMc)
  {
    float centrality = -999.f;
    if constexpr (useCentrality) {
      const auto& collision = candidate.template collision_as<MyEventTableType>();
      centrality = o2::hf_centrality::getCentralityColl(collision);
    }

    auto trackV0PosDau = lfTracks.rawIteratorAt(candidate.posTrackId());
    auto trackV0NegDau = lfTracks.rawIteratorAt(candidate.negTrackId());
    auto trackPiFromCasc = lfTracks.rawIteratorAt(candidate.bachelorId());
    auto trackPiFromCharm = tracks.rawIteratorAt(candidate.prong1Id());

    auto trackPiFromLam = trackV0NegDau;
    auto trackPrFromLam = trackV0PosDau;
    if (trackPiFromCharm.sign() < 0) {
      trackPiFromLam = trackV0PosDau;
      trackPrFromLam = trackV0NegDau;
    }

    if constexpr (svReco == DCAFITTER) {
      if constexpr (tableSize == LITE) {
        rowCandidateLite(candidate.posX(),
                         candidate.posY(),
                         candidate.posZ(),
                         centrality,
                         candidate.template collision_as<MyEventTableType>().numContrib(),
                         candidate.template collision_as<MyEventTableType>().chi2(),
                         candidate.xSecondaryVertex(),
                         candidate.ySecondaryVertex(),
                         candidate.zSecondaryVertex(),
                         candidate.xCascadeVertex(),
                         candidate.yCascadeVertex(),
                         candidate.zCascadeVertex(),
                         candidate.xV0Vertex(),
                         candidate.yV0Vertex(),
                         candidate.zV0Vertex(),
                         trackPiFromCasc.sign(),
                         candidate.pt(),
                         RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1()),
                         candidate.px(),
                         candidate.py(),
                         candidate.pz(),
                         candidate.pxProng1(),
                         candidate.pyProng1(),
                         candidate.pzProng1(),
                         trackPiFromCasc.px(),
                         trackPiFromCasc.py(),
                         trackPiFromCasc.pz(),
                         trackV0PosDau.px(),
                         trackV0PosDau.py(),
                         trackV0PosDau.pz(),
                         trackV0NegDau.px(),
                         trackV0NegDau.py(),
                         trackV0NegDau.pz(),
                         candidate.impactParameter0(),
                         candidate.impactParameter1(),
                         candidate.errorImpactParameter0(),
                         candidate.errorImpactParameter1(),
                         candidate.massLambda(),
                         candidate.massCascade(),
                         HfHelper::invMassXiczeroToXiPi(candidate),
                         candidate.cpaV0(),
                         candidate.cpa(),
                         candidate.cpaCascade(),
                         trackV0PosDau.eta(),
                         trackV0NegDau.eta(),
                         trackPiFromCasc.eta(),
                         trackPiFromCharm.eta(),
                         candidate.eta(),
                         candidate.impactParameterXYV0PosDau(),
                         candidate.impactParameterXYV0NegDau(),
                         candidate.impactParameterXYCascadeBach(),
                         candidate.dcaCascadeDaughters(),
                         candidate.dcaV0Daughters(),
                         candidate.dcaDaughters(),
                         candidate.errorDecayLength(),
                         candidate.impactParameterNormalised0(),
                         candidate.impactParameterNormalised1(),
                         trackPiFromCharm.isGlobalTrackWoDCA(),
                         trackPiFromCharm.itsNCls(),
                         trackPiFromCharm.tpcNClsCrossedRows(),
                         trackPiFromCasc.tpcNClsCrossedRows(),
                         trackV0PosDau.tpcNClsCrossedRows(),
                         trackV0NegDau.tpcNClsCrossedRows(),
                         trackPiFromCharm.tpcNSigmaPi(),
                         trackPiFromCasc.tpcNSigmaPi(),
                         trackPiFromLam.tpcNSigmaPi(),
                         trackPrFromLam.tpcNSigmaPr(),
                         trackPiFromCharm.tofNSigmaPi(),
                         trackPiFromCasc.tofNSigmaPi(),
                         trackPiFromLam.tofNSigmaPi(),
                         trackPrFromLam.tofNSigmaPr(),
                         flagMc,
                         originMc);
      } else {
        rowCandidateFull(candidate.posX(),
                         candidate.posY(),
                         candidate.posZ(),
                         centrality,
                         candidate.template collision_as<MyEventTableType>().numContrib(),
                         candidate.template collision_as<MyEventTableType>().chi2(),
                         candidate.xSecondaryVertex(),
                         candidate.ySecondaryVertex(),
                         candidate.zSecondaryVertex(),
                         candidate.xCascadeVertex(),
                         candidate.yCascadeVertex(),
                         candidate.zCascadeVertex(),
                         candidate.xV0Vertex(),
                         candidate.yV0Vertex(),
                         candidate.zV0Vertex(),
                         trackPiFromCasc.sign(),
                         candidate.pt(),
                         RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1()),
                         candidate.px(),
                         candidate.py(),
                         candidate.pz(),
                         candidate.pxProng0(),
                         candidate.pyProng0(),
                         candidate.pzProng0(),
                         candidate.pxProng1(),
                         candidate.pyProng1(),
                         candidate.pzProng1(),
                         trackPiFromCasc.px(),
                         trackPiFromCasc.py(),
                         trackPiFromCasc.pz(),
                         trackV0PosDau.px(),
                         trackV0PosDau.py(),
                         trackV0PosDau.pz(),
                         trackV0NegDau.px(),
                         trackV0NegDau.py(),
                         trackV0NegDau.pz(),
                         candidate.impactParameter0(),
                         candidate.impactParameter1(),
                         candidate.impactParameterZ0(),
                         candidate.impactParameterZ1(),
                         candidate.errorImpactParameter0(),
                         candidate.errorImpactParameter1(),
                         candidate.massLambda(),
                         candidate.massCascade(),
                         HfHelper::invMassXiczeroToXiPi(candidate),
                         candidate.cpaV0(),
                         candidate.cpa(),
                         candidate.cpaCascade(),
                         candidate.cpaXYV0(),
                         candidate.cpaXY(),
                         candidate.cpaXYCascade(),
                         HfHelper::ctXiczero(candidate),
                         trackV0PosDau.eta(),
                         trackV0NegDau.eta(),
                         trackPiFromCasc.eta(),
                         trackPiFromCharm.eta(),
                         candidate.eta(),
                         candidate.etaCascade(),
                         candidate.etaV0(),
                         candidate.impactParameterXYV0PosDau(),
                         candidate.impactParameterXYV0NegDau(),
                         candidate.impactParameterXYCascadeBach(),
                         candidate.impactParameterZV0PosDau(),
                         candidate.impactParameterZV0NegDau(),
                         candidate.impactParameterZCascadeBach(),
                         candidate.dcaCascadeDaughters(),
                         candidate.dcaV0Daughters(),
                         candidate.dcaDaughters(),
                         candidate.decayLength(),
                         candidate.decayLengthCascade(),
                         candidate.decayLengthV0(),
                         candidate.errorDecayLength(),
                         candidate.errorDecayLengthXY(),
                         candidate.impactParameterNormalised0(),
                         candidate.impactParameterNormalised1(),
                         candidate.decayLengthNormalised(),
                         trackPiFromCharm.isGlobalTrackWoDCA(),
                         trackPiFromCharm.itsNCls(),
                         trackPiFromCharm.tpcNClsCrossedRows(),
                         trackPiFromCasc.tpcNClsCrossedRows(),
                         trackV0PosDau.tpcNClsCrossedRows(),
                         trackV0NegDau.tpcNClsCrossedRows(),
                         trackPiFromCharm.tpcNSigmaPi(),
                         trackPiFromCasc.tpcNSigmaPi(),
                         trackPiFromLam.tpcNSigmaPi(),
                         trackPrFromLam.tpcNSigmaPr(),
                         trackPiFromCharm.tofNSigmaPi(),
                         trackPiFromCasc.tofNSigmaPi(),
                         trackPiFromLam.tofNSigmaPi(),
                         trackPrFromLam.tofNSigmaPr(),
                         flagMc,
                         originMc);
      }
    } else {
      if constexpr (tableSize == LITE) {
        // currently, no lite sized table for KFParticle
      } else {
        rowKfCandidate(centrality,
                       trackPiFromCharm.tpcNSigmaPi(),
                       trackPiFromCharm.tofNSigmaPi(),
                       trackPiFromCasc.tpcNSigmaPi(),
                       trackPiFromCasc.tofNSigmaPi(),
                       trackPiFromLam.tpcNSigmaPi(),
                       trackPiFromLam.tofNSigmaPi(),
                       trackPrFromLam.tpcNSigmaPr(),
                       trackPrFromLam.tofNSigmaPr(),
                       candidate.impactParameter1(),
                       candidate.dcaCascadeDaughters(),
                       candidate.dcaV0Daughters(),
                       candidate.dcaDaughters(),
                       candidate.impactParameter0(),
                       candidate.impactParameterXYV0PosDau(),
                       candidate.impactParameterXYV0NegDau(),
                       candidate.impactParameterXYCascadeBach(),
                       candidate.chi2GeoV0(),
                       candidate.chi2GeoCascade(),
                       candidate.chi2Geo(),
                       candidate.decayLengthNormKfV0(),
                       candidate.decayLengthNormKfCascade(),
                       candidate.chi2TopoOverNdfV0ToPv(),
                       candidate.chi2TopoOverNdfCascadeToPv(),
                       candidate.chi2TopoOverNdfToPv(),
                       candidate.chi2TopoOverNdfV0ToCascade(),
                       candidate.chi2TopoOverNdfCascadeToCharmHadron(),
                       candidate.decayLengthXYV0(),
                       candidate.decayLengthXYCascade(),
                       candidate.decayLengthXY(),
                       candidate.cpaV0ToCascade(),
                       candidate.cpaV0(),
                       candidate.cpaCascadeToCharmHadron(),
                       candidate.cpaCascade(),
                       candidate.massLambda(),
                       candidate.massCascade(),
                       HfHelper::invMassXiczeroToXiPi(candidate),
                       HfHelper::yXiczero(candidate),
                       RecoDecay::sqrtSumOfSquares(candidate.pxProng1(), candidate.pyProng1()),
                       candidate.pt(),
                       HfHelper::cosThetaStarXiczeroToXiPi(candidate),
                       HfHelper::ctXiczero(candidate),
                       candidate.eta(),
                       candidate.ndfGeoV0(),
                       candidate.ndfGeoCascade(),
                       candidate.ndfGeo(),
                       candidate.chi2GeoOverNdfV0(),
                       candidate.chi2GeoOverNdfCascade(),
                       candidate.chi2GeoOverNdf(),
                       flagMc,
                       originMc);
      }
    }
  }

  template <typename McCollType, typename CollType, typename CandType>
  void fillParticle(const CandType& mcParticles, const CollType& collisions)
  {
    for (const auto& particle : mcParticles) {

      auto ptGen = particle.pt();
      auto yGen = particle.y();

      int nContribMax = 0;
      bool recoCollPassedTvx = false;
      bool recoCollPassedSel8 = false;
      bool recoCollPassedSel8Zvtx = false;

      auto mcCollision = particle.template mcCollision_as<McCollType>();
      const auto& recoCollsPerMcColl = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());

      for (const auto& recoCol : recoCollsPerMcColl) {
        nContribMax = recoCol.numContrib() > nContribMax ? recoCol.numContrib() : nContribMax;
        if (recoCol.selection_bit(aod::evsel::kIsTriggerTVX)) {
          recoCollPassedTvx = true;
          if (recoCol.sel8()) {
            recoCollPassedSel8 = true;
            if (std::abs(recoCol.posZ()) < zPvCut) {
              recoCollPassedSel8Zvtx = true;
              break;
            }
          }
        }
      }

      float ptGenBhad = (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) ? mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt() : -999.f;

      rowCandidateParticles(RecoDecay::m(particle.pVector(), particle.e()),
                            ptGen,
                            particle.eta(),
                            particle.phi(),
                            yGen,
                            particle.flagMcMatchGen(),
                            particle.originMcGen(),
                            particle.pdgCode(),
                            ptGenBhad,
                            nContribMax,
                            recoCollsPerMcColl.size(),
                            recoCollPassedTvx,
                            recoCollPassedSel8,
                            recoCollPassedSel8Zvtx);
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

  void processDataFull(MyEventTable const& collisions,
                       MyTrackTable const& tracks,
                       MyTrackTableLf const& lfTracks,
                       MyCandTableDca const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTable>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processDataLite(MyEventTable const& collisions,
                       MyTrackTable const& tracks,
                       MyTrackTableLf const& lfTracks,
                       MyCandTableDca const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTable>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processDataLiteWithFT0M(MyEventTableWithFT0M const& collisions,
                               MyTrackTable const& tracks,
                               MyTrackTableLf const& lfTracks,
                               MyCandTableDca const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0M>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processDataLiteWithFT0C(MyEventTableWithFT0C const& collisions,
                               MyTrackTable const& tracks,
                               MyTrackTableLf const& lfTracks,
                               MyCandTableDca const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithFT0C>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processDataLiteWithNTracksPV(MyEventTableWithNTracksPV const& collisions,
                                    MyTrackTable const& tracks,
                                    MyTrackTableLf const& lfTracks,
                                    MyCandTableDca const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithNTracksPV>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
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
  void processKfData(MyEventTable const& collisions,
                     MyTrackTable const& tracks,
                     MyTrackTableLf const& lfTracks,
                     MyCandTableKf const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, false, MyEventTable>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processKfDataWithFT0M(MyEventTableWithFT0M const& collisions,
                             MyTrackTable const& tracks,
                             MyTrackTableLf const& lfTracks,
                             MyCandTableKf const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithFT0M>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processKfDataWithFT0C(MyEventTableWithFT0C const& collisions,
                             MyTrackTable const& tracks,
                             MyTrackTableLf const& lfTracks,
                             MyCandTableKf const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithFT0C>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processKfDataWithNTracksPV(MyEventTableWithNTracksPV const& collisions,
                                  MyTrackTable const& tracks,
                                  MyTrackTableLf const& lfTracks,
                                  MyCandTableKf const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithNTracksPV>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfData, "Process KF data, no cent", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithFT0M, "Process KF data, with FT0M", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithFT0C, "Process KF data, with FT0C", false);
  PROCESS_SWITCH(HfTreeCreatorToXiPiQa, processKfDataWithNTracksPV, "Process KF data, with NTracksPV", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~MC with DCAFitter~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processMcFullXic0(MyEventTableWithMcLabels const& collsWithMcLable,
                         aod::McCollisions const&,
                         MyTrackTable const& tracks,
                         MyTrackTableLf const& lfTracks,
                         MyCandTableXic0McDca const& candidates,
                         MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTableWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
    }
  }

  void processMcFullOmegac0(MyEventTableWithMcLabels const& collsWithMcLable,
                            aod::McCollisions const&,
                            MyTrackTable const& tracks,
                            MyTrackTableLf const& lfTracks,
                            MyCandTableOmegac0McDca const& candidates,
                            MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, FULL, false, MyEventTableWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
    }
  }

  void processMcLiteXic0(MyEventTableWithMcLabels const& collsWithMcLable,
                         aod::McCollisions const&,
                         MyTrackTable const& tracks,
                         MyTrackTableLf const& lfTracks,
                         MyCandTableXic0McDca const& candidates,
                         MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTableWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
    }
  }

  void processMcLiteXic0WithFT0C(MyEventTableWithMcLabelsFT0C const& collsWithMcLable,
                                 McCollisionsCentFT0Cs const&,
                                 MyTrackTable const& tracks,
                                 MyTrackTableLf const& lfTracks,
                                 MyCandTableXic0McDca const& candidates,
                                 MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithMcLabelsFT0C>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<McCollisionsCentFT0Cs>(mcParticles, collsWithMcLable);
    }
  }

  void processMcLiteXic0WithFT0M(MyEventTableWithMcLabelsFT0M const& collsWithMcLable,
                                 McCollisionsCentFT0Ms const&,
                                 MyTrackTable const& tracks,
                                 MyTrackTableLf const& lfTracks,
                                 MyCandTableXic0McDca const& candidates,
                                 MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithMcLabelsFT0M>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<McCollisionsCentFT0Ms>(mcParticles, collsWithMcLable);
    }
  }

  void processMcLiteXic0WithNTracksPV(MyEventTableWithMcLabelsNTracksPV const& collsWithMcLable,
                                      aod::McCollisions const&,
                                      MyTrackTable const& tracks,
                                      MyTrackTableLf const& lfTracks,
                                      MyCandTableXic0McDca const& candidates,
                                      MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, true, MyEventTableWithMcLabelsNTracksPV>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
    }
  }

  void processMcLiteOmegac0(MyEventTableWithMcLabels const& collsWithMcLable,
                            aod::McCollisions const&,
                            MyTrackTable const& tracks,
                            MyTrackTableLf const& lfTracks,
                            MyCandTableOmegac0McDca const& candidates,
                            MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, LITE, false, MyEventTableWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
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
  void processKfMcXic0(MyEventTableWithMcLabels const& collsWithMcLable,
                       aod::McCollisions const&,
                       MyTrackTable const& tracks,
                       MyTrackTableLf const& lfTracks,
                       MyCandTableXic0McKf const& candidates,
                       MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, false, MyEventTableWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
    }
  }

  void processKfMcXic0WithFT0C(MyEventTableWithMcLabelsFT0C const& collsWithMcLable,
                               McCollisionsCentFT0Cs const&,
                               MyTrackTable const& tracks,
                               MyTrackTableLf const& lfTracks,
                               MyCandTableXic0McKf const& candidates,
                               MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithMcLabelsFT0C>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<McCollisionsCentFT0Cs>(mcParticles, collsWithMcLable);
    }
  }

  void processKfMcXic0WithFT0M(MyEventTableWithMcLabelsFT0M const& collsWithMcLable,
                               McCollisionsCentFT0Ms const&,
                               MyTrackTable const& tracks,
                               MyTrackTableLf const& lfTracks,
                               MyCandTableXic0McKf const& candidates,
                               MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithMcLabelsFT0M>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<McCollisionsCentFT0Ms>(mcParticles, collsWithMcLable);
    }
  }

  void processKfMcXic0WithNTracksPV(MyEventTableWithMcLabelsNTracksPV const& collsWithMcLable,
                                    aod::McCollisions const&,
                                    MyTrackTable const& tracks,
                                    MyTrackTableLf const& lfTracks,
                                    MyCandTableXic0McKf const& candidates,
                                    MatchedGenXiPi const& mcParticles)
  {
    // Filling event properties
    rowEv.reserve(collsWithMcLable.size());
    for (const auto& collision : collsWithMcLable) {
      fillEvent<true>(collision);
    }

    // Filling candidate table
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, FULL, true, MyEventTableWithMcLabelsNTracksPV>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }

    // Filling particle properties if requested
    if (fillGenTable) {
      rowCandidateParticles.reserve(mcParticles.size());
      fillParticle<aod::McCollisions>(mcParticles, collsWithMcLable);
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
