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

/// \file treeCreatorOmegacToOmegaPi.cxx
/// \brief Writer of the omegac0 to Omega Pi candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences

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
DECLARE_SOA_COLUMN(PxPiFromCharmBaryon, pxPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyPiFromCharmBaryon, pyPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzPiFromCharmBaryon, pzPiFromCharmBaryon, float);
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
DECLARE_SOA_COLUMN(ImpactParPiFromCharmBaryonXY, impactParPiFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(ErrImpactParCascXY, errImpactParCascXY, float);
DECLARE_SOA_COLUMN(ErrImpactParPiFromCharmBaryonXY, errImpactParPiFromCharmBaryonXY, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaKaFromCasc, etaKaFromCasc, float);
DECLARE_SOA_COLUMN(EtaPiFromCharmBaryon, etaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(NormImpParCascade, normImpParCascade, double);
DECLARE_SOA_COLUMN(NormImpParPiFromCharmBar, normImpParPiFromCharmBar, double);
DECLARE_SOA_COLUMN(IsPionGlbTrkWoDca, isPionGlbTrkWoDca, bool);
DECLARE_SOA_COLUMN(PionItsNCls, pionItsNCls, uint8_t);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
// from selector
DECLARE_SOA_COLUMN(PidTpcInfoStored, pidTpcInfoStored, int);
DECLARE_SOA_COLUMN(PidTofInfoStored, pidTofInfoStored, int);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCasc, tpcNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCasc, tofNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
// from creator KF
DECLARE_SOA_COLUMN(NSigmaTPCPiFromOmegac, nSigmaTPCPiFromOmegac, float);
DECLARE_SOA_COLUMN(NSigmaTOFPiFromOmegac, nSigmaTOFPiFromOmegac, float);
DECLARE_SOA_COLUMN(NSigmaTPCKaFromCasc, nSigmaTPCKaFromCasc, float);
DECLARE_SOA_COLUMN(NSigmaTOFKaFromCasc, nSigmaTOFKaFromCasc, float);
DECLARE_SOA_COLUMN(NSigmaTPCPiFromV0, nSigmaTPCPiFromV0, float);
DECLARE_SOA_COLUMN(NSigmaTPCPrFromV0, nSigmaTPCPrFromV0, float);
DECLARE_SOA_COLUMN(KfDcaXYPiFromOmegac, kfDcaXYPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfDcaXYCascToPv, kfDcaXYCascToPv, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCasc, chi2GeoCasc, float);
DECLARE_SOA_COLUMN(Chi2GeoOmegac, chi2GeoOmegac, float);
DECLARE_SOA_COLUMN(Chi2MassV0, chi2MassV0, float);
DECLARE_SOA_COLUMN(Chi2MassCasc, chi2MassCasc, float);
DECLARE_SOA_COLUMN(V0ldl, v0ldl, float);
DECLARE_SOA_COLUMN(Cascldl, cascldl, float);
DECLARE_SOA_COLUMN(Omegacldl, omegacldl, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToPv, chi2TopoCascToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoPiFromOmegacToPv, chi2TopoPiFromOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoOmegacToPv, chi2TopoOmegacToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCasc, chi2TopoV0ToCasc, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToOmegac, chi2TopoCascToOmegac, float);
DECLARE_SOA_COLUMN(DecayLenXYLambda, decayLenXYLambda, float);
DECLARE_SOA_COLUMN(DecayLenXYCasc, decayLenXYCasc, float);
DECLARE_SOA_COLUMN(DecayLenXYOmegac, decayLenXYOmegac, float);
DECLARE_SOA_COLUMN(CosPaV0ToCasc, cosPaV0ToCasc, float);
DECLARE_SOA_COLUMN(CosPaV0ToPv, cosPaV0ToPv, float);
DECLARE_SOA_COLUMN(CosPaCascToOmegac, cosPaCascToOmegac, float);
DECLARE_SOA_COLUMN(CosPaCascToPv, cosPaCascToPv, float);
DECLARE_SOA_COLUMN(CosPaOmegacToPv, cosPaOmegacToPv, float);
DECLARE_SOA_COLUMN(KfRapOmegac, kfRapOmegac, float);
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromOmegac, cosThetaStarPiFromOmegac, float);
DECLARE_SOA_COLUMN(CtOmegac, ctOmegac, float);
DECLARE_SOA_COLUMN(EtaOmegac, etaOmegac, float);
DECLARE_SOA_COLUMN(V0Ndf, v0Ndf, float);
DECLARE_SOA_COLUMN(CascNdf, cascNdf, float);
DECLARE_SOA_COLUMN(OmegacNdf, omegacNdf, float);
DECLARE_SOA_COLUMN(MassV0Ndf, massV0Ndf, float);
DECLARE_SOA_COLUMN(MassCascNdf, massCascNdf, float);
DECLARE_SOA_COLUMN(V0Chi2OverNdf, v0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(CascChi2OverNdf, cascChi2OverNdf, float);
DECLARE_SOA_COLUMN(OmegacChi2OverNdf, omegacChi2OverNdf, float);
DECLARE_SOA_COLUMN(MassV0Chi2OverNdf, massV0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(MassCascChi2OverNdf, massCascChi2OverNdf, float);
DECLARE_SOA_COLUMN(CascRejectInvmass, cascRejectInvmass, float);
} // namespace full

DECLARE_SOA_TABLE(HfToOmegaPiEvs, "AOD", "HFTOOMEPIEV",
                  full::IsEventSel8, full::IsEventSelZ);

DECLARE_SOA_TABLE(HfOmegac0ToOmegaPiLites, "AOD", "HFTOOMEPILITE",
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::SignDecay,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxKaFromCasc, full::PyKaFromCasc, full::PzKaFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaKaFromCasc, full::EtaPiFromCharmBaryon,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::ErrorDecayLengthCharmBaryon, full::NormImpParCascade, full::NormImpParPiFromCharmBar,
                  full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::PidTpcInfoStored, full::PidTofInfoStored,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaKaFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaKaFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginRec, full::CollisionMatched);

DECLARE_SOA_TABLE(HfKfOmegacFulls, "AOD", "HFKFOMEGACFULL",
                  full::NSigmaTPCPiFromOmegac, full::NSigmaTOFPiFromOmegac, full::NSigmaTPCKaFromCasc, full::NSigmaTOFKaFromCasc,
                  full::NSigmaTPCPiFromV0, full::NSigmaTPCPrFromV0,
                  full::KfDcaXYPiFromOmegac, full::DcaCascDau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoOmegac,
                  full::Chi2MassV0, full::Chi2MassCasc,
                  full::V0ldl, full::Cascldl, full::Omegacldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv, full::Chi2TopoPiFromOmegacToPv, full::Chi2TopoOmegacToPv,
                  full::Chi2TopoV0ToCasc, full::Chi2TopoCascToOmegac,
                  full::DecayLenXYLambda, full::DecayLenXYCasc, full::DecayLenXYOmegac,
                  full::CosPaV0ToCasc, full::CosPaV0ToPv, full::CosPaCascToOmegac, full::CosPaCascToPv, full::CosPaOmegacToPv,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfRapOmegac, full::KfptPiFromOmegac, full::KfptOmegac,
                  full::CosThetaStarPiFromOmegac, full::CtOmegac, full::EtaOmegac,
                  full::V0Ndf, full::CascNdf, full::OmegacNdf,
                  full::MassV0Ndf, full::MassCascNdf,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::OmegacChi2OverNdf,
                  full::MassV0Chi2OverNdf, full::MassCascChi2OverNdf, full::CascRejectInvmass,
                  full::FlagMcMatchRec, full::OriginRec, full::CollisionMatched);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorOmegac0ToOmegaPi {

  Produces<o2::aod::HfOmegac0ToOmegaPiLites> rowCandidateLite;
  Produces<o2::aod::HfKfOmegacFulls> rowKfCandidateFull;
  Produces<o2::aod::HfToOmegaPiEvs> rowEv;

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
        candidate.tpcNSigmaPiFromCharmBaryon(),
        candidate.tpcNSigmaKaFromCasc(),
        candidate.tpcNSigmaPiFromLambda(),
        candidate.tpcNSigmaPrFromLambda(),
        candidate.tofNSigmaPiFromCharmBaryon(),
        candidate.tofNSigmaKaFromCasc(),
        candidate.tofNSigmaPiFromLambda(),
        candidate.tofNSigmaPrFromLambda(),
        flagMc,
        originMc,
        collisionMatched);
    }
  }

  template <typename T>
  void fillKfCandidate(const T& candidate, int8_t flagMc, int8_t originMc, bool collisionMatched)
  {
    if (candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {

      rowKfCandidateFull(
        candidate.tpcNSigmaPiFromCharmBaryon(),
        candidate.tofNSigmaPiFromCharmBaryon(),
        candidate.tpcNSigmaKaFromCasc(),
        candidate.tofNSigmaKaFromCasc(),
        candidate.tpcNSigmaPiFromLambda(),
        candidate.tpcNSigmaPrFromLambda(),
        candidate.kfDcaXYPiFromOmegac(),
        candidate.dcaCascDau(),
        candidate.dcaCharmBaryonDau(),
        candidate.kfDcaXYCascToPv(),
        candidate.dcaXYToPvV0Dau0(),
        candidate.dcaXYToPvV0Dau1(),
        candidate.dcaXYToPvCascDau(),
        candidate.chi2GeoV0(),
        candidate.chi2GeoCasc(),
        candidate.chi2GeoOmegac(),
        candidate.chi2MassV0(),
        candidate.chi2MassCasc(),
        candidate.v0ldl(),
        candidate.cascldl(),
        candidate.omegacldl(),
        candidate.chi2TopoV0ToPv(),
        candidate.chi2TopoCascToPv(),
        candidate.chi2TopoPiFromOmegacToPv(),
        candidate.chi2TopoOmegacToPv(),
        candidate.chi2TopoV0ToCasc(),
        candidate.chi2TopoCascToOmegac(),
        candidate.decayLenXYLambda(),
        candidate.decayLenXYCasc(),
        candidate.decayLenXYOmegac(),
        candidate.cosPaV0ToCasc(),
        candidate.cosPAV0(),
        candidate.cosPaCascToOmegac(),
        candidate.cosPACasc(),
        candidate.cosPACharmBaryon(),
        candidate.invMassLambda(),
        candidate.invMassCascade(),
        candidate.invMassCharmBaryon(),
        candidate.kfRapOmegac(),
        candidate.kfptPiFromOmegac(),
        candidate.kfptOmegac(),
        candidate.cosThetaStarPiFromOmegac(),
        candidate.ctauOmegac(),
        candidate.etaCharmBaryon(),
        candidate.v0Ndf(),
        candidate.cascNdf(),
        candidate.omegacNdf(),
        candidate.massV0Ndf(),
        candidate.massCascNdf(),
        candidate.v0Chi2OverNdf(),
        candidate.cascChi2OverNdf(),
        candidate.omegacChi2OverNdf(),
        candidate.massV0Chi2OverNdf(),
        candidate.massCascChi2OverNdf(),
        candidate.cascRejectInvmass(),
        flagMc,
        originMc,
        collisionMatched);
    }
  }

  void processDataLite(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidateLite(candidate, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processDataLite, "Process data", true);

  void processKfDataFull(MyEventTable const& collisions, MyTrackTable const&,
                         soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf, aod::HfSelToOmegaPi> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate(candidate, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKfDataFull, "Process KF data", false);

  void processMcLite(MyEventTable const& collisions, MyTrackTable const&,
                     soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfToOmegaPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidateLite(candidate, candidate.flagMcMatchRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processMcLite, "Process MC", false);

  void processKFMcFull(MyEventTable const& collisions, MyTrackTable const&,
                       soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf, aod::HfSelToOmegaPi, aod::HfToOmegaPiMCRec> const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate(candidate, candidate.flagMcMatchRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKFMcFull, "Process KF MC", false);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorOmegac0ToOmegaPi>(cfgc)};
}
