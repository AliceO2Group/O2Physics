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

/// \file treeCreatorOmegac0ToOmegaPi.cxx
/// \brief Writer of the omegac0 to Omega Pi candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Federica Zanone <federica.zanone@cern.ch>, Heidelberg University
/// \author Yunfan Liu <yunfan.liu@cern.ch>, China University of Geosciences
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston
/// \author Ruiqi Yin <ruiqi.yin@cern.ch>, Fudan University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

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
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
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
DECLARE_SOA_COLUMN(DeviationPiFromOmegacToPv, deviationPiFromOmegacToPv, float);
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
                  full::XPv, full::YPv, full::ZPv, collision::NumContrib, collision::Chi2, cent::CentFT0M,
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
                  full::FlagMcMatchRec, full::OriginMcRec, full::CollisionMatched, hf_track_index::HFflag);

DECLARE_SOA_TABLE(HfKfOmegacFulls, "AOD", "HFKFOMEGACFULL",
                  full::NSigmaTPCPiFromOmegac, full::NSigmaTOFPiFromOmegac, full::NSigmaTPCKaFromCasc, full::NSigmaTOFKaFromCasc,
                  full::NSigmaTPCPiFromV0, full::NSigmaTPCPrFromV0,
                  full::KfDcaXYPiFromOmegac, full::DcaCascDau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoOmegac,
                  full::Chi2MassV0, full::Chi2MassCasc,
                  full::V0ldl, full::Cascldl, full::Omegacldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv, full::Chi2TopoPiFromOmegacToPv, full::Chi2TopoOmegacToPv, full::DeviationPiFromOmegacToPv,
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
                  full::FlagMcMatchRec, full::OriginMcRec, full::CollisionMatched, hf_track_index::HFflag, collision::NumContrib, cent::CentFT0M);

DECLARE_SOA_TABLE(HfKfOmegacLites, "AOD", "HFKFOMEGACLITE",
                  full::NSigmaTPCPiFromOmegac, full::NSigmaTOFPiFromOmegac, full::NSigmaTPCKaFromCasc, full::NSigmaTOFKaFromCasc,
                  full::NSigmaTPCPiFromV0, full::NSigmaTPCPrFromV0,
                  full::KfDcaXYPiFromOmegac, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv, full::DcaCascDau,
                  full::V0ldl, full::Cascldl, full::Omegacldl, full::Chi2TopoPiFromOmegacToPv, full::Chi2TopoOmegacToPv, full::DeviationPiFromOmegacToPv,
                  full::DecayLenXYOmegac,
                  full::CosPaCascToPv, full::CosPaOmegacToPv,
                  full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfptPiFromOmegac, full::KfptOmegac,
                  full::CosThetaStarPiFromOmegac, full::CtOmegac, full::EtaOmegac,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::OmegacChi2OverNdf,
                  full::CascRejectInvmass,
                  full::FlagMcMatchRec, full::OriginMcRec, full::CollisionMatched, hf_track_index::HFflag, collision::NumContrib, cent::CentFT0M);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorOmegac0ToOmegaPi {

  Produces<o2::aod::HfOmegac0ToOmegaPiLites> rowCandidateLite;
  Produces<o2::aod::HfKfOmegacFulls> rowKfCandidateFull;
  Produces<o2::aod::HfKfOmegacLites> rowKfCandidateLite;
  Produces<o2::aod::HfToOmegaPiEvs> rowEv;

  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};
  Configurable<bool> keepOnlyMcSignal{"keepOnlyMcSignal", true, "Fill MC tree only with signal candidates"};

  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>;

  using CandSel = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi>>;
  using CandKfSel = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf, aod::HfSelToOmegaPi>>;
  using CandMcSel = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfSelToOmegaPi, aod::HfToOmegaPiMCRec>>;
  using CandKfMcSel = soa::Filtered<soa::Join<aod::HfCandToOmegaPi, aod::HfOmegacKf, aod::HfSelToOmegaPi, aod::HfToOmegaPiMCRec>>;

  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using CollsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using McCollsWithFT0M = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterOmegaCToOmegaPiFlag = (o2::aod::hf_track_index::hfflag & static_cast<uint8_t>(BIT(aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi))) != static_cast<uint8_t>(0);

  void init(InitContext const&)
  {
  }

  template <typename C>
  void fillEvent(const C& collision, float cutZPv)
  {
    rowEv(collision.sel8(), std::abs(collision.posZ()) < cutZPv);
  }

  template <typename P, typename C>
  void fillCandidateLite(const P& candidate, const C&, int8_t flagMc, int8_t originMc, bool collisionMatched, float centFt0m)
  {
    if (candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {

      rowCandidateLite(
        candidate.xPv(),
        candidate.yPv(),
        candidate.zPv(),
        candidate.template collision_as<C>().numContrib(),
        candidate.template collision_as<C>().chi2(),
        centFt0m,
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
        candidate.template bachelorFromCharmBaryon_as<Tracks>().isGlobalTrackWoDCA(),
        candidate.template bachelorFromCharmBaryon_as<Tracks>().itsNCls(),
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
        collisionMatched,
        candidate.hfflag());
    }
  }

  template <typename P, typename C>
  void fillKfCandidate(const P& candidate, const C&, int8_t flagMc, int8_t originMc, bool collisionMatched, float centFt0m)
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
        candidate.deviationPiFromOmegacToPv(),
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
        candidate.cTauOmegac(),
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
        collisionMatched,
        candidate.hfflag(),
        candidate.template collision_as<C>().numContrib(),
        centFt0m);
    }
  }

  template <typename P, typename C>
  void fillKfCandidateLite(const P& candidate, const C&, int8_t flagMc, int8_t originMc, bool collisionMatched, float centFt0m)
  {
    if (candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {

      rowKfCandidateLite(
        candidate.tpcNSigmaPiFromCharmBaryon(),
        candidate.tofNSigmaPiFromCharmBaryon(),
        candidate.tpcNSigmaKaFromCasc(),
        candidate.tofNSigmaKaFromCasc(),
        candidate.tpcNSigmaPiFromLambda(),
        candidate.tpcNSigmaPrFromLambda(),
        candidate.kfDcaXYPiFromOmegac(),
        candidate.dcaCharmBaryonDau(),
        candidate.kfDcaXYCascToPv(),
        candidate.dcaCascDau(),
        candidate.v0ldl(),
        candidate.cascldl(),
        candidate.omegacldl(),
        candidate.chi2TopoPiFromOmegacToPv(),
        candidate.chi2TopoOmegacToPv(),
        candidate.deviationPiFromOmegacToPv(),
        candidate.decayLenXYOmegac(),
        candidate.cosPACasc(),
        candidate.cosPACharmBaryon(),
        candidate.invMassCascade(),
        candidate.invMassCharmBaryon(),
        candidate.kfptPiFromOmegac(),
        candidate.kfptOmegac(),
        candidate.cosThetaStarPiFromOmegac(),
        candidate.cTauOmegac(),
        candidate.etaCharmBaryon(),
        candidate.v0Chi2OverNdf(),
        candidate.cascChi2OverNdf(),
        candidate.omegacChi2OverNdf(),
        candidate.cascRejectInvmass(),
        flagMc,
        originMc,
        collisionMatched,
        candidate.hfflag(),
        candidate.template collision_as<C>().numContrib(),
        centFt0m);
    }
  } // fillKfCandidateLite end

  void processData(Colls const& collisions,
                   CandSel const& candidates,
                   Tracks const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidateLite(candidate, collisions, -7, RecoDecay::OriginType::None, false, -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processData, "Process data", false);

  void processKfData(Colls const& collisions,
                     CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate(candidate, collisions, -7, RecoDecay::OriginType::None, false, -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKfData, "Process KF data", false);

  void processKfDataLite(Colls const& collisions,
                         CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidateLite(candidate, collisions, -7, RecoDecay::OriginType::None, false, -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKfDataLite, "Process KF data lite", false);

  void processDataCent(CollsWithFT0M const& collisions,
                       CandSel const& candidates,
                       Tracks const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      auto collision = candidate.collision_as<CollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(collision);
      fillCandidateLite(candidate, collisions, -7, RecoDecay::OriginType::None, false, centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processDataCent, "Process data with FT0M info", false);

  void processKfDataCent(CollsWithFT0M const& collisions,
                         CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      auto collision = candidate.collision_as<CollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(collision);
      fillKfCandidate(candidate, collisions, -7, RecoDecay::OriginType::None, false, centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKfDataCent, "Process KF data with FT0M info", false);

  void processKfDataCentLite(CollsWithFT0M const& collisions,
                             CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowKfCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      auto collision = candidate.collision_as<CollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(collision);
      fillKfCandidateLite(candidate, collisions, -7, RecoDecay::OriginType::None, false, centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKfDataCentLite, "Process KF data lite with FT0M info", false);

  void processMc(Colls const& collisions,
                 CandMcSel const& candidates,
                 Tracks const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      fillCandidateLite(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processMc, "Process MC", false);

  void processKFMc(Colls const& collisions,
                   CandKfMcSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      fillKfCandidate(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKFMc, "Process KF MC", false);

  void processKFMcLite(Colls const& collisions,
                       CandKfMcSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      fillKfCandidateLite(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), -1.);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKFMcLite, "Process KF MC Lite", false);

  void processMcCent(CollsWithMcLabels const& collisions,
                     CandMcSel const& candidates,
                     Tracks const&,
                     McCollsWithFT0M const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      auto mcCollision = candidate.collision_as<CollsWithMcLabels>().mcCollision_as<McCollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(mcCollision);
      fillCandidateLite(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processMcCent, "Process MC with FT0M info", false);

  void processKFMcCent(CollsWithMcLabels const& collisions,
                       CandKfMcSel const& candidates,
                       McCollsWithFT0M const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      auto mcCollision = candidate.collision_as<CollsWithMcLabels>().mcCollision_as<McCollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(mcCollision);
      fillKfCandidate(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKFMcCent, "Process KF MC with FT0M info", false);

  void processKFMcCentLite(CollsWithMcLabels const& collisions,
                           CandKfMcSel const& candidates,
                           McCollsWithFT0M const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent(collision, zPvCut);
    }

    // Filling candidate properties
    rowCandidateLite.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      auto mcCollision = candidate.collision_as<CollsWithMcLabels>().mcCollision_as<McCollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(mcCollision);
      fillKfCandidateLite(candidate, collisions, candidate.flagMcMatchRec(), candidate.originMcRec(), candidate.collisionMatched(), centFt0m);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegac0ToOmegaPi, processKFMcCentLite, "Process KF MC Lite with FT0M info", false);

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorOmegac0ToOmegaPi>(cfgc)};
}
