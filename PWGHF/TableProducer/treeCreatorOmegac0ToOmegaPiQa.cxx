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

/// \file treeCreatorOmegac0ToOmegaPiQa.cxx
/// \brief Writer of the omegac0 to Omega Pi candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston
/// \author Maria Fernanda Torres Cabrera <maria.fernanda.torres.cabrera@cern.ch>, University of Houston

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTablesQa.h"
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

#include <array>
#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// SV Reco method
enum {
  DCAFITTER = 0,
  KFPARTICLE
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
DECLARE_SOA_COLUMN(PxCharmBaryon, pxCharmBaryon, float);
DECLARE_SOA_COLUMN(PyCharmBaryon, pyCharmBaryon, float);
DECLARE_SOA_COLUMN(PzCharmBaryon, pzCharmBaryon, float);
DECLARE_SOA_COLUMN(PxPiFromCharmBaryon, pxPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PyPiFromCharmBaryon, pyPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PzPiFromCharmBaryon, pzPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(PtCharmBaryon, ptCharmBaryon, float);
DECLARE_SOA_COLUMN(PtPiFromCharmBaryon, ptPiFromCharmBaryon, float);
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
DECLARE_SOA_COLUMN(CosPAV0, cosPAV0, float);
DECLARE_SOA_COLUMN(CosPACharmBaryon, cosPACharmBaryon, float);
DECLARE_SOA_COLUMN(CosPACasc, cosPACasc, float);
DECLARE_SOA_COLUMN(EtaV0PosDau, etaV0PosDau, float);
DECLARE_SOA_COLUMN(EtaV0NegDau, etaV0NegDau, float);
DECLARE_SOA_COLUMN(EtaKaFromCasc, etaKaFromCasc, float);
DECLARE_SOA_COLUMN(EtaPiFromCharmBaryon, etaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(EtaCharmBaryon, etaCharmBaryon, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaV0Dau, dcaV0Dau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
DECLARE_SOA_COLUMN(DecLenCharmBaryon, decLenCharmBaryon, float);
DECLARE_SOA_COLUMN(ErrorDecayLengthCharmBaryon, errorDecayLengthCharmBaryon, float);
DECLARE_SOA_COLUMN(NormImpParCascade, normImpParCascade, double);
DECLARE_SOA_COLUMN(NormImpParPiFromCharmBar, normImpParPiFromCharmBar, double);
DECLARE_SOA_COLUMN(NormDecayLenCharmBar, normDecayLenCharmBar, double);
DECLARE_SOA_COLUMN(IsPionGlbTrkWoDca, isPionGlbTrkWoDca, bool);
DECLARE_SOA_COLUMN(PionItsNCls, pionItsNCls, uint8_t);
DECLARE_SOA_COLUMN(NTpcRowsPion, nTpcRowsPion, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsKaFromCasc, nTpcRowsKaFromCasc, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsPosV0Dau, nTpcRowsPosV0Dau, int16_t);
DECLARE_SOA_COLUMN(NTpcRowsNegV0Dau, nTpcRowsNegV0Dau, int16_t);
// from creator KF
DECLARE_SOA_COLUMN(KfDcaXYPiFromOmegac, kfDcaXYPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfDcaXYCascToPv, kfDcaXYCascToPv, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(Chi2GeoCasc, chi2GeoCasc, float);
DECLARE_SOA_COLUMN(Chi2GeoOmegac, chi2GeoOmegac, float);
DECLARE_SOA_COLUMN(V0ldl, v0ldl, float);
DECLARE_SOA_COLUMN(Cascldl, cascldl, float);
DECLARE_SOA_COLUMN(Omegacldl, omegacldl, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascToPv, chi2TopoCascToPv, float);
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
DECLARE_SOA_COLUMN(KfptPiFromOmegac, kfptPiFromOmegac, float);
DECLARE_SOA_COLUMN(KfptOmegac, kfptOmegac, float);
DECLARE_SOA_COLUMN(CosThetaStarPiFromOmegac, cosThetaStarPiFromOmegac, float);
DECLARE_SOA_COLUMN(EtaOmegac, etaOmegac, float);
DECLARE_SOA_COLUMN(V0Ndf, v0Ndf, float);
DECLARE_SOA_COLUMN(CascNdf, cascNdf, float);
DECLARE_SOA_COLUMN(OmegacNdf, omegacNdf, float);
DECLARE_SOA_COLUMN(V0Chi2OverNdf, v0Chi2OverNdf, float);
DECLARE_SOA_COLUMN(CascChi2OverNdf, cascChi2OverNdf, float);
DECLARE_SOA_COLUMN(OmegacChi2OverNdf, omegacChi2OverNdf, float);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
// particle PID
DECLARE_SOA_COLUMN(TpcNSigmaPiFromCharmBaryon, tpcNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TpcNSigmaKaFromCasc, tpcNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TpcNSigmaPiFromLambda, tpcNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TpcNSigmaPrFromLambda, tpcNSigmaPrFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromCharmBaryon, tofNSigmaPiFromCharmBaryon, float);
DECLARE_SOA_COLUMN(TofNSigmaKaFromCasc, tofNSigmaKaFromCasc, float);
DECLARE_SOA_COLUMN(TofNSigmaPiFromLambda, tofNSigmaPiFromLambda, float);
DECLARE_SOA_COLUMN(TofNSigmaPrFromLambda, tofNSigmaPrFromLambda, float);
} // namespace full

DECLARE_SOA_TABLE(HfToOmegaPiEvs, "AOD", "HFTOOMEPIEV",
                  full::IsEventSel8, full::IsEventSelZ);

DECLARE_SOA_TABLE(HfOmegac0ToOmegaPi, "AOD", "HFTOOMEGAPI",
                  full::XPv, full::YPv, full::ZPv, full::Centrality, collision::NumContrib, collision::Chi2,
                  full::XDecayVtxCharmBaryon, full::YDecayVtxCharmBaryon, full::ZDecayVtxCharmBaryon,
                  full::XDecayVtxCascade, full::YDecayVtxCascade, full::ZDecayVtxCascade,
                  full::XDecayVtxV0, full::YDecayVtxV0, full::ZDecayVtxV0,
                  full::PtCharmBaryon,
                  full::PtPiFromCharmBaryon,
                  full::PxCharmBaryon, full::PyCharmBaryon, full::PzCharmBaryon,
                  full::PxPiFromCharmBaryon, full::PyPiFromCharmBaryon, full::PzPiFromCharmBaryon,
                  full::PxKaFromCasc, full::PyKaFromCasc, full::PzKaFromCasc,
                  full::PxPosV0Dau, full::PyPosV0Dau, full::PzPosV0Dau,
                  full::PxNegV0Dau, full::PyNegV0Dau, full::PzNegV0Dau,
                  full::ImpactParCascXY, full::ImpactParPiFromCharmBaryonXY,
                  full::ErrImpactParCascXY, full::ErrImpactParPiFromCharmBaryonXY,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::CosPAV0, full::CosPACharmBaryon, full::CosPACasc,
                  full::EtaV0PosDau, full::EtaV0NegDau, full::EtaKaFromCasc,
                  full::EtaPiFromCharmBaryon,
                  full::EtaCharmBaryon,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau,
                  full::DecLenCharmBaryon,
                  full::ErrorDecayLengthCharmBaryon, full::NormImpParCascade, full::NormImpParPiFromCharmBar,
                  full::NormDecayLenCharmBar,
                  full::IsPionGlbTrkWoDca, full::PionItsNCls,
                  full::NTpcRowsPion, full::NTpcRowsKaFromCasc, full::NTpcRowsPosV0Dau, full::NTpcRowsNegV0Dau,
                  full::TpcNSigmaPiFromCharmBaryon, full::TpcNSigmaKaFromCasc, full::TpcNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda,
                  full::TofNSigmaPiFromCharmBaryon, full::TofNSigmaKaFromCasc, full::TofNSigmaPiFromLambda, full::TofNSigmaPrFromLambda,
                  full::FlagMcMatchRec, full::OriginMcRec);

DECLARE_SOA_TABLE(HfKfOmegac, "AOD", "HFKFOMEGAC",
                  full::Centrality,
                  collision::NumContrib,
                  full::TpcNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCharmBaryon, full::TpcNSigmaKaFromCasc, full::TofNSigmaKaFromCasc,
                  full::TpcNSigmaPiFromLambda, full::TofNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda, full::TofNSigmaPrFromLambda,
                  full::KfDcaXYPiFromOmegac, full::DcaCascDau, full::DcaV0Dau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoOmegac,
                  full::V0ldl, full::Cascldl, full::Omegacldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv,
                  full::Chi2TopoOmegacToPv,
                  full::Chi2TopoV0ToCasc, full::Chi2TopoCascToOmegac,
                  full::DecayLenXYLambda, full::DecayLenXYCasc, full::DecayLenXYOmegac,
                  full::CosPaV0ToCasc, full::CosPaV0ToPv, full::CosPaCascToOmegac, full::CosPaCascToPv, full::CosPaOmegacToPv,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfptPiFromOmegac, full::KfptOmegac, full::CosThetaStarPiFromOmegac,
                  full::EtaOmegac,
                  full::V0Ndf, full::CascNdf, full::OmegacNdf,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::OmegacChi2OverNdf,
                  full::FlagMcMatchRec, full::OriginMcRec);
} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorToOmegaPiQa {

  Produces<o2::aod::HfOmegac0ToOmegaPi> rowCandidate;
  Produces<o2::aod::HfKfOmegac> rowKfCandidate;
  Produces<o2::aod::HfToOmegaPiEvs> rowEv;

  Configurable<int> selectionFlagOmegac0{"selectionFlagOmegac0", 7, "Selection Flag for Omegac0 candidates. 1 -> skimmed, 3 -> topological, 7 -> PID, 15 -> ML"};
  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};
  Configurable<bool> keepOnlyMcSignal{"keepOnlyMcSignal", true, "Fill MC tree only with signal candidates"};

  using TracksSel = soa::Join<aod::TracksWExtra, aod::TrackSelection, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;
  using TracksSelLf = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksPidPi, aod::TracksPidPr, aod::TracksPidKa>;

  using CandSel = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfSelOmegaC0ToOmegaPi>>;
  using CandKfSel = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfCandCascBachKfQa, aod::HfSelOmegaC0ToOmegaPi>>;
  using CandMcSel = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfSelOmegaC0ToOmegaPi, aod::HfCandCascBachMCRecQa>>;
  using CandKfMcSel = soa::Filtered<soa::Join<aod::HfCandCascBachQa, aod::HfCandCascBachKfQa, aod::HfSelOmegaC0ToOmegaPi, aod::HfCandCascBachMCRecQa>>;

  using CollsSel = soa::Join<aod::Collisions, aod::EvSels>;
  using CollsSelWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollsWithMcLabels = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using McCollsWithFT0M = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Filter filterOmegaCMatchedRec = nabs(aod::hf_cand_mc_flag::flagMcMatchRec) == static_cast<int8_t>(hf_decay::hf_cand_casc_bach_qa::DecayChannelMain::Omegac0ToOmegaPi);
  Filter filterSelectCandidates = aod::hf_sel_candidate_omegac0_qa::isSelOmegaC0ToOmegaPi >= selectionFlagOmegac0;

  void init(InitContext const&)
  {
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

  template <int svReco, bool useCentrality, typename CollTable, typename T>
  void fillCandidate(const T& candidate, TracksSel const& tracks, TracksSelLf const& lfTracks, int8_t flagMc, int8_t originMc, float centrality = -999.f)
  {

    // save all candidate information
    if constexpr (useCentrality) {
      if (centrality < -998.f) {
        auto const& collision = candidate.template collision_as<CollTable>();
        centrality = o2::hf_centrality::getCentralityColl(collision);
      }
    }

    // retrive PID information
    auto trackV0PosDau = lfTracks.rawIteratorAt(candidate.posTrackId());
    auto trackV0NegDau = lfTracks.rawIteratorAt(candidate.negTrackId());
    auto trackKaFromCasc = lfTracks.rawIteratorAt(candidate.bachelorId());
    auto trackPiFromCharm = tracks.rawIteratorAt(candidate.prong1Id());

    auto trackPiFromLam = trackV0NegDau;
    auto trackPrFromLam = trackV0PosDau;
    if (trackPiFromCharm.sign() < 0) {
      trackPiFromLam = trackV0PosDau;
      trackPrFromLam = trackV0NegDau;
    }

    if constexpr (svReco == DCAFITTER) {
      auto const& collision = candidate.template collision_as<CollTable>();
      rowCandidate(collision.posX(),
                   collision.posY(),
                   collision.posZ(),
                   centrality,
                   collision.numContrib(),
                   collision.chi2(),
                   candidate.xSecondaryVertex(),
                   candidate.ySecondaryVertex(),
                   candidate.zSecondaryVertex(),
                   candidate.xCascadeVertex(),
                   candidate.yCascadeVertex(),
                   candidate.zCascadeVertex(),
                   candidate.xV0Vertex(),
                   candidate.yV0Vertex(),
                   candidate.zV0Vertex(),
                   candidate.pt(),
                   RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1()),
                   candidate.px(),
                   candidate.py(),
                   candidate.pz(),
                   candidate.pxProng1(),
                   candidate.pyProng1(),
                   candidate.pzProng1(),
                   trackKaFromCasc.px(),
                   trackKaFromCasc.py(),
                   trackKaFromCasc.pz(),
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
                   HfHelperQa::invMassOmegaczeroToOmegaPi(candidate),
                   candidate.cpaV0(),
                   candidate.cpa(),
                   candidate.cpaCascade(),
                   trackV0PosDau.eta(),
                   trackV0NegDau.eta(),
                   trackKaFromCasc.eta(),
                   RecoDecay::eta(std::array{candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()}),
                   candidate.eta(),
                   candidate.impactParameterXYV0PosDau(),
                   candidate.impactParameterXYV0NegDau(),
                   candidate.impactParameterXYCascadeBach(),
                   candidate.dcaCascadeDaughters(),
                   candidate.dcaV0Daughters(),
                   candidate.dcaDaughters(),
                   candidate.decayLength(),
                   candidate.errorDecayLength(),
                   candidate.impactParameterNormalised0(),
                   candidate.impactParameterNormalised1(),
                   candidate.decayLengthNormalised(),
                   trackPiFromCharm.isGlobalTrackWoDCA(),
                   trackPiFromCharm.itsNCls(),
                   trackPiFromCharm.tpcNClsCrossedRows(),
                   trackKaFromCasc.tpcNClsCrossedRows(),
                   trackV0PosDau.tpcNClsCrossedRows(),
                   trackV0NegDau.tpcNClsCrossedRows(),
                   trackPiFromCharm.tpcNSigmaPi(),
                   trackKaFromCasc.tpcNSigmaKa(),
                   trackPiFromLam.tpcNSigmaPi(),
                   trackPrFromLam.tpcNSigmaPr(),
                   trackPiFromCharm.tofNSigmaPi(),
                   trackKaFromCasc.tofNSigmaKa(),
                   trackPiFromLam.tofNSigmaPi(),
                   trackPrFromLam.tofNSigmaPr(),
                   flagMc,
                   originMc);
    } else {
      rowKfCandidate(centrality,
                     candidate.template collision_as<CollTable>().numContrib(),
                     trackPiFromCharm.tpcNSigmaPi(),
                     trackPiFromCharm.tofNSigmaPi(),
                     trackKaFromCasc.tpcNSigmaKa(),
                     trackKaFromCasc.tofNSigmaKa(),
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
                     candidate.decayLengthNormKf(),
                     candidate.chi2TopoV0ToPv(),
                     candidate.chi2TopoCascadeToPv(),
                     candidate.chi2TopoToPv(),
                     candidate.chi2TopoV0ToCascade(),
                     candidate.chi2TopoCascadeToCharmHadron(),
                     candidate.decayLengthXYV0(),
                     candidate.decayLengthXYCascade(),
                     candidate.decayLengthXY(),
                     candidate.cpaV0ToCascade(),
                     candidate.cpaV0(),
                     candidate.cpaCascadeToCharmHadron(),
                     candidate.cpaCascade(),
                     candidate.cpa(),
                     candidate.massLambda(),
                     candidate.massCascade(),
                     HfHelperQa::invMassOmegaczeroToOmegaPi(candidate),
                     RecoDecay::pt(candidate.pxProng1(), candidate.pyProng1()),
                     candidate.pt(),
                     HfHelperQa::cosThetaStarOmegaczeroToOmegaPi(candidate),
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

  ////////////////////////////////////
  //                                //
  //       Process functions        //
  //                                //
  ////////////////////////////////////

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~Data with DCAFitter~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processData(CollsSel const& collisions,
                   TracksSel const& tracks,
                   TracksSelLf const& lfTracks,
                   CandSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, false, CollsSel>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processDataCent(CollsSelWithFT0M const& collisions,
                       TracksSel const& tracks,
                       TracksSelLf const& lfTracks,
                       CandSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<DCAFITTER, true, CollsSelWithFT0M>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processData, "Process DCAFitter data, no cent", true);
  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processDataCent, "Process DCAFitter data, with FT0M", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~Data with KFParticle~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processKfData(CollsSel const& collisions,
                     TracksSel const& tracks,
                     TracksSelLf const& lfTracks,
                     CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, false, CollsSel>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  void processKfDataCent(CollsSelWithFT0M const& collisions,
                         TracksSel const& tracks,
                         TracksSelLf const& lfTracks,
                         CandKfSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<true>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillCandidate<KFPARTICLE, true, CollsSelWithFT0M>(candidate, tracks, lfTracks, -7, RecoDecay::OriginType::None);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processKfData, "Process KF data, no cent", false);
  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processKfDataCent, "Process KF data, with FT0M", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~MC with DCAFitter~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processMc(CollsSel const& collisions,
                 TracksSel const& tracks,
                 TracksSelLf const& lfTracks,
                 CandMcSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      fillCandidate<DCAFITTER, false, CollsSel>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }
  }

  void processMcCent(CollsWithMcLabels const& collisions,
                     TracksSel const& tracks,
                     TracksSelLf const& lfTracks,
                     CandMcSel const& candidates,
                     McCollsWithFT0M const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      auto mcCollision = candidate.collision_as<CollsWithMcLabels>().mcCollision_as<McCollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(mcCollision);
      fillCandidate<DCAFITTER, true, CollsWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec(), centFt0m);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processMc, "Process DCAFitter MC for Omegac0", false);
  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processMcCent, "Process DCAFitter MC for Omegac0 with FT0M", false);

  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//
  //*~~~~~~~MC with KFParticle~~~~~~~~*//
  //*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*//

  void processKfMc(CollsSel const& collisions,
                   TracksSel const& tracks,
                   TracksSelLf const& lfTracks,
                   CandKfMcSel const& candidates)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      fillCandidate<KFPARTICLE, false, CollsSel>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec());
    }
  }

  void processKfMcCent(CollsWithMcLabels const& collisions,
                       TracksSel const& tracks,
                       TracksSelLf const& lfTracks,
                       CandKfMcSel const& candidates,
                       McCollsWithFT0M const&)
  {
    // Filling event properties
    rowEv.reserve(collisions.size());
    for (const auto& collision : collisions) {
      fillEvent<false>(collision);
    }

    // Filling candidate properties
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      if (keepOnlyMcSignal && candidate.originMcRec() == RecoDecay::OriginType::None) {
        continue;
      }
      auto mcCollision = candidate.collision_as<CollsWithMcLabels>().mcCollision_as<McCollsWithFT0M>();
      float centFt0m = o2::hf_centrality::getCentralityColl(mcCollision);
      fillCandidate<KFPARTICLE, true, CollsWithMcLabels>(candidate, tracks, lfTracks, candidate.flagMcMatchRec(), candidate.originMcRec(), centFt0m);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processKfMc, "Process KF MC for Omegac0", false);
  PROCESS_SWITCH(HfTreeCreatorToOmegaPiQa, processKfMcCent, "Process KF MC for Omegac0 with FT0M", false);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorToOmegaPiQa>(cfgc)};
}
