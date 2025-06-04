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

/// \file treeCreatorXic0ToXiPiKf.cxx
/// \brief Writer of the xic0 to Xi Pi candidates in the form of flat tables to be stored in TTrees.
///        In this file are defined and filled the output tables
///
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University
/// \author Tao Fang <tao.fang@cern.ch>, Central China Normal University

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(InvMassLambda, invMassLambda, float);
DECLARE_SOA_COLUMN(InvMassCascade, invMassCascade, float);
DECLARE_SOA_COLUMN(InvMassCharmBaryon, invMassCharmBaryon, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau0, dcaXYToPvV0Dau0, float);
DECLARE_SOA_COLUMN(DcaXYToPvV0Dau1, dcaXYToPvV0Dau1, float);
DECLARE_SOA_COLUMN(DcaXYToPvCascDau, dcaXYToPvCascDau, float);
DECLARE_SOA_COLUMN(DcaCascDau, dcaCascDau, float);
DECLARE_SOA_COLUMN(DcaCharmBaryonDau, dcaCharmBaryonDau, float);
// from creator - MC
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t); // reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);         // debug flag for mis-association reconstruction level
DECLARE_SOA_COLUMN(OriginRec, originRec, int8_t);
DECLARE_SOA_COLUMN(CollisionMatched, collisionMatched, bool);
// from selector
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

} // namespace full

DECLARE_SOA_TABLE(HfKfXicFulls, "AOD", "HFKFXICFULL",
                  full::Centrality,
                  full::TpcNSigmaPiFromCharmBaryon, full::TofNSigmaPiFromCharmBaryon, full::TpcNSigmaPiFromCasc, full::TofNSigmaPiFromCasc,
                  full::TpcNSigmaPiFromLambda, full::TofNSigmaPiFromLambda, full::TpcNSigmaPrFromLambda, full::TofNSigmaPrFromLambda,
                  full::KfDcaXYPiFromXic, full::DcaCascDau, full::DcaCharmBaryonDau, full::KfDcaXYCascToPv,
                  full::DcaXYToPvV0Dau0, full::DcaXYToPvV0Dau1, full::DcaXYToPvCascDau,
                  full::Chi2GeoV0, full::Chi2GeoCasc, full::Chi2GeoXic,
                  full::Chi2MassV0, full::Chi2MassCasc,
                  full::V0ldl, full::Cascldl, full::Xicldl,
                  full::Chi2TopoV0ToPv, full::Chi2TopoCascToPv, full::Chi2TopoPiFromXicToPv, full::Chi2TopoXicToPv,
                  full::Chi2TopoV0ToCasc, full::Chi2TopoCascToXic,
                  full::DecayLenXYLambda, full::DecayLenXYCasc, full::DecayLenXYXic,
                  full::CosPaV0ToCasc, full::CosPaV0ToPv, full::CosPaCascToXic, full::CosPaCascToPv, full::CosPaXicToPv,
                  full::InvMassLambda, full::InvMassCascade, full::InvMassCharmBaryon,
                  full::KfRapXic, full::KfptPiFromXic, full::KfptXic,
                  full::CosThetaStarPiFromXic, full::CtXic, full::EtaXic,
                  full::V0Ndf, full::CascNdf, full::XicNdf,
                  full::MassV0Ndf, full::MassCascNdf,
                  full::V0Chi2OverNdf, full::CascChi2OverNdf, full::XicChi2OverNdf,
                  full::MassV0Chi2OverNdf, full::MassCascChi2OverNdf,
                  full::FlagMcMatchRec, full::DebugMcRec, full::OriginRec, full::CollisionMatched);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorXic0ToXiPiKf {

  Produces<o2::aod::HfKfXicFulls> rowKfCandidate;

  Configurable<float> zPvCut{"zPvCut", 10., "Cut on absolute value of primary vertex z coordinate"};

  using MyTrackTable = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>;
  using MyEventTable = soa::Join<aod::Collisions, aod::EvSels>;
  using MyEventTableWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using MyEventTableWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using MyEventTableWithNTracksPV = soa::Join<aod::Collisions, aod::EvSels, aod::CentNTPVs>;

  void init(InitContext const&)
  {
  }

  template <bool useCentrality, typename MyEventTableType, typename T>
  void fillKfCandidate(const T& candidate, int8_t flagMc, int8_t debugMc, int8_t originMc, bool collisionMatched)
  {

    if (candidate.resultSelections() && candidate.statusPidCharmBaryon() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {

      float centrality = -999.f;
      if constexpr (useCentrality) {
        auto const& collision = candidate.template collision_as<MyEventTableType>();
        centrality = o2::hf_centrality::getCentralityColl(collision);
      }

      rowKfCandidate(
        centrality,
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
        candidate.xicldl(),
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
        candidate.cosPACharmBaryon(),
        candidate.invMassLambda(),
        candidate.invMassCascade(),
        candidate.invMassCharmBaryon(),
        candidate.kfRapXic(),
        candidate.kfptPiFromXic(),
        candidate.kfptXic(),
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

  void processKfData(MyTrackTable const&, MyEventTable const&,
                     soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<false, MyEventTable>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfData, "Process KF data", false);

  void processKfDataWithFT0C(MyTrackTable const&, MyEventTableWithFT0C const&,
                             soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithFT0C>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfDataWithFT0C, "Process KF data with FT0C", false);

  void processKfDataWithFT0M(MyTrackTable const&, MyEventTableWithFT0M const&,
                             soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithFT0M>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfDataWithFT0M, "Process KF data with FT0M", false);

  void processDataLiteWithNTracksPV(MyTrackTable const&,
                                    soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithNTracksPV>(candidate, -7, -7, RecoDecay::OriginType::None, false);
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processDataLiteWithNTracksPV, "Process KF data with Ntracks", false);

  void processKfMcXic0(MyTrackTable const&,
                       soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<false, MyEventTable>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfMcXic0, "Process MC with information for xic0", false);

  void processKfMCWithFT0C(MyTrackTable const&,
                           soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithFT0C>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfMCWithFT0C, "Process MC with information for xic0 at FT0C", false);

  void processKfMCWithFT0M(MyTrackTable const&,
                           soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithFT0M>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processKfMCWithFT0M, "Process MC with information for xic0 at FT0M", false);

  void processMCLiteWithNTracksPV(MyTrackTable const&,
                                  soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfXicToXiPiMCRec> const& candidates)
  {
    rowKfCandidate.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      fillKfCandidate<true, MyEventTableWithNTracksPV>(candidate, candidate.flagMcMatchRec(), candidate.debugMcRec(), candidate.originRec(), candidate.collisionMatched());
    }
  }
  PROCESS_SWITCH(HfTreeCreatorXic0ToXiPiKf, processMCLiteWithNTracksPV, "Process MC with information for xic0 at Ntrack", false);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorXic0ToXiPiKf>(cfgc)};
}
