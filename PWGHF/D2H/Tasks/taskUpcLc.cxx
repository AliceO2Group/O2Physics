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

/// \file taskUpcLc.cxx
/// \brief Λc± → p± K∓ π± analysis task
/// \note Extended from taskLc
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University
/// \author Ran Tu <ran.tu@cern.ch>, Fudan University

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsUpcHf.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCHelpers.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <vector> // std::vector

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::analysis::hf_upc;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(BkgScore, bkgScore, float);
DECLARE_SOA_COLUMN(PromptScore, promptScore, float);
DECLARE_SOA_COLUMN(FdScore, fdScore, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(Cpa, cpa, float);
DECLARE_SOA_COLUMN(PvContributors, pvContributors, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);
DECLARE_SOA_COLUMN(Vtz, vtz, float);
DECLARE_SOA_COLUMN(AmpFV0A, ampFV0A, float);
DECLARE_SOA_COLUMN(AmpFT0A, ampFT0A, float);
DECLARE_SOA_COLUMN(AmpFT0C, ampFT0C, float);
DECLARE_SOA_COLUMN(ZdcTimeZNA, zdcTimeZNA, float);
DECLARE_SOA_COLUMN(ZdcTimeZNC, zdcTimeZNC, float);
} // namespace full
DECLARE_SOA_TABLE(HfUpcQa, "AOD", "HFUPCQA",
                  full::PvContributors,
                  full::Multiplicity,
                  full::Vtz,
                  full::AmpFV0A,
                  full::AmpFT0A,
                  full::AmpFT0C,
                  full::ZdcTimeZNA,
                  full::ZdcTimeZNC);

DECLARE_SOA_TABLE(HfUpcLcBdtInfos, "AOD", "HFUPCLCBDTINFOS",
                  full::M,
                  full::Pt,
                  full::BkgScore,
                  full::PromptScore,
                  full::FdScore,
                  full::AmpFV0A,
                  full::AmpFT0A,
                  full::AmpFT0C,
                  full::ZdcTimeZNA,
                  full::ZdcTimeZNC);

DECLARE_SOA_TABLE(HfUpcLcInfos, "AOD", "HFUPCLCINFOS",
                  full::M,
                  full::Pt,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  full::Chi2PCA,
                  full::DecayLength,
                  full::Cpa,
                  full::AmpFV0A,
                  full::AmpFT0A,
                  full::AmpFT0C,
                  full::ZdcTimeZNA,
                  full::ZdcTimeZNC);
} // namespace o2::aod

/// Λc± → p± K∓ π± analysis task
struct HfTaskUpcLc {
  Produces<o2::aod::HfUpcLcBdtInfos> rowCandUpcBdt;
  Produces<o2::aod::HfUpcLcInfos> rowCandUpc;
  Produces<o2::aod::HfUpcQa> rowUpcQa;

  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fillTreeOnlySingleGap{"fillTreeOnlySingleGap", false, "Only fill the tree for candidates that pass the single-gap UPC events"};
  Configurable<bool> fillTreeUpcQa{"fillTreeUpcQa", false, "Fill Tree for UPC QA"};
  Configurable<bool> verticesWithUpc{"verticesWithUpc", false, "Consider vertices with UPC settings"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfEventSelection hfEvSel;         // event selection and monitoring
  HfUpcGapThresholds upcThresholds; // UPC gap determination thresholds
  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;

  using LcCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using LcCandidatesMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;
  Preslice<aod::HfCand3Prong> candLcPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

  HistogramRegistry registry{"registry", {}};

  enum MlClasses : int {
    MlClassBackground = 0,
    MlClassPrompt,
    MlClassNonPrompt,
    NumberOfMlClasses
  };

  void init(InitContext&)
  {
    const std::array<bool, 2> doprocess{doprocessDataWithMlWithUpc, doprocessDataStdWithUpc};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

    auto vbins = (std::vector<double>)binsPt;
    registry.add("Data/fitInfo/ampFT0A_vs_ampFT0C", "FT0-A vs FT0-C amplitude;FT0-A amplitude (a.u.);FT0-C amplitude (a.u.)", {HistType::kTH2F, {{500, 0., 500}, {500, 0., 500}}});
    registry.add("Data/zdc/energyZNA_vs_energyZNC", "ZNA vs ZNC common energy;E_{ZNA}^{common} (a.u.);E_{ZNC}^{common} (a.u.)", {HistType::kTH2F, {{100, 0., 10}, {100, 0., 10}}});
    registry.add("Data/zdc/timeZNA_vs_timeZNC", "ZNA vs ZNC time;ZNA Time;ZNC time", {HistType::kTH2F, {{200, -10., 10}, {200, -10., 10}}});
    registry.add("Data/hUpcGapAfterSelection", "UPC gap type after selection;Gap side;Counts", {HistType::kTH1F, {{7, -1.5, 5.5}}});
    registry.add("Data/hUpcMulti", "Multiplicity of UPC events;Multiplicity;Counts", {HistType::kTH1F, {{200, -0.5, 199.5}}});
    registry.add("Data/hUpcVtz", "Vertex Z position of UPC events;Vz (cm);Counts", {HistType::kTH1F, {{200, -10., 10.}}});

    hfEvSel.addHistograms(registry);
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param collision is collision
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  template <bool FillMl, typename CollType, typename CandType, typename BCsType>
  void runAnalysisPerCollisionDataWithUpc(CollType const& collisions,
                                          CandType const& candidates,
                                          BCsType const& bcs,
                                          aod::FT0s const& ft0s,
                                          aod::FV0As const& fv0as,
                                          aod::FDDs const& fdds

  )
  {
    for (const auto& collision : collisions) {
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMaskWithUpc<true, CentralityEstimator::None, BCsType>(collision, centrality, ccdb, registry, bcs);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }
      const auto thisCollId = collision.globalIndex();
      const auto& groupedLcCandidates = candidates.sliceBy(candLcPerCollision, thisCollId);
      const auto numPvContributors = collision.numContrib();
      const auto& bc = collision.template bc_as<BCsType>();

      // Determine gap type using SGSelector with BC range checking
      const auto gapResult = hf_upc::determineGapType(collision, bcs, upcThresholds);
      const int gap = gapResult.value;
      const int upcFlag = (collision.flags() & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) ? 1 : 0;

      // Use the BC with FIT activity if available from SGSelector
      auto bcForUPC = bc;
      if (gapResult.bc) {
        bcForUPC = *(gapResult.bc);
      }

      // Get FIT information from the UPC BC
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, bcForUPC, bcs, ft0s, fv0as, fdds);

      // Get ZDC energies if available (extract once and reuse)
      const bool hasZdc = bcForUPC.has_zdc();
      float zdcEnergyZNA = -1.f;
      float zdcEnergyZNC = -1.f;
      float zdcTimeZNA = -1.f;
      float zdcTimeZNC = -1.f;
      if (verticesWithUpc && !upcFlag) {
        continue;
      }
      if (hasZdc) {
        const auto zdc = bcForUPC.zdc();
        zdcEnergyZNA = zdc.energyCommonZNA();
        zdcEnergyZNC = zdc.energyCommonZNC();
        zdcTimeZNA = zdc.timeZNA();
        zdcTimeZNC = zdc.timeZNC();
        registry.fill(HIST("Data/fitInfo/ampFT0A_vs_ampFT0C"), fitInfo.ampFT0A, fitInfo.ampFT0C);
        registry.fill(HIST("Data/zdc/energyZNA_vs_energyZNC"), zdcEnergyZNA, zdcEnergyZNC);
        registry.fill(HIST("Data/zdc/timeZNA_vs_timeZNC"), zdcTimeZNA, zdcTimeZNC);
        registry.fill(HIST("Data/hUpcGapAfterSelection"), static_cast<int>(gap));
      }
      if (gap == o2::aod::sgselector::TrueGap::SingleGapA || gap == o2::aod::sgselector::TrueGap::SingleGapC) {
        registry.fill(HIST("Data/hUpcMulti"), collision.multNTracksPV());
        registry.fill(HIST("Data/hUpcVtz"), collision.posZ());
      }
      if (fillTreeUpcQa) {
        rowUpcQa(numPvContributors, collision.multNTracksPV(), collision.posZ(), fitInfo.ampFV0A, fitInfo.ampFT0A, fitInfo.ampFT0C, zdcTimeZNA, zdcTimeZNC);
      }

      for (const auto& candidate : groupedLcCandidates) {
        if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
          continue;
        }
        if (yCandRecoMax >= 0. && std::abs(HfHelper::yLc(candidate)) > yCandRecoMax) {
          continue;
        }
        const auto pt = candidate.pt();
        const auto ptProng0 = candidate.ptProng0();
        const auto ptProng1 = candidate.ptProng1();
        const auto ptProng2 = candidate.ptProng2();
        const auto decayLength = candidate.decayLength();
        const auto chi2PCA = candidate.chi2PCA();
        const auto cpa = candidate.cpa();

        double outputBkg(-1), outputPrompt(-1), outputFD(-1);

        auto fillTHnData = [&](bool isPKPi) {
          const auto massLc = isPKPi ? HfHelper::invMassLcToPKPi(candidate) : HfHelper::invMassLcToPiKP(candidate);

          if constexpr (FillMl) {
            const auto& mlProb = isPKPi ? candidate.mlProbLcToPKPi() : candidate.mlProbLcToPiKP();
            if (mlProb.size() == NumberOfMlClasses) {
              outputBkg = mlProb[MlClassBackground]; /// bkg score
              outputPrompt = mlProb[MlClassPrompt];  /// prompt score
              outputFD = mlProb[MlClassNonPrompt];   /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            if (fillTreeOnlySingleGap) {
              if (gap == o2::aod::sgselector::TrueGap::SingleGapA || gap == o2::aod::sgselector::TrueGap::SingleGapC) {
                rowCandUpcBdt(massLc, pt, outputBkg, outputPrompt, outputFD, fitInfo.ampFV0A, fitInfo.ampFT0A, fitInfo.ampFT0C, zdcTimeZNA, zdcTimeZNC);
              }
            } else {
              rowCandUpcBdt(massLc, pt, outputBkg, outputPrompt, outputFD, fitInfo.ampFV0A, fitInfo.ampFT0A, fitInfo.ampFT0C, zdcTimeZNA, zdcTimeZNC);
            }

          } else {
            if (fillTreeOnlySingleGap) {
              if (gap == o2::aod::sgselector::TrueGap::SingleGapA || gap == o2::aod::sgselector::TrueGap::SingleGapC) {
                rowCandUpc(massLc, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, fitInfo.ampFV0A, fitInfo.ampFT0A, fitInfo.ampFT0C, zdcTimeZNA, zdcTimeZNC);
              }
            } else {
              rowCandUpc(massLc, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, fitInfo.ampFV0A, fitInfo.ampFT0A, fitInfo.ampFT0C, zdcTimeZNA, zdcTimeZNC);
            }
          }
        };

        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          fillTHnData(true);
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          fillTHnData(false);
        }
      }
    }
  }

  void processDataWithMlWithUpc(soa::Join<aod::Collisions, aod::EvSels, aod::Mults> const& collisions,
                                aod::BcFullInfos const& bcs,
                                LcCandidatesMl const& selectedLcCandidatesMl,
                                aod::Tracks const&,
                                aod::FT0s const& ft0s,
                                aod::FV0As const& fv0as,
                                aod::FDDs const& fdds,
                                aod::Zdcs const& /*zdcs*/)
  {
    runAnalysisPerCollisionDataWithUpc<true>(collisions, selectedLcCandidatesMl, bcs, ft0s, fv0as, fdds);
  }
  PROCESS_SWITCH(HfTaskUpcLc, processDataWithMlWithUpc, "Process real data with the ML method with UPC", false);

  void processDataStdWithUpc(soa::Join<aod::Collisions, aod::EvSels, aod::Mults> const& collisions,
                             aod::BcFullInfos const& bcs,
                             LcCandidates const& selectedLcCandidates,
                             aod::Tracks const&,
                             aod::FT0s const& ft0s,
                             aod::FV0As const& fv0as,
                             aod::FDDs const& fdds,
                             aod::Zdcs const& /*zdcs*/)
  {
    runAnalysisPerCollisionDataWithUpc<false>(collisions, selectedLcCandidates, bcs, ft0s, fv0as, fdds);
  }
  PROCESS_SWITCH(HfTaskUpcLc, processDataStdWithUpc, "Process real data with the standard method with UPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskUpcLc>(cfgc)};
}
