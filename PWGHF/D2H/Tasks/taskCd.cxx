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

/// \file taskCd.cxx
/// \brief Cd± → d± K∓ π±  analysis task
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg Universiity

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>

#include <Rtypes.h>

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
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;

namespace o2::aod
{
namespace full
{
// Candidate kinematics
DECLARE_SOA_COLUMN(M, m, float);                                //! Invariant mass of candidate (GeV/c^2)
DECLARE_SOA_COLUMN(Pt, pt, float);                              //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                  //! Transverse momentum of prong 0 (GeV/c)
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                  //! Transverse momentum of prong 1 (GeV/c)
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);                  //! Transverse momentum of prong 2 (GeV/c)
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);  //! Impact parameter (DCA to PV) of prong 0 (cm)
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);  //! Impact parameter (DCA to PV) of prong 1 (cm)
DECLARE_SOA_COLUMN(ImpactParameter2, impactParameter2, float);  //! Impact parameter (DCA to PV) of prong 2 (cm)
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);            //! Decay length (3D) of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);        //! Decay length in transverse plane (cm)
DECLARE_SOA_COLUMN(Cpa, cpa, float);                            //! Cosine of pointing angle (3D)
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                        //! Cosine of pointing angle in XY plane
DECLARE_SOA_COLUMN(NSigmaTpcDe, nSigmaTpcDe, float);            //! TPC nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcKa, nSigmaTpcKa, float);            //! TPC nσ for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTpcPi, nSigmaTpcPi, float);            //! TPC nσ for pion hypothesis
DECLARE_SOA_COLUMN(NSigmaItsDe, nSigmaItsDe, float);            //! ITS nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTofDe, nSigmaTofDe, float);            //! TOF nσ for deuteron hypothesis
DECLARE_SOA_COLUMN(NSigmaTofKa, nSigmaTofKa, float);            //! TOF nσ for kaon hypothesis
DECLARE_SOA_COLUMN(NSigmaTofPi, nSigmaTofPi, float);            //! TOF nσ for pion hypothesis
DECLARE_SOA_COLUMN(NItsClusters, nItsClusters, int8_t);         //! Number of ITS clusters used in the track fit
DECLARE_SOA_COLUMN(NItsNClusterSize, nItsNClusterSize, int8_t); //! Number of ITS clusters size used in the track fit
DECLARE_SOA_COLUMN(NTpcClusters, nTpcClusters, int8_t);         //! Number of TPC clusters used in the track fit
DECLARE_SOA_COLUMN(NTpcSignalsDe, nTpcSignalsDe, int8_t);       //! Number of TPC signas
DECLARE_SOA_COLUMN(NItsSignalsDe, nItsSignalsDe, int8_t);       //! Number of ITS signas
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t); //! Candidates falg
DECLARE_SOA_COLUMN(Cent, cent, float);                          //! Centrality
DECLARE_SOA_COLUMN(GIndexCol, gIndexCol, int);                  //! Global index for the collisionAdd commentMore actions
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, int64_t);              //! Timestamp for the collision
} // namespace full

// Full table: include ALL columns declared above
DECLARE_SOA_TABLE(HfCandCd, "AOD", "HFCANDCD",
                  full::M,
                  full::Pt,
                  full::PtProng0,
                  full::PtProng1,
                  full::PtProng2,
                  full::ImpactParameter0,
                  full::ImpactParameter1,
                  full::ImpactParameter2,
                  full::DecayLength,
                  full::Cpa,
                  full::NSigmaTpcDe,
                  full::NSigmaItsDe,
                  full::NSigmaTofDe,
                  full::NItsClusters,
                  full::NItsNClusterSize,
                  full::NTpcClusters,
                  full::NTpcSignalsDe,
                  full::NItsSignalsDe,
                  full::CandidateSelFlag,
                  full::Cent,
                  full::GIndexCol,
                  full::TimeStamp);
} // namespace o2::aod

struct HfTaskCd {

  Produces<o2::aod::HfCandCd> rowCandCd;
  Configurable<int> selectionFlagCd{"selectionFlagCd", 1, "Selection Flag for Cd"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_cd_to_de_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fillTHn{"fillTHn", false, "fill THn"};
  Configurable<bool> fillTree{"fillTree", false, "Flag to fill candiates tree"};

  SliceCache cache;

  using CollisionsWEvSel = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsWithEvSelFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithEvSelFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;

  using CdCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelCd, aod::HfCand3ProngWPidPiKaDe>>;
  using HFTracks = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullDe>;

  Filter filterSelectCandidates = aod::hf_sel_candidate_cd::isSelCdToDeKPi >= selectionFlagCd;
  Preslice<aod::HfCand3Prong> candCdPerCollision = aod::hf_cand::collisionId;

  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {72, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {400, 2.4, 4.4}, ""};
  ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisCentrality{"thnConfigAxisCentrality", {100, 0, 100}, ""};

  HistogramRegistry registry{
    "registry",
    {/// mass candidate
     {"Data/hMass", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{400, 2.4, 4.4}}}},
     /// pT
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// cosine of pointing angle
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// cosine of pointing angle xy
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// Chi 2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    std::array<bool, 3> doprocess{doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }
    /// mass candidate
    registry.add("Data/hMassVsPtVsNPvContributors", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}; Number of PV contributors", {HistType::kTH3F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}, {500, 0., 5000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (de K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{400, 2.4, 4.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    registry.add("Data/hImpParErrProng0", "3-prong candidates;prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng1", "3-prong candidates;prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng2", "3-prong candidates;prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hNsigmaTPCDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{TPC}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{TOF}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaITSDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{ITS}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hTPCSignalDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{ITS}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {2000, 0, 2000}}});
    registry.add("Data/hITSSignalDeVsP", "deuteron;#it{p} (GeV/#it{c}); n#sigma^{ITS}_{d}", {HistType::kTH2F, {{200, -10.f, 10.f}, {20, 0, 20}}});
    registry.add("Data/hNsigmaTPCPiVsP", "Pion;#it{p} (GeV/#it{c});n#sigma^{TPC}_{pi};", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFPiVsP", "Pion;#it{p} (GeV/#it{c});n#sigma^{TOF}_{pi};", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTPCKaVsP", "Kaon;#it{p} (GeV/#it{c}); n#sigma^{TPC}_{Kaon}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    registry.add("Data/hNsigmaTOFKaVsP", "Kaon;#it{p} (GeV/#it{c}); n#sigma^{TOF}_{Kaon}", {HistType::kTH2F, {{200, -10.f, 10.f}, {200, -6.f, 6.f}}});
    if (fillTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (de K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(C_{d}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng0{thnConfigAxisPtProng, "#it{p}_{T}(prong0) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng1{thnConfigAxisPtProng, "#it{p}_{T}(prong1) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng2{thnConfigAxisPtProng, "#it{p}_{T}(prong2) (GeV/#it{c})"};
      const AxisSpec thnAxisChi2PCA{thnConfigAxisChi2PCA, "Chi2PCA to sec. vertex (cm)"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle"};
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, "centrality (FT0C)"};

      std::vector axesStd{thnAxisMass, thnAxisPt, thnAxisPtProng0, thnAxisPtProng1, thnAxisPtProng2, thnAxisChi2PCA, thnAxisDecLength, thnAxisCPA, thnAxisCentrality};
      registry.add("hnCdVars", "THn for Reconstructed Cd candidates for data", HistType::kTHnSparseF, axesStd);
    }
  }

  // taken from: https://github.com/AliceO2Group/O2Physics/blob/master/EventFiltering/PWGCF/CFFilterAll.cxx
  template <typename T>
  float itsSignal(T const& track)
  {
    uint32_t clsizeflag = track.itsClusterSizes();
    auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
    auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
    auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
    auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
    auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
    auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
    auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
    int numLayers = 7;
    int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
    float cosLamnda = 1. / std::cosh(track.eta());
    return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
  };

  /// Fill histograms for real data
  template <typename CollType, typename CandType, typename TrackType, typename TrackWithItsType, typename BcType>
  void fillHistosData(CollType const& collision, CandType const& candidates, TrackType const& /*tracks*/, TrackWithItsType const& tracksWithItsPid, BcType const& /*bcs*/)
  {
    auto thisCollId = collision.globalIndex();
    auto groupedCdCandidates = candidates.sliceBy(candCdPerCollision, thisCollId);
    auto numPvContributors = collision.numContrib();
    auto bc = collision.template bc_as<BcType>();
    int64_t timeStamp = bc.timestamp();

    for (const auto& candidate : groupedCdCandidates) {
      if (!TESTBIT(candidate.hfflag(), aod::hf_cand_3prong::DecayType::CdToDeKPi)) {
        continue;
      }

      const auto pt = candidate.pt();
      const auto ptProng0 = candidate.ptProng0();
      const auto ptProng1 = candidate.ptProng1();
      const auto ptProng2 = candidate.ptProng2();
      const auto decayLength = candidate.decayLength();
      const auto decayLengthXY = candidate.decayLengthXY();
      const auto chi2PCA = candidate.chi2PCA();
      const auto cpa = candidate.cpa();
      const auto cpaXY = candidate.cpaXY();
      float invMassCd = 0.f;
      if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
        invMassCd = HfHelper::invMassCdToDeKPi(candidate);
      }
      if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
        invMassCd = HfHelper::invMassCdToPiKDe(candidate);
      }

      if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToDeKPi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToDeKPi(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToDeKPi(candidate), pt);
      }
      if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
        registry.fill(HIST("Data/hMass"), HfHelper::invMassCdToPiKDe(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), HfHelper::invMassCdToPiKDe(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), HfHelper::invMassCdToPiKDe(candidate), pt);
      }
      registry.fill(HIST("Data/hPt"), pt);
      registry.fill(HIST("Data/hPtProng0"), ptProng0);
      registry.fill(HIST("Data/hPtProng1"), ptProng1);
      registry.fill(HIST("Data/hPtProng2"), ptProng2);
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), pt);
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), pt);
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), pt);
      registry.fill(HIST("Data/hDecLength"), decayLength);
      registry.fill(HIST("Data/hDecLengthVsPt"), decayLength, pt);
      registry.fill(HIST("Data/hDecLengthxy"), decayLengthXY);
      registry.fill(HIST("Data/hDecLengthxyVsPt"), decayLengthXY, pt);
      registry.fill(HIST("Data/hCPA"), cpa);
      registry.fill(HIST("Data/hCPAVsPt"), cpa, pt);
      registry.fill(HIST("Data/hCPAxy"), cpaXY);
      registry.fill(HIST("Data/hCPAxyVsPt"), cpaXY, pt);
      registry.fill(HIST("Data/hDca2"), chi2PCA);
      registry.fill(HIST("Data/hDca2VsPt"), chi2PCA, pt);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), pt);
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToDeKPi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelCdToPiKDe(), pt);
      registry.fill(HIST("Data/hImpParErrProng0"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2"), candidate.errorImpactParameter2(), pt);

      float const cent = o2::hf_centrality::getCentralityColl(collision);

      if (fillTHn) {
        double massCd(-1);
        if (candidate.isSelCdToDeKPi() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToDeKPi(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
        if (candidate.isSelCdToPiKDe() >= selectionFlagCd) {
          massCd = HfHelper::invMassCdToPiKDe(candidate);
          std::vector<double> valuesToFill{massCd, pt, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, cent};
          registry.get<THnSparse>(HIST("hnCdVars"))->Fill(valuesToFill.data());
        }
      }

      if (fillTree) {
        int candFlag = -999;

        float nSigmaTpcDe = 0.f, nSigmaTpcKa = 0.f, nSigmaTpcPi = 0.f;
        float nSigmaItsDe = 0.f;
        float nSigmaTofDe = 0.f, nSigmaTofKa = 0.f, nSigmaTofPi = 0.f;

        int itsNClusterDe = 0;
        int itsNClusterSizeDe = 0;
        int tpcNClusterDe = 0;

        float tpcSignalsDe = 0.f;
        float itsSignalsDe = 0.f;

        float pSignedDe = -999.f;
        float pSignedPi = -999.f;

        nSigmaTpcKa = candidate.nSigTpcKa1();
        nSigmaTofKa = candidate.nSigTofKa1();

        const bool selDeKPi = (candidate.isSelCdToDeKPi() >= 1);
        const bool selPiKDe = (candidate.isSelCdToPiKDe() >= 1);

        auto prong0 = candidate.template prong0_as<TrackType>();
        auto prong1 = candidate.template prong1_as<TrackType>();
        auto prong2 = candidate.template prong2_as<TrackType>();

        auto prong0Its = tracksWithItsPid.iteratorAt(candidate.prong0Id() - tracksWithItsPid.offset());
        auto prong2Its = tracksWithItsPid.iteratorAt(candidate.prong2Id() - tracksWithItsPid.offset());

        if (selDeKPi) {
          candFlag = 1;
          pSignedDe = prong0.p() * prong0.sign();
          pSignedPi = prong2.p() * prong2.sign();
          nSigmaTpcDe = candidate.nSigTpcDe0();
          nSigmaTofDe = candidate.nSigTofDe0();
          nSigmaTpcPi = candidate.nSigTpcPi2();
          nSigmaTofPi = candidate.nSigTofPi2();
          nSigmaItsDe = prong0Its.itsNSigmaDe();
          itsNClusterDe = prong0.itsNCls();
          itsNClusterSizeDe = prong0.itsClusterSizes();
          tpcNClusterDe = prong0.tpcNClsCrossedRows();
          tpcSignalsDe = prong0.tpcSignal();
          itsSignalsDe = itsSignal(prong0);
        } else if (selPiKDe) {
          candFlag = -1;
          pSignedDe = prong2.p() * prong2.sign();
          pSignedPi = prong0.p() * prong0.sign();
          nSigmaTpcDe = candidate.nSigTpcDe2();
          nSigmaTofDe = candidate.nSigTofDe2();
          nSigmaTpcPi = candidate.nSigTpcPi0();
          nSigmaTofPi = candidate.nSigTofPi0();
          nSigmaItsDe = prong2Its.itsNSigmaDe();
          itsNClusterDe = prong2.itsNCls();
          itsNClusterSizeDe = prong2.itsClusterSizes();
          tpcNClusterDe = prong2.tpcNClsCrossedRows();
          tpcSignalsDe = prong2.tpcSignal();
          itsSignalsDe = itsSignal(prong2);
        }

        //  PID QA
        registry.fill(HIST("Data/hNsigmaTPCDeVsP"), pSignedDe, nSigmaTpcDe);
        registry.fill(HIST("Data/hNsigmaTOFDeVsP"), pSignedDe, nSigmaTofDe);
        registry.fill(HIST("Data/hNsigmaITSDeVsP"), pSignedDe, nSigmaItsDe);
        registry.fill(HIST("Data/hTPCSignalDeVsP"), pSignedDe, tpcSignalsDe);
        registry.fill(HIST("Data/hITSSignalDeVsP"), pSignedDe, itsSignalsDe);
        registry.fill(HIST("Data/hNsigmaTPCPiVsP"), pSignedPi, nSigmaTpcPi);
        registry.fill(HIST("Data/hNsigmaTOFPiVsP"), pSignedPi, nSigmaTofPi);
        registry.fill(HIST("Data/hNsigmaTPCKaVsP"), prong1.p() * prong1.sign(), nSigmaTpcKa);
        registry.fill(HIST("Data/hNsigmaTOFKaVsP"), prong1.p() * prong1.sign(), nSigmaTofKa);

        rowCandCd(
          invMassCd,
          pt,
          ptProng0,
          ptProng1,
          ptProng2,
          candidate.impactParameter0(),
          candidate.impactParameter1(),
          candidate.impactParameter2(),
          decayLength,
          cpa,
          nSigmaTpcDe,
          nSigmaItsDe,
          nSigmaTofDe,
          itsNClusterDe,
          itsNClusterSizeDe,
          tpcNClusterDe,
          tpcSignalsDe,
          itsSignalsDe,
          candFlag,
          cent,
          collision.globalIndex(),
          timeStamp);
      }
    }
  }
  /// Run the analysis on real data
  template <typename CollType, typename CandType, typename TrackType, typename TrackWithItsType, typename BcType>
  void runAnalysisPerCollisionData(CollType const& collisions,
                                   CandType const& candidates,
                                   TrackType const& tracks,
                                   TrackWithItsType const& tracksWithItsPid,
                                   BcType const& bcs)
  {

    for (const auto& collision : collisions) {
      fillHistosData(collision, candidates, tracks, tracksWithItsPid, bcs);
    }
  }

  void processDataStd(CollisionsWEvSel const& collisions,
                      CdCandidates const& selectedCdCandidates,
                      HFTracks const& tracks,
                      aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStd, "Process Data with the standard method", true);

  void processDataStdWithFT0C(CollisionsWithEvSelFT0C const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              HFTracks const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0C, "Process real data with the standard method and with FT0C centrality", false);

  void processDataStdWithFT0M(CollisionsWithEvSelFT0M const& collisions,
                              CdCandidates const& selectedCdCandidates,
                              HFTracks const& tracks,
                              aod::BCsWithTimestamps const& bcWithTimeStamps)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<HFTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);
    runAnalysisPerCollisionData(collisions, selectedCdCandidates, tracks, tracksWithItsPid, bcWithTimeStamps);
  }
  PROCESS_SWITCH(HfTaskCd, processDataStdWithFT0M, "Process real data with the standard method and with FT0M centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCd>(cfgc)};
}
