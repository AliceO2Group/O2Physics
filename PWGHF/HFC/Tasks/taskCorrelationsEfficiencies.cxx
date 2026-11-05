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

/// \file taskCorrelationsEfficiencies.cxx
/// \brief Charm hadrons and associated tracks efficiency evaluation for correlations studies
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"
#include "PWGHF/HFC/Utils/utilsFlow.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
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
#include <Framework/runDataProcessing.h>

#include <TF1.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::hfc_utils_flow;

namespace o2::aod
{
namespace hf_correlations_eff
{
// Charm-hadron candidate columns
DECLARE_SOA_COLUMN(Pt, pt, float);             //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PvContrib, pvContrib, int); //! Number of primary vertex contributors
DECLARE_SOA_COLUMN(Eta, eta, float);           //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);           //! Azimuthal angle of candidate
DECLARE_SOA_COLUMN(Origin, origin, int);       //! Origin of the candidate (prompt or non-prompt)
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float); //! ML score of the first configured index
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float); //! ML score of the second configured index

// Associated track columns
DECLARE_SOA_COLUMN(PtAssoc, ptAssoc, float);         //! Transverse momentum of associated track (GeV/c)
DECLARE_SOA_COLUMN(EtaAssoc, etaAssoc, float);       //! Pseudorapidity of associated track
DECLARE_SOA_COLUMN(PhiAssoc, phiAssoc, float);       //! Azimuthal angle of associated track
DECLARE_SOA_COLUMN(ZVtx, zVtx, float);               //! Primary vertex z position (cm)
DECLARE_SOA_COLUMN(AssTrkSpecie, assTrkSpecie, int); //! Associated track species (0=pi, 1=K, 2=p, 3=el, 4=mu)
DECLARE_SOA_COLUMN(AssTrkOrigin, assTrkOrigin, int); //! Associated track origin (0=primary, 1=secondary)
} // namespace hf_correlations_eff

DECLARE_SOA_TABLE(HfcMcrigInfo, "AOD", "HFCMCTRIGINFO",
                  hf_correlations_eff::Pt,
                  hf_correlations_eff::PvContrib,
                  hf_correlations_eff::Eta,
                  hf_correlations_eff::Phi,
                  hf_correlations_eff::Origin,
                  hf_correlations_eff::MlScore0,
                  hf_correlations_eff::MlScore1);

DECLARE_SOA_TABLE(HfcMcAssocInfo, "AOD", "HFCMCASSOCINFO",
                  hf_correlations_eff::PtAssoc,
                  hf_correlations_eff::EtaAssoc,
                  hf_correlations_eff::PhiAssoc,
                  hf_correlations_eff::ZVtx,
                  hf_correlations_eff::AssTrkSpecie,
                  hf_correlations_eff::AssTrkOrigin);
} // namespace o2::aod

/// Code to obtain efficiency correction factors for (Delta Eta, Delta Phi) correlation maps
struct HfTaskCorrelationsEfficiencies {

  Produces<o2::aod::HfcMcrigInfo> rowTrigMc;
  Produces<o2::aod::HfcMcAssocInfo> rowAssocMc;
  Produces<o2::aod::HfcRedAssTracks> rowAssocTrkSels;

  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for hadron (ML score tables are required to run the task)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 70, "min. TPC crossed rows for associated tracks"};
  Configurable<float> etaTrkMax{"etaTrkMax", 1., "max. track eta"};
  Configurable<float> ptTrkMin{"ptTrkMin", 0.3, "min. track pT"};
  Configurable<float> ptTrkMax{"ptTrkMax", 3., "max. track pT"};
  Configurable<float> dcaXYTrkMax{"dcaXYTrkMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrkMax{"dcaZTrkMax", 1., "max. track DCA Z"};
  Configurable<bool> usePtDiffDcaXYCut{"usePtDiffDcaXYCut", false, "Use pt-differential DCAxy cut for associated tracks"};
  Configurable<float> dcaXYTrkNSigmaMax{"dcaXYTrkNSigmaMax", 7, "Cut on number of sigma deviations from expected DCA in the transverse direction"};
  Configurable<std::string> dcaXYPtPrimTrkFunc{"dcaXYPtPrimTrkFunc", "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut"};

  // Configurables for tree filling
  Configurable<bool> fillRecoTrees{"fillRecoTrees", false, "Fill trees of reconstructed charmed hadrons and associated tracks"};
  Configurable<float> downSampleTracks{"downSampleTracks", 1., "Fraction of tracks to keep"};
  Configurable<float> ptMaxDownSampleTracks{"ptMaxDownSampleTracks", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> downSampleCharm{"downSampleCharm", 1., "Fraction of pairs to keep"};
  Configurable<float> ptMaxDownSampleCharm{"ptMaxDownSampleCharm", 10., "Maximum pt for the application of the downsampling factor"};

  HfcFlowHelper hfcFlowHelper; // flow utilities
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  SliceCache cache;

  std::unique_ptr<TF1> funcDcaXYPtCutPrimTrk{nullptr};

  using RecoCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::PVMults>;
  using McCandsDs = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi, aod::HfCand3ProngMcRec>>;
  using McCandsDplus = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi, aod::HfCand3ProngMcRec>>;
  using McCandsD0 = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0, aod::HfCand2ProngMcRec>>;
  using McCandsLc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi, aod::HfCand3ProngMcRec>>;
  using McTracks = soa::Filtered<soa::Join<aod::McTrackLabels, aod::TracksWDca, aod::TrackSelection, aod::TracksExtra>>;

  Filter filterSelectDsCands = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCands = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Cands = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectLcCands = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Filter filterSelectTrk = (nabs(aod::track::eta) < etaTrkMax) && (aod::track::pt > ptTrkMin) && (aod::track::pt < ptTrkMax) && (nabs(aod::track::dcaXY) < dcaXYTrkMax) && (nabs(aod::track::dcaZ) < dcaZTrkMax);

  Preslice<McCandsDs> candsDsPerColl = aod::hf_cand::collisionId;
  Preslice<McCandsDplus> candsDplusPerColl = aod::hf_cand::collisionId;
  Preslice<McCandsD0> candsD0PerColl = aod::hf_cand::collisionId;
  Preslice<McCandsLc> candsLcPerColl = aod::hf_cand::collisionId;
  Preslice<McTracks> trackIndicesPerColl = aod::track::collisionId;

  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  ConfigurableAxis binsPtCharm{"binsPtCharm", {VARIABLE_WIDTH, 0., 2., 4., 6., 8., 24.}, ""};
  ConfigurableAxis binsPtAssoc{"binsPtAssoc", {20, 0, 10}, ""};
  ConfigurableAxis binsZVtx{"binsZVtx", {120, -11, 11}, ""};
  ConfigurableAxis binsPvContrib{"binsPvContrib", {5001, -0.5, 5000.5}, ""};
  ConfigurableAxis binsEta{"binsEta", {50, -2., 2.}, ""};
  ConfigurableAxis binsPhi{"binsPhi", {64, 0, 6.28}, ""};
  ConfigurableAxis binsMlOne{"binsMlOne", {100, 0., 1.}, ""};
  ConfigurableAxis binsMlTwo{"binsMlTwo", {100, 0., 1.}, ""};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    std::vector<unsigned int> classesMl;
    for (const auto scoreClass : classMl.value) {
      if (scoreClass < 0) {
        LOG(fatal) << "ML class indexes must be non-negative";
        return;
      }
      classesMl.push_back(static_cast<unsigned int>(scoreClass));
    }

    if (classesMl.size() > 2) {
      LOG(fatal) << "At most two ML classes can be selected";
      return;
    }

    hfcFlowHelper.init(ccdb, registry, classesMl, yCandMax.value);

    // Setup pt-dependent DCAxy cut function
    if (usePtDiffDcaXYCut) {
      funcDcaXYPtCutPrimTrk = std::make_unique<TF1>("funcDcaXYPtCutPrimTrk", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()), 0.001, 100.);
      funcDcaXYPtCutPrimTrk->SetParameter(0, dcaXYTrkNSigmaMax);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()));
    }

    // Add efficiency sparses
    // Axes: mass, pt, pv contributors, pseudo-rapidity, ml_score_0, ml_score_1
    const AxisSpec axisPtTrig = {binsPtCharm, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPvContrib = {binsPvContrib, "PV contributors"};
    const AxisSpec axisEtaTrig = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhiTrig = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisMlScore0 = {binsMlOne, "ML score 0"};
    const AxisSpec axisMlScore1 = {binsMlTwo, "ML score 1"};
    std::vector<AxisSpec> axesCharmReco = {axisPtTrig, axisPvContrib, axisEtaTrig, axisPhiTrig, axisMlScore0, axisMlScore1};

    registry.add("hCharmPromptGen", "Generated prompt charm-hadron candidates;;", HistType::kTHnSparseF, {{axisPtTrig}, {axisEtaTrig}, {axisPhiTrig}});
    registry.add("hCharmFdGen", "Generated non-prompt charm-hadron candidates;;", HistType::kTHnSparseF, {{axisPtTrig}, {axisEtaTrig}, {axisPhiTrig}});
    registry.add("hCharmPromptReco", "Reconstructed prompt charm-hadron candidates;;", HistType::kTHnSparseF, axesCharmReco);
    registry.add("hCharmFdReco", "Reconstructed non-prompt charm-hadron candidates;;", HistType::kTHnSparseF, axesCharmReco);

    // Axes: pt, eta, phi, z_vtx, specie
    const AxisSpec axisPtAssoc = {binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisEtaAssoc = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhiAssoc = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisZVtx = {binsZVtx, "#it{z}_{vtx}"};
    const AxisSpec axisAssTrkSpecie = {static_cast<int>(NAssTracksTypes), -0.5, static_cast<int>(NAssTracksTypes) - 0.5, "Particle species"};
    const AxisSpec axisTrackOrigin = {2, -0.5, 1.5, "Particle origin"};
    std::vector<AxisSpec> axesAssocTracks = {axisPtAssoc, axisEtaAssoc, axisPhiAssoc, axisZVtx, axisAssTrkSpecie, axisTrackOrigin};
    registry.add("hAssTrackReco", "Reconstructed associated tracks;;", HistType::kTHnSparseF, axesAssocTracks);
    registry.add("hAssTrackGen", "Generated associated tracks;;", HistType::kTHnSparseF, axesAssocTracks);
  }

  /// Match charm-hadron based on its reconstruction flag
  /// \param assTrk is the associated track
  /// \return the associated track specie index
  template <bool IsReco, TriggerChannels CandType, typename TCharmHad>
  bool isMcMatched(TCharmHad const& cand)
  {
    int flagMc{-1};
    if constexpr (IsReco) {
      flagMc = std::abs(cand.flagMcMatchRec());
    } else {
      flagMc = std::abs(cand.flagMcMatchGen());
    }

    if constexpr (CandType == TriggerChannels::D0ToPiK) {
      return flagMc == static_cast<int8_t>(hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK);
    } else if constexpr (CandType == TriggerChannels::DplusToPiKPi) {
      return flagMc == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi);
    } else if constexpr (CandType == TriggerChannels::LcToPKPi) {
      return flagMc == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi);
    } else { // Ds
      int8_t reso{0};
      if constexpr (IsReco) {
        reso = std::abs(cand.flagMcDecayChanRec());
      } else {
        reso = std::abs(cand.flagMcDecayChanGen());
      }
      return flagMc == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) &&
             reso == static_cast<int8_t>(hf_decay::hf_cand_3prong::DecayChannelResonant::DsToPhiPi);
    }
  }

  /// Fill info for charm-hadron efficiencies
  /// \param trigCands are the trigger candidates
  /// \param assTrks are the associated tracks
  /// \param cent is the collision centrality
  /// \param mult is the collision multiplicity
  template <TriggerChannels CandType, typename TCharmCands>
  void fillCharmEfficiencies(TCharmCands const& charmCands, int numContrib)
  {
    for (const auto& cand : charmCands) {
      if (!hfcFlowHelper.cutCandRapidity<CandType>(cand)) {
        continue;
      }

      if (!isMcMatched<true, CandType>(cand)) {
        continue;
      }

      // Fill the sparse
      std::array<float, 2> outputMl = hfcFlowHelper.getCandMlScores<CandType>(cand);
      if (cand.originMcRec() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hCharmPromptReco"), cand.pt(), numContrib, cand.eta(), cand.phi(), outputMl[0], outputMl[1]);
      }
      if (cand.originMcRec() == RecoDecay::OriginType::NonPrompt) {
        registry.fill(HIST("hCharmFdReco"), cand.pt(), numContrib, cand.eta(), cand.phi(), outputMl[0], outputMl[1]);
      }

      if (fillRecoTrees) {
        // Downsample charm candidates based on pt
        if (downSampleCharm < 1.) {
          float const pseudoRndm = cand.pt() * 1000. - static_cast<int64_t>(cand.pt() * 1000);
          if (cand.pt() < ptMaxDownSampleCharm && pseudoRndm >= downSampleCharm) {
            continue;
          }
        }
        rowTrigMc(cand.pt(), numContrib, cand.eta(), cand.phi(), cand.originMcRec(), outputMl[0], outputMl[1]);
      }
    }
  }

  // Dplus efficiency
  void processEffDplus(RecoCollisions const& recoCollisions,
                       McCandsDplus const& candsDplus,
                       soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& genParticles,
                       aod::McCollisions const&)
  {
    for (const auto& coll : recoCollisions) {
      auto candsInColl = candsDplus.sliceBy(candsDplusPerColl, coll.globalIndex());
      if (candsInColl.size() < 1) {
        continue;
      }
      float cent{-1.}, mult{-1.};
      if (!hfcFlowHelper.checkCollision<o2::hf_centrality::CentralityEstimator::None>(coll, cent, mult)) {
        continue;
      }
      fillCharmEfficiencies<TriggerChannels::DplusToPiKPi>(candsInColl, coll.numContrib());
    }

    for (const auto& genCand : genParticles) {
      if (!genCand.has_mcCollision()) {
        continue;
      }
      if (std::abs(genCand.y()) > yCandMax) {
        continue;
      }
      if (!isMcMatched<false, TriggerChannels::DplusToPiKPi>(genCand)) {
        continue;
      }

      if (genCand.originMcGen() == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hCharmPromptGen"), genCand.pt(), genCand.eta(), genCand.phi());
      }
      if (genCand.originMcGen() == RecoDecay::OriginType::NonPrompt) {
        registry.fill(HIST("hCharmFdGen"), genCand.pt(), genCand.eta(), genCand.phi());
      }
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationsEfficiencies, processEffDplus, "Process efficiencies for Dplus candidates", false);

  /// Match associated tracks based on their PDG code
  /// \param assTrk is the associated track
  /// \return the associated track specie index
  template <typename TAssocTrk>
  int getAssocTrkSpecie(TAssocTrk const& assTrk)
  {
    switch (std::abs(assTrk.pdgCode())) {
      case kPiPlus:
        return static_cast<int>(AssTracksTypes::Pi);
      case kKPlus:
        return static_cast<int>(AssTracksTypes::Ka);
      case kProton:
        return static_cast<int>(AssTracksTypes::Pr);
      case kElectron:
        return static_cast<int>(AssTracksTypes::El);
      case kMuonMinus:
        return static_cast<int>(AssTracksTypes::Mu);
      default:
        return -1;
    }
  }

  // Associated Tracks
  void processEffAssTracks(RecoCollisions const& recoCollisions,
                           McTracks const& recoTracks,
                           aod::McParticles const& genParticles,
                           aod::McCollisions const&)
  {
    for (const auto& coll : recoCollisions) {
      float cent{-1.}, mult{-1.};
      if (!hfcFlowHelper.checkCollision<o2::hf_centrality::CentralityEstimator::None>(coll, cent, mult)) {
        continue;
      }
      if (!coll.has_mcCollision()) {
        continue;
      }
      auto mcCollision = coll.template mcCollision_as<aod::McCollisions>();

      auto tracksInColl = recoTracks.sliceBy(trackIndicesPerColl, coll.globalIndex());
      for (const auto& assTrk : tracksInColl) {
        if (!assTrk.has_mcParticle()) {
          continue;
        }
        auto mcParticle = assTrk.template mcParticle_as<aod::McParticles>();
        int matchedSpecie = getAssocTrkSpecie(mcParticle);
        if (matchedSpecie < 0) {
          continue;
        }

        // DCA selection
        double assTrkPt = assTrk.pt();
        if (usePtDiffDcaXYCut) {
          float const dcaXYTrkCut = funcDcaXYPtCutPrimTrk->Eval(assTrkPt);
          if (std::fabs(assTrk.dcaXY()) > dcaXYTrkCut) {
            continue;
          }
        }
        int assTrkOrigin = mcParticle.isPhysicalPrimary() ? 0 : 1; // 0 for primary, 1 for secondary

        std::vector<double> values = {assTrk.pt(), assTrk.eta(), assTrk.phi(), mcCollision.posZ(), static_cast<double>(matchedSpecie), static_cast<double>(assTrkOrigin)};

        // Fill histograms for associated tracks
        registry.get<THnSparse>(HIST("hAssTrackReco"))->Fill(values.data());

        if (fillRecoTrees) {
          // Downsample tracks based on pt
          if (downSampleTracks < 1.) {
            float const pseudoRndm = assTrk.pt() * 1000. - static_cast<int64_t>(assTrk.pt() * 1000);
            if (assTrk.pt() < ptMaxDownSampleTracks && pseudoRndm >= downSampleTracks) {
              continue;
            }
          }
          rowAssocMc(assTrk.pt(), assTrk.eta(), assTrk.phi(), mcCollision.posZ(), matchedSpecie, assTrkOrigin);
          rowAssocTrkSels(assTrk.tpcNClsCrossedRows(), assTrk.itsClusterMap(), assTrk.itsNCls(), assTrk.dcaXY(), assTrk.dcaZ());
        }
      }
    }

    // Fill histograms for generated associated tracks
    for (const auto& genTrk : genParticles) {
      int matchedSpecie = getAssocTrkSpecie(genTrk);
      if (matchedSpecie < 0) {
        continue;
      }
      if (!genTrk.has_mcCollision()) {
        continue;
      }

      // Pt and eta selection
      if (genTrk.pt() < ptTrkMin || genTrk.pt() > ptTrkMax || std::fabs(genTrk.eta()) > etaTrkMax) {
        continue;
      }

      int genTrkOrigin = genTrk.isPhysicalPrimary() ? 0 : 1; // 0 for primary, 1 for secondary
      auto genCollision = genTrk.template mcCollision_as<aod::McCollisions>();
      registry.fill(HIST("hAssTrackGen"), genTrk.pt(), genTrk.eta(), genTrk.phi(), genCollision.posZ(), matchedSpecie, genTrkOrigin);
    }
  }
  PROCESS_SWITCH(HfTaskCorrelationsEfficiencies, processEffAssTracks, "Process efficiencies for associated tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCorrelationsEfficiencies>(cfgc)};
}
