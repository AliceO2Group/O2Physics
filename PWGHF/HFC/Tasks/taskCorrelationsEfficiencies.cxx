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
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/HFC/Utils/utilsFlow.h"

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
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_evsel;
using namespace o2::hfc_utils_flow;

enum axesCharmReco {
  kAxisPtTrig,
  kAxisPvContrib,
  kAxisEtaTrig,
  kAxisPhiTrig,
  kAxisMlScore0,
  kAxisMlScore1
};

enum SparseAxesCharmGen {
  kAxisPtTrigGen,
  kAxisEtaTrigGen,
  kAxisPhiTrigGen
};

enum axesAssocTracks {
  kAxisPtAssoc,
  kAxisEtaAssoc,
  kAxisPhiAssoc,
  kAxisZVtx,
  kAxisAssTrkSpecie,
  kAxisAssTrkOrigin
};

/// Code to select collisions with at least one Ds meson
struct HfTaskCorrelationsEfficiencies {

  Configurable<int> selectionFlag{"selectionFlag", 15, "Selection Flag for hadron (ML score tables are required to run the task)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};
  Configurable<float> yCandMax{"yCandMax", 0.8, "max. cand. rapidity"};
  Configurable<float> ptCandMin{"ptCandMin", 0., "min. cand. pT"};
  Configurable<float> ptCandMax{"ptCandMax", 24., "max. cand. pT"};
  Configurable<int> tpcNClsCrossedRowsMin{"tpcNClsCrossedRowsMin", 70, "min. TPC crossed rows for associated tracks"};
  Configurable<float> etaTrkMax{"etaTrkMax", 1., "max. track eta"};
  Configurable<float> ptTrkMin{"ptTrkMin", 0.3, "min. track pT"};
  Configurable<float> ptTrkMax{"ptTrkMax", 3., "max. track pT"};
  Configurable<float> dcaXYTrkMax{"dcaXYTrkMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrkMax{"dcaZTrkMax", 1., "max. track DCA Z"};
  Configurable<bool> usePtDiffDcaXYCut{"usePtDiffDcaXYCut", false, "Use pt-differential DCAxy cut for associated tracks"};
  Configurable<float> dcaXYTrkNSigmaMax{"dcaXYTrkNSigmaMax", 7, "Cut on number of sigma deviations from expected DCA in the transverse direction"};
  Configurable<std::string> dcaXYPtPrimTrkFunc{"dcaXYPtPrimTrkFunc", "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut"};
  Configurable<std::vector<double>> bdtScore0PtMaxs{"bdtScore0PtMaxs", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential maximum score 0 for charm candidates"};
  Configurable<std::vector<double>> bdtScore1PtMins{"bdtScore1PtMins", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential minimum score 1 for charm candidates"};

  // Configurables for efficiency validation
  Configurable<bool> runEfficiencyValidation{"runEfficiencyValidation", false, "Flag for running efficiency validation"};
  Configurable<std::vector<double>> binsPtEffCharm{"binsPtEffCharm", std::vector<double>{0., 2., 4., 6., 8., 10., 24.}, "pT bin limits for efficiency"};
  Configurable<std::vector<double>> binsPtEffAssTrk{"binsPtEffAssTrk", std::vector<double>{0.3, 2., 4., 8., 12., 50.}, "pT bin limits for associated particle efficiency"};
  Configurable<std::vector<float>> effPromptCharm{"effPromptCharm", {1., 1., 1., 1., 1., 1.}, "efficiency values for prompt D+ meson"};
  Configurable<std::vector<float>> effFdCharm{"effFdCharm", {1., 1., 1., 1., 1., 1.}, "efficiency values for beauty feed-down D+ meson"};
  Configurable<std::vector<float>> effAssTrk{"effAssTrk", {1., 1., 1., 1., 1.}, "efficiency values for associated particles"};
  Configurable<bool> loadAccXEffFromCCDB{"loadAccXEffFromCCDB", false, "Flag for loading efficiency distributions from CCDB"};
  Configurable<std::string> assTrkEffPath{"assTrkEffPath", "", "CCDB path for ass. tracks efficiency"};
  Configurable<std::string> promptEffPath{"promptEffPath", "", "CCDB path for trigger prompt efficiency"};
  Configurable<std::string> fdEffPath{"fdEffPath", "", "CCDB path for trigger fd efficiency"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "timestamp of the efficiency files used to query in CCDB"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HfcFlowHelper hfcFlowHelper; // flow utilities
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb{};
  std::shared_ptr<THnSparse> mRecoPrompt = nullptr;
  std::shared_ptr<THnSparse> mGenPrompt = nullptr;
  std::shared_ptr<THnSparse> mRecoFD = nullptr;
  std::shared_ptr<THnSparse> mGenFD = nullptr;
  std::shared_ptr<THnSparse> mRecoAssTrk = nullptr;
  std::shared_ptr<THnSparse> mGenAssTrk = nullptr;
  SliceCache cache;

  TF1* funcDcaXYPtCutPrimTrk = nullptr;

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

  ConfigurableAxis binsInvMass{"binsInvMass", {600, 1.6, 2.2}, ""};
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
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);

    // Fetch CCDB correction histograms for efficiency validation
    if (runEfficiencyValidation && loadAccXEffFromCCDB) {

      mRecoPrompt = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/PromptReco.root", promptEffPath.value.c_str()), timestampCcdb));
      if (mRecoPrompt == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffPath.value.c_str());
      }
      mGenPrompt = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/PromptGen.root", promptEffPath.value.c_str()), timestampCcdb));
      if (mGenPrompt == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", promptEffPath.value.c_str());
      }
      LOGF(info, "Loaded trigger efficiency (prompt D) histogram from %s", promptEffPath.value.c_str());

      mRecoFD = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/FDReco.root", fdEffPath.value.c_str()), timestampCcdb));
      if (mRecoFD == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", fdEffPath.value.c_str());
      }
      mGenFD = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/FDGen.root", fdEffPath.value.c_str()), timestampCcdb));
      if (mGenFD == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", fdEffPath.value.c_str());
      }
      LOGF(info, "Loaded feed-down D meson efficiency histogram from %s", fdEffPath.value.c_str());

      mRecoAssTrk = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/TracksReco.root", assTrkEffPath.value.c_str()), timestampCcdb));
      if (mRecoAssTrk == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", assTrkEffPath.value.c_str());
      }
      mGenAssTrk = std::shared_ptr<THnSparse>(ccdb->getForTimeStamp<THnSparse>(Form("%s/TracksGen.root", assTrkEffPath.value.c_str()), timestampCcdb));
      if (mGenAssTrk == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", assTrkEffPath.value.c_str());
      }
      LOGF(info, "Loaded associated efficiency histogram from %s", assTrkEffPath.value.c_str());
    }

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
      funcDcaXYPtCutPrimTrk = new TF1("funcDcaXYPtCutPrimTrk", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()), 0.001, 100);
      funcDcaXYPtCutPrimTrk->SetParameter(0, dcaXYTrkNSigmaMax);
      LOGF(info, "DCAxy pt-dependence function: %s", Form("[0]*%s", dcaXYPtPrimTrkFunc.value.data()));
    }

    if (((binsPtCharm.value.size() - 1) != (bdtScore0PtMaxs.value.size() + 1) ||
         (binsPtCharm.value.size() - 1) != (bdtScore1PtMins.value.size() + 1))) {
      LOG(info) << "binsPtCharm.value.size(): " << binsPtCharm.value.size();
      LOG(info) << "bdtScore0PtMaxs.value.size(): " << bdtScore0PtMaxs.value.size();
      LOG(info) << "bdtScore1PtMins.value.size(): " << bdtScore1PtMins.value.size();
      LOGP(fatal, "The size of bdtScore0PtMaxs and bdtScore1PtMins must be the one of binsPtCharm minus one!");
    }

    // Add efficiency sparses
    // Axes: mass, pt, pv contributors, pseudo-rapidity, ml_score_0, ml_score_1
    const AxisSpec axisInvMass = {binsInvMass, "Mass (GeV/c^{2})"};
    const AxisSpec axisPtTrig = {binsPtCharm, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPvContrib = {binsPvContrib, "PV contributors"};
    const AxisSpec axisEtaTrig = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhiTrig = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisMlScore0 = {binsMlOne, "ML score 0"};
    const AxisSpec axisMlScore1 = {binsMlTwo, "ML score 1"};
    std::vector<AxisSpec> axesCharmReco = {axisPtTrig, axisPvContrib, axisEtaTrig, axisPhiTrig, axisMlScore0, axisMlScore1};

    registry.add("hCharmPromptGen", "Generated prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, {{axisPtTrig}, {axisEtaTrig}, {axisPhiTrig}});
    registry.add("hCharmFdGen", "Generated non-prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, {{axisPtTrig}, {axisEtaTrig}, {axisPhiTrig}});
    registry.add("hCharmPromptReco", "Reconstructed prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, axesCharmReco);
    registry.add("hCharmFdReco", "Reconstructed non-prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, axesCharmReco);
    if (runEfficiencyValidation) {
      registry.add("hCharmPromptRecoPtRew", "Reconstructed prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, axesCharmReco);
      registry.add("hCharmFdRecoPtRew", "Reconstructed non-prompt charm-hadron candidates;;", {HistType::kTHnSparseF}, axesCharmReco);
      registry.get<THnSparse>(HIST("hCharmPromptRecoPtRew"))->Sumw2();
      registry.get<THnSparse>(HIST("hCharmFdRecoPtRew"))->Sumw2();
    }

    // Axes: pt, eta, phi, z_vtx, specie
    const AxisSpec axisPtAssoc = {binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisEtaAssoc = {binsEta, "#it{#eta}"};
    const AxisSpec axisPhiAssoc = {binsPhi, "#it{#varphi}"};
    const AxisSpec axisZVtx = {binsZVtx, "#it{z}_{vtx}"};
    const AxisSpec axisAssTrkSpecie = {static_cast<int>(NAssTracksTypes),
                                       -0.5, static_cast<int>(NAssTracksTypes) - 0.5,
                                       "Particle species"};
    const AxisSpec axisTrackOrigin = {2, -0.5, 1.5, "Particle origin"};
    std::vector<AxisSpec> axesAssocTracks = {axisPtAssoc, axisEtaAssoc, axisPhiAssoc, axisZVtx, axisAssTrkSpecie, axisTrackOrigin};
    registry.add("hAssTrackReco", "Reconstructed associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);
    if (runEfficiencyValidation) {
      registry.add("hAssTrackRecoPtRew", "Reconstructed associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);
      registry.get<THnSparse>(HIST("hAssTrackRecoPtRew"))->Sumw2();
      if (loadAccXEffFromCCDB) {
        registry.add("hAssTrackRecoPtEtaRew", "Reconstructed associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);
        registry.add("hAssTrackRecoPtEtaPhiRew", "Reconstructed associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);
        registry.add("hAssTrackRecoPtEtaPhiZVtxRew", "Reconstructed associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);
        registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaRew"))->Sumw2();
        registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaPhiRew"))->Sumw2();
        registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaPhiZVtxRew"))->Sumw2();
      }
    }
    registry.add("hAssTrackGen", "Generated associated tracks;;", {HistType::kTHnSparseF}, axesAssocTracks);

  } // end init

  /// Apply pT-differential ML BDT bkg score cut
  /// \param ptTrig is the pT of the charm candidate
  template <TriggerChannels CandType, typename TCand>
  bool isSelBdtScoreCut(TCand const& cand,
                        double ptTrig)
  {
    std::array<float, 2> outputMl = hfcFlowHelper.getCandMlScores<CandType>(cand);
    for (size_t iPt = 0; iPt < binsPtCharm.value.size() - 1; iPt++) {
      if (ptTrig >= binsPtCharm.value[iPt] && ptTrig < binsPtCharm.value[iPt + 1]) {
        return (outputMl[0] < bdtScore0PtMaxs.value[iPt]) && (outputMl[1] > bdtScore1PtMins.value[iPt]);
      }
    }
    return false;
  }

  /// Match associated tracks based on their PDG code
  /// \param assTrk is the associated track
  /// \return the associated track specie index
  template <bool IsReco, TriggerChannels CandType, typename TCharmHad>
  bool isMcMatched(TCharmHad const& cand)
  {
    int flagMc{-1}, resoChannel{-1};
    if constexpr (IsReco) {
      flagMc = std::abs(cand.flagMcMatchRec());
      resoChannel = std::abs(cand.flagMcDecayChanRec());
    } else {
      flagMc = std::abs(cand.flagMcMatchGen());
      resoChannel = std::abs(cand.flagMcDecayChanGen());
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

  /// Fill histograms and tables for same-event correlations
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

      int channel = isMcMatched<true, CandType>(cand);
      if (channel < 0) {
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
    }
  }

  template <typename TCharmCand>
  double getCharmEffWeight(const TCharmCand& charmCand, bool isPrompt)
  {
    double eff = 1.;
    double pt = charmCand.pt();
    if (loadAccXEffFromCCDB) {

      int ptBin = o2::analysis::findBin(binsPtCharm, pt);
      double mlScore0Max = bdtScore0PtMaxs.value[ptBin];
      double mlScore1Min = bdtScore1PtMins.value[ptBin];

      if (isPrompt) {

        int const sparsePtBin = mRecoPrompt->GetAxis(kAxisPtTrig)->FindBin(pt);
        mRecoPrompt->GetAxis(kAxisPtTrig)->SetRange(sparsePtBin, sparsePtBin);
        mRecoPrompt->GetAxis(kAxisMlScore0)->SetRangeUser(0., mlScore0Max);
        mRecoPrompt->GetAxis(kAxisMlScore1)->SetRangeUser(mlScore1Min, 1.);
        mGenPrompt->GetAxis(kAxisPtTrig)->SetRange(sparsePtBin, sparsePtBin);
        eff = mRecoPrompt->Integral(kTRUE) / mGenPrompt->Integral(kTRUE);

        // Restore axes to full range
        for (int i = 0; i < mRecoPrompt->GetNdimensions(); ++i) mRecoPrompt->GetAxis(i)->SetRange(0, 0);
        for (int i = 0; i < mGenPrompt->GetNdimensions(); ++i) mGenPrompt->GetAxis(i)->SetRange(0, 0);

      } else {

        int const sparsePtBin = mRecoFD->GetAxis(kAxisPtTrig)->FindBin(pt);
        mRecoFD->GetAxis(kAxisPtTrig)->SetRange(sparsePtBin, sparsePtBin);
        mRecoFD->GetAxis(kAxisMlScore0)->SetRangeUser(0., mlScore0Max);
        mRecoFD->GetAxis(kAxisMlScore1)->SetRangeUser(mlScore1Min, 1.);
        mGenFD->GetAxis(kAxisPtTrig)->SetRange(sparsePtBin, sparsePtBin);
        eff = mRecoFD->Integral(kTRUE) / mGenFD->Integral(kTRUE);

        // Restore axes to full range
        for (int i = 0; i < mRecoFD->GetNdimensions(); ++i) mRecoFD->GetAxis(i)->SetRange(0, 0);
        for (int i = 0; i < mGenFD->GetNdimensions(); ++i) mGenFD->GetAxis(i)->SetRange(0, 0);

      }

    } else {
      int const bin = o2::analysis::findBin(binsPtEffCharm, pt);
      eff = isPrompt ? effPromptCharm->at(bin) : effFdCharm->at(bin);
    }

    return 1 / eff;
  }

  /// Fill histograms and tables for same-event correlations
  /// \param trigCands are the trigger candidates
  /// \param assTrks are the associated tracks
  /// \param cent is the collision centrality
  /// \param mult is the collision multiplicity
  template <TriggerChannels CandType, typename TCharmCands>
  void validateCharmEfficiencies(TCharmCands const& charmCands)
  {
    for (const auto& charmCand : charmCands) {
      if (!hfcFlowHelper.cutCandRapidity<CandType>(charmCand)) {
        continue;
      }

      int channel = isMcMatched<true, CandType>(charmCand);
      if (channel < 0) {
        continue;
      }
      if (!isSelBdtScoreCut<CandType>(charmCand, charmCand.pt())) {
        continue;
      }

      // Fill the sparse
      std::array<float, 2> outputMl = hfcFlowHelper.getCandMlScores<CandType>(charmCand);
      bool isPrompt = (charmCand.originMcRec() == RecoDecay::OriginType::Prompt);
      double effWeight = getCharmEffWeight(charmCand, isPrompt);
      if (isPrompt) {
        registry.fill(HIST("hCharmPromptRecoPtRew"), charmCand.pt(), charmCand.template collision_as<RecoCollisions>().numContrib(),
                                                     charmCand.eta(), charmCand.phi(), outputMl[0], outputMl[1], effWeight);
      } else {
        registry.fill(HIST("hCharmFdRecoPtRew"), charmCand.pt(), charmCand.template collision_as<RecoCollisions>().numContrib(),
                                                 charmCand.eta(), charmCand.phi(), outputMl[0], outputMl[1], effWeight);
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
      if (runEfficiencyValidation) {
        validateCharmEfficiencies<TriggerChannels::DplusToPiKPi>(candsInColl);
      }
    }

    for (const auto& genCand : genParticles) {
      if (!genCand.has_mcCollision()) {
        continue;
      }
      if (std::abs(genCand.y()) > yCandMax) {
        continue;
      }
      int channel = isMcMatched<false, TriggerChannels::DplusToPiKPi>(genCand);
      if (channel < 0) {
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
  PROCESS_SWITCH(HfTaskCorrelationsEfficiencies, processEffDplus, "Process efficiencies for Dplus candidates", true);

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

        std::vector<double> values = {mcParticle.pt(), mcParticle.eta(), mcParticle.phi(), mcCollision.posZ(), static_cast<double>(matchedSpecie), static_cast<double>(assTrkOrigin)};

        // Fill histograms for associated tracks
        registry.get<THnSparse>(HIST("hAssTrackReco"))->Fill(values.data());

        // Fill histograms for efficiency validation if enabled
        if (runEfficiencyValidation) {

          float effWeightPt{1.f};
          if (!loadAccXEffFromCCDB) {
            int const bin = o2::analysis::findBin(binsPtEffAssTrk, assTrkPt);
            effWeightPt = 1. / effAssTrk->at(bin);
          } else {

            float effWeightPtEta{1.f}, effWeightPtEtaPhi{1.f}, effWeightPtEtaPhiZVtx{1.f};

            // Pt-reweighting
            mRecoAssTrk->GetAxis(kAxisPtAssoc)->SetRangeUser(assTrkPt, assTrkPt);
            mGenAssTrk->GetAxis(kAxisPtAssoc)->SetRangeUser(assTrkPt, assTrkPt);
            effWeightPt = mRecoAssTrk->Integral(kTRUE) > 0 ? mGenAssTrk->Integral(kTRUE) / mRecoAssTrk->Integral(kTRUE) : 0;

            // Pt-eta-reweighting
            mRecoAssTrk->GetAxis(kAxisEtaAssoc)->SetRangeUser(mcParticle.eta(), mcParticle.eta());
            mGenAssTrk->GetAxis(kAxisEtaAssoc)->SetRangeUser(mcParticle.eta(), mcParticle.eta());
            effWeightPtEta = mRecoAssTrk->Integral(kTRUE) > 0 ? mGenAssTrk->Integral(kTRUE) / mRecoAssTrk->Integral(kTRUE) : 0;

            // Pt-eta-phi-reweighting
            mRecoAssTrk->GetAxis(kAxisPhiAssoc)->SetRangeUser(mcParticle.phi(), mcParticle.phi());
            mGenAssTrk->GetAxis(kAxisPhiAssoc)->SetRangeUser(mcParticle.phi(), mcParticle.phi());
            effWeightPtEtaPhi = mRecoAssTrk->Integral(kTRUE) > 0 ? mGenAssTrk->Integral(kTRUE) / mRecoAssTrk->Integral(kTRUE) : 0;

            // Pt-eta-phi-zvtx-reweighting
            mRecoAssTrk->GetAxis(kAxisZVtx)->SetRangeUser(mcCollision.posZ(), mcCollision.posZ());
            mGenAssTrk->GetAxis(kAxisZVtx)->SetRangeUser(mcCollision.posZ(), mcCollision.posZ());
            effWeightPtEtaPhiZVtx = mRecoAssTrk->Integral(kTRUE) > 0 ? mGenAssTrk->Integral(kTRUE) / mRecoAssTrk->Integral(kTRUE) : 0;

            // Restore axes to full range
            for (int i = 0; i < mRecoAssTrk->GetNdimensions(); ++i) mRecoAssTrk->GetAxis(i)->SetRange(0, 0);
            for (int i = 0; i < mGenAssTrk->GetNdimensions(); ++i) mGenAssTrk->GetAxis(i)->SetRange(0, 0);

            registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaRew"))->Fill(values.data(), effWeightPtEta);
            registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaPhiRew"))->Fill(values.data(), effWeightPtEtaPhi);
            registry.get<THnSparse>(HIST("hAssTrackRecoPtEtaPhiZVtxRew"))->Fill(values.data(), effWeightPtEtaPhiZVtx);
          }

          registry.get<THnSparse>(HIST("hAssTrackRecoPtRew"))->Fill(values.data(), effWeightPt);

        }
      }
    }

    // Fill histograms for generated associated tracks
    for (const auto& genTrk : genParticles) {
      int matchedSpecie = getAssocTrkSpecie(genTrk);
      if (matchedSpecie < 0) {
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
