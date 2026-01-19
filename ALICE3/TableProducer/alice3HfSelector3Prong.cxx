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

/// \file alice3HfSelector3prong.cxx
/// \brief 3-prong candidates selection task
///
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN Turin

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/ML/HfMlResponse3Prong.h"
#include "ALICE3/Utils/utilsHfAlice3.h"
#include "ALICE3/Utils/utilsSelectionsAlice3.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <array>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::aod::a3_hf_sel_3prong;

/// Struct for applying Lc selection cuts
struct Alice3HfSelector3Prong {
  Produces<aod::Alice3Sel3Ps> candSelFlags; // flags for isSelLc
  Produces<aod::Alice3Ml3Ps> candMlScores;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of cand pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of cand pT"};
  // TRK PID
  Configurable<double> ptPidTrkMin{"ptPidTrkMin", 0.1, "Lower bound of track pT for Trk PID"};
  Configurable<double> ptPidTrkMax{"ptPidTrkMax", 1., "Upper bound of track pT for Trk PID"};
  Configurable<double> nSigmaTrkMax{"nSigmaTrkMax", 3., "Nsigma cut on Trk only"};
  // RICH PID
  Configurable<double> ptPidRichMin{"ptPidRichMin", 0.1, "Lower bound of track pT for Rich PID"};
  Configurable<double> ptPidRichMax{"ptPidRichMax", 1., "Upper bound of track pT for Rich PID"};
  Configurable<double> nSigmaRichMax{"nSigmaRichMax", 3., "Nsigma cut on Rich only"};
  // INNER TOF PID
  Configurable<double> ptPidInnTofMin{"ptPidInnTofMin", 0.5, "Lower bound of track pT for InTOF PID"};
  Configurable<double> ptPidInnTofMax{"ptPidInnTofMax", 2.5, "Upper bound of track pT for InTOF PID"};
  Configurable<double> nSigmaInnTofMax{"nSigmaInnTofMax", 3., "Nsigma cut on InTOF only"};
  // OUTER TOF PID
  Configurable<double> ptPidOutTofMin{"ptPidOutTofMin", 0.5, "Lower bound of track pT for OutTOF PID"};
  Configurable<double> ptPidOutTofMax{"ptPidOutTofMax", 2.5, "Upper bound of track pT for OutTOF PID"};
  Configurable<double> nSigmaOutTofMax{"nSigmaOutTofMax", 3., "Nsigma cut on OutTOF only"};
  // DCA track cuts
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCA XY/Z pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsSingleTrack{"cutsSingleTrack", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_3prongs_alice3::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_3prongs_alice3::Cuts[0], hf_cuts_3prongs_alice3::NBinsPt, hf_cuts_3prongs_alice3::NCutVars, hf_cuts_3prongs_alice3::labelsPt, hf_cuts_3prongs_alice3::labelsCutVar}, "Lc cand selection per pT bin"};
  // QA switch
  // ML inference
  Configurable<bool> applyMl{"applyMl", false, "Flag to apply ML selections"};
  Configurable<std::vector<double>> binsPtMl{"binsPtMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirMl{"cutDirMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsMl{"cutsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesMl{"nClassesMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"EventFiltering/PWGHF/BDTLc"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_MassHypo0.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelperAlice3 hfHelper;
  o2::analysis::HfMlResponse3Prong<float> mlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using CandsLc = soa::Join<aod::Alice3Cand3Ps, aod::Alice3PidLcs>;
  using CandsLcWMcTruth = soa::Join<aod::Alice3Cand3Ps, aod::Alice3PidLcs, aod::Alice3McRecFlags>;

  float MassReference{-1.f};

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    constexpr int kNBinsSelections = 1 + aod::SelectionStep::NSelectionSteps;
    std::string labels[kNBinsSelections];
    labels[0] = "No selection";
    labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
    labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
    labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
    labels[1 + aod::SelectionStep::RecoMl] = "ML selection";
    static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
    registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
      registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    auto h = registry.add<TH1>("hSelectionsTopology", "hSelectionsTopology", {HistType::kTH1D, {{11, -0.5, 10.5, "Selection step"}}});
    h->GetXaxis()->SetBinLabel(1, "All candidates");
    h->GetXaxis()->SetBinLabel(2, "pT cand");
    h->GetXaxis()->SetBinLabel(3, "pT prong cuts");
    h->GetXaxis()->SetBinLabel(4, "cos pointing angle");
    h->GetXaxis()->SetBinLabel(5, "chi2PCA");
    h->GetXaxis()->SetBinLabel(6, "decay length");
    h->GetXaxis()->SetBinLabel(7, "decay length XY");
    h->GetXaxis()->SetBinLabel(8, "norm decay length XY");
    h->GetXaxis()->SetBinLabel(9, "impPar XY");
    h->GetXaxis()->SetBinLabel(10, "prong DCA");
    h->GetXaxis()->SetBinLabel(11, "finally accepted");

    registry.add("Tried/hChi2PCA", "Chi2PCA;Chi2PCA;entries", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("Tried/hDecayLength", "Decay Length;Decay Length;entries", {HistType::kTH1F, {{100, 0., 200.}}});
    registry.addClone("Tried/", "Accepted/");

    if (applyMl) {
      mlResponse.configure(binsPtMl, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        mlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        mlResponse.setModelPathsLocal(onnxFileNames);
      }
      mlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      mlResponse.init();
    }

    if (doprocessLc) {
      MassReference = o2::constants::physics::MassLambdaCPlus;
    }
  }

  /// Conjugate-independent topological cuts
  /// \tparam T is candidate type
  /// \param cand is cand
  /// \param candPt is candidate pT
  /// \return true if cand passes all cuts
  template <typename T>
  bool selectionTopol(const T& cand, float candPt)
  {
    int const ptBin = findBin(binsPt, candPt);
    registry.fill(HIST("hSelectionsTopology"), 0.f);
    // check that the cand pT is within the analysis range
    if (candPt < ptCandMin || candPt >= ptCandMax) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 1.f);

    // cut on daughter pT
    if (cand.ptProng0() < cuts->get(ptBin, "pT prong 0") ||
        cand.ptProng1() < cuts->get(ptBin, "pT prong 1") ||
        cand.ptProng2() < cuts->get(ptBin, "pT prong 2")) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 2.f);

    // cosine of pointing angle
    if (cand.cpa() <= cuts->get(ptBin, "cos pointing angle")) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 3.f);

    // cand chi2PCA
    registry.fill(HIST("Tried/hChi2PCA"), cand.chi2PCA());
    if (cand.chi2PCA() > cuts->get(ptBin, "Chi2PCA")) {
      return false;
    }
    registry.fill(HIST("Accepted/hChi2PCA"), cand.chi2PCA());
    registry.fill(HIST("hSelectionsTopology"), 4.f);

    // cand decay length
    registry.fill(HIST("Tried/hDecayLength"), cand.decayLength());
    if (cand.decayLength() <= cuts->get(ptBin, "decay length")) {
      return false;
    }
    registry.fill(HIST("Accepted/hDecayLength"), cand.decayLength());
    registry.fill(HIST("hSelectionsTopology"), 5.f);

    // cand decay length XY
    if (cand.decayLengthXY() <= cuts->get(ptBin, "decLengthXY")) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 6.f);

    // cand normalized decay length XY
    if (cand.decayLengthXYNormalised() < cuts->get(ptBin, "normDecLXY")) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 7.f);

    // cand impact parameter XY
    if (std::abs(cand.impactParameterXY()) > cuts->get(ptBin, "impParXY")) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 8.f);

    // cand daughter prong DCA
    if (!isSelectedCandidateProngDca(cand)) {
      return false;
    }
    registry.fill(HIST("hSelectionsTopology"), 9.f);

    registry.fill(HIST("hSelectionsTopology"), 10.f);
    return true;
  }

  /// Candidate mass selection
  /// \tparam CharmHad is the charm hadron type
  /// \tparam SwapHypo indicates whether to swap mass hypothesis or not
  /// \tparam TCandidate is candidate type
  /// \param ptBin is candidate pT bin
  /// \param cand is candidate
  /// \return true if candidate passes mass selection
  template <CharmHadAlice3 CharmHad, bool SwapHypo, typename TCandidate>
  bool selectionCandidateMass(int const ptBin, const TCandidate& cand)
  {
    const float massCand = hfHelper.getCandMass<CharmHad, SwapHypo>(cand);
    // cut on mass window
    if (std::abs(massCand - MassReference) > cuts->get(ptBin, "m")) {
      return false;
    }

    return true;
  }

  /// Single-track dca_xy and dca_z cuts
  /// \tparam T1 is candidate type
  /// \param cand is the Lc cand
  /// \return true if all the prongs pass the selections
  template <typename T1>
  bool isSelectedCandidateProngDca(const T1& cand)
  {
    return (isSelectedTrackDca(binsPtTrack, cutsSingleTrack, cand.ptProng0(), cand.impactParameterY0(), cand.impactParameterZ0()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, cand.ptProng1(), cand.impactParameterY1(), cand.impactParameterZ1()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, cand.ptProng2(), cand.impactParameterY2(), cand.impactParameterZ2()));
  }

  /// Apply PID selection
  /// \tparam CharmHad is charm hadron type
  /// \tparam TCand is candidate type
  /// \param cand is candidate
  /// \param pidMask is bitmask to be configured
  template <CharmHadAlice3 CharmHad, typename TCand>
  void configurePidMask(const TCand& cand, uint32_t& pidMask)
  {

    auto isSelPid = [&](int selCut, float nsigma, float pt, float nSigmaMax, float ptMin, float ptMax) {
      bool isSelected = !(pt >= ptMin && pt < ptMax && std::abs(nsigma) > nSigmaMax);
      if (isSelected)
        SETBIT(pidMask, selCut);
      return isSelected;
    };

    // prong 0
    float ptProng0{cand.ptProng0()};
    if constexpr (CharmHad == CharmHadAlice3::Lc) {
      isSelPid(cand.nSigTrkPr0(), PidSels::TrkProng0, ptProng0, nSigmaTrkMax, ptPidTrkMin, ptPidTrkMax);
      isSelPid(cand.nSigRichPr0(), PidSels::RichProng0, ptProng0, nSigmaRichMax, ptPidRichMin, ptPidRichMax);
      isSelPid(cand.nSigInnTofPr0(), PidSels::InnTofProng0, ptProng0, nSigmaInnTofMax, ptPidInnTofMin, ptPidInnTofMax);
      isSelPid(cand.nSigOutTofPr0(), PidSels::OutTofProng0, ptProng0, nSigmaOutTofMax, ptPidOutTofMin, ptPidOutTofMax);
    }

    // prong 1
    float ptProng1{cand.ptProng1()};
    if constexpr (CharmHad == CharmHadAlice3::Lc) {
      isSelPid(cand.nSigTrkKa1(), PidSels::TrkProng1, ptProng1, nSigmaTrkMax, ptPidTrkMin, ptPidTrkMax);
      isSelPid(cand.nSigRichKa1(), PidSels::RichProng1, ptProng1, nSigmaRichMax, ptPidRichMin, ptPidRichMax);
      isSelPid(cand.nSigInnTofKa1(), PidSels::InnTofProng1, ptProng1, nSigmaInnTofMax, ptPidInnTofMin, ptPidInnTofMax);
      isSelPid(cand.nSigOutTofKa1(), PidSels::OutTofProng1, ptProng1, nSigmaOutTofMax, ptPidOutTofMin, ptPidOutTofMax);
    }

    // prong 2
    float ptProng2{cand.ptProng2()};
    if constexpr (CharmHad == CharmHadAlice3::Lc) {
      isSelPid(cand.nSigTrkPi2(), PidSels::TrkProng2, ptProng2, nSigmaTrkMax, ptPidTrkMin, ptPidTrkMax);
      isSelPid(cand.nSigRichPi2(), PidSels::RichProng2, ptProng2, nSigmaRichMax, ptPidRichMin, ptPidRichMax);
      isSelPid(cand.nSigInnTofPi2(), PidSels::InnTofProng2, ptProng2, nSigmaInnTofMax, ptPidInnTofMin, ptPidInnTofMax);
      isSelPid(cand.nSigOutTofPi2(), PidSels::OutTofProng2, ptProng2, nSigmaOutTofMax, ptPidOutTofMin, ptPidOutTofMax);
    }

    return;
  }

  /// \brief function to apply Lc selections
  /// \tparam CharmHad is charm hadron type
  /// \tparam CandType is candidate type
  /// \param cands are 3-prong candidates
  template <CharmHadAlice3 CharmHad, typename CandType>
  void runSelect3Prong(CandType const& cands)
  {
    std::vector<float> outputMl{-1.f, -1.f, -1.f};
    uint32_t pidMask = 0;

    // looping over 3-prong cands
    for (const auto& cand : cands) {
      registry.fill(HIST("hSelections"), 1, cand.pt());
      outputMl = {-1.f, -1.f, -1.f};
      pidMask = 0;

      const float ptCand = cand.pt();
      const int ptBin = findBin(binsPt, ptCand);
      if (ptBin == -1) {
        candSelFlags(false, false, pidMask);
        if (applyMl) {
          candMlScores(outputMl[0], outputMl[1], outputMl[2]);
        }
        continue;
      }

      // Here all cands pass the cut on the mass selection
      const bool selMassHypo0 = selectionCandidateMass<CharmHad, false>(ptBin, cand);
      const bool selMassHypo1 = selectionCandidateMass<CharmHad, true>(ptBin, cand);
      if (!selMassHypo0 && !selMassHypo1) {
        candSelFlags(false, false, pidMask);
        if (applyMl) {
          candMlScores(outputMl[0], outputMl[1], outputMl[2]);
        }
        continue;
      }
      registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);

      // Topological selection (TODO: track quality selection)
      if (!selectionTopol(cand, ptCand)) {
        candSelFlags(false, false, pidMask);
        if (applyMl) {
          candMlScores(outputMl[0], outputMl[1], outputMl[2]);
        }
        continue;
      }
      registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, ptCand);

      // PID selection
      configurePidMask<CharmHad>(cand, pidMask);
      if (pidMask == 0) {
        candSelFlags(false, false, pidMask);
        if (applyMl) {
          candMlScores(outputMl[0], outputMl[1], outputMl[2]);
        }
        continue;
      }
      registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, ptCand);

      bool isSelectedMl = true;
      // ML selections
      if (applyMl) {

        std::vector<float> inputFeaturesMassHypo0 = mlResponse.getInputFeatures(cand);
        isSelectedMl = mlResponse.isSelectedMl(inputFeaturesMassHypo0, ptCand, outputMl);
        candMlScores(outputMl[0], outputMl[1], outputMl[2]);
        if (!isSelectedMl) {
          candSelFlags(false, false, pidMask);
          continue;
        }

        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoMl, ptCand);
      }

      candSelFlags(selMassHypo0, selMassHypo1, pidMask);
    }
  }

  /// \brief process function for cand selection
  /// \param cands Lc cand table
  void processLc(CandsLcWMcTruth const& cands)
  {
    runSelect3Prong<CharmHadAlice3::Lc>(cands);
  }
  PROCESS_SWITCH(Alice3HfSelector3Prong, processLc, "Process 3 prong selection for Lc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3HfSelector3Prong>(cfgc)};
}
