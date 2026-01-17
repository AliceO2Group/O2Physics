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

/// \file candidateSelectorCd.cxx
/// \brief Cd± → d± K∓ π± selection task
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg Universiity

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/TrackSelectorPID.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH2.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Struct for applying Cd selection cuts
struct HfCandidateSelectorCd {
  Produces<aod::HfSelCd> hfSelCdCandidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> usePid{"usePid", true, "Bool to use or not the PID based on nSigma cut at filtering level"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.1, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 1., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.5, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 2.5, "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 3., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // Combined PID options
  Configurable<bool> usePidTpcAndTof{"usePidTpcAndTof", false, "Bool to decide how to combine TPC and TOF PID: true = both (if present, only one otherwise); false = one is enough"};
  // TPC quality track cuts
  Configurable<int> tpcNClustersFoundMin{"tpcNClustersFoundMin", 0, "min number of found TPC clusters"};
  Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 0, "min number of crossed rows in TPC"};
  Configurable<float> tpcNCrossedRowsOverFindableClustersMin{"tpcNCrossedRowsOverFindableClustersMin", 0., "min ratio crossed rows / findable clusters"};
  Configurable<float> tpcChi2PerClusterMax{"tpcChi2PerClusterMax", 1e10f, "max tpc fit chi2 per TPC cluster"};
  // ITS quality track cuts
  Configurable<int> itsNClustersFoundMin{"itsNClustersFoundMin", 0, "min. number of found ITS clusters"};
  Configurable<float> itsChi2PerClusterMax{"itsChi2PerClusterMax", 1e10f, "max its fit chi2 per ITS cluster"};
  // DCA track cuts
  Configurable<std::vector<double>> binsPtTrack{"binsPtTrack", std::vector<double>{hf_cuts_single_track::vecBinsPtTrack}, "track pT bin limits for DCA XY/Z pT-dependent cut"};
  Configurable<LabeledArray<double>> cutsSingleTrack{"cutsSingleTrack", {hf_cuts_single_track::CutsTrack[0], hf_cuts_single_track::NBinsPtTrack, hf_cuts_single_track::NCutVarsTrack, hf_cuts_single_track::labelsPtTrack, hf_cuts_single_track::labelsCutVarTrack}, "Single-track selections"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_cd_to_de_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_cd_to_de_k_pi::Cuts[0], hf_cuts_cd_to_de_k_pi::NBinsPt, hf_cuts_cd_to_de_k_pi::NCutVars, hf_cuts_cd_to_de_k_pi::labelsPt, hf_cuts_cd_to_de_k_pi::labelsCutVar}, "Cd candidate selection per pT bin"};
  // QA switch
  Configurable<bool> activateQA{"activateQA", false, "Flag to enable QA histogram"};

  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;
  TrackSelectorDe selectorDeuteron;

  const float massCharmDeuteron = 3.23; // possible mass

  using TracksSel = soa::Join<aod::TracksWExtra,
                              aod::TracksPidPi, aod::PidTpcTofFullPi, aod::TracksPidKa, aod::PidTpcTofFullKa, aod::TracksPidDe, aod::PidTpcTofFullDe>;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {

    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
    selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
    selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
    selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
    selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    selectorKaon = selectorPion;
    selectorDeuteron = selectorPion;

    if (activateQA) {
      constexpr int kNBinsSelections = aod::SelectionStep::NSelectionSteps;
      std::string labels[kNBinsSelections];
      labels[0] = "No selection";
      labels[1 + aod::SelectionStep::RecoSkims] = "Skims selection";
      labels[1 + aod::SelectionStep::RecoTopol] = "Skims & Topological selections";
      labels[1 + aod::SelectionStep::RecoPID] = "Skims & Topological & PID selections";
      static const AxisSpec axisSelections = {kNBinsSelections, 0.5, kNBinsSelections + 0.5, ""};
      registry.add("hSelections", "Selections;;#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {axisSelections, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
      for (int iBin = 0; iBin < kNBinsSelections; ++iBin) {
        registry.get<TH2>(HIST("hSelections"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
      }
    }
  }

  /// Single track quality cuts
  /// \param track is track
  /// \return true if track passes all cuts
  template <typename T>
  bool isSelectedCandidateProngQuality(const T& trackPos1, const T& trackNeg, const T& trackPos2)
  {
    if (!isSelectedTrackTpcQuality(trackPos1, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value) ||
        !isSelectedTrackTpcQuality(trackNeg, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value) ||
        !isSelectedTrackTpcQuality(trackPos2, tpcNClustersFoundMin.value, tpcNCrossedRowsMin.value, tpcNCrossedRowsOverFindableClustersMin.value, tpcChi2PerClusterMax.value)) {
      return false;
    }
    if (!isSelectedTrackItsQuality(trackPos1, itsNClustersFoundMin.value, itsChi2PerClusterMax.value) ||
        !isSelectedTrackItsQuality(trackNeg, itsNClustersFoundMin.value, itsChi2PerClusterMax.value) ||
        !isSelectedTrackItsQuality(trackPos2, itsNClustersFoundMin.value, itsChi2PerClusterMax.value)) {
      return false;
    }
    return true;
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <aod::hf_cand::VertexerType ReconstructionType, typename T>
  bool selectionTopol(const T& candidate)
  {
    auto ptCand = candidate.pt();

    int const binPt = findBin(binsPt, ptCand);
    if (binPt == -1) {
      return false;
    }

    // check that the candidate pT is within the analysis range
    if (ptCand < ptCandMin || ptCand >= ptCandMax) {
      return false;
    }

    // cosine of pointing angle
    if (candidate.cpa() <= cuts->get(binPt, "cos pointing angle")) {
      return false;
    }

    // candidate chi2PCA
    if (candidate.chi2PCA() > cuts->get(binPt, "Chi2PCA")) {
      return false;
    }

    if (candidate.decayLength() <= cuts->get(binPt, "decay length")) {
      return false;
    }

    // candidate decay length XY
    if (candidate.decayLengthXY() <= cuts->get(binPt, "decLengthXY")) {
      return false;
    }

    // candidate normalized decay length XY
    if (candidate.decayLengthXYNormalised() < cuts->get(binPt, "normDecLXY")) {
      return false;
    }

    // candidate impact parameter XY
    if (std::abs(candidate.impactParameterXY()) > cuts->get(binPt, "impParXY")) {
      return false;
    }

    if (!isSelectedCandidateProngDca(candidate)) {
      return false;
    }

    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackDeuteron is the track with the deuteron hypothesis
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \return true if candidate passes all cuts for the given Conjugate
  template <aod::hf_cand::VertexerType ReconstructionType, typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackDeuteron, const T2& trackKaon, const T2& trackPion)
  {

    auto ptCand = candidate.pt();
    int const binPt = findBin(binsPt, ptCand);
    if (binPt == -1) {
      return false;
    }

    // cut on daughter pT
    if (trackDeuteron.pt() < cuts->get(binPt, "pT De") || trackKaon.pt() < cuts->get(binPt, "pT K") || trackPion.pt() < cuts->get(binPt, "pT Pi")) {
      return false;
    }

    float massCd{0.f};
    if (trackDeuteron.globalIndex() == candidate.prong0Id()) {
      massCd = HfHelper::invMassCdToDeKPi(candidate);
    } else {
      massCd = HfHelper::invMassCdToPiKDe(candidate);
    }

    // cut on Cd->deKpi, piKde mass values
    if (std::abs(massCd - massCharmDeuteron) > cuts->get(binPt, "m")) {
      return false;
    }

    return true;
  }

  /// Single-track dca_xy and dca_z cuts
  /// \param candidate is the Cd candidate
  /// \return true if all the prongs pass the selections
  template <typename T1>
  bool isSelectedCandidateProngDca(const T1& candidate)
  {
    return (isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng0(), candidate.impactParameter0(), candidate.impactParameterZ0()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng1(), candidate.impactParameter1(), candidate.impactParameterZ1()) &&
            isSelectedTrackDca(binsPtTrack, cutsSingleTrack, candidate.ptProng2(), candidate.impactParameter2(), candidate.impactParameterZ2()));
  }

  /// Apply PID selection
  /// \param pidTrackDeuteron is the PID status of deuteron candidate track
  /// \param pidTrackKaon is the PID status of kaon candidate track
  /// \param pidTrackPion is the PID status of pion candidate track
  /// \return true if prongs pass all selections
  bool isSelectedPID(const TrackSelectorPID::Status pidTrackDeuteron, const TrackSelectorPID::Status pidTrackKaon, const TrackSelectorPID::Status pidTrackPion)
  {
    return pidTrackDeuteron != TrackSelectorPID::Rejected &&
           pidTrackKaon != TrackSelectorPID::Rejected &&
           pidTrackPion != TrackSelectorPID::Rejected;
  }

  /// \brief function to apply Cd selections
  /// \param reconstructionType is the reconstruction type (DCAFitterN )
  /// \param candidates Cd candidate table
  /// \param tracks track table
  template <aod::hf_cand::VertexerType ReconstructionType, typename CandType, typename TTracks>
  void runSelectCd(CandType const& candidates, TTracks const&)
  {
    // looping over 3-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag
      auto statusCdToDeKPi = 0;
      auto statusCdToPiKDe = 0;

      auto ptCand = candidate.pt();

      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::CdToDeKPi)) {
        hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
        if (activateQA) {
          registry.fill(HIST("hSelections"), 1, ptCand);
        }
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoSkims, ptCand);
      }

      auto trackPos1 = candidate.template prong0_as<TTracks>(); // positive daughter (negative for the antiparticles)
      auto trackNeg = candidate.template prong1_as<TTracks>();  // negative daughter (positive for the antiparticles)
      auto trackPos2 = candidate.template prong2_as<TTracks>(); // positive daughter (negative for the antiparticles)

      // implement filter bit 4 cut - should be done before this task at the track selection level

      // track quality selection
      bool const trackQualitySel = isSelectedCandidateProngQuality(trackPos1, trackNeg, trackPos2);
      if (!trackQualitySel) {
        hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
        continue;
      }

      // conjugate-independent topological selection
      if (!selectionTopol<ReconstructionType>(candidate)) {
        hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
        continue;
      }

      // conjugate-dependent topological selection for Cd
      bool const topolCdToDeKPi = selectionTopolConjugate<ReconstructionType>(candidate, trackPos1, trackNeg, trackPos2);
      bool const topolCdToPiKDe = selectionTopolConjugate<ReconstructionType>(candidate, trackPos2, trackNeg, trackPos1);

      if (!topolCdToDeKPi && !topolCdToPiKDe) {
        hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoTopol, candidate.pt());
      }

      // PID not applied, accepted by default
      auto pidCdToDeKPi = 1;
      auto pidCdToPiKDe = 1;

      if (usePid) {
        // track-level PID selection
        TrackSelectorPID::Status pidTrackPos1Deuteron;
        TrackSelectorPID::Status pidTrackPos2Deuteron;
        TrackSelectorPID::Status pidTrackPos1Pion;
        TrackSelectorPID::Status pidTrackPos2Pion;
        TrackSelectorPID::Status pidTrackNegKaon;
        if (usePidTpcAndTof) {
          pidTrackPos1Deuteron = selectorDeuteron.statusTpcAndTof(trackPos1, candidate.nSigTpcDe0(), candidate.nSigTofDe0());
          pidTrackPos2Deuteron = selectorDeuteron.statusTpcAndTof(trackPos2, candidate.nSigTpcDe2(), candidate.nSigTofDe2());
          pidTrackPos1Pion = selectorPion.statusTpcAndTof(trackPos1, candidate.nSigTpcPi0(), candidate.nSigTofPi0());
          pidTrackPos2Pion = selectorPion.statusTpcAndTof(trackPos2, candidate.nSigTpcPi2(), candidate.nSigTofPi2());
          pidTrackNegKaon = selectorKaon.statusTpcAndTof(trackNeg, candidate.nSigTpcKa1(), candidate.nSigTofKa1());
        } else {
          pidTrackPos1Deuteron = selectorDeuteron.statusTpcOrTof(trackPos1, candidate.nSigTpcDe0(), candidate.nSigTofDe0());
          pidTrackPos2Deuteron = selectorDeuteron.statusTpcOrTof(trackPos2, candidate.nSigTpcDe2(), candidate.nSigTofDe2());
          pidTrackPos1Pion = selectorPion.statusTpcOrTof(trackPos1, candidate.nSigTpcPi0(), candidate.nSigTofPi0());
          pidTrackPos2Pion = selectorPion.statusTpcOrTof(trackPos2, candidate.nSigTpcPi2(), candidate.nSigTofPi2());
          pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg, candidate.nSigTpcKa1(), candidate.nSigTofKa1());
        }

        if (!isSelectedPID(pidTrackPos1Deuteron, pidTrackNegKaon, pidTrackPos2Pion)) {
          pidCdToDeKPi = 0; // reject CdToDeKPi
        }
        if (!isSelectedPID(pidTrackPos2Deuteron, pidTrackNegKaon, pidTrackPos1Pion)) {
          pidCdToPiKDe = 0; // accept CdToPiKDe
        }
      }

      if (pidCdToDeKPi == 0 && pidCdToPiKDe == 0) {
        hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
        continue;
      }

      if (activateQA) {
        registry.fill(HIST("hSelections"), 2 + aod::SelectionStep::RecoPID, candidate.pt());
      }

      if (pidCdToDeKPi == 1 && topolCdToDeKPi && trackQualitySel) {
        statusCdToDeKPi = 1; // identified as CdToDeKPi
      }
      if (pidCdToPiKDe == 1 && topolCdToPiKDe && trackQualitySel) {
        statusCdToPiKDe = 1; // identified as CdToPiKDe
      }

      hfSelCdCandidate(statusCdToDeKPi, statusCdToPiKDe);
    }
  }

  /// \brief process function  with DCAFitterN
  /// \param candidates Cd candidate table
  /// \param tracks track table
  void processCdWithDCAFitterN(aod::HfCand3ProngWPidPiKaDe const& candidates,
                               TracksSel const& tracks)
  {
    runSelectCd<aod::hf_cand::VertexerType::DCAFitter>(candidates, tracks);
  }
  PROCESS_SWITCH(HfCandidateSelectorCd, processCdWithDCAFitterN, "Process Cd selection  with DCAFitterN", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCandidateSelectorCd>(cfgc)};
}
