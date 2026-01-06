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

/// \file tableMakerJpsiHf.cxx
/// \brief Task for the production of the derived data of J/psi - open HF associate production
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::aod::hf_cand_2prong;

// bit maps used for the Fill functions of the VarManager
constexpr static uint32_t gkTrackFillMapWithColl = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID | VarManager::ObjTypes::ReducedTrackCollInfo;
constexpr static uint32_t gkMuonFillMapWithColl = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCollInfo;

namespace
{
struct CandidateDilepton {
  float fMass;
  float fPt;
  float fPhi;
  float fEta;
  float fTauz;
  float fLz;
  float fLxy;
  int fSign;
  uint32_t fMcDecision;

  float pt() const { return fPt; }
  float mass() const { return fMass; }
  float phi() const { return fPhi; }
  float eta() const { return fEta; }
  int sign() const { return fSign; }
  uint32_t mcDecision() const { return fMcDecision; }
  float tauz() const { return fTauz; }
  float lz() const { return fLz; }
  float lxy() const { return fLxy; }
  float px() const { return fPt * std::cos(fPhi); }
  float py() const { return fPt * std::sin(fPhi); }
  float pz() const { return fPt * std::sinh(fEta); }
  float p() const { return fPt * std::cosh(fEta); }
  float rap() const { return std::log((std::sqrt(fMass * fMass + fPt * fPt * std::cosh(fEta) * std::cosh(fEta)) + fPt * std::sinh(fEta)) / std::sqrt(fMass * fMass + fPt * fPt)); }
};
} // namespace

// Declarations of various short names
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyD0CandidatesSelected = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0>;
using MyD0CandidatesSelectedWithBdt = soa::Join<aod::HfCand2ProngWPid, aod::HfSelD0, aod::HfMlD0>;
using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

using MyBarrelTracksSelectedWithColl = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::ReducedTracksBarrelInfo>;
using MyMuonTracksSelectedWithColl = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsInfo>;

HfHelper hfHelper;

struct tableMakerJpsiHf {
  //
  // This task combines dilepton candidates with a open charm hadron
  //

  // Produce derived tables
  Produces<RedJpDmColls> redCollisions;
  Produces<RedJpDmDmesons> redDmesons;
  Produces<RedJpDmDmDau0s> redDmesDau0;
  Produces<RedJpDmDmDau1s> redDmesDau1;
  Produces<RedJpDmDmesBdts> redDmesBdts;
  Produces<RedJpDmD0Masss> redD0Masses;
  Produces<RedJpDmDileptons> redDileptons;
  Produces<RedJpDmColCounts> redCollCounter;

  // HF configurables
  // cuts on BDT output scores to be applied only for the histograms
  Configurable<double> yCandDmesonMax{"yCandDmesonMax", -1., "max. cand. rapidity"};
  // DQ configurables
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "jpsiO2MCdebugCuts2", "Comma separated list of barrel track cuts"};
  // Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Comma separated list of pair cuts"}; It seems to me that they are not used
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<double> massDileptonCandMin{"massDileptonCandMin", 1, "minimum dilepton mass"};
  Configurable<double> massDileptonCandMax{"massDileptonCandMax", 5, "maximum dilepton mass"};
  // General configurables
  Configurable<bool> configDebug{"configDebug", true, "If true, fill D0 - J/psi histograms separately"};
  Configurable<bool> storeTableForNorm{"storeTableForNorm", true, "If true, store a table with number of processed collisions for normalisation"};

  Preslice<MyD0CandidatesSelected> perCollisionDmeson = aod::hf_cand::collisionId;
  Preslice<MyD0CandidatesSelectedWithBdt> perCollisionDmesonWithBdt = aod::hf_cand::collisionId;
  PresliceUnsorted<MyMuonTracksSelectedWithColl> perCollisionMuons = aod::reducedmuon::collisionId;
  PresliceUnsorted<MyBarrelTracksSelectedWithColl> perCollisionElectrons = aod::reducedtrack::collisionId;

  SliceCache cache;
  Filter filterD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= 1 || aod::hf_sel_candidate_d0::isSelD0bar >= 1;

  // Define histograms manager
  float* fValuesDileptonCharmHadron{};
  HistogramManager* fHistMan{};
  OutputObj<THashList> fOutputList{"output"};

  std::vector<AnalysisCompositeCut> fTrackCuts;
  std::vector<AnalysisCompositeCut> fMuonCuts;

  void init(o2::framework::InitContext&)
  {
    std::array<bool, 4> doprocess{doprocessJspiToMuMuD0, doprocessJspiToMuMuD0WithBdt, doprocessJspiToEED0, doprocessJspiToEED0WithBdt};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "Only one process function should be enabled! Please check your configuration!");
    }

    fValuesDileptonCharmHadron = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
    fHistMan->AddHistClass("Dmeson");
    fHistMan->AddHistClass("JPsi");
    fHistMan->AddHistClass("JPsiDmeson");
    dqhistograms::DefineHistograms(fHistMan, "Dmeson", "dilepton-charmhadron", "dmeson");
    if (doprocessJspiToMuMuD0 || doprocessJspiToMuMuD0WithBdt) {
      dqhistograms::DefineHistograms(fHistMan, "JPsi", "dilepton-charmhadron", "jpsitomumu");
      dqhistograms::DefineHistograms(fHistMan, "JPsiDmeson", "dilepton-charmhadron", "jpsitomumudmeson");
    }
    if (doprocessJspiToEED0 || doprocessJspiToEED0WithBdt) {
      dqhistograms::DefineHistograms(fHistMan, "JPsi", "dilepton-charmhadron", "jpsitoee");
      dqhistograms::DefineHistograms(fHistMan, "JPsiDmeson", "dilepton-charmhadron", "jpsitoeedmeson");
    }
    VarManager::SetUseVars(fHistMan->GetUsedVars());
    fOutputList.setObject(fHistMan->GetMainHistogramList());

    // cut strings
    if (doprocessJspiToMuMuD0 || doprocessJspiToMuMuD0WithBdt) {
      TString cutNamesMuon = fConfigMuonCuts.value;
      if (!cutNamesMuon.IsNull()) {
        std::unique_ptr<TObjArray> objArray(cutNamesMuon.Tokenize(","));
        for (int iCut{0}; iCut < objArray->GetEntries(); ++iCut) {
          fMuonCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(iCut)->GetName()));
        }
      }
    }
    if (doprocessJspiToEED0 || doprocessJspiToEED0WithBdt) {
      TString cutNamesElectron = fConfigTrackCuts.value;
      if (!cutNamesElectron.IsNull()) { // if track cuts
        std::unique_ptr<TObjArray> objArray(cutNamesElectron.Tokenize(","));
        for (int iCut{0}; iCut < objArray->GetEntries(); ++iCut) {
          fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(iCut)->GetName()));
        }
      }
    }

    // It seems they are not applied when filling the table
    // TString cutNamesPair = fConfigPairCuts.value;
    // if (!cutNamesPair.IsNull()) {
    //   std::unique_ptr<TObjArray> objArray(cutNamesPair.Tokenize(","));
    //   for (int iCut{0}; iCut < objArray->GetEntries(); ++iCut) {
    //     fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(iCut)->GetName()));
    //   }
    // }
  }

  // template function for lepton selection
  template <int TPairType, uint32_t TTrackFillMap /* gkMuonFillMapWithColl or gkTrackFillMapWithColl*/, typename TTrack>
  bool isLeptonSelected(TTrack const& lepton)
  {

    if constexpr (TPairType == VarManager::kDecayToEE) {
      VarManager::FillTrack<TTrackFillMap>(lepton);
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++) {
        if (!(*cut).IsSelected(VarManager::fgValues)) {
          return false;
        }
      }
    } else if constexpr (TPairType == VarManager::kDecayToMuMu) {
      VarManager::FillTrack<TTrackFillMap>(lepton);
      for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++) {
        if (!(*cut).IsSelected(VarManager::fgValues)) {
          return false;
        }
      }
    }

    return true;
  }

  // template function to compute the dilepton combinatorial
  template <int TPairType, uint32_t TTrackFillMap /* gkMuonFillMapWithColl or gkTrackFillMapWithColl*/, typename TTracks>
  std::vector<CandidateDilepton> computeDileptonCombinatorial(TTracks const& tracks)
  {
    std::vector<CandidateDilepton> dileptons{};
    uint32_t dileptonMcDecision{0u}; // placeholder, copy of the dqEfficiency.cxx one

    for (auto& [trackFirst, trackSecond] : combinations(tracks, tracks)) {
      bool isFirstSelected = isLeptonSelected<TPairType, TTrackFillMap>(trackFirst);
      bool isSecondSelected = isLeptonSelected<TPairType, TTrackFillMap>(trackSecond);

      if (!isFirstSelected || !isSecondSelected || (trackFirst.sign() + trackSecond.sign()) != 0) {
        continue;
      }

      // we fill a struct for the dileptons, which has the same signatures of the DQ tables to be used in the VarManager
      VarManager::FillPair<TPairType, TTrackFillMap>(trackFirst, trackSecond);

      if (VarManager::fgValues[VarManager::kMass] < massDileptonCandMin || VarManager::fgValues[VarManager::kMass] > massDileptonCandMax) {
        continue;
      }

      CandidateDilepton dilepton{};
      dilepton.fMass = VarManager::fgValues[VarManager::kMass];
      dilepton.fPt = VarManager::fgValues[VarManager::kPt];
      dilepton.fEta = VarManager::fgValues[VarManager::kEta];
      dilepton.fPhi = VarManager::fgValues[VarManager::kPhi];
      dilepton.fSign = trackFirst.sign() + trackSecond.sign();
      dilepton.fMcDecision = dileptonMcDecision;
      dileptons.push_back(dilepton);
    }

    return dileptons;
  }

  // Template function to run pair - hadron combinations
  // TODO: generalise to all charm-hadron species
  template <bool withDca, bool withBdt, int TPairType, uint32_t TTrackFillMap /* gkMuonFillMapWithColl or gkTrackFillMapWithColl*/, typename TDqTrack, typename THfTrack>
  void runDileptonDmeson(TDqTrack const& leptons, THfTrack const& dmesons, MyEvents::iterator const& collision, TracksWithExtra const&)
  {
    VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);

    bool isCollSel{false};
    if (configDebug) {
      for (auto const& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), DecayType::D0ToPiK)) {
          continue;
        }

        std::array<float, 6> scores = {999., -999., -999., 999., -999., -999.}; // D0 + D0bar
        if constexpr (withBdt) {
          if (dmeson.mlProbD0().size() == 3) {
            for (auto iScore{0u}; iScore < dmeson.mlProbD0().size(); ++iScore) {
              scores[iScore] = dmeson.mlProbD0()[iScore];
            }
          }
          if (dmeson.mlProbD0bar().size() == 3) {
            for (auto iScore{0u}; iScore < dmeson.mlProbD0bar().size(); ++iScore) {
              scores[iScore + 3] = dmeson.mlProbD0bar()[iScore];
            }
          }
        }

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(rapD0) > yCandDmesonMax) {
          continue;
        }

        if (dmeson.isSelD0() >= 1) {
          VarManager::FillSingleDileptonCharmHadron<VarManager::kD0ToPiK>(dmeson, hfHelper, scores[0], fValuesDileptonCharmHadron);
          fHistMan->FillHistClass("Dmeson", fValuesDileptonCharmHadron);
          VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
        }
        if (dmeson.isSelD0bar() >= 1) {
          VarManager::FillSingleDileptonCharmHadron<VarManager::kD0barToKPi>(dmeson, hfHelper, scores[3], fValuesDileptonCharmHadron);
          fHistMan->FillHistClass("Dmeson", fValuesDileptonCharmHadron);
          VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
        }
      }
    }

    // we compute the dileptons on the fly
    auto dileptons = computeDileptonCombinatorial<TPairType, TTrackFillMap>(leptons);

    // loop over dileptons
    std::vector<int> filledDmesonIds{}; // keep track of the D-mesons filled in the table to avoid repetitions
    for (auto dilepton : dileptons) {
      auto massJPsi = dilepton.mass();
      bool isJPsiFilled{false};

      if (massJPsi < massDileptonCandMin || massJPsi > massDileptonCandMax) {
        continue;
      }

      if (configDebug) {
        float dummyScore{-999.f};
        VarManager::FillSingleDileptonCharmHadron<VarManager::kJPsi>(dilepton, hfHelper, dummyScore, fValuesDileptonCharmHadron);
        fHistMan->FillHistClass("JPsi", fValuesDileptonCharmHadron);
        VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
      }

      // loop over D mesons
      for (auto const& dmeson : dmesons) {
        if (!TESTBIT(dmeson.hfflag(), DecayType::D0ToPiK)) {
          continue;
        }

        int dmesonIdx = dmeson.globalIndex();
        bool isDmesonFilled{false};
        if (std::find(filledDmesonIds.begin(), filledDmesonIds.end(), dmesonIdx) != filledDmesonIds.end()) { // we already included this D meson in the table, skip it
          isDmesonFilled = true;
        }

        auto rapD0 = hfHelper.yD0(dmeson);

        if (yCandDmesonMax >= 0. && std::abs(rapD0) > yCandDmesonMax) {
          continue;
        }

        auto massD0 = -1.;
        auto massD0bar = -1.;

        if (dmeson.isSelD0() >= 1 || dmeson.isSelD0bar() >= 1) {
          if (!isCollSel) {
            redCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib());
            isCollSel = true;
          }
          auto indexRed = redCollisions.lastIndex();
          if (!isJPsiFilled) {
            if constexpr (withDca) {
              redDileptons(indexRed, dilepton.px(), dilepton.py(), dilepton.pz(), dilepton.mass(), dilepton.sign(), dilepton.mcDecision(), dilepton.tauz(), dilepton.lz(), dilepton.lxy());
            } else {
              redDileptons(indexRed, dilepton.px(), dilepton.py(), dilepton.pz(), dilepton.mass(), dilepton.sign(), dilepton.mcDecision(), 0, 0, 0);
            }
            isJPsiFilled = true;
          }
          if (!isDmesonFilled) {
            redDmesons(indexRed, dmeson.px(), dmeson.py(), dmeson.pz(), dmeson.xSecondaryVertex(), dmeson.ySecondaryVertex(), dmeson.zSecondaryVertex(), 0, 0);
            auto trackProng0 = dmeson.template prong0_as<TracksWithExtra>();
            auto trackProng1 = dmeson.template prong1_as<TracksWithExtra>();
            // one table for each daughter with single track variables
            redDmesDau0(trackProng0.pt(), trackProng0.eta(), trackProng0.itsNCls(), trackProng0.tpcNClsCrossedRows(), dmeson.nSigTpcPi0(), dmeson.nSigTofPi0(), dmeson.nSigTpcKa0(), dmeson.nSigTofKa0());
            redDmesDau1(trackProng1.pt(), trackProng1.eta(), trackProng1.itsNCls(), trackProng1.tpcNClsCrossedRows(), dmeson.nSigTpcPi1(), dmeson.nSigTofPi1(), dmeson.nSigTpcKa1(), dmeson.nSigTofKa1());
            filledDmesonIds.push_back(dmesonIdx);
          }
          std::array<float, 6> scores = {999., -999., -999., 999., -999., -999.}; // D0 + D0bar
          if constexpr (withBdt) {
            if (!isDmesonFilled) {
              if (dmeson.mlProbD0().size() == 3) {
                for (auto iScore{0u}; iScore < dmeson.mlProbD0().size(); ++iScore) {
                  scores[iScore] = dmeson.mlProbD0()[iScore];
                }
              }
              if (dmeson.mlProbD0bar().size() == 3) {
                for (auto iScore{0u}; iScore < dmeson.mlProbD0bar().size(); ++iScore) {
                  scores[iScore + 3] = dmeson.mlProbD0bar()[iScore];
                }
              }
              redDmesBdts(scores[0], scores[1], scores[2], scores[3], scores[4], scores[5]);
            }
          }

          if (dmeson.isSelD0() >= 1) {
            massD0 = hfHelper.invMassD0ToPiK(dmeson);
            VarManager::FillDileptonCharmHadron<VarManager::kD0ToPiK>(dilepton, dmeson, hfHelper, scores[0], fValuesDileptonCharmHadron);
            fHistMan->FillHistClass("JPsiDmeson", fValuesDileptonCharmHadron);
            VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
          }
          if (dmeson.isSelD0bar() >= 1) {
            massD0bar = hfHelper.invMassD0barToKPi(dmeson);
            VarManager::FillDileptonCharmHadron<VarManager::kD0barToKPi>(dilepton, dmeson, hfHelper, scores[3], fValuesDileptonCharmHadron);
            fHistMan->FillHistClass("JPsiDmeson", fValuesDileptonCharmHadron);
            VarManager::ResetValues(0, VarManager::kNVars, fValuesDileptonCharmHadron);
          }
          if (!isDmesonFilled) {
            redD0Masses(massD0, massD0bar);
          }
        }
      }
    }
  }

  // process J/psi(->mumu) - D0
  void processJspiToMuMuD0(MyEvents const& collisions, MyMuonTracksSelectedWithColl const& muonCandidates, soa::Filtered<MyD0CandidatesSelected> const& selectedD0Candidates, TracksWithExtra const& barrelTracks)
  {
    if (storeTableForNorm) {
      redCollCounter(collisions.size());
    }
    for (auto const& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0Candidates.sliceBy(perCollisionDmeson, collision.globalIndex());
      auto groupedLeptonCandidates = muonCandidates.sliceBy(perCollisionMuons, collision.globalIndex());
      runDileptonDmeson<false, false, VarManager::kDecayToMuMu, gkMuonFillMapWithColl>(groupedLeptonCandidates, groupedDmesonCandidates, collision, barrelTracks);
    }
  }

  // process J/psi(->ee) - D0
  void processJspiToEED0(MyEvents const& collisions, MyBarrelTracksSelectedWithColl const& electronCandidates, soa::Filtered<MyD0CandidatesSelected> const& selectedD0Candidates, TracksWithExtra const& barrelTracks)
  {
    if (storeTableForNorm) {
      redCollCounter(collisions.size());
    }
    for (auto const& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0Candidates.sliceBy(perCollisionDmeson, collision.globalIndex());
      auto groupedLeptonCandidates = electronCandidates.sliceBy(perCollisionElectrons, collision.globalIndex());
      runDileptonDmeson<false, false, VarManager::kDecayToEE, gkTrackFillMapWithColl>(groupedLeptonCandidates, groupedDmesonCandidates, collision, barrelTracks);
    }
  }

  // process J/psi(->mumu) - D0 adding the BDT output scores to the D0 table
  void processJspiToMuMuD0WithBdt(MyEvents const& collisions, MyMuonTracksSelectedWithColl const& muonCandidates, soa::Filtered<MyD0CandidatesSelectedWithBdt> const& selectedD0CandidatesWithBdt, TracksWithExtra const& barrelTracks)
  {
    if (storeTableForNorm) {
      redCollCounter(collisions.size());
    }
    for (auto const& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0CandidatesWithBdt.sliceBy(perCollisionDmesonWithBdt, collision.globalIndex());
      auto groupedLeptonCandidates = muonCandidates.sliceBy(perCollisionMuons, collision.globalIndex());
      runDileptonDmeson<false, true, VarManager::kDecayToMuMu, gkMuonFillMapWithColl>(groupedLeptonCandidates, groupedDmesonCandidates, collision, barrelTracks);
    }
  }

  // process J/psi(->ee) - D0 adding the BDT output scores to the D0 table
  void processJspiToEED0WithBdt(MyEvents const& collisions, MyBarrelTracksSelectedWithColl const& electronCandidates, soa::Filtered<MyD0CandidatesSelectedWithBdt> const& selectedD0CandidatesWithBdt, TracksWithExtra const& barrelTracks)
  {
    if (storeTableForNorm) {
      redCollCounter(collisions.size());
    }
    for (auto const& collision : collisions) {
      auto groupedDmesonCandidates = selectedD0CandidatesWithBdt.sliceBy(perCollisionDmesonWithBdt, collision.globalIndex());
      auto groupedLeptonCandidates = electronCandidates.sliceBy(perCollisionElectrons, collision.globalIndex());
      runDileptonDmeson<false, true, VarManager::kDecayToEE, gkTrackFillMapWithColl>(groupedLeptonCandidates, groupedDmesonCandidates, collision, barrelTracks);
    }
  }

  PROCESS_SWITCH(tableMakerJpsiHf, processJspiToMuMuD0, "Process J/psi(->mumu) - D0", false);
  PROCESS_SWITCH(tableMakerJpsiHf, processJspiToEED0, "Process J/psi(->ee) - D0", false);
  PROCESS_SWITCH(tableMakerJpsiHf, processJspiToMuMuD0WithBdt, "Process J/psi(->mumu) - D0 with BDT", false);
  PROCESS_SWITCH(tableMakerJpsiHf, processJspiToEED0WithBdt, "Process J/psi(->ee) - D0 with BDT", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tableMakerJpsiHf>(cfgc)};
}
