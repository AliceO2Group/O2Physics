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
//
/// \file taskDeuteronFromLb.cxx
/// \brief A filter task for non prompt deuterons from beauty-hadron decays
/// \author Marta Razza <marta.razza@cern.ch>
/// \author Francesca Ercolessi <francesca.ercolessi@cern.ch>
/// \since May 25, 2026

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <TH1.h>

#include <array>
#include <chrono>
#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskDeuteronFromLb {

  Zorro zorro;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> applySkimming{"applySkimming", false, "Skimmed dataset processing"};
  Configurable<std::string> cfgSkimming{"cfgSkimming", "fH2fromLb", "Configurable for skimming"};
  Configurable<bool> sel8{"sel8", true, "sel8 event selection"};
  Configurable<bool> separateAntideuterons{"separateAntideuterons", true, "Fill FromLb histos for primary deuterons whose mother is Lb, Primary histos for all others"};
  Configurable<float> cfgEta{"cfgEta", 0.8f, "Track eta selection"};
  Configurable<float> cfgTPCNclsFound{"cfgTPCNclsFound", 100, "Minimum TPC clusters found"};
  Configurable<float> cfgTPCChi2Ncl{"cfgTPCChi2Ncl", 4.0f, "Maximum TPC chi2 per N clusters"};
  Configurable<float> cfgITSChi2Ncl{"cfgITSChi2Ncl", 36.0f, "Maximum ITS chi2 per N clusters"};
  Configurable<float> cfgITScls{"cfgITScls", 2, "Minimum ITS clusters"};
  Configurable<float> cfgMaxPt{"cfgMaxPt", 5.0f, "Maximum pT cut"};
  Configurable<float> cfgMinPt{"cfgMinPt", 0.5f, "Minimum pT cut"};
  Configurable<float> cfgTPCNsigma{"cfgTPCNsigma", 4.0f, "TPC n sigma for deuteron PID"};
  Configurable<float> cfgTofNsigmaMin{"cfgTofNsigmaMin", 3.0f, "TOF n sigma min for deuteron PID"};
  Configurable<float> cfgTofNsigmaMax{"cfgTofNsigmaMax", 4.0f, "TOF n sigma max for deuteron PID"};
  Configurable<float> ptThresholdPid{"ptThresholdPid", 0.5f, "pT threshold to switch between 4 and 3 sigmas for TOF PID"};
  Configurable<float> cfgDCAmin{"cfgDCAmin", 0.05f, "Minimum DCA for deuteron PID"};
  Configurable<float> cfgDCAmax{"cfgDCAmax", 1000.0f, "Maximum DCA for deuteron PID"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5f, "Rapidity cut"};
  // PDG codes
  Configurable<int> pdgCodeMother{"pdgCodeMother", -5122, "PDG code of the mother particle (default: anti-Lambda_b)"};
  Configurable<int> pdgCodeDaughter{"pdgCodeDaughter", -1000010020, "PDG code of the daughter particle (default: anti-deuteron)"};

  int mRunNumber = 0;
  float d_bz = 0.f;
  int mCurrentRun = -1;

  framework::Service<ccdb::BasicCCDBManager> ccdb;

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using MCTrackCandidates = o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::McTrackLabels>;
  using MCCollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe>;

  Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;

  ConfigurableAxis ptAxis{"ptAxis", {100, 0., 10.f}, "p_{T} GeV/c"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.f, 10.f}, "nSigma"};
  ConfigurableAxis dcaXyAxis{"dcaXyAxis", {1000, -0.2f, 0.2f}, "DCA xy (cm)"};
  ConfigurableAxis dcaZAxis{"dcaZAxis", {1000, -0.2f, 0.2f}, "DCA z (cm)"};

  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Event filtered;; Number of events", 4, 0., 4.)};

  void init(framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    if (applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    QAHistos.add("MC/ptGeneratedLb", "ptGeneratedLb", HistType::kTH1F, {ptAxis});
    QAHistos.add("MC/ptAntiDeuteronPrimary", "ptAntiDeuteronPrimaryReco", HistType::kTH1F, {ptAxis});
    QAHistos.add("MC/ptAntiDeuteronFromLb", "ptAntiDeuteronFromLbReco", HistType::kTH1F, {ptAxis});
    QAHistos.add("MC/hDCAxy-Primary", "DCAxy-Primary", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    QAHistos.add("MC/hDCAxy-FromLb", "DCAxy-FromLb", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    QAHistos.add("Data/hDCAxyVsPt", "DCAxy #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    QAHistos.add("Data/hDCAzVsPt", "DCAz #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    QAHistos.add("Data/hnSigmaTPCVsPt", "n#sigma TPC vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    QAHistos.add("Data/hnSigmaTOFVsPt", "n#sigma TOF vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    QAHistos.add("Data/ptAntiDeuteron", "ptAntiDeuteron", {HistType::kTH1F, {ptAxis}});
    QAHistos.add("Data/etaAntideuteron", "etaAntideuteron", {HistType::kTH1F, {{100, -1.0f, 1.0f, "eta #bar{d}"}}});
    QAHistos.add("Data/hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "ZORRO");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "sel8");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "z-vertex");
  }

  void initCCDB(o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    d_bz = 0.f;

    if (applySkimming) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfgSkimming.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
  }

  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
  {
    if (std::abs(track.eta()) > cfgEta)
      return false;
    if (std::abs(track.dcaXY()) < cfgDCAmin || std::abs(track.dcaXY()) > cfgDCAmax)
      return false;
    if (!track.hasITS())
      return false;
    if (!track.hasTPC())
      return false;
    if (!track.hasTOF())
      return false;
    if (track.tpcNClsFound() < cfgTPCNclsFound)
      return false;
    if (track.tpcChi2NCl() > cfgTPCChi2Ncl)
      return false;
    if (track.itsChi2NCl() > cfgITSChi2Ncl)
      return false;
    if (track.itsNCls() < cfgITScls)
      return false;
    if (track.pt() > cfgMaxPt)
      return false;
    if (track.pt() < cfgMinPt)
      return false;
    if (track.sign() > 0)
      return false;

    return true;
  }

  void processData(CollisionCandidates const& collisions,
                   o2::aod::TrackAssoc const& trackIndices,
                   TrackCandidates const& tracks,
                   o2::aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      if (mCurrentRun != collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", collision.bc_as<o2::aod::BCsWithTimestamps>().timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        mCurrentRun = collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber();
      }

      const auto& bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc);
      hProcessedEvents->Fill(0.5);
      if (applySkimming) {
        if (!zorro.isSelected(bc.globalBC())) {
          continue;
        }
      }
      hProcessedEvents->Fill(1.5);
      if (sel8 && !collision.sel8()) {
        continue;
      }
      hProcessedEvents->Fill(2.5);
      if (std::abs(collision.posZ()) > cutzvertex) {
        continue;
      }
      hProcessedEvents->Fill(3.5);
      QAHistos.fill(HIST("Data/hVtxZ"), collision.posZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = tracks.rawIteratorAt(trackId.trackId());
        std::array<float, 2> dca{track.dcaXY(), track.dcaZ()};

        if (track.collisionId() != collision.globalIndex()) {
          auto trackPar = getTrackParCov(track);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dca);
        }

        if (!passedSingleTrackSelection(track)) {
          continue;
        }

        const bool isTPCDe = std::abs(track.tpcNSigmaDe()) < cfgTPCNsigma;
        const bool isTOFDe_min = std::abs(track.tofNSigmaDe()) > cfgTofNsigmaMin;
        const bool isTOFDe_max = std::abs(track.tofNSigmaDe()) < cfgTofNsigmaMax;

        if (track.pt() < ptThresholdPid) {
          if (isTPCDe && isTOFDe_max) {
            QAHistos.fill(HIST("Data/ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("Data/etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("Data/hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("Data/hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("Data/hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("Data/hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        } else {
          if (isTPCDe && isTOFDe_min && isTOFDe_max) {
            QAHistos.fill(HIST("Data/ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("Data/etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("Data/hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("Data/hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("Data/hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("Data/hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processData, "processData", true);

  void processMC(MCCollisionCandidates::iterator const&, MCTrackCandidates const& tracks, o2::aod::McParticles const&)
  {
    for (const auto& track : tracks) {
      if (!passedSingleTrackSelection(track)) {
        continue;
      }

      if (!track.has_mcParticle()) {
        continue;
      }

      auto mcParticle = track.mcParticle();
      if (mcParticle.pdgCode() == pdgCodeDaughter) {
        if (std::abs(mcParticle.y()) > rapidityCut) {
          continue;
        }
        if (mcParticle.isPhysicalPrimary()) {
          bool isFromLb = false;
          if (separateAntideuterons) {
            for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
              if (mom.pdgCode() == pdgCodeMother) {
                isFromLb = true;
                break;
              }
            }
          }
          if (isFromLb) {
            QAHistos.fill(HIST("MC/hDCAxy-FromLb"), track.dcaXY());
            QAHistos.fill(HIST("MC/ptAntiDeuteronFromLb"), track.pt());
          } else {
            QAHistos.fill(HIST("MC/hDCAxy-Primary"), track.dcaXY());
            QAHistos.fill(HIST("MC/ptAntiDeuteronPrimary"), track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processMC, "processMC", false);

  void processGen(o2::aod::McCollision const&, o2::aod::McParticles const& mcParticles)
  {
    hProcessedEvents->Fill(0.5);
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == pdgCodeMother) {
        if (std::abs(mcParticle.y()) <= rapidityCut) {
          QAHistos.fill(HIST("MC/ptGeneratedLb"), mcParticle.pt());
        }
      }

      if (mcParticle.pdgCode() == pdgCodeDaughter) {
        if (std::abs(mcParticle.y()) > rapidityCut)
          continue;

        bool isFromLb = false;
        if (mcParticle.has_mothers()) {
          for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
            if (mom.pdgCode() == pdgCodeMother) {
              isFromLb = true;
              break;
            }
          }
        }

        if (isFromLb) {
          QAHistos.fill(HIST("MC/ptAntiDeuteronFromLb"), mcParticle.pt());
        } else if (mcParticle.isPhysicalPrimary()) {
          QAHistos.fill(HIST("MC/ptAntiDeuteronPrimary"), mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processGen, "processGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskDeuteronFromLb>(cfgc)};
}
