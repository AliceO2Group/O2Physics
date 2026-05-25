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
/// \brief A filter task for non prompt deuterons
/// \author Marta Razza marta.razza@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \since Dec 17, 2025

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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct H2fromLb {

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Event filtered;; Number of events", 4, 0., 4.)};

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
  Configurable<float> cfgTOFNsigma_min{"cfgTOFNsigma_min", 3.0f, "TOF n sigma min for deuteron PID"};
  Configurable<float> cfgTOFNsigma_max{"cfgTOFNsigma_max", 4.0f, "TOF n sigma max for deuteron PID"};
  Configurable<float> ptThresholdPid{"ptThresholdPid", 1.0f, "pT threshold to switch between 4 and 3 sigmas for TOF PID"};
  // PDG codes
  Configurable<int> pdgCodeMother{"pdgCodeMother", -5122, "PDG code of the mother particle (default: anti-Lambda_b)"};
  Configurable<int> pdgCodeDaughter{"pdgCodeDaughter", -1000010020, "PDG code of the daughter particle (default: anti-deuteron)"};

  int mRunNumber = 0;
  float d_bz = 0.f;
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  framework::Service<ccdb::BasicCCDBManager> ccdb;
  void init(framework::InitContext&)
  {

    ccdb->setURL("http://alice-ccdb.cern.ch"); // Set CCDB URL to get magnetic field
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    if (applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    ConfigurableAxis ptAxis{"ptAxis", {100, 0., 10.f}, "p_{T} GeV/c"};
    ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.f, 10.f}, "nSigma"};
    ConfigurableAxis DCAxyAxis{"DCAxyAxis", {1000, -0.2f, 0.2f}, "DCA xy (cm)"};
    ConfigurableAxis DCAzAxis{"DCAzAxis", {1000, -0.2f, 0.2f}, "DCA z (cm)"};

    // general QA histograms
    QAHistos.add("hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});
    QAHistos.add("ptGeneratedLb", "ptGeneratedLb", HistType::kTH1F, {ptAxis});
    QAHistos.add("ptAntiDeuteronPrimary", "ptAntiDeuteronPrimaryReco", HistType::kTH1F, {ptAxis});
    QAHistos.add("ptAntiDeuteronFromLb", "ptAntiDeuteronFromLbReco", HistType::kTH1F, {ptAxis});
    QAHistos.add("hDCAxy-Primary", "DCAxy-Primary", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    QAHistos.add("hDCAxy-FromLb", "DCAxy-FromLb", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    QAHistos.add("hDCAxyVsPt", "DCAxy #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, DCAxyAxis}});
    QAHistos.add("hDCAzVsPt", "DCAz #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, DCAzAxis}});
    QAHistos.add("hnSigmaTPCVsPt", "n#sigma TPC vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    QAHistos.add("hnSigmaTOFVsPt", "n#sigma TOF vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    QAHistos.add("ptAntiDeuteron", "ptAntiDeuteron", {HistType::kTH1F, {ptAxis}});
    QAHistos.add("etaAntideuteron", "etaAntideuteron", {HistType::kTH1F, {{100, -1.0f, 1.0f, "eta #bar{d}"}}});

    // processed events
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "ZORRO");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "sel8");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "z-vertex");
  }
  void initCCDB(o2::aod::BCsWithTimestamps::iterator const& bc) // inspired by PWGLF/TableProducer/lambdakzerobuilder.cxx
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
  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  // Tables for MC processing
  using MCTrackCandidates = o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::McTrackLabels>;
  using MCCollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabels>;
  // Tables for Data processing
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe>;

  // Single-Track Selection
  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
  {
    // Single-Track Selections
    if (std::abs(track.eta()) > cfgEta)
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
  Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;
  int mCurrentRun = -1;
  void processData(CollisionCandidates const& collisions,
                   o2::aod::TrackAssoc const& trackIndices,
                   TrackCandidates const& tracks,
                   o2::aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) // start loop over collisions
    {
      if (mCurrentRun != collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber()) { // If the run is new then we need to initialize the propagator field
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
      QAHistos.fill(HIST("hVtxZ"), collision.posZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks

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
        const bool isTOFDe_min = std::abs(track.tofNSigmaDe()) > cfgTOFNsigma_min;
        const bool isTOFDe_max = std::abs(track.tofNSigmaDe()) < cfgTOFNsigma_max;

        if (track.pt() < ptThresholdPid) {
          if (isTPCDe && isTOFDe_max) {
            QAHistos.fill(HIST("ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        } else {
          if (isTPCDe && isTOFDe_min && isTOFDe_max) {
            QAHistos.fill(HIST("ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(H2fromLb, processData, "processData", false);

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
        if (std::abs(mcParticle.y()) > 0.5) {
          continue;
        }
        if (mcParticle.isPhysicalPrimary()) {
          bool isFromLb = false;
          if (separateAntideuterons) {
            for (auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
              if (mom.pdgCode() == pdgCodeMother) { // Lambda_b
                isFromLb = true;

                break;
              }
            }
          }
          if (isFromLb) {
            QAHistos.fill(HIST("hDCAxy-FromLb"), track.dcaXY());
            QAHistos.fill(HIST("ptAntiDeuteronFromLb"), track.pt());
          } else {
            QAHistos.fill(HIST("hDCAxy-Primary"), track.dcaXY());
            QAHistos.fill(HIST("ptAntiDeuteronPrimary"), track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(H2fromLb, processMC, "processMC", true);

  void processGen(o2::aod::McCollision const&, o2::aod::McParticles const& mcParticles)
  {
    hProcessedEvents->Fill(0.5);
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == pdgCodeMother) {
        if (std::abs(mcParticle.y()) <= 0.5) {
          QAHistos.fill(HIST("ptGeneratedLb"), mcParticle.pt()); // rinominato
        }
      }

      if (mcParticle.pdgCode() == pdgCodeDaughter) {
        if (std::abs(mcParticle.y()) > 0.5)
          continue;

        bool isFromLb = false;
        if (mcParticle.has_mothers()) {
          for (auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
            if (mom.pdgCode() == pdgCodeMother) {
              isFromLb = true;
              break;
            }
          }
        }

        if (isFromLb) {
          QAHistos.fill(HIST("ptAntiDeuteronFromLb"), mcParticle.pt());
        } else if (mcParticle.isPhysicalPrimary()) {
          QAHistos.fill(HIST("ptAntiDeuteronPrimary"), mcParticle.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(H2fromLb, processGen, "processGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<H2fromLb>(cfgc, TaskName{"task-h2-from-Lb"})};
}
