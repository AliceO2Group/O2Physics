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
#include <CommonConstants/PhysicsConstants.h>
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
  Configurable<float> cfgTPCNsigma{"cfgTPCNsigma", 3.0f, "TPC n sigma for deuteron PID"};
  Configurable<float> cfgTofNsigma{"cfgTofNsigma", 3.0f, "TOF n sigma for deuteron PID"};
  Configurable<float> ptThresholdPid{"ptThresholdPid", 0.5f, "pT threshold to switch between TPC and TPC+TOF PID"};
  Configurable<float> cfgDCAmin{"cfgDCAmin", 0.05f, "Minimum DCA for deuteron PID"};
  Configurable<float> cfgDCAmax{"cfgDCAmax", 1000.0f, "Maximum DCA for deuteron PID"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5f, "Rapidity cut"};
  // PDG codes
  Configurable<int> pdgCodeBeautyMeson{"pdgCodeBeautyMeson", -521, "PDG code of the beauty meson mother particle (default: B-)"};
  Configurable<int> pdgCodeBeautyBaryon{"pdgCodeBeautyBaryon", -5122, "PDG code of the beauty baryon mother particle (default: anti-Lambda_b)"};
  Configurable<int> pdgCodeDaughter{"pdgCodeDaughter", -1000010020, "PDG code of the daughter particle (default: anti-deuteron)"};

  int mRunNumber = 0;
  float dBz = 0.f;
  int mCurrentRun = -1;

  framework::Service<ccdb::BasicCCDBManager> ccdb{};

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using MCTrackCandidates = o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::McTrackLabels>;
  using MCCollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe>;

  Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;

  ConfigurableAxis ptAxis{"ptAxis", {100, 0.f, 10.f}, "p_{T} GeV/c"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.f, 10.f}, "nSigma"};
  ConfigurableAxis dcaXyAxis{"dcaXyAxis", {1000, -0.2f, 0.2f}, "DCA xy (cm)"};
  ConfigurableAxis dcaZAxis{"dcaZAxis", {1000, -0.2f, 0.2f}, "DCA z (cm)"};

  HistogramRegistry qaHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Event filtered;; Number of events", 4, 0., 4.)};

  bool isSelectedBeautyHadron(int pdgCode)
  {
    return pdgCode == pdgCodeBeautyMeson || pdgCode == pdgCodeBeautyBaryon;
  }

  double getBeautyHadronMass(int pdgCode)
  {
    if (pdgCode == pdgCodeBeautyMeson) {
      return o2::constants::physics::MassBPlus;
    }

    if (pdgCode == pdgCodeBeautyBaryon) {
      return o2::constants::physics::MassLambdaB0;
    }
    return -1.;
  }

  void init(framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    if (applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    qaHistos.add("Data/hDCAxyVsPt", "DCAxy #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("Data/hDCAzVsPt", "DCAz #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("Data/hnSigmaTPCVsPt", "n#sigma TPC vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/hnSigmaTOFVsPt", "n#sigma TOF vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/ptAntiDeuteron", "ptAntiDeuteron", {HistType::kTH1F, {ptAxis}});
    qaHistos.add("Data/etaAntideuteron", "etaAntideuteron", {HistType::kTH1F, {{100, -1.0f, 1.0f, "eta #bar{d}"}}});
    qaHistos.add("Data/hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});
    // MC generated-level histograms
    qaHistos.add("MCGen/ptGeneratedBminus", "p_{T} generated B^{-};p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptGeneratedAntiLambdaB", "p_{T} generated #bar{#Lambda}_{b};p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});

    qaHistos.add("MCGen/ptAntiDeuteronPrimary", "p_{T} #bar{d} primary gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronFromBeautyHadron", "p_{T} #bar{d} from beauty gen hadrons;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});

    qaHistos.add("MCGen/hMotherPdgCode", "PDG code of mother, gen level;PDG code;Counts", HistType::kTH1I, {{12000, -6000, 6000}});
    qaHistos.add("MCGen/ctauBminus", "ctau of B^{-}, gen level;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCGen/ctauAntiLambdaB", "ctau of #bar{#Lambda}_{b}, gen level;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    // MC reco/MC-anchored histograms
    qaHistos.add("MCReco/ptAntiDeuteronFromBminus",
                 "p_{T} #bar{d} from B^{-} reco/MC anchored;p_{T} (GeV/c);Counts",
                 HistType::kTH1F, {ptAxis});

    qaHistos.add("MCReco/ptAntiDeuteronFromAntiLambdaB",
                 "p_{T} #bar{d} from #bar{#Lambda}_{b} reco/MC anchored;p_{T} (GeV/c);Counts",
                 HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronFromBminus", "p_{T} #bar{d} from B^{-} gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronFromAntiLambdaB", "p_{T} #bar{d} from #bar{#Lambda}_{b} gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/hDCAxy-Primary", "DCAxy primary reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    qaHistos.add("MCReco/hDCAxy-FromBeautyHadron", "DCAxy from beauty reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {{400, -0.2f, 0.2f, "DCA xy (cm)"}}});
    qaHistos.add("MCReco/hMotherPdgCode", "PDG code of mother, reco/MC anchored;PDG code;Counts", HistType::kTH1I, {{12000, -6000, 6000}});
    qaHistos.add("MCReco/ctauBminus", "ctau of B^{-}, reco/MC anchored;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCReco/ctauAntiLambdaB", "ctau of #bar{#Lambda}_{b}, reco/MC anchored;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});

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
    dBz = 0.f;

    if (applySkimming) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfgSkimming.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
  }

  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
  {
    if (std::abs(track.eta()) > cfgEta) {
      return false;
    }
    if (!track.hasITS()) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (!track.hasTOF()) {
      return false;
    }
    if (track.tpcNClsFound() < cfgTPCNclsFound) {
      return false;
    }
    if (track.tpcChi2NCl() > cfgTPCChi2Ncl) {
      return false;
    }
    if (track.itsChi2NCl() > cfgITSChi2Ncl) {
      return false;
    }
    if (track.itsNCls() < cfgITScls) {
      return false;
    }
    if (track.pt() > cfgMaxPt) {
      return false;
    }
    if (track.pt() < cfgMinPt) {
      return false;
    }
    if (track.sign() > 0) {
      return false;
    }

    return true;
  }

  void processData(CollisionCandidates const& collisions,
                   o2::aod::TrackAssoc const& trackIndices,
                   TrackCandidates const& tracks,
                   o2::aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      if (mCurrentRun != collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber()) {
        auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", collision.bc_as<o2::aod::BCsWithTimestamps>().timestamp());
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
      qaHistos.fill(HIST("Data/hVtxZ"), collision.posZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = tracks.rawIteratorAt(trackId.trackId());
        std::array<float, 2> dca{track.dcaXY(), track.dcaZ()};

        if (track.collisionId() != collision.globalIndex()) {
          auto trackPar = getTrackParCov(track);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dca);
        }

        if (std::abs(dca[0]) < cfgDCAmin || std::abs(dca[0]) > cfgDCAmax) {
          continue;
        }

        if (!passedSingleTrackSelection(track)) {
          continue;
        }

        const bool isTPCDe = std::abs(track.tpcNSigmaDe()) < cfgTPCNsigma;
        const bool isTOFDe = std::abs(track.tofNSigmaDe()) < cfgTofNsigma;

        if (track.pt() < ptThresholdPid) {
          if (isTPCDe) {
            qaHistos.fill(HIST("Data/ptAntiDeuteron"), track.pt());
            qaHistos.fill(HIST("Data/etaAntideuteron"), track.eta());
            qaHistos.fill(HIST("Data/hDCAxyVsPt"), track.pt(), dca[0]);
            qaHistos.fill(HIST("Data/hDCAzVsPt"), track.pt(), dca[1]);
            qaHistos.fill(HIST("Data/hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            qaHistos.fill(HIST("Data/hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        } else {
          if (isTPCDe && isTOFDe) {
            qaHistos.fill(HIST("Data/ptAntiDeuteron"), track.pt());
            qaHistos.fill(HIST("Data/etaAntideuteron"), track.eta());
            qaHistos.fill(HIST("Data/hDCAxyVsPt"), track.pt(), dca[0]);
            qaHistos.fill(HIST("Data/hDCAzVsPt"), track.pt(), dca[1]);
            qaHistos.fill(HIST("Data/hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            qaHistos.fill(HIST("Data/hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processData, "processData", true);

  void processMC(MCCollisionCandidates::iterator const& collision, MCTrackCandidates const& tracks, o2::aod::McParticles const&)
  {
    for (const auto& track : tracks) {
      if (std::abs(track.dcaXY()) < cfgDCAmin || std::abs(track.dcaXY()) > cfgDCAmax) {
        continue;
      }
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
          int motherPdg = 0;
          if (separateAntideuterons) {
            for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
              qaHistos.fill(HIST("MCReco/hMotherPdgCode"), mom.pdgCode());
              if (mom.pdgCode() == pdgCodeBeautyMeson || mom.pdgCode() == pdgCodeBeautyBaryon) {
                motherPdg = mom.pdgCode();
                double dx = mcParticle.vx() - collision.posX();
                double dy = mcParticle.vy() - collision.posY();
                double dz = mcParticle.vz() - collision.posZ();

                double flightDistance = std::sqrt(dx * dx + dy * dy + dz * dz);
                double massBeauty = getBeautyHadronMass(mom.pdgCode());

                if (massBeauty > 0.) {
                  double ctauBeauty = flightDistance / mom.p() * massBeauty;
                  ctauBeauty *= 1.e4;

                  if (motherPdg == pdgCodeBeautyMeson) {
                    qaHistos.fill(HIST("MCReco/ctauBminus"), ctauBeauty);
                  } else if (motherPdg == pdgCodeBeautyBaryon) {
                    qaHistos.fill(HIST("MCReco/ctauAntiLambdaB"), ctauBeauty);
                  }
                }
                break;
              }
            }
          }
          if (motherPdg == pdgCodeBeautyMeson) {
            qaHistos.fill(HIST("MCReco/hDCAxy-FromBeautyHadron"), track.dcaXY());
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronFromBminus"), track.pt());
          } else if (motherPdg == pdgCodeBeautyBaryon) {
            qaHistos.fill(HIST("MCReco/hDCAxy-FromBeautyHadron"), track.dcaXY());
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronFromAntiLambdaB"), track.pt());
          } else {
            qaHistos.fill(HIST("MCReco/hDCAxy-Primary"), track.dcaXY());
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronPrimary"), track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processMC, "processMC", false);

  void processGen(o2::aod::McCollision const& collision, o2::aod::McParticles const& mcParticles)
  {
    hProcessedEvents->Fill(0.5);

    for (const auto& mcParticle : mcParticles) {

      // Beauty hadron mother
      if (isSelectedBeautyHadron(mcParticle.pdgCode())) {
        if (std::abs(mcParticle.y()) <= rapidityCut) {
          if (mcParticle.pdgCode() == pdgCodeBeautyMeson) {
            qaHistos.fill(HIST("MCGen/ptGeneratedBminus"), mcParticle.pt());
          } else if (mcParticle.pdgCode() == pdgCodeBeautyBaryon) {
            qaHistos.fill(HIST("MCGen/ptGeneratedAntiLambdaB"), mcParticle.pt());
          }
        }

        if (mcParticle.has_daughters()) {

          for (const auto& daughter : mcParticle.daughters_as<o2::aod::McParticles>()) {

            double dx = daughter.vx() - collision.posX();
            double dy = daughter.vy() - collision.posY();
            double dz = daughter.vz() - collision.posZ();

            double flightDistance = std::sqrt(dx * dx + dy * dy + dz * dz);
            double massBeauty = getBeautyHadronMass(mcParticle.pdgCode());
            if (massBeauty > 0.) {
              double ctauBeauty = flightDistance / mcParticle.p() * massBeauty;
              ctauBeauty *= 1.e4; // cm -> micrometers

              if (mcParticle.pdgCode() == pdgCodeBeautyMeson) {
                qaHistos.fill(HIST("MCGen/ctauBminus"), ctauBeauty);
              } else if (mcParticle.pdgCode() == pdgCodeBeautyBaryon) {
                qaHistos.fill(HIST("MCGen/ctauAntiLambdaB"), ctauBeauty);
              }
            }
            break;
          }
        }
      }

      if (mcParticle.pdgCode() == pdgCodeDaughter) {

        if (std::abs(mcParticle.y()) > rapidityCut) {
          continue;
        }

        int motherPdg = 0;
        if (mcParticle.has_mothers()) {
          for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {

            qaHistos.fill(HIST("MCGen/hMotherPdgCode"), mom.pdgCode());

            if (mom.pdgCode() == pdgCodeBeautyMeson) {
              motherPdg = mom.pdgCode();
              break;
            }

            if (mom.pdgCode() == pdgCodeBeautyBaryon) {
              motherPdg = mom.pdgCode();
              break;
            }
          }
        }

        if (motherPdg == pdgCodeBeautyMeson) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronFromBminus"), mcParticle.pt());
        } else if (motherPdg == pdgCodeBeautyBaryon) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronFromAntiLambdaB"), mcParticle.pt());
        } else if (mcParticle.isPhysicalPrimary()) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronPrimary"), mcParticle.pt());
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
