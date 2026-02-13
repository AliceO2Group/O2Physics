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

#include "../filterTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "TVector3.h"

#include <cmath>

struct H2fromLbFilter {

  // Recall the output table
  o2::framework::Produces<o2::aod::H2fromLbFilters> table;

  o2::framework::Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Define a histograms and registries
  o2::framework::HistogramRegistry QAHistos{"QAHistos", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, true};
  o2::framework::OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Event filtered;; Number of events", 6, 0., 6.)};

  o2::framework::ConfigurableAxis pTAxis{"pTAxis", {200, -10.f, 10.f}, "p_{T} GeV/c"};
  o2::framework::ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.f, 10.f}, "p_{T} GeV/c"};
  o2::framework::ConfigurableAxis DCAxyAxis{"DCAxyAxis", {1000, -0.2f, 0.2f}, "DCA xy (cm)"};
  o2::framework::ConfigurableAxis DCAzAxis{"DCAzAxis", {1000, -0.2f, 0.2f}, "DCA z (cm)"};

  o2::framework::Configurable<float> cutzvertex{"cutzvertex", 20.0f, "Accepted z-vertex range"}; // 20 cm
  o2::framework::Configurable<bool> isTVX{"isTVX", true, "isTVX event selection"};
  o2::framework::Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", true, "isNoTimeFrameBorder event selection"};
  o2::framework::Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", true, "isNoITSROFrameBorder event selection"};
  o2::framework::Configurable<float> cfgTPCNsigma{"cfgTPCNsigma", 4.0f, "TPC n sigma for deuteron PID"};
  o2::framework::Configurable<float> cfgTOFNsigma{"cfgTOFNsigma", 4.0f, "TOF n sigma for deuteron PID"};
  o2::framework::Configurable<float> cfgITSNsigma{"cfgITSNsigma", -2.0f, "ITS n sigma for deuteron PID"};
  o2::framework::Configurable<float> cfgEta{"cfgEta", 0.8f, "Track eta selection"};
  o2::framework::Configurable<int> cfgTPCNclsFound{"cfgTPCNclsFound", 100, "Minimum TPC clusters found"};
  o2::framework::Configurable<float> cfgTPCChi2Ncl{"cfgTPCChi2Ncl", 4.0f, "Maximum TPC chi2 per N clusters"};
  o2::framework::Configurable<float> cfgITSChi2Ncl{"cfgITSChi2Ncl", 36.0f, "Maximum ITS chi2 per N clusters"};
  o2::framework::Configurable<int> cfgITScls{"cfgITScls", 2, "Minimum ITS clusters"};
  o2::framework::Configurable<float> cfgMaxPt{"cfgMaxPt", 5.0f, "Maximum pT cut"};
  o2::framework::Configurable<float> cfgMinPt{"cfgMinPt", 0.5f, "Minimum pT cut"};
  o2::framework::Configurable<float> cfgDCAcut{"cfgDCAcut", 0.003f, "DCA cut for non prompt deuteron"};
  o2::framework::Configurable<float> ptThresholdPid{"ptThresholdPid", 1.0f, "pT threshold to switch between ITS and TOF PID"};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::framework::Service<o2::pid::tof::TOFResponse> tofResponse;

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch"); // Set CCDB URL to get magnetic field
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    tofResponse->initSetup(ccdb, initContext);

    o2::framework::AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};

    // general QA histograms
    QAHistos.add("hVtxZ", "Z-Vertex distribution after selection;Z (cm)", o2::framework::HistType::kTH1F, {{100, -50, 50}});
    QAHistos.add("hDCAxyVsPt", "DCAxy #bar{d} vs p_{T}", {o2::framework::HistType::kTH2D, {pTAxis, DCAxyAxis}});
    QAHistos.add("hDCAzVsPt", "DCAz #bar{d} vs p_{T}", {o2::framework::HistType::kTH2D, {pTAxis, DCAzAxis}});
    QAHistos.add("hnSigmaTPCVsPt", "n#sigma TPC vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TPC", {o2::framework::HistType::kTH2D, {pTAxis, nSigmaAxis}});
    QAHistos.add("hnSigmaTOFVsPt", "n#sigma TOF vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TOF", {o2::framework::HistType::kTH2D, {pTAxis, nSigmaAxis}});
    QAHistos.add("hnSigmaITSVsPt", "n#sigma ITS vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma ITS", {o2::framework::HistType::kTH2D, {pTAxis, nSigmaAxis}});
    QAHistos.add("ptAntiDeuteron", "ptAntiDeuteron", {o2::framework::HistType::kTH1F, {ptAxis}});
    QAHistos.add("etaAntideuteron", "etaAntideuteron", {o2::framework::HistType::kTH1F, {{100, -1.0f, 1.0f, "eta #bar{d}"}}});
    QAHistos.add("hDCAxyVsPt-pre_selection", "DCAxy #bar{d} vs p_{T}", {o2::framework::HistType::kTH2D, {pTAxis, DCAxyAxis}});
    QAHistos.add("hDCAzVsPt-pre-selection", "DCAz #bar{d} vs p_{T}", {o2::framework::HistType::kTH2D, {pTAxis, DCAzAxis}});

    // processed events
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "TVX");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "TF border");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "ITS-RO border");
    hProcessedEvents->GetXaxis()->SetBinLabel(5, "z-vertex");
    hProcessedEvents->GetXaxis()->SetBinLabel(6, o2::aod::filtering::H2fromLb::columnLabel());
  }

  // Tables
  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra,
                                        o2::aod::TracksDCA, o2::aod::TrackSelection,
                                        o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe,
                                        o2::aod::TOFSignal, o2::aod::TOFEvTime>;

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

  void fillTriggerTable(bool keepEvent[])
  {
    table(keepEvent[0]);
  }

  o2::framework::Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;
  int mCurrentRun = -1;

  void process(CollisionCandidates const& collisions,
               o2::aod::TrackAssoc const& trackIndices,
               TrackCandidates const& tracks,
               o2::aod::BCsWithTimestamps const& bcs)
  {
    tofResponse->processSetup(bcs.iteratorAt(0)); // Update the response parameters
    for (const auto& collision : collisions)      // start loop over collisions
    {
      if (mCurrentRun != collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber()) { // If the run is new then we need to initialize the propagator field
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", collision.bc_as<o2::aod::BCsWithTimestamps>().timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        mCurrentRun = collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber();
      }

      // Is event good? keepEvent[0] = non promp deuteron
      bool keepEvent[1]{}; // explicitly zero-initialised
      hProcessedEvents->Fill(0.5);
      // Event selection cuts
      if (isTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fillTriggerTable(keepEvent);
        continue;
      }
      hProcessedEvents->Fill(1.5);
      if (isNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fillTriggerTable(keepEvent);
        continue;
      }
      hProcessedEvents->Fill(2.5);
      if (isNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fillTriggerTable(keepEvent);
        continue;
      }
      hProcessedEvents->Fill(3.5);
      if (std::abs(collision.posZ()) > cutzvertex) {
        fillTriggerTable(keepEvent);
        continue;
      }
      hProcessedEvents->Fill(4.5);
      QAHistos.fill(HIST("hVtxZ"), collision.posZ());

      // Loop over tracks
      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      auto tracksWithItsPid = o2::soa::Attach<TrackCandidates, o2::aod::pidits::ITSNSigmaDe>(tracks);

      float tofEventTime = 0.f;
      float tofEventTimeErr = 0.f;

      for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks
        const auto& track = tracksWithItsPid.rawIteratorAt(trackId.trackId());
        if (!passedSingleTrackSelection(track)) {
          continue;
        }
        if (track.collisionId() != collision.globalIndex()) {
          continue;
        }
        tofEventTime = track.tofEvTime();
        tofEventTimeErr = track.tofEvTimeErr();
        break;
      }

      for (const auto& trackId : trackIdsThisCollision) { // start loop over tracks

        const auto& track = tracksWithItsPid.rawIteratorAt(trackId.trackId());

        std::array<float, 2> dca{track.dcaXY(), track.dcaZ()};
        std::array<float, 3> pVec = track.pVector();

        if (track.collisionId() != collision.globalIndex()) {
          auto trackPar = getTrackParCov(track);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dca);
          getPxPyPz(trackPar, pVec);
        }

        if (!passedSingleTrackSelection(track)) {
          continue;
        }
        float recalculatedNSigmaTOFDe = 0.f;
        if (track.collisionId() != collision.globalIndex()) {
          recalculatedNSigmaTOFDe = tofResponse->nSigma<o2::track::PID::Deuteron>(track.tofSignalInAnotherBC(track.collision_as<CollisionCandidates>().bc_as<o2::aod::BCsWithTimestamps>().globalBC(), collision.bc_as<o2::aod::BCsWithTimestamps>().globalBC()),
                                                                                  track.tofExpMom(), track.length(), track.p(), track.eta(), tofEventTime, tofEventTimeErr);
        }
        const bool isTOFDe = std::abs(track.tofNSigmaDe()) < cfgTOFNsigma;
        const bool isTPCDe = std::abs(track.tpcNSigmaDe()) < cfgTPCNsigma;
        const bool isITSDe = track.itsNSigmaDe() > cfgITSNsigma;

        if (track.pt() < ptThresholdPid) {
          if (isTPCDe && isITSDe) {
            QAHistos.fill(HIST("hDCAxyVsPt-pre_selection"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt-pre-selection"), track.pt(), dca[1]);
            if (std::abs(dca[0]) < cfgDCAcut) {
              continue;
            }
            keepEvent[0] = true;
            QAHistos.fill(HIST("ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
            QAHistos.fill(HIST("hnSigmaITSVsPt"), track.pt(), track.itsNSigmaDe());
          }
        } else {
          if (isTPCDe && isTOFDe) {

            QAHistos.fill(HIST("hDCAxyVsPt-pre_selection"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt-pre-selection"), track.pt(), dca[1]);
            if (std::abs(dca[0]) < cfgDCAcut) {
              continue;
            }
            keepEvent[0] = true;
            QAHistos.fill(HIST("ptAntiDeuteron"), track.pt());
            QAHistos.fill(HIST("etaAntideuteron"), track.eta());
            QAHistos.fill(HIST("hDCAxyVsPt"), track.pt(), dca[0]);
            QAHistos.fill(HIST("hDCAzVsPt"), track.pt(), dca[1]);
            QAHistos.fill(HIST("hnSigmaTPCVsPt"), track.pt(), track.tpcNSigmaDe());
            QAHistos.fill(HIST("hnSigmaTOFVsPt"), track.pt(), track.tofNSigmaDe());
            QAHistos.fill(HIST("hnSigmaITSVsPt"), track.pt(), track.itsNSigmaDe());
          }
        }
      } // end track loop
      if (keepEvent[0]) {
        hProcessedEvents->Fill(5.5);
      }

      fillTriggerTable(keepEvent);
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  return o2::framework::WorkflowSpec{o2::framework::adaptAnalysisTask<H2fromLbFilter>(cfgc)};
}
