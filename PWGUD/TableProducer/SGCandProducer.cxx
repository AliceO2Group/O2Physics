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

#include <cmath>
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/SGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SGCandProducer {
  // get an SGCutparHolder
  SGCutParHolder sameCuts = SGCutParHolder(); // SGCutparHolder
  Configurable<SGCutParHolder> SGCuts{"SGCuts", {}, "SG event cuts"};
  Configurable<bool> saveAllTracks{"saveAllTracks", true, "save only PV contributors or all tracks associated to a collision"};
  Configurable<bool> savenonPVCITSOnlyTracks{"savenonPVCITSOnlyTracks", false, "save non PV contributors with ITS only information"};
  Configurable<bool> rejectAtTFBoundary{"rejectAtTFBoundary", true, "reject collisions at a TF boundary"};
  Configurable<bool> noITSROFrameBorder{"noITSROFrameBorder", true, "reject ITS RO Frame Border"};
  Configurable<bool> noSameBunchPileUp{"noSameBunchPileUp", true, "reject SameBunchPileUp"};
  Configurable<bool> IsGoodVertex{"IsGoodVertex", false, "Select FT0 PV vertex matching"};
  Configurable<bool> ITSTPCVertex{"ITSTPCVertex", true, "reject ITS-only vertex"}; // if one wants to look at Single Gap pp events
  //  SG selector
  SGSelector sgSelector;

  // data tables
  Produces<aod::SGCollisions> outputSGCollisions;
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<o2::aod::UDZdcsReduced> udZdcsReduced;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDFwdTracks> outputFwdTracks;
  Produces<aod::UDFwdTracksExtra> outputFwdTracksExtra;
  Produces<aod::UDTracksLabels> outputTracksLabel;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::TOFSignal, aod::pidTOFbeta,
                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;

  // function to update UDFwdTracks, UDFwdTracksExtra
  template <typename TFwdTrack>
  void updateUDFwdTrackTables(TFwdTrack const& fwdtrack, uint64_t const& bcnum)
  {
    outputFwdTracks(outputCollisions.lastIndex(),
                    fwdtrack.px(), fwdtrack.py(), fwdtrack.pz(), fwdtrack.sign(),
                    bcnum, fwdtrack.trackTime(), fwdtrack.trackTimeRes());
    outputFwdTracksExtra(fwdtrack.nClusters(),
                         fwdtrack.pDca(),
                         fwdtrack.rAtAbsorberEnd(),
                         fwdtrack.chi2(),
                         fwdtrack.chi2MatchMCHMID(),
                         fwdtrack.mchBitMap(),
                         fwdtrack.midBitMap(),
                         fwdtrack.midBoards());
  }

  // function to update UDTracks, UDTracksCov, UDTracksDCA, UDTracksPID, UDTracksExtra, UDTracksFlag,
  // and UDTrackCollisionIDs
  template <typename TTrack>
  void updateUDTrackTables(int64_t lastIndex, TTrack const& track, uint64_t const& bcnum)
  {
    outputTracks(lastIndex,
                 track.px(), track.py(), track.pz(), track.sign(),
                 bcnum, track.trackTime(), track.trackTimeRes());

    // float sigmaY = track.sigmaY();
    // float sigmaZ = track.sigmaZ();
    float sigmaY = -1.;
    float sigmaZ = -1.;
    outputTracksCov(track.x(), track.y(), track.z(), sigmaY, sigmaZ);

    outputTracksDCA(track.dcaZ(), track.dcaXY());
    outputTracksPID(track.tpcNSigmaEl(),
                    track.tpcNSigmaMu(),
                    track.tpcNSigmaPi(),
                    track.tpcNSigmaKa(),
                    track.tpcNSigmaPr(),
                    track.beta(),
                    track.betaerror(),
                    track.tofNSigmaEl(),
                    track.tofNSigmaMu(),
                    track.tofNSigmaPi(),
                    track.tofNSigmaKa(),
                    track.tofNSigmaPr());
    outputTracksExtra(track.tpcInnerParam(),
                      track.itsClusterSizes(),
                      track.tpcNClsFindable(),
                      track.tpcNClsFindableMinusFound(),
                      track.tpcNClsFindableMinusCrossedRows(),
                      track.tpcNClsShared(),
                      track.trdPattern(),
                      track.itsChi2NCl(),
                      track.tpcChi2NCl(),
                      track.trdChi2(),
                      track.tofChi2(),
                      track.tpcSignal(),
                      track.tofSignal(),
                      track.trdSignal(),
                      track.length(),
                      track.tofExpMom(),
                      track.detectorMap());
    outputTracksFlag(track.has_collision(),
                     track.isPVContributor());
    outputTracksLabel(track.globalIndex());
  }

  void init(InitContext&)
  {
    sameCuts = (SGCutParHolder)SGCuts;
    registry.add("reco/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}});
  }

  // process function for real data
  void process(CC const& collision, BCs const& bcs, TCs& tracks, FWs& fwdtracks,
               aod::Zdcs& zdcs, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    LOGF(debug, "<SGCandProducer>  collision %d", collision.globalIndex());
    registry.get<TH1>(HIST("reco/Stat"))->Fill(0., 1.);
    // reject collisions at TF boundaries
    if (rejectAtTFBoundary && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    // reject collisions at ITS RO TF boundaries
    if (noITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    // registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // reject Same Bunch PileUp
    if (noSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    // registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // check vertex matching to FT0
    if (IsGoodVertex && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    // registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // reject ITS Only vertices
    if (ITSTPCVertex && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    // registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(2., 1.);
    auto bc = collision.foundBC_as<BCs>();
    auto newbc = bc;

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, sameCuts.NDtcoll(), bcs, sameCuts.minNBCs());
    auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, bc);
    // auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, tracks);
    int issgevent = isSGEvent.value;
    if (isSGEvent.bc) {
      newbc = *(isSGEvent.bc);
    } else {
      LOGF(info, "No Newbc %i", bc.globalBC());
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(issgevent + 3, 1.);
    if (issgevent <= 2) {
      //    LOGF(info, "Current BC: %i, %i, %i", bc.globalBC(), newbc.globalBC(), issgevent);
      if (sameCuts.minRgtrwTOF()) {
        if (udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib()) < sameCuts.minRgtrwTOF())
          return;
      }
      upchelpers::FITInfo fitInfo{};
      udhelpers::getFITinfo(fitInfo, newbc.globalBC(), bcs, ft0s, fv0as, fdds);
      // update SG candidates tables
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       collision.numContrib(), udhelpers::netCharge<true>(tracks),
                       1.); // rtrwTOF); //omit the calculation to speed up the things while skimming
      outputSGCollisions(issgevent);
      outputCollisionsSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C,
                           fitInfo.triggerMaskFT0,
                           fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC,
                           fitInfo.triggerMaskFDD,
                           fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                           fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                           fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                           fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      outputCollsLabels(collision.globalIndex());
      if (newbc.has_zdc()) {
        auto zdc = newbc.zdc();
        udZdcsReduced(outputCollisions.lastIndex(), zdc.timeZNA(), zdc.timeZNC(), zdc.energyCommonZNA(), zdc.energyCommonZNC());
      } else {
        udZdcsReduced(outputCollisions.lastIndex(), -999, -999, -999, -999);
      }
      // update SGTracks tables
      for (auto& track : tracks) {
        if (track.pt() > sameCuts.minPt() && track.eta() > sameCuts.minEta() && track.eta() < sameCuts.maxEta()) {
          if (track.isPVContributor()) {
            updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
          } else if (saveAllTracks) {
            if (track.itsClusterSizes() && track.itsChi2NCl() > 0 && ((track.tpcNClsFindable() == 0 && savenonPVCITSOnlyTracks) || track.tpcNClsFindable() > 50))
              updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
            // if (track.isPVContributor())  updateUDTrackTables(outputCollisions.lastIndex(), track, bc.globalBC());
          }
        }
      }
      // update SGFwdTracks tables
      if (sameCuts.withFwdTracks()) {
        for (auto& fwdtrack : fwdtracks) {
          if (!sgSelector.FwdTrkSelector(fwdtrack))
            updateUDFwdTrackTables(fwdtrack, bc.globalBC());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<SGCandProducer>(cfgc, TaskName{"sgcandproducer"})};

  return workflow;
}
