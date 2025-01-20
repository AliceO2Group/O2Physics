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
#include <vector>
#include <map>
#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/SGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct SGCandProducer {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // data inputs
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BC = BCs::iterator;
  using TCs = soa::Join<aod::Tracks, /*aod::TracksCov,*/ aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                        aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                        aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                        aod::TOFSignal, aod::pidTOFbeta,
                        aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                        aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using FWs = aod::FwdTracks;
  // get an SGCutparHolder
  SGCutParHolder sameCuts = SGCutParHolder(); // SGCutparHolder
  Configurable<SGCutParHolder> SGCuts{"SGCuts", {}, "SG event cuts"};
  Configurable<bool> verboseInfo{"verboseInfo", false, "Print general info to terminal; default it false."};
  Configurable<bool> saveAllTracks{"saveAllTracks", true, "save only PV contributors or all tracks associated to a collision"};
  Configurable<bool> savenonPVCITSOnlyTracks{"savenonPVCITSOnlyTracks", false, "save non PV contributors with ITS only information"};
  Configurable<bool> rejectAtTFBoundary{"rejectAtTFBoundary", true, "reject collisions at a TF boundary"};
  Configurable<bool> noITSROFrameBorder{"noITSROFrameBorder", true, "reject ITS RO Frame Border"};
  Configurable<bool> noSameBunchPileUp{"noSameBunchPileUp", true, "reject SameBunchPileUp"};
  Configurable<bool> IsGoodVertex{"IsGoodVertex", false, "Select FT0 PV vertex matching"};
  Configurable<bool> ITSTPCVertex{"ITSTPCVertex", true, "reject ITS-only vertex"}; // if one wants to look at Single Gap pp events

  // Configurables to decide which tables are filled
  Configurable<bool> fillTrackTables{"fillTrackTables", true, "Fill track tables"};
  Configurable<bool> fillFwdTrackTables{"fillFwdTrackTables", true, "Fill forward track tables"};
  //  SG selector
  SGSelector sgSelector;
  ctpRateFetcher mRateFetcher;

  // data tables
  Produces<aod::SGCollisions> outputSGCollisions;
  Produces<aod::UDCollisions> outputCollisions;
  Produces<aod::UDCollisionsSels> outputCollisionsSels;
  Produces<aod::UDCollisionSelExtras> outputCollisionSelExtras;
  Produces<aod::UDCollsLabels> outputCollsLabels;
  Produces<aod::UDZdcs> outputZdcs;
  Produces<aod::UDZdcsReduced> udZdcsReduced;
  Produces<aod::UDTracks> outputTracks;
  Produces<aod::UDTracksCov> outputTracksCov;
  Produces<aod::UDTracksDCA> outputTracksDCA;
  Produces<aod::UDTracksPID> outputTracksPID;
  Produces<aod::UDTracksPIDExtra> outputTracksPIDExtra;
  Produces<aod::UDTracksExtra> outputTracksExtra;
  Produces<aod::UDTracksFlags> outputTracksFlag;
  Produces<aod::UDFwdTracks> outputFwdTracks;
  Produces<aod::UDFwdTracksExtra> outputFwdTracksExtra;
  Produces<aod::UDTracksLabels> outputTracksLabel;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  // function to update UDFwdTracks, UDFwdTracksExtra
  template <typename TFwdTrack>
  void updateUDFwdTrackTables(TFwdTrack const& fwdtrack, uint64_t const& bcnum)
  {
    outputFwdTracks(outputCollisions.lastIndex(),
                    fwdtrack.px(), fwdtrack.py(), fwdtrack.pz(), fwdtrack.sign(),
                    bcnum, fwdtrack.trackTime(), fwdtrack.trackTimeRes());
    outputFwdTracksExtra(fwdtrack.trackType(),
                         fwdtrack.nClusters(),
                         fwdtrack.pDca(),
                         fwdtrack.rAtAbsorberEnd(),
                         fwdtrack.chi2(),
                         fwdtrack.chi2MatchMCHMID(),
                         fwdtrack.chi2MatchMCHMFT(),
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
    outputTracksPIDExtra(track.tpcNSigmaDe(),
                         track.tpcNSigmaTr(),
                         track.tpcNSigmaHe(),
                         track.tpcNSigmaAl(),
                         track.tofNSigmaDe(),
                         track.tofNSigmaTr(),
                         track.tofNSigmaHe(),
                         track.tofNSigmaAl());
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
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
    sameCuts = (SGCutParHolder)SGCuts;
    registry.add("reco/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{14, -0.5, 13.5}}});
  }

  // process function for real data
  void process(CC const& collision, BCs const& bcs, TCs& tracks, FWs& fwdtracks,
               aod::Zdcs& /*zdcs*/, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    if (verboseInfo)
      LOGF(debug, "<SGCandProducer>  collision %d", collision.globalIndex());
    registry.get<TH1>(HIST("reco/Stat"))->Fill(0., 1.);
    // reject collisions at TF boundaries
    if (rejectAtTFBoundary && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(1., 1.);
    // reject collisions at ITS RO TF boundaries
    if (noITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(2., 1.);
    // reject Same Bunch PileUp
    if (noSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(3., 1.);
    // check vertex matching to FT0
    if (IsGoodVertex && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(4., 1.);
    // reject ITS Only vertices
    if (ITSTPCVertex && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(5., 1.);
    // nominal BC
    if (!collision.has_foundBC()) {
      return;
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(6., 1.);
    int trs = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      trs = 1;
    }
    int trofs = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      trofs = 1;
    }
    int hmpr = 0;
    if (collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      hmpr = 1;
    }
    auto bc = collision.foundBC_as<BCs>();
    double ir = 0.;
    const uint64_t ts = bc.timestamp();
    const int runnumber = bc.runNumber();
    if (bc.has_zdc()) {
      ir = mRateFetcher.fetch(ccdb.service, ts, runnumber, "ZNC hadronic") * 1.e-3;
    }
    auto newbc = bc;

    // obtain slice of compatible BCs
    auto bcRange = udhelpers::compatibleBCs(collision, sameCuts.NDtcoll(), bcs, sameCuts.minNBCs());
    auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, bc);
    // auto isSGEvent = sgSelector.IsSelected(sameCuts, collision, bcRange, tracks);
    int issgevent = isSGEvent.value;
    if (isSGEvent.bc && issgevent < 2) {
      newbc = *(isSGEvent.bc);
    } else {
      if (verboseInfo)
        LOGF(info, "No Newbc %i", bc.globalBC());
    }
    registry.get<TH1>(HIST("reco/Stat"))->Fill(issgevent + 8, 1.);
    if (issgevent <= 2) {
      if (verboseInfo)
        LOGF(info, "Current BC: %i, %i, %i", bc.globalBC(), newbc.globalBC(), issgevent);
      if (sameCuts.minRgtrwTOF()) {
        if (udhelpers::rPVtrwTOF<true>(tracks, collision.numContrib()) < sameCuts.minRgtrwTOF())
          return;
      }
      upchelpers::FITInfo fitInfo{};
      uint8_t chFT0A = 0;
      uint8_t chFT0C = 0;
      uint8_t chFDDA = 0;
      uint8_t chFDDC = 0;
      uint8_t chFV0A = 0;
      int occ = 0;
      occ = collision.trackOccupancyInTimeRange();
      udhelpers::getFITinfo(fitInfo, newbc, bcs, ft0s, fv0as, fdds);
      // update SG candidates tables
      int upc_flag = 0;
      ushort flags = collision.flags();
      if (flags & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode)
        upc_flag = 1;
      outputCollisions(bc.globalBC(), bc.runNumber(),
                       collision.posX(), collision.posY(), collision.posZ(), upc_flag,
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
      outputCollisionSelExtras(chFT0A, chFT0C, chFDDA, chFDDC, chFV0A, occ, ir, trs, trofs, hmpr);
      outputCollsLabels(collision.globalIndex());
      if (newbc.has_zdc()) {
        auto zdc = newbc.zdc();
        udZdcsReduced(outputCollisions.lastIndex(), zdc.timeZNA(), zdc.timeZNC(), zdc.energyCommonZNA(), zdc.energyCommonZNC());
      } else {
        udZdcsReduced(outputCollisions.lastIndex(), -999, -999, -999, -999);
      }
      // update SGTracks tables
      if (fillTrackTables) {
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
      }
      // update SGFwdTracks tables
      if (fillFwdTrackTables) {
        if (sameCuts.withFwdTracks()) {
          for (auto& fwdtrack : fwdtracks) {
            if (!sgSelector.FwdTrkSelector(fwdtrack))
              updateUDFwdTrackTables(fwdtrack, bc.globalBC());
          }
        }
      }
    }
  }
};

struct McSGCandProducer {
  // MC tables
  Configurable<bool> verboseInfoMC{"verboseInfoMC", false, "Print general info to terminal; default it false."};
  Produces<aod::UDMcCollisions> outputMcCollisions;
  Produces<aod::UDMcParticles> outputMcParticles;
  Produces<aod::UDMcCollsLabels> outputMcCollsLabels;
  Produces<aod::UDMcTrackLabels> outputMcTrackLabels;

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::McTrackLabels>;
  using UDCCs = soa::Join<aod::UDCollisions, aod::UDCollsLabels, aod::SGCollisions>;
  using UDTCs = soa::Join<aod::UDTracks, aod::UDTracksLabels>;

  // prepare slices
  SliceCache cache;
  PresliceUnsorted<aod::McParticles> mcPartsPerMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<UDTCs> udtracksPerUDCollision = aod::udtrack::udCollisionId;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  template <typename TMcCollision>
  void updateUDMcCollisions(TMcCollision const& mccol, uint64_t globBC)
  {
    // save mccol
    outputMcCollisions(globBC,
                       mccol.generatorsID(),
                       mccol.posX(),
                       mccol.posY(),
                       mccol.posZ(),
                       mccol.t(),
                       mccol.weight(),
                       mccol.impactParameter());
  }

  template <typename TMcParticle>
  void updateUDMcParticle(TMcParticle const& McPart, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // save McPart
    // mother and daughter indices are set to -1
    // ATTENTION: this can be improved to also include mother and daughter indices
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};

    // update UDMcParticles
    if (mcPartIsSaved.find(McPart.globalIndex()) == mcPartIsSaved.end()) {
      outputMcParticles(McCollisionId,
                        McPart.pdgCode(),
                        McPart.statusCode(),
                        McPart.flags(),
                        newmids,
                        newdids,
                        McPart.weight(),
                        McPart.px(),
                        McPart.py(),
                        McPart.pz(),
                        McPart.e());
      mcPartIsSaved[McPart.globalIndex()] = outputMcParticles.lastIndex();
    }
  }

  template <typename TMcParticles>
  void updateUDMcParticles(TMcParticles const& McParts, int64_t McCollisionId, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // save McParts
    // new mother and daughter ids
    std::vector<int32_t> newmids;
    int32_t newdids[2] = {-1, -1};
    int64_t newval = -1;

    // Determine the particle indices within the UDMcParticles table
    // before filling the table
    // This is needed to be able to assign the new daughter indices
    std::map<int64_t, int64_t> oldnew;
    auto lastId = outputMcParticles.lastIndex();
    for (auto mcpart : McParts) {
      auto oldId = mcpart.globalIndex();
      if (mcPartIsSaved.find(oldId) != mcPartIsSaved.end()) {
        oldnew[oldId] = mcPartIsSaved[oldId];
      } else {
        lastId++;
        oldnew[oldId] = lastId;
      }
    }

    // all particles of the McCollision are saved
    for (auto mcpart : McParts) {
      if (mcPartIsSaved.find(mcpart.globalIndex()) == mcPartIsSaved.end()) {
        // mothers
        newmids.clear();
        auto oldmids = mcpart.mothersIds();
        for (auto oldmid : oldmids) {
          auto m = McParts.rawIteratorAt(oldmid);
          if (verboseInfoMC)
            LOGF(debug, "    m %d", m.globalIndex());
          if (mcPartIsSaved.find(oldmid) != mcPartIsSaved.end()) {
            newval = mcPartIsSaved[oldmid];
          } else {
            newval = -1;
          }
          newmids.push_back(newval);
        }
        // daughters
        auto olddids = mcpart.daughtersIds();
        for (uint ii = 0; ii < olddids.size(); ii++) {
          if (oldnew.find(olddids[ii]) != oldnew.end()) {
            newval = oldnew[olddids[ii]];
          } else {
            newval = -1;
          }
          newdids[ii] = newval;
        }
        if (verboseInfoMC)
          LOGF(debug, " ms %i ds %i", oldmids.size(), olddids.size());

        // update UDMcParticles
        outputMcParticles(McCollisionId,
                          mcpart.pdgCode(),
                          mcpart.statusCode(),
                          mcpart.flags(),
                          newmids,
                          newdids,
                          mcpart.weight(),
                          mcpart.px(),
                          mcpart.py(),
                          mcpart.pz(),
                          mcpart.e());
        mcPartIsSaved[mcpart.globalIndex()] = outputMcParticles.lastIndex();
      }
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabel(TTrack const& udtrack, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
    auto trackId = udtrack.trackId();
    if (trackId >= 0) {
      auto track = udtrack.template track_as<TCs>();
      auto mcTrackId = track.mcParticleId();
      if (mcTrackId >= 0) {
        if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
          outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, track.mcMask());
      }
    } else {
      outputMcTrackLabels(-1, -1);
    }
  }

  template <typename TTrack>
  void updateUDMcTrackLabels(TTrack const& udtracks, std::map<int64_t, int64_t>& mcPartIsSaved)
  {
    // loop over all tracks
    for (auto udtrack : udtracks) {
      // udtrack (UDTCs) -> track (TCs) -> mcTrack (McParticles) -> udMcTrack (UDMcParticles)
      auto trackId = udtrack.trackId();
      if (trackId >= 0) {
        auto track = udtrack.template track_as<TCs>();
        auto mcTrackId = track.mcParticleId();
        if (mcTrackId >= 0) {
          if (mcPartIsSaved.find(mcTrackId) != mcPartIsSaved.end()) {
            outputMcTrackLabels(mcPartIsSaved[mcTrackId], track.mcMask());
          } else {
            outputMcTrackLabels(-1, track.mcMask());
          }
        } else {
          outputMcTrackLabels(-1, track.mcMask());
        }
      } else {
        outputMcTrackLabels(-1, -1);
      }
    }
  }

  void init(InitContext& context)
  {
    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMC")) {
      registry.add("mcTruth/collisions", "Number of associated collisions", {HistType::kTH1F, {{11, -0.5, 10.5}}});
      registry.add("mcTruth/collType", "Collision type", {HistType::kTH1F, {{5, -0.5, 4.5}}});
      registry.add("mcTruth/IVMpt", "Invariant mass versus p_{T}", {HistType::kTH2F, {{150, 0.0, 3.0}, {150, 0.0, 3.0}}});
    }
  }

  // process function for MC data
  // save the MC truth of all events of interest and of the DG events
  void processMC(aod::McCollisions const& mccols, aod::McParticles const& mcparts,
                 UDCCs const& sgcands, UDTCs const& udtracks,
                 CCs const& /*collisions*/, BCs const& /*bcs*/, TCs const& /*tracks*/)
  {
    if (verboseInfoMC) {
      LOGF(info, "Number of McCollisions %d", mccols.size());
      LOGF(info, "Number of SG candidates %d", sgcands.size());
      LOGF(info, "Number of UD tracks %d", udtracks.size());
    }
    if (sgcands.size() <= 0) {
      if (verboseInfoMC)
        LOGF(info, "No DG candidates to save!");
      return;
    }

    // use a hash table to keep track of the McCollisions which have been added to the UDMcCollision table
    // {McCollisionId : udMcCollisionId}
    // similar for the McParticles which have been added to the UDMcParticle table
    // {McParticleId : udMcParticleId}
    std::map<int64_t, int64_t> mcColIsSaved;
    std::map<int64_t, int64_t> mcPartIsSaved;

    // loop over McCollisions and UDCCs simultaneously
    auto mccol = mccols.iteratorAt(0);
    auto sgcand = sgcands.iteratorAt(0);
    auto lastmccol = mccols.iteratorAt(mccols.size() - 1);
    auto lastsgcand = sgcands.iteratorAt(sgcands.size() - 1);

    // advance dgcand and mccol until both are AtEnd
    int64_t mccolId = mccol.globalIndex();
    int64_t mcsgId = -1;
    // int64_t colId = -1;
    auto sgcandAtEnd = sgcand == lastsgcand;
    auto mccolAtEnd = mccol == lastmccol;
    bool goon = !sgcandAtEnd || !mccolAtEnd;
    while (goon) {
      auto bcIter = mccol.bc_as<BCs>();
      uint64_t globBC = bcIter.globalBC();
      // uint64_t globBC = 0;
      // check if dgcand has an associated McCollision
      if (sgcand.has_collision()) {
        auto sgcandCol = sgcand.collision_as<CCs>();
        // colId = sgcandCol.globalIndex();
        if (sgcandCol.has_mcCollision()) {
          mcsgId = sgcandCol.mcCollision().globalIndex();
        } else {
          mcsgId = -1;
        }
      } else {
        //  colId = -1;
        mcsgId = -1;
      }
      if (verboseInfoMC)
        LOGF(info, "\nStart of loop mcsgId %d mccolId %d", mcsgId, mccolId);

      // two cases to consider
      // 1. mcdgId <= mccolId: the event to process is a dgcand. In this case the Mc tables as well as the McLabel tables are updated
      // 2. mccolId < mcdgId: the event to process is an MC event of interest without reconstructed dgcand. In this case only the Mc tables are updated
      if ((!sgcandAtEnd && !mccolAtEnd && (mcsgId <= mccolId)) || mccolAtEnd) {
        // this is case 1.
        if (verboseInfoMC)
          LOGF(info, "Doing case 1 with mcsgId %d", mcsgId);

        // update UDMcCollisions and UDMcColsLabels (for each UDCollision -> UDMcCollisions)
        // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
        // get dgcand tracks
        auto sgTracks = udtracks.sliceByCached(aod::udtrack::udCollisionId, sgcand.globalIndex(), cache);

        // If the sgcand has an associated McCollision then the McCollision and all associated
        // McParticles are saved
        if (mcsgId >= 0) {
          if (mcColIsSaved.find(mcsgId) == mcColIsSaved.end()) {
            if (verboseInfoMC)
              LOGF(info, "  Saving McCollision %d", mcsgId);
            // update UDMcCollisions
            auto sgcandMcCol = sgcand.collision_as<CCs>().mcCollision();
            updateUDMcCollisions(sgcandMcCol, globBC);
            mcColIsSaved[mcsgId] = outputMcCollisions.lastIndex();
          }

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(mcColIsSaved[mcsgId]);

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mcsgId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mcsgId], mcPartIsSaved);

          // update UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          updateUDMcTrackLabels(sgTracks, mcPartIsSaved);

        } else {
          // If the sgcand has no associated McCollision then only the McParticles which are associated
          // with the tracks of the sgcand are saved
          if (verboseInfoMC)
            LOGF(info, "  Saving McCollision %d", -1);

          // update UDMcColsLabels (for each UDCollision -> UDMcCollisions)
          outputMcCollsLabels(-1);

          // update UDMcParticles and UDMcTrackLabels (for each UDTrack -> UDMcParticles)
          // loop over tracks of dgcand
          for (auto sgtrack : sgTracks) {
            if (sgtrack.has_track()) {
              auto track = sgtrack.track_as<TCs>();
              if (track.has_mcParticle()) {
                auto mcPart = track.mcParticle();
                updateUDMcParticle(mcPart, -1, mcPartIsSaved);
                updateUDMcTrackLabel(sgtrack, mcPartIsSaved);
              } else {
                outputMcTrackLabels(-1, track.mcMask());
              }
            } else {
              outputMcTrackLabels(-1, -1);
            }
          }
        }
        // advance sgcand
        if (sgcand != lastsgcand) {
          sgcand++;
        } else {
          sgcandAtEnd = true;
        }
      } else {
        // this is case 2.
        if (verboseInfoMC)
          LOGF(info, "Doing case 2");

        // update UDMcCollisions and UDMcParticles
        if (mcColIsSaved.find(mccolId) == mcColIsSaved.end()) {
          if (verboseInfoMC)
            LOGF(info, "  Saving McCollision %d", mccolId);
          // update UDMcCollisions
          updateUDMcCollisions(mccol, globBC);
          mcColIsSaved[mccolId] = outputMcCollisions.lastIndex();

          // update UDMcParticles
          auto mcPartsSlice = mcparts.sliceBy(mcPartsPerMcCollision, mccolId);
          updateUDMcParticles(mcPartsSlice, mcColIsSaved[mccolId], mcPartIsSaved);
        }

        // advance mccol
        if (mccol != lastmccol) {
          mccol++;
          mccolId = mccol.globalIndex();
        } else {
          mccolAtEnd = true;
        }
      }

      goon = !sgcandAtEnd || !mccolAtEnd;
      if (verboseInfoMC)
        LOGF(info, "End of loop mcsgId %d mccolId %d", mcsgId, mccolId);
    }
  }
  PROCESS_SWITCH(McSGCandProducer, processMC, "Produce MC tables", false);
  void processDummy(aod::Collisions const& /*collisions*/)
  {
    // do nothing
    if (verboseInfoMC)
      LOGF(info, "Running dummy process function!");
  }
  PROCESS_SWITCH(McSGCandProducer, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<SGCandProducer>(cfgc, TaskName{"sgcandproducer"}),
    adaptAnalysisTask<McSGCandProducer>(cfgc, TaskName{"mcsgcandproducer"})};
  return workflow;
}
