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
// Integration tester for quick and dirty cross checks that
// the framework is working reasonably
//
// Includes further QA if option enabled

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

#include <RtypesCore.h>

#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include <Framework/runDataProcessing.h>

using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

struct integrationTest {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsTracks{"nBinsTracks", 500, "number of bins in track histo"};
  Configurable<int> MaxNTrack{"MaxNTrack", 500, "NTrack max"};
  Configurable<int> nBinsTracks2D{"nBinsTracks2D", 200, "number of bins in 2D Ntrack corerelation plots"};
  Configurable<int> MaxNTrack2D{"MaxNTrack2D", 200, "NTrack max"};
  Configurable<int> nBinsCollisions{"nBinsCollisions", 50, "number of bins in collision histo"};
  Configurable<bool> do2DNTrackCorr{"do2DNTrackCorr", true, "Do 2D Ntrack correlation plots"};

  Configurable<bool> doBasicQA{"doBasicQA", true, "Do basic QA"};
  Configurable<bool> checkThinning{"checkThinning", true, "Check thinning-related quantities (retaining of TPC-only tracks)"};
  ConfigurableAxis axisHasDetector{"axisHasDetector", {16, -0.5f, 15.5f}, ""};
  ConfigurableAxis axisEta{"axisEta", {200, -2.0f, 2.0f}, ""};
  ConfigurableAxis axisPhi{"axisPhi", {200, 0.0f, +2 * TMath::Pi()}, ""};
  ConfigurableAxis axisNclu{"axisNclu", {10, -0.5f, 9.5f}, ""};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for analysis"};

  enum kTrackType { kTrackUsedV0 = 0,
                    kTrackUsedCascade, 
                    kTrackUsedDecay3Body, 
                    kTrackUsedTrackedV0,
                    kTrackUsedTrackedCascade,
                    kTrackUsedQA}; 

  enum kTable { kBC = 0,
                kBCFlag,
                kCollision,
                kOrigin,
                kTrack,
                kTrackCov,
                kTrackExtra,
                kMftTrack,
                kMftTrackCov,
                kFwdTrack,
                kFwdTrackCov,
                kFwdTrackCl,
                kAmbiguousTrack,
                kAmbiguousMftTrack,
                kAmbiguousFwdTrack,
                kTracked3Body, 
                kTrackedCascade,
                kTrackedV0,
                kTrackQA,
                kV0,
                kCascade,
                kDecay3Body,
                kCalo,
                kCaloTrigger,
                kCpvCluster,
                kFDD,
                kFDDExtra,
                kFT0,
                kFT0Extra,
                kFV0A,
                kFV0AExtra,
                kZDC,
                kHMPID,
                kMcCollision,
                kMcCollisionLabel,
                kMcParticle,
                kMcTrackLabel,
                kMcMftTrackLabel,
                kMcFwdTrackLabel,
                kFrameCounter,
                kNTables
  };

  void init(InitContext&)
  {
    TString lTableNames[] =
      {
        "bc",
        "bcflag",
        "collision",
        "origin",
        "track_iu",
        "trackcov_iu",
        "trackextra",
        "mfttrack",
        "mfttrackcov",
        "fwdtrack",
        "fwdtrackcov",
        "fwdtrkcl",
        "ambiguoustrack",
        "ambiguousmfttrack",
        "ambiguousfwdtrack",
        "tracked3body", 
        "trackedcascade", 
        "trackedv0",
        "trackqa",
        "v0",
        "cascade",
        "decay3body",
        "calo",
        "calotrigger",
        "cpvcluster",
        "fdd",
        "fddextra",
        "ft0",
        "ft0extra", 
        "fv0a",
        "fv0aextra",
        "zdc",
        "hmpid",
        "mccollision",
        "mccollisionlabel",
        "mcparticle",
        "mctracklabel",
        "mcmfttracklabel",
        "mcfwdtracklabel",
        "Total count",
        "" // empty (last)
      };
    const AxisSpec axisTables{40, 0.0f, 40.0f, ""};
    const AxisSpec axisTracks{nBinsTracks, -0.5f, MaxNTrack - 0.5f, "N_{tracks}"};
    const AxisSpec axisCollisions{nBinsCollisions, -0.5f, nBinsCollisions - 0.5f, "N_{collisions}"};
    // Label correctly to avoid confusion
    const AxisSpec axisTracks2Dtof{nBinsTracks2D, -0.5f, MaxNTrack2D - 0.5f, "N_{tracks}^{TOF}"};
    const AxisSpec axisTracks2Dtrd{nBinsTracks2D, -0.5f, MaxNTrack2D - 0.5f, "N_{tracks}^{TRD}"};
    const AxisSpec axisTracks2Dtpc{nBinsTracks2D, -0.5f, MaxNTrack2D - 0.5f, "N_{tracks}^{TPC}"};
    const AxisSpec axisTracks2Dits{nBinsTracks2D, -0.5f, MaxNTrack2D - 0.5f, "N_{tracks}^{ITS}"};

    // Table size bookkeeping / informational purposes
    auto hs = histos.add<TH1>("hTableSizes", "hTableSizes", HistType::kTH1D, {axisTables});
    for (Int_t ii = 0; ii < kNTables; ii++)
      hs->GetXaxis()->SetBinLabel(ii + 1, lTableNames[ii].Data());

    // Histograms with related indices: per BC
    histos.add<TH1>("hCollisionsPerBC", "hCollisionsPerBC", HistType::kTH1D, {axisCollisions});

    // Histograms with related indices: per collision
    histos.add<TH1>("hTracks", "hTracks", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hTracksITS", "hTracksITS", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hTracksTPC", "hTracksTPC", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hTracksNoTPCOnly", "hTracksNoTPCOnly", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hTracksTRD", "hTracksTRD", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hTracksTOF", "hTracksTOF", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hMFTTracks", "hMFTTracks", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hFWDTracks", "hFWDTracks", HistType::kTH1D, {axisTracks});

    histos.add<TH1>("hV0s", "hV0s", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hCascades", "hCascades", HistType::kTH1D, {axisTracks});
    histos.add<TH1>("hHasDetector", "hHasDetector", HistType::kTH1D, {axisHasDetector});

    if (do2DNTrackCorr) {
      histos.add<TH2>("hTOFvsTRD", "hTOFvsTRD", HistType::kTH2D, {axisTracks2Dtrd, axisTracks2Dtof});
      histos.add<TH2>("hTOFvsTPC", "hTOFvsTPC", HistType::kTH2D, {axisTracks2Dtpc, axisTracks2Dtof});
      histos.add<TH2>("hTRDvsTPC", "hTRDvsTPC", HistType::kTH2D, {axisTracks2Dtpc, axisTracks2Dtrd});
      histos.add<TH2>("hITSvsTPC", "hITSvsTPC", HistType::kTH2D, {axisTracks2Dtpc, axisTracks2Dits});
    }

    // pT histograms
    histos.add<TH1>("hPt", "hPt", HistType::kTH1D, {axisPt});
    histos.add<TH1>("hPtITS", "hPtITS", HistType::kTH1D, {axisPt});
    histos.add<TH1>("hPtTPC", "hPtTPC", HistType::kTH1D, {axisPt});
    histos.add<TH1>("hPtNoTPCOnly", "hPtNoTPCOnly", HistType::kTH1D, {axisPt});
    histos.add<TH1>("hPtTOF", "hPtTOF", HistType::kTH1D, {axisPt});
    histos.add<TH1>("hPtTRD", "hPtTRD", HistType::kTH1D, {axisPt});

    if (doBasicQA) {
      // general QA
      histos.add<TH2>("h2dPhiVsEtaAll", "h2dPhiVsEtaAll", HistType::kTH2D, {axisEta, axisPhi});
      histos.add<TH2>("h2dPhiVsEtaNoTPCOnly", "h2dPhiVsEtaNoTPCOnly", HistType::kTH2D, {axisEta, axisPhi});
      histos.add<TH1>("hNCluAll", "hNCluAll", HistType::kTH1D, {axisNclu});
      histos.add<TH1>("hNCluNoTPCOnly", "hNCluNoTPCOnly", HistType::kTH1D, {axisNclu});
    }
    if (checkThinning) {
      auto hThinningQA = histos.add<TH2>("hThinningQA", "hThinningQA", HistType::kTH2D, {{64, -0.5f, 63.5f}, {2, -0.5f, 1.5f}});

      TString kTrackTypeNames[] =
      {
        "V",
        "tV",
        "3b", 
        "C",
        "tC",
        "QA"
      };

      // construct labels to make this plot easier to understand
      for(uint8_t i=0; i<64; i++){ 
         TString trackTypeString = ""; 
         for(uint8_t j=0; j<kTrackUsedQA+1; j++){
          if( ((i) & (1<<(j))) ){
            trackTypeString.Append(Form("%s ", kTrackTypeNames[j].Data()));
          }
         }
         hThinningQA->GetXaxis()->SetBinLabel(i+1, trackTypeString.Data());
      } 
      hThinningQA->GetYaxis()->SetBinLabel(0, "Not TPConly");
      hThinningQA->GetYaxis()->SetBinLabel(1, "TPConly");
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Standard sizes (uncorrelated)
  void processDataModel(
    aod::BCs const& bcs,
    aod::Collisions const& collisions,
    aod::TracksIU const& tracks,
    aod::TracksCovIU const& trackcovs,
    aod::TracksExtra const& trackextras,
    aod::MFTTracks const& mfttracks,
    aod::FwdTracks const& fwdtracks,
    aod::FwdTracksCov const& fwdtrackcovs,
    aod::AmbiguousTracks const& ambitracks,
    aod::AmbiguousMFTTracks const& ambimfttracks,
    aod::AmbiguousFwdTracks const& ambifwdtracks,
    aod::V0s const& v0s,
    aod::Cascades const& cascades,
    aod::Calos const& calos,
    aod::CaloTriggers const& calotriggers,
    aod::FDDs const& fdds,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    aod::Zdcs const& zdcs)
  {
    histos.fill(HIST("hTableSizes"), static_cast<float>(kBC) + 0.5f, bcs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCollision) + 0.5f, collisions.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrack) + 0.5f, tracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackCov) + 0.5f, trackcovs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackExtra) + 0.5f, trackextras.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMftTrack) + 0.5f, mfttracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFwdTrack) + 0.5f, fwdtracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFwdTrackCov) + 0.5f, fwdtrackcovs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousTrack) + 0.5f, ambitracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousMftTrack) + 0.5f, ambimfttracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousFwdTrack) + 0.5f, ambifwdtracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kV0) + 0.5f, v0s.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCascade) + 0.5f, cascades.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCalo) + 0.5f, calos.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCaloTrigger) + 0.5f, calotriggers.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFDD) + 0.5f, fdds.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFT0) + 0.5f, ft0s.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFV0A) + 0.5f, fv0as.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kZDC) + 0.5f, zdcs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFrameCounter) + 0.5f);
  }
  PROCESS_SWITCH(integrationTest, processDataModel, "Check data model", true);

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  // Revised tables for newer data 
  // Warning: will not work for legacy AO2Ds
  void processDataModel2026(
    aod::BCs const& bcs,
    aod::BCFlags const& bcflags,
    aod::Collisions const& collisions,
    aod::Origins const& origins,
    aod::TracksIU const& tracks,
    aod::TracksCovIU const& trackcovs,
    aod::TracksExtra const& trackextras,
    aod::MFTTracks const& mfttracks,
    aod::MFTTracksCov const& mfttrackcovs,
    aod::FwdTracks const& fwdtracks,
    aod::FwdTracksCov const& fwdtrackcovs,
    aod::FwdTrkCls const& fwdtrkcls,
    aod::AmbiguousTracks const& ambitracks,
    aod::AmbiguousMFTTracks const& ambimfttracks,
    aod::AmbiguousFwdTracks const& ambifwdtracks,
    aod::Tracked3Bodys const& tracked3bodys, 
    aod::TrackedCascades const& trackedcascades,
    aod::TrackedV0s const& trackedv0s,
    aod::TracksQAVersion const& tracksQA, // does not comply with normal iterator naming
    aod::V0s const& v0s,
    aod::Cascades const& cascades,
    aod::Decay3Bodys const& decay3bodys, 
    aod::Calos const& calos,
    aod::CaloTriggers const& calotriggers,
    aod::CPVClusters const& cpvclusters, 
    aod::FDDs const& fdds,
    aod::FDDsExtra const& fddsextra,
    aod::FT0s const& ft0s,
    aod::FT0sExtra const& ft0sextra,
    aod::FV0As const& fv0as,
    aod::FV0AsExtra const& fv0asextra,
    aod::Zdcs const& zdcs,
    aod::HMPIDs const& hmpids)
  {
    histos.fill(HIST("hTableSizes"), static_cast<float>(kBC) + 0.5f, bcs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kBCFlag) + 0.5f, bcflags.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCollision) + 0.5f, collisions.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kOrigin) + 0.5f, origins.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrack) + 0.5f, tracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackCov) + 0.5f, trackcovs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackExtra) + 0.5f, trackextras.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMftTrack) + 0.5f, mfttracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMftTrackCov) + 0.5f, mfttrackcovs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFwdTrack) + 0.5f, fwdtracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFwdTrackCov) + 0.5f, fwdtrackcovs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFwdTrackCl) + 0.5f, fwdtrkcls.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousTrack) + 0.5f, ambitracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousMftTrack) + 0.5f, ambimfttracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kAmbiguousFwdTrack) + 0.5f, ambifwdtracks.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTracked3Body) + 0.5f, tracked3bodys.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackedCascade) + 0.5f, trackedcascades.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackedV0) + 0.5f, trackedv0s.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kTrackQA) + 0.5f, tracksQA.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kV0) + 0.5f, v0s.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCascade) + 0.5f, cascades.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kDecay3Body) + 0.5f, decay3bodys.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCalo) + 0.5f, calos.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCaloTrigger) + 0.5f, calotriggers.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kCpvCluster) + 0.5f, cpvclusters.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFDD) + 0.5f, fdds.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFDDExtra) + 0.5f, fddsextra.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFT0) + 0.5f, ft0s.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFT0Extra) + 0.5f, ft0sextra.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFV0A) + 0.5f, fv0as.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFV0AExtra) + 0.5f, fv0asextra.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kZDC) + 0.5f, zdcs.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kHMPID) + 0.5f, hmpids.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kFrameCounter) + 0.5f);

    if (checkThinning) {
      // perform checks regarding thinning
      std::vector<uint8_t> trackBits(tracks.size(), 0); 

      if(tracks.size() != trackextras.size()){ 
        LOGF(fatal, "Check thinning requested and Tracks.size() != TracksExtra.size(). Very bad. Quitting...");
      }

      // mark tracks that are used by each kind of composite object
      for (const auto& v0 : v0s) {
        trackBits[v0.negTrackId()] = trackBits[v0.negTrackId()] | (1 << kTrackUsedV0);
        trackBits[v0.posTrackId()] = trackBits[v0.posTrackId()] | (1 << kTrackUsedV0);
      }
      for (const auto& cascade : cascades) {
        auto v0 = cascade.v0(); // simple de-ref, no table joining
        trackBits[v0.negTrackId()] = trackBits[v0.negTrackId()] | (1 << kTrackUsedCascade);
        trackBits[v0.posTrackId()] = trackBits[v0.posTrackId()] | (1 << kTrackUsedCascade);
        trackBits[cascade.bachelorId()] = trackBits[cascade.bachelorId()] | (1 << kTrackUsedCascade);
      }
      for (const auto& trackedv0 : trackedv0s) {
        auto v0 = trackedv0.v0(); // simple de-ref, no table joining
        trackBits[v0.negTrackId()] = trackBits[v0.negTrackId()] | (1 << kTrackUsedTrackedV0);
        trackBits[v0.posTrackId()] = trackBits[v0.posTrackId()] | (1 << kTrackUsedTrackedV0);
      }
      for (const auto& trackedcascade : trackedcascades) {
        auto cascade = trackedcascade.cascade(); // simple de-ref, no table joining
        auto v0 = cascade.v0(); // simple de-ref, no table joining
        trackBits[v0.negTrackId()] = trackBits[v0.negTrackId()] | (1 << kTrackUsedTrackedCascade);
        trackBits[v0.posTrackId()] = trackBits[v0.posTrackId()] | (1 << kTrackUsedTrackedCascade);
        trackBits[cascade.bachelorId()] = trackBits[cascade.bachelorId()] | (1 << kTrackUsedTrackedCascade);
      }
      for (const auto& trackQA : tracksQA) {
        trackBits[trackQA.trackId()] = trackBits[trackQA.trackId()] | (1 << kTrackUsedQA);
      }
      for (uint32_t iTrack = 0; iTrack < trackextras.size(); iTrack++) { 
        auto track = trackextras.rawIteratorAt(iTrack); 
        bool isTPConly = track.hasTPC() && !track.hasTOF() && !track.hasTRD() && !track.hasITS();
        histos.fill(HIST("hThinningQA"), static_cast<float>(trackBits[iTrack]), static_cast<float>(isTPConly));
      }
    }    
  }
  PROCESS_SWITCH(integrationTest, processDataModel2026, "Check data model, 2026 version", false);

  void processDataModelMC(
    aod::McCollisions const& mccollisions,
    aod::McCollisionLabels const& mccollisionlabels,
    aod::McParticles const& mcparticles,
    aod::McTrackLabels const& mctracklabels,
    aod::McMFTTrackLabels const& mcmfttracklabels,
    aod::McFwdTrackLabels const& mcfwdtracklabels)
  {
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcCollision) + 0.5f, mccollisions.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcCollisionLabel) + 0.5f, mccollisionlabels.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcParticle) + 0.5f, mcparticles.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcTrackLabel) + 0.5f, mctracklabels.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcMftTrackLabel) + 0.5f, mcmfttracklabels.size());
    histos.fill(HIST("hTableSizes"), static_cast<float>(kMcFwdTrackLabel) + 0.5f, mcfwdtracklabels.size());
  }
  PROCESS_SWITCH(integrationTest, processDataModelMC, "Check data model for MC tables", false);

  void processBCs(aod::BC const&, aod::Collisions const& collisions)
  {
    histos.fill(HIST("hCollisionsPerBC"), static_cast<float>(collisions.size()));
  }
  PROCESS_SWITCH(integrationTest, processBCs, "Check collisions per BC", true);

  void processCollisions(aod::Collision const&, FullTracksIU const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades)
  {
    Int_t lHasITS = 0, lHasTPC = 0, lHasTRD = 0, lHasTOF = 0, lNotTPCOnly = 0;
    for (const auto& track : tracks) {
      // TPC only bool
      bool isTPConly = track.hasTPC() && !track.hasTOF() && !track.hasTRD() && !track.hasITS();
      histos.fill(HIST("hPt"), track.pt());
      if (track.hasITS()) {
        lHasITS++;
        histos.fill(HIST("hPtITS"), track.pt());
      }
      if (track.hasTPC()) {
        lHasTPC++;
        histos.fill(HIST("hPtTPC"), track.pt());
      }
      if (!isTPConly) {
        lNotTPCOnly++;
        histos.fill(HIST("hPtNoTPCOnly"), track.pt());
      }
      if (track.hasTRD()) {
        lHasTRD++;
        histos.fill(HIST("hPtTRD"), track.pt());
      }
      if (track.hasTOF()) {
        lHasTOF++;
        histos.fill(HIST("hPtTOF"), track.pt());
      }

      if (doBasicQA) {
        histos.fill(HIST("hNCluAll"), track.itsNCls());
        histos.fill(HIST("h2dPhiVsEtaAll"), track.eta(), track.phi());

        if (!isTPConly) {
          histos.fill(HIST("hNCluNoTPCOnly"), track.itsNCls());
          histos.fill(HIST("h2dPhiVsEtaNoTPCOnly"), track.eta(), track.phi());
        }

        // encode particle has??? properties here
        uint8_t encodedHasInfo = (track.hasITS() << 0) | (track.hasTPC() << 1) | (track.hasTRD() << 2) | (track.hasTOF() << 3);
        histos.fill(HIST("hHasDetector"), encodedHasInfo);
      }
    }
    histos.fill(HIST("hTracks"), tracks.size());
    histos.fill(HIST("hTracksITS"), lHasITS);
    histos.fill(HIST("hTracksTPC"), lHasTPC);
    histos.fill(HIST("hTracksNoTPCOnly"), lNotTPCOnly);
    histos.fill(HIST("hTracksTRD"), lHasTRD);
    histos.fill(HIST("hTracksTOF"), lHasTOF);
    if (do2DNTrackCorr) {
      histos.fill(HIST("hTOFvsTRD"), lHasTRD, lHasTOF);
      histos.fill(HIST("hTOFvsTPC"), lHasTPC, lHasTOF);
      histos.fill(HIST("hTRDvsTPC"), lHasTPC, lHasTRD);
      histos.fill(HIST("hITSvsTPC"), lHasTPC, lHasITS);
    }
    histos.fill(HIST("hV0s"), v0s.size());
    histos.fill(HIST("hCascades"), cascades.size());
  }
  PROCESS_SWITCH(integrationTest, processCollisions, "Check collision-associated basics", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<integrationTest>(cfgc)};
}
