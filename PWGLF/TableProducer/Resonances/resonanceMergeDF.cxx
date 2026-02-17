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
///
/// \file resonanceMergeDF.cxx
/// \brief Merges multiple dataframes into a single dataframe
///
///
///  In typical dataframes (DF), we usually observe a range of 200 to 300 collisions.
/// This limited number of collisions often results in most events having very few currentwindowneighbors() for event mixing.
///  However, for resonances analysis, a minimum of 10 currentwindowneighbors() is required.
///  To address this limitation, this script is designed to aggregate information from multiple dataframes into a single dataframe,
///  thereby increasing the number of available current window neighbors for analysis. Here, nDF refers to the number of events or collisions.
///  For instance, if the total number of collisions across all dataframes is, say, 10,836, setting nDF to 10,836 will result in the creation of a single table.
///  Conversely, if nDF is set to 2,709, it will generate four tables, each containing approximately 2,709 collisions.
///  If nDF is set to 1, it will generate tables equal to the number of parent tables.
///
/// ///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///         Nasir Mehdi Malik <nasir.mehdi.malik@cern.ch>
///         Min-jae Kim <minjae.kim@cern.ch>
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers

struct ResonanceMergeDF {
  //  SliceCache cache;
  Configurable<int> nDF{"nDF", 1, "no of combination of collision"};
  Configurable<bool> isLoggingEnabled{"isLoggingEnabled", 0, "print log"};
  Configurable<bool> cpidCut{"cpidCut", 0, "pid cut"};
  Configurable<bool> crejtpc{"crejtpc", 0, "reject electron pion"};
  Configurable<bool> crejtof{"crejtof", 0, "reject electron pion tof"};
  Configurable<bool> isPrimary{"isPrimary", 0, "is Primary only"};
  Configurable<bool> isGlobal{"isGlobal", 0, "Global tracks only"};
  Configurable<float> cDCAXY{"cDCAXY", 1., "value of dcaxy"};
  Configurable<float> cDCAZ{"cDCAZ", 1., "value of dcaz"};
  Configurable<float> nsigmaPr{"nsigmaPr", 6., "nsigma value for proton"};
  Configurable<float> nsigmaKa{"nsigmaKa", 6., "nsigma value for kaon"};
  Configurable<float> nsigmatofPr{"nsigmatofPr", 6., "nsigma value for tof prot"};
  Configurable<float> nsigmatofKa{"nsigmatofKa", 6., "nsigma value for tof kaon"};

  // Xi1530 candidate cuts
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoDCA"};
  Configurable<bool> requireTOF{"requireTOF", false, "Require TOF"};
  Configurable<float> applyTOFveto{"applyTOFveto", 999, "Apply TOF veto with value, 999 for passing all"};
  Configurable<float> nsigmaPi{"nsigmaPi", 5., "nsigma value for pion"};
  Configurable<float> minCent{"minCent", 0., "Minimum centrality"};
  Configurable<float> maxCent{"maxCent", 100., "Maximum centrality"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    const AxisSpec axisCent(110, 0, 110, "FT0 (%)");
    histos.add("Event/h1d_ft0_mult_percentile", "FT0 (%)", kTH1F, {axisCent});
    histos.add("Event/h1d_ft0_mult_percentile_CASC", "FT0 (%)", kTH1F, {axisCent});
  }
  Produces<aod::ResoCollisionDFs> resoCollisionsdf;
  Produces<aod::ResoTrackDFs> reso2trksdf;
  Produces<aod::ResoCascadeDFs> reso2cascadesdf;
  int df = 0;

  std::vector<std::tuple<float, float, float, float, float, float, int>> vecOfTuples;
  std::vector<std::vector<std::tuple<float, float, float, float,
                                     unsigned char, unsigned char,
                                     int16_t, int16_t, int8_t, int8_t, int8_t,
                                     int8_t, int8_t, int8_t, float,
                                     uint8_t>>>
    vecOfVecOfTuples;
  std::vector<std::vector<std::tuple<float, float, float, float,
                                     int*,
                                     int8_t, int8_t, int8_t, // TPC Pos Trk
                                     int8_t, int8_t, int8_t, // TPC Neg Trk
                                     int8_t, int8_t, int8_t, // TPC Bach Trk
                                     int8_t, int8_t, int8_t, // TOF Pos Trk
                                     int8_t, int8_t, int8_t, // TOF Neg Trk
                                     int8_t, int8_t, int8_t, // TOF Bach Trk
                                     float, float, float, float,
                                     float, float, float, float,
                                     float, float, int, float,
                                     float, float, float,
                                     float, float, float>>>
    vecOfVecOfTuplesCasc;
  void processTrackDataDF(aod::ResoCollisions::iterator const& collision, aod::ResoTracks const& tracks)
  {

    int nCollisions = nDF;
    vecOfTuples.push_back(std::make_tuple(collision.posX(), collision.posY(), collision.posZ(), collision.cent(), 0, 0, 0));
    std::vector<std::tuple<float, float, float, float,
                           unsigned char, unsigned char,
                           int16_t, int16_t, int8_t, int8_t, int8_t,
                           int8_t, int8_t, int8_t, float,
                           uint8_t>>
      innerVector;
    for (const auto& track : tracks) {
      if (cpidCut) {
        if (!track.hasTOF()) {
          if (std::abs(track.tpcNSigmaPr()) > nsigmaPr && std::abs(track.tpcNSigmaKa()) > nsigmaKa)
            continue;

          if (crejtpc && (std::abs(track.tpcNSigmaPr()) > std::abs(track.tpcNSigmaPi()) && std::abs(track.tpcNSigmaKa()) > std::abs(track.tpcNSigmaPi())))
            continue;

        } else {
          if (std::abs(track.tofNSigmaPr()) > nsigmatofPr && std::abs(track.tofNSigmaKa()) > nsigmatofKa)
            continue;

          if (crejtof && (std::abs(track.tofNSigmaPr()) > std::abs(track.tofNSigmaPi()) && std::abs(track.tofNSigmaKa()) > std::abs(track.tofNSigmaPi())))
            continue;
        }

        if (std::abs(track.dcaXY()) > cDCAXY)
          continue;
        if (std::abs(track.dcaZ()) > cDCAZ)
          continue;
      }

      innerVector.push_back(std::make_tuple(
        //  track.trackId(),
        track.pt(),
        track.px(),
        track.py(),
        track.pz(),
        (uint8_t)track.tpcNClsCrossedRows(),
        (uint8_t)track.tpcNClsFound(),
        static_cast<int16_t>(track.dcaXY() * 10000),
        static_cast<int16_t>(track.dcaZ() * 10000),
        (int8_t)(track.tpcNSigmaPi() * 10),
        (int8_t)(track.tpcNSigmaKa() * 10),
        (int8_t)(track.tpcNSigmaPr() * 10),
        (int8_t)(track.tofNSigmaPi() * 10),
        (int8_t)(track.tofNSigmaKa() * 10),
        (int8_t)(track.tofNSigmaPr() * 10),
        (int8_t)(track.tpcSignal() * 10),
        track.trackFlags()));
    }

    vecOfVecOfTuples.push_back(innerVector);
    innerVector.clear();
    df++;
    if (isLoggingEnabled) 
      LOGF(info, "collisions: df = %i", df);
    if (df < nCollisions)
      return;
    df = 0;

    for (size_t i = 0; i < vecOfTuples.size(); ++i) {
      const auto& tuple = vecOfTuples[i];
      const auto& innerVector = vecOfVecOfTuples[i];

      histos.fill(HIST("Event/h1d_ft0_mult_percentile"), std::get<3>(tuple));
      resoCollisionsdf(0, std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple), std::get<5>(tuple), 0., 0., 0., 0., 0, std::get<6>(tuple));
      //  LOGF(info, "collisions: Index = %d ) %f - %f - %f %f %d -- %d", std::get<0>(tuple).globalIndex(),std::get<1>(tuple),std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple), std::get<5>(tuple).size(),resoCollisionsdf.lastIndex());

      for (const auto& tuple : innerVector) {
        reso2trksdf(resoCollisionsdf.lastIndex(),
                    std::get<0>(tuple),
                    std::get<1>(tuple),
                    std::get<2>(tuple),
                    std::get<3>(tuple),
                    std::get<4>(tuple),
                    std::get<5>(tuple),
                    std::get<6>(tuple),
                    std::get<7>(tuple),
                    std::get<8>(tuple),
                    std::get<9>(tuple),
                    std::get<10>(tuple),
                    std::get<11>(tuple),
                    std::get<12>(tuple),
                    std::get<13>(tuple),
                    std::get<14>(tuple),
                    std::get<15>(tuple));
      }
    }

    vecOfTuples.clear();
    vecOfVecOfTuples.clear();
  }

  PROCESS_SWITCH(ResonanceMergeDF, processTrackDataDF, "Process for data merged DF", true);

  void processTrackDataDFCasc(aod::ResoCollisions::iterator const& collision, aod::ResoTracks const& tracks, aod::ResoCascades const& trackCascs)
  {

    int nCollisions = nDF;
    vecOfTuples.push_back(std::make_tuple(collision.posX(), collision.posY(), collision.posZ(), collision.cent(), 0, 0, 0));
    std::vector<std::tuple<float, float, float, float,
                           unsigned char, unsigned char,
                           int16_t, int16_t, int8_t, int8_t, int8_t,
                           int8_t, int8_t, int8_t, float,
                           uint8_t>>
      innerVector;
    std::vector<std::tuple<float, float, float, float,
                           int*,
                           int8_t, int8_t, int8_t, // TPC Pos Trk
                           int8_t, int8_t, int8_t, // TPC Neg Trk
                           int8_t, int8_t, int8_t, // TPC Bach Trk
                           int8_t, int8_t, int8_t, // TOF Pos Trk
                           int8_t, int8_t, int8_t, // TOF Neg Trk
                           int8_t, int8_t, int8_t, // TOF Bach Trk
                           float, float, float, float,
                           float, float, float, float,
                           float, float, int, float,
                           float, float, float,
                           float, float, float>>
      innerVectorCasc;
    for (const auto& track : tracks) {
      if (cpidCut) {
        if (!track.hasTOF()) {
          if (std::abs(track.tpcNSigmaPr()) > nsigmaPr && std::abs(track.tpcNSigmaKa()) > nsigmaKa)
            continue;

          if (crejtpc && (std::abs(track.tpcNSigmaPr()) > std::abs(track.tpcNSigmaPi()) && std::abs(track.tpcNSigmaKa()) > std::abs(track.tpcNSigmaPi())))
            continue;

        } else {
          if (std::abs(track.tofNSigmaPr()) > nsigmatofPr && std::abs(track.tofNSigmaKa()) > nsigmatofKa)
            continue;

          if (crejtof && (std::abs(track.tofNSigmaPr()) > std::abs(track.tofNSigmaPi()) && std::abs(track.tofNSigmaKa()) > std::abs(track.tofNSigmaPi())))
            continue;
        }

        if (std::abs(track.dcaXY()) > cDCAXY)
          continue;
        if (std::abs(track.dcaZ()) > cDCAZ)
          continue;
      }

      innerVector.push_back(std::make_tuple(
        //  track.trackId(),
        track.pt(),
        track.px(),
        track.py(),
        track.pz(),
        (uint8_t)track.tpcNClsCrossedRows(),
        (uint8_t)track.tpcNClsFound(),
        static_cast<int16_t>(track.dcaXY() * 10000),
        static_cast<int16_t>(track.dcaZ() * 10000),
        (int8_t)(track.tpcNSigmaPi() * 10),
        (int8_t)(track.tpcNSigmaKa() * 10),
        (int8_t)(track.tpcNSigmaPr() * 10),
        (int8_t)(track.tofNSigmaPi() * 10),
        (int8_t)(track.tofNSigmaKa() * 10),
        (int8_t)(track.tofNSigmaPr() * 10),
        (int8_t)(track.tpcSignal() * 10),
        track.trackFlags()));
    }

    for (const auto& trackCasc : trackCascs) {
      innerVectorCasc.push_back(std::make_tuple(
        trackCasc.pt(),
        trackCasc.px(),
        trackCasc.py(),
        trackCasc.pz(),
        const_cast<int*>(trackCasc.cascadeIndices()),
        (int8_t)(trackCasc.daughterTPCNSigmaPosPi() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaPosKa() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaPosPr() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaNegPi() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaNegKa() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaNegPr() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaBachPi() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaBachKa() * 10),
        (int8_t)(trackCasc.daughterTPCNSigmaBachPr() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaPosPi() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaPosKa() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaPosPr() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaNegPi() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaNegKa() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaNegPr() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaBachPi() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaBachKa() * 10),
        (int8_t)(trackCasc.daughterTOFNSigmaBachPr() * 10),
        trackCasc.v0CosPA(),
        trackCasc.cascCosPA(),
        trackCasc.daughDCA(),
        trackCasc.cascDaughDCA(),
        trackCasc.dcapostopv(),
        trackCasc.dcanegtopv(),
        trackCasc.dcabachtopv(),
        trackCasc.dcav0topv(),
        trackCasc.dcaXYCascToPV(),
        trackCasc.dcaZCascToPV(),
        trackCasc.sign(),
        trackCasc.mLambda(),
        trackCasc.mXi(),
        trackCasc.transRadius(), trackCasc.cascTransRadius(), trackCasc.decayVtxX(), trackCasc.decayVtxY(), trackCasc.decayVtxZ()));
    }

    vecOfVecOfTuples.push_back(innerVector);
    vecOfVecOfTuplesCasc.push_back(innerVectorCasc);
    innerVector.clear();
    innerVectorCasc.clear();

    df++;
    LOGF(info, "collisions: df = %i", df);
    if (df < nCollisions)
      return;
    df = 0;

    for (size_t i = 0; i < vecOfTuples.size(); ++i) {
      const auto& tuple = vecOfTuples[i];
      const auto& innerVector = vecOfVecOfTuples[i];
      const auto& innerVectorCasc = vecOfVecOfTuplesCasc[i];

      histos.fill(HIST("Event/h1d_ft0_mult_percentile"), std::get<3>(tuple));
      resoCollisionsdf(0, std::get<0>(tuple), std::get<1>(tuple), std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple), std::get<5>(tuple), 0., 0., 0., 0., 0, std::get<6>(tuple));
      //  LOGF(info, "collisions: Index = %d ) %f - %f - %f %f %d -- %d", std::get<0>(tuple).globalIndex(),std::get<1>(tuple),std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple), std::get<5>(tuple).size(),resoCollisionsdf.lastIndex());

      for (const auto& tuple : innerVector) {
        reso2trksdf(resoCollisionsdf.lastIndex(),
                    std::get<0>(tuple),
                    std::get<1>(tuple),
                    std::get<2>(tuple),
                    std::get<3>(tuple),
                    std::get<4>(tuple),
                    std::get<5>(tuple),
                    std::get<6>(tuple),
                    std::get<7>(tuple),
                    std::get<8>(tuple),
                    std::get<9>(tuple),
                    std::get<10>(tuple),
                    std::get<11>(tuple),
                    std::get<12>(tuple),
                    std::get<13>(tuple),
                    std::get<14>(tuple),
                    std::get<15>(tuple));
      }

      for (const auto& tuple : innerVectorCasc) {
        reso2cascadesdf(resoCollisionsdf.lastIndex(),
                        std::get<0>(tuple),
                        std::get<1>(tuple),
                        std::get<2>(tuple),
                        std::get<3>(tuple),
                        std::get<4>(tuple),
                        std::get<5>(tuple),
                        std::get<6>(tuple),
                        std::get<7>(tuple),
                        std::get<8>(tuple),
                        std::get<9>(tuple),
                        std::get<10>(tuple),
                        std::get<11>(tuple),
                        std::get<12>(tuple),
                        std::get<13>(tuple),
                        std::get<14>(tuple),
                        std::get<15>(tuple),
                        std::get<16>(tuple),
                        std::get<17>(tuple),
                        std::get<18>(tuple),
                        std::get<19>(tuple),
                        std::get<20>(tuple),
                        std::get<21>(tuple),
                        std::get<22>(tuple),
                        std::get<23>(tuple),
                        std::get<24>(tuple),
                        std::get<25>(tuple),
                        std::get<26>(tuple),
                        std::get<27>(tuple),
                        std::get<28>(tuple),
                        std::get<29>(tuple),
                        std::get<30>(tuple),
                        std::get<31>(tuple),
                        std::get<32>(tuple),
                        std::get<33>(tuple),
                        std::get<34>(tuple),
                        std::get<35>(tuple),
                        std::get<36>(tuple),
                        std::get<37>(tuple),
                        std::get<38>(tuple),
                        std::get<39>(tuple),
                        std::get<40>(tuple));
      }
    }

    vecOfTuples.clear();
    vecOfVecOfTuples.clear();
    vecOfVecOfTuplesCasc.clear(); //
  }

  PROCESS_SWITCH(ResonanceMergeDF, processTrackDataDFCasc, "Process for data merged DF for cascade", false);

  void processLambdaStarCandidate(aod::ResoCollisions::iterator const& collision, aod::ResoTracks const& tracks)
  {

    if (doprocessTrackDataDF)
      LOG(fatal) << "Disable processTrackDataDF first!";

    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent());

    resoCollisionsdf(0, collision.posX(), collision.posY(), collision.posZ(), collision.cent(), 0, 0, 0., 0., 0., 0., 0, 0);

    for (const auto& track : tracks) {
      if (isPrimary && !track.isPrimaryTrack())
        continue;
      if (isGlobal && !track.isGlobalTrack())
        continue;
      if (!track.hasTOF()) {
        if (std::abs(track.tpcNSigmaPr()) > nsigmaPr && std::abs(track.tpcNSigmaKa()) > nsigmaKa)
          continue;

        if (crejtpc && (std::abs(track.tpcNSigmaPr()) > std::abs(track.tpcNSigmaPi()) && std::abs(track.tpcNSigmaKa()) > std::abs(track.tpcNSigmaPi())))
          continue;

      } else {
        if (std::abs(track.tofNSigmaPr()) > nsigmatofPr && std::abs(track.tofNSigmaKa()) > nsigmatofKa)
          continue;

        if (crejtof && (std::abs(track.tofNSigmaPr()) > std::abs(track.tofNSigmaPi()) && std::abs(track.tofNSigmaKa()) > std::abs(track.tofNSigmaPi())))
          continue;
      }
      if (std::abs(track.dcaXY()) > cDCAXY)
        continue;
      if (std::abs(track.dcaZ()) > cDCAZ)
        continue;
      reso2trksdf(resoCollisionsdf.lastIndex(),
                  // track.trackId(),
                  track.pt(),
                  track.px(),
                  track.py(),
                  track.pz(),
                  (uint8_t)track.tpcNClsCrossedRows(),
                  (uint8_t)track.tpcNClsFound(),
                  static_cast<int16_t>(track.dcaXY() * 10000),
                  static_cast<int16_t>(track.dcaZ() * 10000),
                  (int8_t)(track.tpcNSigmaPi() * 10),
                  (int8_t)(track.tpcNSigmaKa() * 10),
                  (int8_t)(track.tpcNSigmaPr() * 10),
                  (int8_t)(track.tofNSigmaPi() * 10),
                  (int8_t)(track.tofNSigmaKa() * 10),
                  (int8_t)(track.tofNSigmaPr() * 10),
                  (int8_t)(track.tpcSignal() * 10),
                  track.trackFlags());
    }
  }
  PROCESS_SWITCH(ResonanceMergeDF, processLambdaStarCandidate, "Process for lambda star candidate", false);

  void processXi1530Candidate(aod::ResoCollisions::iterator const& collision, aod::ResoTracks const& tracks, aod::ResoCascades const& resocasctracks)
  {
    if (doprocessTrackDataDF)
      LOG(fatal) << "Disable processTrackDataDF first!";
    if (doprocessLambdaStarCandidate)
      LOG(fatal) << "Disable processLambdaStarCandidate first!";

    if (collision.cent() < minCent || collision.cent() > maxCent)
      return;

    resoCollisionsdf(0, collision.posX(), collision.posY(), collision.posZ(), collision.cent(), 0, 0, 0., 0., 0., 0., 0, 0);
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent());

    for (const auto& track : tracks) {
      if (trackSelection == 1) {
        if (!track.isGlobalTrack())
          continue;
      } else if (trackSelection == 2) {
        if (!track.isGlobalTrackWoDCA())
          continue;
      }
      if (!track.hasTOF()) {
        if (requireTOF) {
          continue;
        }
        // TPC selection
        if (std::abs(track.tpcNSigmaPi()) > nsigmaPi)
          continue;
      } else {
        if (applyTOFveto > 998 && std::abs(track.tofNSigmaPi()) > applyTOFveto)
          continue;
        // TPC selection
        if (std::abs(track.tpcNSigmaPi()) > nsigmaPi)
          continue;
      }

      if (std::abs(track.dcaXY()) > cDCAXY)
        continue;
      if (std::abs(track.dcaZ()) > cDCAZ)
        continue;

      reso2trksdf(resoCollisionsdf.lastIndex(),
                  // track.trackId(),
                  track.pt(),
                  track.px(),
                  track.py(),
                  track.pz(),
                  (uint8_t)track.tpcNClsCrossedRows(),
                  (uint8_t)track.tpcNClsFound(),
                  static_cast<int16_t>(track.dcaXY() * 10000),
                  static_cast<int16_t>(track.dcaZ() * 10000),
                  (int8_t)(track.tpcNSigmaPi() * 10),
                  (int8_t)(track.tpcNSigmaKa() * 10),
                  (int8_t)(track.tpcNSigmaPr() * 10),
                  (int8_t)(track.tofNSigmaPi() * 10),
                  (int8_t)(track.tofNSigmaKa() * 10),
                  (int8_t)(track.tofNSigmaPr() * 10),
                  (int8_t)(track.tpcSignal() * 10),
                  track.trackFlags());
    }
    // Cascade candidate
    for (const auto& track : resocasctracks) {
      // TODO: add cascade cuts
      reso2cascadesdf(resoCollisionsdf.lastIndex(),
                      // casc.globalIndex(),
                      track.pt(),
                      track.px(),
                      track.py(),
                      track.pz(),
                      const_cast<int*>(track.cascadeIndices()),
                      (int8_t)(track.daughterTPCNSigmaPosPi() * 10),
                      (int8_t)(track.daughterTPCNSigmaPosKa() * 10),
                      (int8_t)(track.daughterTPCNSigmaPosPr() * 10),
                      (int8_t)(track.daughterTPCNSigmaNegPi() * 10),
                      (int8_t)(track.daughterTPCNSigmaNegKa() * 10),
                      (int8_t)(track.daughterTPCNSigmaNegPr() * 10),
                      (int8_t)(track.daughterTPCNSigmaBachPi() * 10),
                      (int8_t)(track.daughterTPCNSigmaBachKa() * 10),
                      (int8_t)(track.daughterTPCNSigmaBachPr() * 10),
                      (int8_t)(track.daughterTOFNSigmaPosPi() * 10),
                      (int8_t)(track.daughterTOFNSigmaPosKa() * 10),
                      (int8_t)(track.daughterTOFNSigmaPosPr() * 10),
                      (int8_t)(track.daughterTOFNSigmaNegPi() * 10),
                      (int8_t)(track.daughterTOFNSigmaNegKa() * 10),
                      (int8_t)(track.daughterTOFNSigmaNegPr() * 10),
                      (int8_t)(track.daughterTOFNSigmaBachPi() * 10),
                      (int8_t)(track.daughterTOFNSigmaBachKa() * 10),
                      (int8_t)(track.daughterTOFNSigmaBachPr() * 10),
                      track.v0CosPA(),
                      track.cascCosPA(),
                      track.daughDCA(),
                      track.cascDaughDCA(),
                      track.dcapostopv(),
                      track.dcanegtopv(),
                      track.dcabachtopv(),
                      track.dcav0topv(),
                      track.dcaXYCascToPV(),
                      track.dcaZCascToPV(),
                      track.sign(),
                      track.mLambda(),
                      track.mXi(),
                      track.transRadius(), track.cascTransRadius(), track.decayVtxX(), track.decayVtxY(), track.decayVtxZ());
    }
  }
  PROCESS_SWITCH(ResonanceMergeDF, processXi1530Candidate, "Process for Xi(1530) candidate", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ResonanceMergeDF>(cfgc)};
}
