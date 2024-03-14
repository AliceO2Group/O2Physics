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

/// \file LFResonanceInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
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
///    Nasir Mehdi Malik

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFResonanceTablesMergeDF.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers

struct reso2dfmerged {
  //  SliceCache cache;

  Configurable<int> nDF{"nDF", 1, "no of combination of collision"};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using resoCols = aod::ResoCollisions;
  using resoTracks = aod::ResoTracks;

  void init(InitContext const&)
  {
    const AxisSpec axisCent(110, 0, 110, "FT0 (%)");
    histos.add("Event/h1d_ft0_mult_percentile", "FT0 (%)", kTH1F, {axisCent});
  }
  Produces<aod::ResoCollisionDFs> resoCollisionsdf;
  Produces<aod::ResoTrackDFs> reso2trksdf;
  int df = 0;

  std::vector<std::tuple<resoCols::iterator, float, float, float, float, resoTracks>> vecOfTuples;
  std::vector<std::vector<std::tuple<float, float, float, float,
                                     float, float, signed char, unsigned char,
                                     float, float, float, float,
                                     bool, float, float, float,
                                     float, float, float, float,
                                     float, float, bool, bool,
                                     bool, float, float, float>>>
    vecOfVecOfTuples;
  void processTrackDataDF(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    int nCollisions = nDF;
    vecOfTuples.push_back(std::make_tuple(collision, collision.posX(), collision.posY(), collision.posZ(), collision.cent(), tracks));
    std::vector<std::tuple<float, float, float, float,
                           float, float, signed char, unsigned char,
                           float, float, float, float,
                           bool, float, float, float,
                           float, float, float, float,
                           float, float, bool, bool,
                           bool, float, float, float>>
      innerVector;
    for (auto& track : tracks) {
      innerVector.push_back(std::make_tuple(
        track.pt(),
        track.px(),
        track.py(),
        track.pz(),
        track.eta(),
        track.phi(),
        track.sign(),
        (uint8_t)track.tpcNClsCrossedRows(),
        track.dcaXY(),
        track.dcaZ(),
        track.x(),
        track.alpha(),
        track.hasTOF(),
        track.tpcNSigmaPi(),
        track.tpcNSigmaKa(),
        track.tpcNSigmaPr(),
        track.tofNSigmaPi(),
        track.tofNSigmaKa(),
        track.tofNSigmaPr(),
        track.tpcSignal(),
        track.passedITSRefit(),
        track.passedTPCRefit(),
        track.isGlobalTrackWoDCA(),
        track.isPrimaryTrack(),
        track.isPVContributor(),
        track.tpcCrossedRowsOverFindableCls(),
        track.itsChi2NCl(),
        track.tpcChi2NCl()));
    }

    vecOfVecOfTuples.push_back(innerVector);
    innerVector.clear();
    df++;
    if (df < nCollisions)
      return;
    df = 0;

    for (size_t i = 0; i < vecOfTuples.size(); ++i) {
      const auto& tuple = vecOfTuples[i];
      const auto& innerVector = vecOfVecOfTuples[i];

      histos.fill(HIST("Event/h1d_ft0_mult_percentile"), std::get<4>(tuple));
      resoCollisionsdf(std::get<1>(tuple), std::get<2>(tuple), std::get<3>(tuple), std::get<4>(tuple), 0, 0., 0., 0., 0., 0, 0);
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
                    std::get<27>(tuple));
      }
    }

    vecOfTuples.clear();
    vecOfVecOfTuples.clear(); //
  }

  PROCESS_SWITCH(reso2dfmerged, processTrackDataDF, "Process for data merged DF", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<reso2dfmerged>(cfgc)};
}
