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

/// \file multiparticle-correlations-ar.cxx
/// \brief multiparticle-correlations-ar - Task belonging to Anton Riedel for computing multiparticle correlations
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <fairlogger/Logger.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;

namespace MultiParticleCorrelationsARTaskGlobalConfig
{
// setup for event variables
enum EventVariable {
  kVX,
  kVY,
  kVZ,
  kVABS,
  kCEN,
  kMULQ,
  kMULW,
  kMULNC,
  kMULTPC,
  kLAST_EventVariable
};
static constexpr std::string_view EventVariableNames[kLAST_EventVariable] = {"VertexX",
                                                                             "VertexY",
                                                                             "VertexZ",
                                                                             "VertexAbs",
                                                                             "Centrality",
                                                                             "MultiplicityQvector",
                                                                             "MultiplicityWeights",
                                                                             "MultiplicityNumContrib",
                                                                             "MultiplicityTPC"};
// setup for track variables
enum TrackVariable {
  kPT,
  kPHI,
  kETA,
  kCHARGE,
  kDCAZ,
  kDCAXY,
  kTPCCLUSTERS,
  kTPCCROSSEDROWS,
  kTPCCHI2,
  kITSCLUSTERS,
  kLAST_TrackVariable
};
static constexpr std::string_view TrackVariableNames[kLAST_TrackVariable] = {"Pt",
                                                                             "Phi",
                                                                             "Eta",
                                                                             "Charge",
                                                                             "DCAZ",
                                                                             "DCAXY",
                                                                             "TPCClusters",
                                                                             "TPCCrossedRows",
                                                                             "TPCChi2",
                                                                             "ITSClusters"};
// prefixes
static constexpr std::string_view ECrecoBefore = std::string_view("reco/EventControl/before/RB_");
static constexpr std::string_view ECrecoAfter = std::string_view("reco/EventControl/after/RA_");
static constexpr std::string_view ECsimBefore = std::string_view("sim/EventControl/before/SB_");
static constexpr std::string_view ECsimAfter = std::string_view("sim/EventControl/after/SA_");
static constexpr std::string_view TCrecoBefore = std::string_view("reco/TrackControl/before/RB_");
static constexpr std::string_view TCrecoAfter = std::string_view("reco/TrackControl/after/RA_");
static constexpr std::string_view TCsimBefore = std::string_view("sim/TrackControl/before/SB_");
static constexpr std::string_view TCsimAfter = std::string_view("sim/TrackControl/after/SA_");

std::string info = std::string(": Hist bins, Hist lower edge, Hist upper edge, lower cut, upper cut, cut ON(1)/OFF(-1)");

}; // namespace MultiParticleCorrelationsARTaskGlobalConfig

namespace AR = MultiParticleCorrelationsARTaskGlobalConfig;

struct MultiParticleCorrelationsARTask {

  // configurables
  // for event variables
  Configurable<std::vector<double>> cfgVX = {std::string(AR::EventVariableNames[AR::kVX]),
                                             {400., -2., 2., -1., 1., 1.},
                                             std::string("Vertex X") + AR::info};
  Configurable<std::vector<double>> cfgVY = {std::string(AR::EventVariableNames[AR::kVY]),
                                             {400., -2., 2., -1., 1., 1.},
                                             std::string("Vertex Y") + AR::info};
  Configurable<std::vector<double>> cfgVZ = {std::string(AR::EventVariableNames[AR::kVZ]),
                                             {2400., -12., 12., -10., 10., 1.},
                                             std::string("Vertex Z") + AR::info};
  Configurable<std::vector<double>> cfgVABS = {std::string(AR::EventVariableNames[AR::kVABS]),
                                               {150., 0., 15, 1.e-6, 15., 1.},
                                               std::string("Vertex distance from origin") + AR::info};
  Configurable<std::vector<double>> cfgCEN = {std::string(AR::EventVariableNames[AR::kCEN]),
                                              {120., 0., 120., 0., 80., 1.},
                                              std::string("Centrality") + AR::info};
  Configurable<std::vector<double>> cfgMULQ = {std::string(AR::EventVariableNames[AR::kMULQ]),
                                               {3000., 0., 3000., 10., 3000., 1.},
                                               std::string("Multiplicity (QVector)") + AR::info};
  Configurable<std::vector<double>> cfgMULW = {std::string(AR::EventVariableNames[AR::kMULW]),
                                               {3000., 0., 3000., 10., 3000., 1.},
                                               std::string("Multiplicity (Weights)") + AR::info};
  Configurable<std::vector<double>> cfgMULNC = {std::string(AR::EventVariableNames[AR::kMULNC]),
                                                {3000., 0., 3000., 10., 3000., 1.},
                                                std::string("Multiplicity (NumContrib)") + AR::info};
  Configurable<std::vector<double>> cfgMULTPC = {std::string(AR::EventVariableNames[AR::kMULTPC]),
                                                 {3000., 0., 3000., 12., 3000., 1.},
                                                 std::string("Multiplicity (TPC)") + AR::info};
  std::vector<Configurable<std::vector<double>>> cfgEvent = {cfgVX, cfgVY, cfgVZ, cfgVABS, cfgCEN, cfgMULQ, cfgMULW, cfgMULNC, cfgMULTPC};

  // for track variables
  Configurable<std::vector<double>> cfgPT = {std::string(AR::TrackVariableNames[AR::kPT]),
                                             {600., 0., 6., 0.2, 5., 1.},
                                             std::string("pt") + AR::info};
  Configurable<std::vector<double>> cfgPHI = {std::string(AR::TrackVariableNames[AR::kPHI]),
                                              {360., 0., 2. * M_PI, 0., 2. * M_PI, 1.},
                                              std::string("phi") + AR::info};
  Configurable<std::vector<double>> cfgETA = {std::string(AR::TrackVariableNames[AR::kETA]),
                                              {1000., -1., 1., -0.8, 0.8, 1.},
                                              std::string("eta") + AR::info};
  Configurable<std::vector<double>> cfgCHARGE = {std::string(AR::TrackVariableNames[AR::kCHARGE]),
                                                 {5., -2.5, 2.5, -1.5, 1.5, 1.},
                                                 std::string("charge") + AR::info};
  Configurable<std::vector<double>> cfgDCAZ = {std::string(AR::TrackVariableNames[AR::kDCAZ]),
                                               {100., -4., 4., -3.2, 3.2, 1.},
                                               std::string("DCA in Z") + AR::info};
  Configurable<std::vector<double>> cfgDCAXY = {std::string(AR::TrackVariableNames[AR::kDCAXY]),
                                                {100., -3., 3., -2.4, 2.4, 1.},
                                                std::string("DCA in XY") + AR::info};
  Configurable<std::vector<double>> cfgTPCCLUSTERS = {std::string(AR::TrackVariableNames[AR::kTPCCLUSTERS]),
                                                      {160., 0., 160., 80., 161., 1.},
                                                      std::string("TPC clusters") + AR::info};
  Configurable<std::vector<double>> cfgTPCCROSSEDROWS = {std::string(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]),
                                                         {160., 0., 160., 80., 161., 1.},
                                                         std::string("TPC crossed rows") + AR::info};
  Configurable<std::vector<double>> cfgTPCCHI2 = {std::string(AR::TrackVariableNames[AR::kTPCCHI2]),
                                                  {500., 0., 5., 0.4, 4., 1.},
                                                  std::string("TPC chi2") + AR::info};
  Configurable<std::vector<double>> cfgITSCLUSTERS = {std::string(AR::TrackVariableNames[AR::kITSCLUSTERS]),
                                                      {6., 0., 6., 0, 7., 1.},
                                                      std::string("ITS clusters") + AR::info};
  std::vector<Configurable<std::vector<double>>> cfgTrack = {cfgPT, cfgPHI, cfgETA, cfgCHARGE, cfgDCAZ, cfgDCAXY, cfgTPCCLUSTERS, cfgTPCCROSSEDROWS, cfgTPCCHI2, cfgITSCLUSTERS};

  // declare histogram registry
  HistogramRegistry registry{"MultiParticleCorrelationsARTask", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {

    // add control histograms for event observables to registry
    for (auto cfg : cfgEvent) {
      registry.add((std::string(AR::ECrecoBefore) + cfg.name).c_str(), "", HistType::kTH1D, {{static_cast<Int_t>(cfg.value.at(0)), cfg.value.at(1), cfg.value.at(2)}});
      registry.addClone((std::string(AR::ECrecoBefore) + cfg.name).c_str(), (std::string(AR::ECrecoAfter) + cfg.name).c_str());
      registry.addClone((std::string(AR::ECrecoBefore) + cfg.name).c_str(), (std::string(AR::ECsimBefore) + cfg.name).c_str());
      registry.addClone((std::string(AR::ECrecoBefore) + cfg.name).c_str(), (std::string(AR::ECsimAfter) + cfg.name).c_str());
    }

    // add control histograms for track observables to registry
    for (auto cfg : cfgTrack) {
      registry.add((std::string(AR::TCrecoBefore) + cfg.name).c_str(), "", HistType::kTH1D, {{static_cast<Int_t>(cfg.value.at(0)), cfg.value.at(1), cfg.value.at(2)}});
      registry.addClone((std::string(AR::TCrecoBefore) + cfg.name).c_str(), (std::string(AR::TCrecoAfter) + cfg.name).c_str());
      registry.addClone((std::string(AR::TCrecoBefore) + cfg.name).c_str(), (std::string(AR::TCsimBefore) + cfg.name).c_str());
      registry.addClone((std::string(AR::TCrecoBefore) + cfg.name).c_str(), (std::string(AR::TCsimAfter) + cfg.name).c_str());
    }
  }

  template <typename CollisionInstance>
  void FillEventControlHistBC(CollisionInstance const& collision, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kVX]), collision.posX());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kVY]), collision.posY());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kVZ]), collision.posZ());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kVABS]), std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)));
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kCEN]), collision.centRun2V0M());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULQ]), collision.size());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULW]), collision.size());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULNC]), collision.numContrib());
    registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULTPC]), collision.multTPC());
  }
  template <typename CollisionInstance>
  void FillEventControlHistAC(CollisionInstance const& collision, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kVX]), collision.posX());
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kVY]), collision.posY());
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kVZ]), collision.posZ());
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kVABS]), std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)));
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kCEN]), collision.centRun2V0M());
    // registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULQ]), collision.size()); not here
    // registry.fill(HIST(AR::ECrecoBefore) + HIST(AR::EventVariableNames[AR::kMULW]), collision.size()); not here
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kMULNC]), collision.numContrib());
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kMULTPC]), collision.multTPC());
  }

  template <typename CollisionInstance>
  bool SurviveCollisionCut(CollisionInstance const& collision)
  {
    bool SurviveCut = true;

    // cut vx
    if (cfgEvent.at(AR::kVX).value.at(5) > 0 &&
        (collision.posX() < cfgEvent.at(AR::kVX).value.at(3) || collision.posX() > cfgEvent.at(AR::kVX).value.at(4))) {
      SurviveCut = false;
    }
    // cut vy
    if (cfgEvent.at(AR::kVY).value.at(5) > 0 &&
        (collision.posY() < cfgEvent.at(AR::kVY).value.at(3) || collision.posY() > cfgEvent.at(AR::kVY).value.at(4))) {
      SurviveCut = false;
    }
    // cut vz
    if (cfgEvent.at(AR::kVZ).value.at(5) > 0 &&
        (collision.posZ() < cfgEvent.at(AR::kVZ).value.at(3) || collision.posZ() > cfgEvent.at(AR::kVZ).value.at(4))) {
      SurviveCut = false;
    }
    // cut vabs
    if (cfgEvent.at(AR::kVABS).value.at(5) > 0 &&
        (std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)) < cfgEvent.at(AR::kVABS).value.at(3) || std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)) > cfgEvent.at(AR::kVABS).value.at(4))) {
      SurviveCut = false;
    }
    // cut centrality
    if (cfgEvent.at(AR::kCEN).value.at(5) > 0 &&
        (collision.centRun2V0M() < cfgEvent.at(AR::kCEN).value.at(3) || collision.centRun2V0M() > cfgEvent.at(AR::kCEN).value.at(4))) {
      SurviveCut = false;
    }
    // cut multiplicity (qvector)
    if (cfgEvent.at(AR::kMULQ).value.at(5) > 0 &&
        (collision.size() < cfgEvent.at(AR::kMULQ).value.at(3) || collision.size() > cfgEvent.at(AR::kMULQ).value.at(4))) {
      SurviveCut = false;
    }
    // cut multiplicity (weights)
    if (cfgEvent.at(AR::kMULW).value.at(5) > 0 &&
        (collision.size() < cfgEvent.at(AR::kMULW).value.at(3) || collision.size() > cfgEvent.at(AR::kMULW).value.at(4))) {
      SurviveCut = false;
    }
    // cut multiplicity (numContrib)
    if (cfgEvent.at(AR::kMULNC).value.at(5) > 0 &&
        (collision.numContrib() < cfgEvent.at(AR::kMULNC).value.at(3) || collision.numContrib() > cfgEvent.at(AR::kMULNC).value.at(4))) {
      SurviveCut = false;
    }
    // cut multiplicity (tpc)
    if (cfgEvent.at(AR::kMULTPC).value.at(5) > 0 &&
        (collision.multTPC() < cfgEvent.at(AR::kMULTPC).value.at(3) || collision.multTPC() > cfgEvent.at(AR::kMULTPC).value.at(4))) {
      SurviveCut = false;
    }

    return SurviveCut;
  }

  template <typename TrackInstance>
  void FillTrackControlHistBC(TrackInstance const& track, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kPT]), track.pt());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kPHI]), track.phi());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kETA]), track.eta());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kCHARGE]), track.sign());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kDCAZ]), track.dcaZ());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kDCAXY]), track.dcaXY());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kTPCCLUSTERS]), track.tpcNClsFound());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]), track.tpcNClsCrossedRows());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kTPCCHI2]), track.tpcChi2NCl());
    registry.fill(HIST(AR::TCrecoBefore) + HIST(AR::TrackVariableNames[AR::kITSCLUSTERS]), track.itsNCls());
  }
  template <typename TrackInstance>
  void FillTrackControlHistAC(TrackInstance const& track, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kPT]), track.pt());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kPHI]), track.phi());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kETA]), track.eta());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kCHARGE]), track.sign());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kDCAZ]), track.dcaZ());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kDCAXY]), track.dcaXY());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kTPCCLUSTERS]), track.tpcNClsFound());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]), track.tpcNClsCrossedRows());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kTPCCHI2]), track.tpcChi2NCl());
    registry.fill(HIST(AR::TCrecoAfter) + HIST(AR::TrackVariableNames[AR::kITSCLUSTERS]), track.itsNCls());
  }
  template <typename TrackInstance>
  bool SurviveTrackCut(TrackInstance const& track)
  {
    bool SurviveCut = true;

    // cut pt
    if (cfgTrack.at(AR::kPT).value.at(5) > 0 &&
        (track.pt() < cfgTrack.at(AR::kPT).value.at(3) || track.pt() > cfgTrack.at(AR::kPT).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kETA).value.at(5) > 0 &&
        (track.eta() < cfgTrack.at(AR::kETA).value.at(3) || track.eta() > cfgTrack.at(AR::kETA).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kPHI).value.at(5) > 0 &&
        (track.phi() < cfgTrack.at(AR::kPHI).value.at(3) || track.phi() > cfgTrack.at(AR::kPHI).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kCHARGE).value.at(5) > 0 &&
        (track.sign() < cfgTrack.at(AR::kCHARGE).value.at(3) || track.sign() > cfgTrack.at(AR::kCHARGE).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kDCAZ).value.at(5) > 0 &&
        (track.dcaZ() < cfgTrack.at(AR::kDCAZ).value.at(3) || track.dcaZ() > cfgTrack.at(AR::kDCAZ).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kDCAXY).value.at(5) > 0 &&
        (track.dcaXY() < cfgTrack.at(AR::kDCAXY).value.at(3) || track.dcaXY() > cfgTrack.at(AR::kDCAXY).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kTPCCLUSTERS).value.at(5) > 0 &&
        (track.tpcNClsFound() < cfgTrack.at(AR::kTPCCLUSTERS).value.at(3) || track.tpcNClsFound() > cfgTrack.at(AR::kTPCCLUSTERS).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kTPCCROSSEDROWS).value.at(5) > 0 &&
        (track.tpcNClsCrossedRows() < cfgTrack.at(AR::kTPCCROSSEDROWS).value.at(3) || track.tpcNClsCrossedRows() > cfgTrack.at(AR::kTPCCROSSEDROWS).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kTPCCHI2).value.at(5) > 0 &&
        (track.tpcChi2NCl() < cfgTrack.at(AR::kTPCCHI2).value.at(3) || track.tpcChi2NCl() > cfgTrack.at(AR::kTPCCHI2).value.at(4))) {
      SurviveCut = false;
    }
    if (cfgTrack.at(AR::kITSCLUSTERS).value.at(5) > 0 &&
        (track.itsNCls() < cfgTrack.at(AR::kITSCLUSTERS).value.at(3) || track.itsNCls() > cfgTrack.at(AR::kITSCLUSTERS).value.at(4))) {
      SurviveCut = false;
    }

    return SurviveCut;
  }

  using CollisionsInstanceIterator = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>::iterator;
  using TracksInstance = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;
  using TracksInstanceIterator = TracksInstance::iterator;

  void process(CollisionsInstanceIterator const& collision, TracksInstance const& tracks)
  {

    // print collision index
    LOGF(info, "Collision Index : %d/%d", collision.index(), collision.size());

    // fill event control histograms before cutting on the collision
    FillEventControlHistBC<CollisionsInstanceIterator>(collision, registry);

    // cut event
    if (!SurviveCollisionCut<CollisionsInstanceIterator>(collision)) {
      LOGF(info, "Cut Collision %d -> Break", collision.index());
      return;
    }

    // fill event control histograms after cutting on the collision
    FillEventControlHistAC<CollisionsInstanceIterator>(collision, registry);

    LOGF(info, "Number of Tracks: %d", tracks.size());
    UInt_t NumberOfTracks = 0;
    // loop over all tracks in the event
    for (auto const& track : tracks) {

      // fill track control histograms before track cut
      FillTrackControlHistBC<TracksInstanceIterator>(track, registry);

      // cut track
      if (!SurviveTrackCut<TracksInstance::iterator>(track)) {
        // LOGF(info, "Cut Track %d -> Continue", track.index());
        continue;
      }

      // fill track control histograms after surviving track cut
      FillTrackControlHistAC<TracksInstanceIterator>(track, registry);

      NumberOfTracks++;
    }
    LOGF(info, "Surviving Tracks: %d ", NumberOfTracks);
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kMULQ]), NumberOfTracks);
    registry.fill(HIST(AR::ECrecoAfter) + HIST(AR::EventVariableNames[AR::kMULW]), NumberOfTracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiParticleCorrelationsARTask>(cfgc)};
}
