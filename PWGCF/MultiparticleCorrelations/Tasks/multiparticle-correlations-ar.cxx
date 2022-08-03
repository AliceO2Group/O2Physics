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

#include "fairlogger/Logger.h"
#include <cstdint>
#include <string_view>
#include <string>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Define useful constant and enums in a seperate name space
namespace MultiParticleCorrelationsARTaskGlobalConfig
{
// for before and after applying cuts
enum BA {
  kBEFORE,
  kAFTER,
  kLAST_BA
};
static constexpr std::string_view BeforeAfter[kLAST_BA] = {"before/", "after/"};

// for reconstructed and simulated data
enum RS {
  kRECO,
  kSIM,
  kLAST_RS
};
static constexpr std::string_view RecoSim[kLAST_RS] = {"reco/", "sim/"};
enum HistConfig {
  kBIN,   // number of bins
  kLEDGE, // lower edge of the histogram
  kUEDGE, // upper edge of the histogram
  kLCUT,  // lower cut
  kUCUT,  // upper cut
  kOCUT,  // option whether to apply cut at all
  kLAST_HistConfig
};

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

// info string for configurables
std::string info = std::string(": Hist bins, Hist lower edge, Hist upper edge, lower cut, upper cut, cut ON(1)/OFF(-1)");

inline float abs(float vx, float vy, float vz)
{
  return std::sqrt(vx * vx + vy * vy * vz * vz);
}
}; // namespace MultiParticleCorrelationsARTaskGlobalConfig

// use an alias for the namespace
namespace AR = MultiParticleCorrelationsARTaskGlobalConfig;

struct MultiParticleCorrelationsARTask {

  // configurables and cuts (where possible)
  // for event variables
  Configurable<std::vector<float>> cfgVX = {std::string(AR::EventVariableNames[AR::kVX]),
                                            {400., -2., 2., -1., 1., 1.},
                                            std::string("Vertex X") + AR::info};
  Filter filterVX = aod::collision::posX > cfgVX.value.at(AR::kLCUT) &&
                    aod::collision::posX<cfgVX.value.at(AR::kUCUT) &&
                                         cfgVX.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgVY = {std::string(AR::EventVariableNames[AR::kVY]),
                                            {400., -2., 2., -1., 1., 1.},
                                            std::string("Vertex Y") + AR::info};
  Filter filterVY = aod::collision::posY > cfgVY.value.at(AR::kLCUT) &&
                    aod::collision::posY<cfgVY.value.at(AR::kUCUT) &&
                                         cfgVY.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgVZ = {std::string(AR::EventVariableNames[AR::kVZ]),
                                            {2400., -12., 12., -10., 10., 1.},
                                            std::string("Vertex Z") + AR::info};
  Filter filterVZ = aod::collision::posZ > cfgVZ.value.at(AR::kLCUT) &&
                    aod::collision::posZ<cfgVZ.value.at(AR::kUCUT) &&
                                         cfgVZ.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgVABS = {std::string(AR::EventVariableNames[AR::kVABS]),
                                              {150., 0., 15, 1.e-6, 15., 1.},
                                              std::string("Vertex distance from origin") + AR::info};
  Configurable<std::vector<float>> cfgCEN = {std::string(AR::EventVariableNames[AR::kCEN]),
                                             {120., 0., 120., 0., 80., 1.},
                                             std::string("Centrality") + AR::info};
  Filter filterCEN = aod::cent::centRun2V0M > cfgCEN.value.at(AR::kLCUT) &&
                     aod::cent::centRun2V0M<cfgCEN.value.at(AR::kUCUT) &&
                                            cfgCEN.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgMULQ = {std::string(AR::EventVariableNames[AR::kMULQ]),
                                              {3000., 0., 3000., 10., 3000., 1.},
                                              std::string("Multiplicity (QVector)") + AR::info};
  Configurable<std::vector<float>> cfgMULW = {std::string(AR::EventVariableNames[AR::kMULW]),
                                              {3000., 0., 3000., 10., 3000., 1.},
                                              std::string("Multiplicity (Weights)") + AR::info};
  Configurable<std::vector<float>> cfgMULNC = {std::string(AR::EventVariableNames[AR::kMULNC]),
                                               {3000., 0., 3000., 10., 3000., 1.},
                                               std::string("Multiplicity (NumContrib)") + AR::info};
  Filter filterMULNC = aod::collision::numContrib > static_cast<uint16_t>(cfgMULNC.value.at(AR::kLCUT)) &&
                       aod::collision::numContrib<static_cast<uint16_t>(cfgMULNC.value.at(AR::kUCUT)) &&
                                                  cfgMULNC.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgMULTPC = {std::string(AR::EventVariableNames[AR::kMULTPC]),
                                                {3000., 0., 3000., 12., 3000., 1.},
                                                std::string("Multiplicity (TPC)") + AR::info};
  Filter filterMULTPC = aod::collision::numContrib > static_cast<uint16_t>(cfgMULTPC.value.at(AR::kLCUT)) &&
                        aod::collision::numContrib<static_cast<uint16_t>(cfgMULTPC.value.at(AR::kUCUT)) &&
                                                   cfgMULTPC.value.at(AR::kOCUT)> 0.;

  // write all event configurables into a vector
  std::vector<Configurable<std::vector<float>>> cfgEvent = {cfgVX,
                                                            cfgVY,
                                                            cfgVZ,
                                                            cfgVABS,
                                                            cfgCEN,
                                                            cfgMULQ,
                                                            cfgMULW,
                                                            cfgMULNC,
                                                            cfgMULTPC};

  // for track variables and cuts (where possible)
  Configurable<std::vector<float>> cfgPT = {std::string(AR::TrackVariableNames[AR::kPT]),
                                            {600., 0., 6., 0.2, 5., 1.},
                                            std::string("pt") + AR::info};
  Filter filterPT = aod::track::pt > cfgPT.value.at(AR::kLCUT) &&
                    aod::track::pt<cfgPT.value.at(AR::kUCUT) &&
                                   cfgPT.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgPHI = {std::string(AR::TrackVariableNames[AR::kPHI]),
                                             {360., 0., 2. * M_PI, 0., 2. * M_PI, 1.},
                                             std::string("phi") + AR::info};
  Filter filterPHI = aod::track::phi > cfgPHI.value.at(AR::kLCUT) &&
                     aod::track::phi<cfgPHI.value.at(AR::kUCUT) &&
                                     cfgPHI.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgETA = {std::string(AR::TrackVariableNames[AR::kETA]),
                                             {1000., -1., 1., -0.8, 0.8, 1.},
                                             std::string("eta") + AR::info};
  Filter filterETA = aod::track::eta > cfgETA.value.at(AR::kLCUT) &&
                     aod::track::eta<cfgETA.value.at(AR::kUCUT) &&
                                     cfgETA.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgCHARGE = {std::string(AR::TrackVariableNames[AR::kCHARGE]),
                                                {5., -2.5, 2.5, -1.5, 1.5, 1.},
                                                std::string("charge") + AR::info};
  Configurable<std::vector<float>> cfgDCAZ = {std::string(AR::TrackVariableNames[AR::kDCAZ]),
                                              {100., -4., 4., -3.2, 3.2, 1.},
                                              std::string("DCA in Z") + AR::info};
  Filter filterDCAZ = aod::track::dcaZ > cfgDCAZ.value.at(AR::kLCUT) &&
                      aod::track::dcaZ<cfgDCAZ.value.at(AR::kUCUT) &&
                                       cfgDCAZ.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgDCAXY = {std::string(AR::TrackVariableNames[AR::kDCAXY]),
                                               {100., -3., 3., -2.4, 2.4, 1.},
                                               std::string("DCA in XY") + AR::info};
  Filter filterDCAXY = aod::track::dcaXY > cfgDCAXY.value.at(AR::kLCUT) &&
                       aod::track::dcaXY<cfgDCAXY.value.at(AR::kUCUT) &&
                                         cfgDCAXY.value.at(AR::kOCUT)> 0.;
  Configurable<std::vector<float>> cfgTPCCLUSTERS = {std::string(AR::TrackVariableNames[AR::kTPCCLUSTERS]),
                                                     {160., 0., 160., 80., 161., 1.},
                                                     std::string("TPC clusters") + AR::info};
  Configurable<std::vector<float>> cfgTPCCROSSEDROWS = {std::string(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]),
                                                        {160., 0., 160., 80., 161., 1.},
                                                        std::string("TPC crossed rows") + AR::info};
  Configurable<std::vector<float>> cfgTPCCHI2 = {std::string(AR::TrackVariableNames[AR::kTPCCHI2]),
                                                 {500., 0., 5., 0.4, 4., 1.},
                                                 std::string("TPC chi2") + AR::info};
  Configurable<std::vector<float>> cfgITSCLUSTERS = {std::string(AR::TrackVariableNames[AR::kITSCLUSTERS]),
                                                     {6., 0., 6., 0, 7., 1.},
                                                     std::string("ITS clusters") + AR::info};
  // write all track configurables into a vector
  std::vector<Configurable<std::vector<float>>> cfgTrack = {cfgPT,
                                                            cfgPHI,
                                                            cfgETA,
                                                            cfgCHARGE,
                                                            cfgDCAZ,
                                                            cfgDCAXY,
                                                            cfgTPCCLUSTERS,
                                                            cfgTPCCROSSEDROWS,
                                                            cfgTPCCHI2,
                                                            cfgITSCLUSTERS};

  // declare histogram registry
  HistogramRegistry registry{"MultiParticleCorrelationsARTask",
                             {},
                             OutputObjHandlingPolicy::AnalysisObject,
                             false,
                             false};

  void init(InitContext&)
  {

    // add control histograms for event/track observables to registry
    for (int rs = 0; rs < AR::kLAST_RS; rs++) {
      for (int ba = 0; ba < AR::kLAST_BA; ba++) {

        // iterate over event configurables
        for (auto cfg : cfgEvent) {
          registry.add((std::string(AR::RecoSim[rs]) +
                        std::string("EventControl/") +
                        std::string(AR::BeforeAfter[ba]) +
                        cfg.name)
                         .c_str(),
                       "",
                       HistType::kTH1D,
                       {{static_cast<Int_t>(cfg.value.at(AR::kBIN)),
                         cfg.value.at(AR::kLEDGE),
                         cfg.value.at(AR::kUEDGE)}});
        }
        // iterate over track configurables
        for (auto cfg : cfgTrack) {
          registry.add((std::string(AR::RecoSim[rs]) +
                        std::string("TrackControl/") +
                        std::string(AR::BeforeAfter[ba]) +
                        cfg.name)
                         .c_str(),
                       "",
                       HistType::kTH1D,
                       {{static_cast<Int_t>(cfg.value.at(AR::kBIN)),
                         cfg.value.at(AR::kLEDGE),
                         cfg.value.at(AR::kUEDGE)}});
        }
      }
    }
  }

  // function for filling control histograms of event variables
  template <AR::RS rs, AR::BA ba, typename CollisionObject>
  void FillEventControlHist(CollisionObject const& collision, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVX]),
                  collision.posX());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVY]),
                  collision.posY());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVZ]),
                  collision.posZ());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kVABS]),
                  std::sqrt(std::pow(collision.posX(), 2) + std::pow(collision.posY(), 2) + std::pow(collision.posZ(), 2)));
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kCEN]),
                  collision.centRun2V0M());
    // registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULQ]), collision.size()); fill separately
    // registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULW]), collision.size()); fill separately
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULNC]),
                  collision.numContrib());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULTPC]),
                  collision.multTPC());
  }
  template <AR::RS rs, AR::BA ba>
  void FillEventControlHistMul(HistogramRegistry& registry, double mulQvector, double mulWeights)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULQ]),
                  mulQvector);
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULW]),
                  mulWeights);
  }

  // function for fill control histograms of track variables
  template <AR::RS rs, AR::BA ba, typename TrackObject>
  void FillTrackControlHist(TrackObject const& track, HistogramRegistry& registry)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kPT]),
                  track.pt());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kPHI]),
                  track.phi());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kETA]),
                  track.eta());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kCHARGE]),
                  track.sign());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kDCAZ]),
                  track.dcaZ());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kDCAXY]),
                  track.dcaXY());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCLUSTERS]),
                  track.tpcNClsFound());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCROSSEDROWS]),
                  track.tpcNClsCrossedRows());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kTPCCHI2]),
                  track.tpcChi2NCl());
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("TrackControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::TrackVariableNames[AR::kITSCLUSTERS]),
                  track.itsNCls());
  }

  // templated function for processing data
  template <AR::RS rs, AR::BA ba, typename CollisionObject, typename TrackObject>
  void ProcessReconstructed(CollisionObject const& collision, TrackObject const& tracks)
  {
    // fill event control histograms before cutting on the collision
    FillEventControlHist<rs, ba, CollisionObject>(collision, registry);
    // fill multiplicity for qvector and weights separately
    FillEventControlHistMul<rs, ba>(registry, collision.size(), collision.size());

    // loop over all tracks in the event
    for (auto const& track : tracks) {
      // fill track control histograms before track cut
      FillTrackControlHist<rs, ba, typename TrackObject::iterator>(track, registry);
    }
  }

  using CollisionsInstance = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>;
  using TracksInstance = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;

  using CollisionsInstanceIterator = CollisionsInstance::iterator;
  // using TracksInstanceIterator = TracksInstance::iterator;

  // process reconstructed data before applying cuts
  void ProcessReconstructedBeforeCut(CollisionsInstanceIterator const& collision, TracksInstance const& tracks)
  {
    LOGF(info, "Process reconstructed data before applying cuts");
    ProcessReconstructed<AR::kRECO, AR::kBEFORE, CollisionsInstanceIterator, TracksInstance>(collision, tracks);
  };
  PROCESS_SWITCH(MultiParticleCorrelationsARTask, ProcessReconstructedBeforeCut, "Process Reco before Cuts", true);

  using FilteredCollisionsInstance = soa::Filtered<CollisionsInstance>;
  using FilteredTracksInstance = soa::Filtered<TracksInstance>;

  using FilteredCollisionsInstanceIterator = FilteredCollisionsInstance::iterator;
  // using FilteredTracksInstanceIterator = FilteredTracksInstance::iterator;

  // process reconstructed data after applying cuts
  void ProcessReconstructedAfterCut(FilteredCollisionsInstanceIterator const& collision, FilteredTracksInstance const& tracks)
  {
    LOGF(info, "Process reconstructed data after applying cuts");
    ProcessReconstructed<AR::kRECO, AR::kAFTER, FilteredCollisionsInstanceIterator, FilteredTracksInstance>(collision, tracks);
  };
  PROCESS_SWITCH(MultiParticleCorrelationsARTask, ProcessReconstructedAfterCut, "Process Reco after Cuts", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiParticleCorrelationsARTask>(cfgc)};
}
