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
#include <vector>
#include <array>
#include <numeric>
#include "TComplex.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Define useful constants and enums in a separate name space
namespace MultiParticleCorrelationsARTaskGlobalVariables
{
// for before and after applying cuts
enum BeforeAfterEnum {
  kBEFORE,
  kAFTER,
  kLAST_BA
};
static constexpr std::string_view BeforeAfter[kLAST_BA] = {"before/", "after/"};

// for reconstructed and simulated data
enum RecoSimEnum {
  kRECO,
  kSIM,
  kLAST_RS
};
static constexpr std::string_view RecoSim[kLAST_RS] = {"reco/", "sim/"};

// for configuring control histograms
enum HistConfig {
  kBIN,   // number of bins
  kLEDGE, // lower edge of the histogram
  kUEDGE, // upper edge of the histogram
  kLCUT,  // lower cut
  kUCUT,  // upper cut
  kOCUT,  // option whether to apply cut at all
  kLAST_HistConfig
};

// event variables
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
// track variables
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

// common info string for all configurables
std::string info = std::string(": Hist bins, Hist lower edge, Hist upper edge, lower cut, upper cut, cut ON(1)/OFF(-1)");

const int MaxHarmonic = 10;
const int MaxPower = 10;

// function for computing the absolute distance of the primary vertex form the origin
inline float abs(float vx, float vy, float vz)
{
  return std::sqrt(vx * vx + vy * vy * vz * vz);
}

// generice function for checking if the value of a variable passes a cut
inline bool SurviveCut(std::vector<float> ConfigValue, float Value)
{
  bool flag = true;
  // check if the cut is configured to be use in the first place
  if (ConfigValue.at(kOCUT) > 0.) {
    // check if the value of the variable is not lower than the lower bound and
    // not not larger than the upper bound
    if (!(Value > ConfigValue.at(kLCUT) && Value < ConfigValue.at(kUCUT))) {
      flag = false;
    }
  }
  return flag;
}
}; // namespace MultiParticleCorrelationsARTaskGlobalVariables

// use an alias for our namespace
namespace AR = MultiParticleCorrelationsARTaskGlobalVariables;

struct MultiParticleCorrelationsARTask {

  // event configurables and cuts
  Configurable<std::vector<float>> cfgVX = {std::string(AR::EventVariableNames[AR::kVX]),
                                            {400., -2., 2., -1., 1., 1.},
                                            std::string("Vertex X") + AR::info};
  Configurable<std::vector<float>> cfgVY = {std::string(AR::EventVariableNames[AR::kVY]),
                                            {400., -2., 2., -1., 1., 1.},
                                            std::string("Vertex Y") + AR::info};
  Configurable<std::vector<float>> cfgVZ = {std::string(AR::EventVariableNames[AR::kVZ]),
                                            {2400., -12., 12., -10., 10., 1.},
                                            std::string("Vertex Z") + AR::info};
  Configurable<std::vector<float>> cfgVABS = {std::string(AR::EventVariableNames[AR::kVABS]),
                                              {150., 0., 15, 1.e-6, 15., 1.},
                                              std::string("Vertex distance from origin") + AR::info};
  Configurable<std::vector<float>> cfgCEN = {std::string(AR::EventVariableNames[AR::kCEN]),
                                             {120., 0., 120., 0., 80., 1.},
                                             std::string("Centrality") + AR::info};
  Configurable<std::vector<float>> cfgMULQ = {std::string(AR::EventVariableNames[AR::kMULQ]),
                                              {3000., 0., 3000., 10., 3000., 1.},
                                              std::string("Multiplicity (QVector)") + AR::info};
  Configurable<std::vector<float>> cfgMULW = {std::string(AR::EventVariableNames[AR::kMULW]),
                                              {3000., 0., 3000., 10., 3000., 1.},
                                              std::string("Multiplicity (Weights)") + AR::info};
  Configurable<std::vector<float>> cfgMULNC = {std::string(AR::EventVariableNames[AR::kMULNC]),
                                               {3000., 0., 3000., 10., 3000., 1.},
                                               std::string("Multiplicity (NumContrib)") + AR::info};
  Configurable<std::vector<float>> cfgMULTPC = {std::string(AR::EventVariableNames[AR::kMULTPC]),
                                                {3000., 0., 3000., 12., 3000., 1.},
                                                std::string("Multiplicity (TPC)") + AR::info};
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

  // track configurables and cuts
  Configurable<std::vector<float>> cfgPT = {std::string(AR::TrackVariableNames[AR::kPT]),
                                            {600., 0., 6., 0.2, 5., 1.},
                                            std::string("pt") + AR::info};
  Configurable<std::vector<float>> cfgPHI = {std::string(AR::TrackVariableNames[AR::kPHI]),
                                             {360., 0., 2. * M_PI, 0., 2. * M_PI, 1.},
                                             std::string("phi") + AR::info};
  Configurable<std::vector<float>> cfgETA = {std::string(AR::TrackVariableNames[AR::kETA]),
                                             {1000., -1., 1., -0.8, 0.8, 1.},
                                             std::string("eta") + AR::info};
  Configurable<std::vector<float>> cfgCHARGE = {std::string(AR::TrackVariableNames[AR::kCHARGE]),
                                                {5., -2.5, 2.5, -1.5, 1.5, 1.},
                                                std::string("charge") + AR::info};
  Configurable<std::vector<float>> cfgDCAZ = {std::string(AR::TrackVariableNames[AR::kDCAZ]),
                                              {100., -4., 4., -3.2, 3.2, 1.},
                                              std::string("DCA in Z") + AR::info};
  Configurable<std::vector<float>> cfgDCAXY = {std::string(AR::TrackVariableNames[AR::kDCAXY]),
                                               {100., -3., 3., -2.4, 2.4, 1.},
                                               std::string("DCA in XY") + AR::info};
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

  // declare 2d array for qvectors
  std::array<std::array<TComplex, AR::MaxHarmonic>, AR::MaxPower> QVectors;

  // declare objects for computing qvectors
  std::vector<double> Angles;
  std::vector<double> Weights;

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

  // function for filling event control histograms
  // exclude MultiplicityQ, i.e. number of tracks in the QVector, and
  // MultiplicityW, i.e. the weighted number of tracks in the QVector, and fill them separably
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba, typename CollisionObject>
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

  // function for filling event control histograms for Multiplicity{Q,W},
  // since they have to be computed separately
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba>
  void FillEventControlHistMul(HistogramRegistry& registry, double mulQvector, double mulWeights)
  {
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULQ]),
                  mulQvector);
    registry.fill(HIST(AR::RecoSim[rs]) + HIST("EventControl/") + HIST(AR::BeforeAfter[ba]) + HIST(AR::EventVariableNames[AR::kMULW]),
                  mulWeights);
  }

  // function for filling track control histograms
  template <AR::RecoSimEnum rs, AR::BeforeAfterEnum ba, typename TrackObject>
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

  // function for checking if collision survives event cuts
  template <typename CollisionObject, typename TrackObject>
  bool SurviveEventCuts(CollisionObject collision, TrackObject tracks)
  {

    // Check if event survives event cuts, where we can get the values for the variables immediately
    if (!(AR::SurviveCut(cfgVX.value, collision.posX()) &&
          AR::SurviveCut(cfgVY.value, collision.posY()) &&
          AR::SurviveCut(cfgVZ.value, collision.posZ()) &&
          AR::SurviveCut(cfgVABS.value, AR::abs(collision.posX(), collision.posY(), collision.posZ())) &&
          AR::SurviveCut(cfgCEN.value, collision.centRun2V0M()) &&
          AR::SurviveCut(cfgMULNC.value, collision.numContrib()) &&
          AR::SurviveCut(cfgMULTPC.value, collision.multTPC()))) {
      return false;
    }

    int MultiplicityQ = 0.;
    float MultiplicityW = 0.;

    // only compute Multiplicity{Q,W} if the other checks pass, i.e. our flag is still set to true
    for (auto& track : tracks) {
      if (SurviveTrackCuts(track) == true) {
        MultiplicityQ += 1.;
        MultiplicityW += 1.;
      }
    }

    // at last, check if event also passes multiplicity cuts
    return AR::SurviveCut(cfgMULQ.value, MultiplicityQ) && AR::SurviveCut(cfgMULW.value, MultiplicityW);
  }

  // function for checking if track survices trach cuts
  template <typename TrackObject>
  bool SurviveTrackCuts(TrackObject track)
  {
    // if all SurviveCut return true, the function will return true
    // if at least one fails, it will return false
    return AR::SurviveCut(cfgPT.value, track.pt()) &&
           AR::SurviveCut(cfgPHI.value, track.phi()) &&
           AR::SurviveCut(cfgETA.value, track.eta()) &&
           AR::SurviveCut(cfgCHARGE.value, track.sign()) &&
           AR::SurviveCut(cfgDCAZ.value, track.dcaZ()) &&
           AR::SurviveCut(cfgDCAXY.value, track.dcaXY()) &&
           AR::SurviveCut(cfgTPCCLUSTERS.value, track.tpcNClsFound()) &&
           AR::SurviveCut(cfgTPCCROSSEDROWS.value, track.tpcNClsCrossedRows()) &&
           AR::SurviveCut(cfgTPCCHI2.value, track.tpcChi2NCl()) &&
           AR::SurviveCut(cfgITSCLUSTERS.value, track.itsNCls());
  }

  // Calculate all Q-vectors
  void CalculateQvectors()
  {
    // Make sure all Q-vectors are initially zero
    for (int h = 0; h < AR::MaxHarmonic; h++) {
      for (int p = 0; p < AR::MaxPower; p++) {
        QVectors[h][p] = TComplex(0., 0.);
      }
    }

    // Calculate Q-vectors for available angles and weights
    double dPhi = 0.;
    double wPhi = 1.;         // particle weight
    double wPhiToPowerP = 1.; // particle weight raised to power p
    for (std::size_t i = 0; i < Angles.size(); i++) {
      dPhi = Angles.at(i);
      wPhi = Weights.at(i);
      for (int h = 0; h < AR::MaxHarmonic; h++) {
        for (int p = 0; p < AR::MaxPower; p++) {
          wPhiToPowerP = TMath::Power(wPhi, p);
          QVectors[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi),
                                     wPhiToPowerP * TMath::Sin(h * dPhi));
        }
      }
    }
  }

  using CollisionsInstance = soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::Mults>;
  using TracksInstance = soa::Join<aod::Tracks, aod::TracksDCA, aod::TracksExtra>;

  using CollisionsInstanceIterator = CollisionsInstance::iterator;
  // using TracksInstanceIterator = TracksInstance::iterator;

  void process(CollisionsInstanceIterator const& collision, TracksInstance const& tracks)
  {

    // clear angles and weights
    Angles.clear();
    Weights.clear();

    LOGF(info, "Process reconstructed event: %d", collision.index());

    FillEventControlHist<AR::kRECO, AR::kBEFORE, CollisionsInstanceIterator>(collision, registry);
    FillEventControlHistMul<AR::kRECO, AR::kBEFORE>(registry, collision.size(), collision.size());

    if (!SurviveEventCuts(collision, tracks)) {
      LOGF(info, "Event was CUT");
      return;
    }

    FillEventControlHist<AR::kRECO, AR::kAFTER, CollisionsInstanceIterator>(collision, registry);

    // loop over all tracks in the event
    for (auto const& track : tracks) {
      // fill track control histograms before track cut
      FillTrackControlHist<AR::kRECO, AR::kBEFORE, TracksInstance::iterator>(track, registry);
      if (!SurviveTrackCuts(track)) {
        continue;
      }
      // fill track control histograms after cut
      FillTrackControlHist<AR::kRECO, AR::kAFTER, TracksInstance::iterator>(track, registry);

      // fill angles into vector for processing
      Angles.push_back(track.phi());
      Weights.push_back(1.);
    }

    FillEventControlHistMul<AR::kRECO, AR::kAFTER>(registry, Angles.size(), std::accumulate(Weights.begin(), Weights.end(), 0.));

    // calculate qvectors from filled angles and weights
    CalculateQvectors();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiParticleCorrelationsARTask>(cfgc)};
}
