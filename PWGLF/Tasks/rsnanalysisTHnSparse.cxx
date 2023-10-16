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
/// \author Veronika Barbasova (veronika.barbasova@cern.ch)
/// \since October 12, 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include <TLorentzVector.h>
#include "ReconstructionDataFormats/PID.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct rsnanalysisTHnSparse {

  SliceCache cache;

  Configurable<bool> test{"test", false, "Running in test mode."};
  Configurable<int> refresh{"print-refresh", 0, "Freqency of print event information."};
  Configurable<int> refresh_index{"print-refresh-index", 1, "Freqency of print event information index."};
  Configurable<float> tpcnSigma1{"tpcnSigma1", 3.0f, "TPC NSigma cut of the first particle."};
  Configurable<float> tpcnSigma2{"tpcnSigma2", 3.0f, "TPC NSigma cut of the second particle."};
  Configurable<int> dauther1{"dauther1", 3, "Particle type of the first dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> dauther2{"dauther", 3, "Particle type of the second dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<float> vertexfilter{"vertex-filter", 10.0f, "Z vertex range."};

  ConfigurableAxis invaxis{"invAxis", {200, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptAxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis maxis{"mAxis", {200, 0., 2000.}, "Multiplicity axis binning."};
  ConfigurableAxis nsigmaaxis{"nsigma-axis", {300, -15., 15.}, "NSigma axis binning."};

  HistogramRegistry registry{"registry",
                             {{"hNsigmaPos", "hNsigmaPos", {HistType::kTH1F, {nsigmaaxis}}},
                              {"hNsigmaNeg", "hNsigmaNeg", {HistType::kTH1F, {nsigmaaxis}}}}};

  // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
  float mass1 = o2::track::PID::getMass(dauther1);
  float mass2 = o2::track::PID::getMass(dauther2);

  float multiplicity;

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < vertexfilter);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullKa>;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma1);
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma2);

  TLorentzVector d1, d2, mother;

  void init(o2::framework::InitContext&)
  {
    if (test)
      LOGF(info, "Running test");

    AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
    AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
    AxisSpec mAxis = {maxis, "N", "m"};
    HistogramConfigSpec Hist({HistType::kTHnSparseF, {invAxis, ptAxis, mAxis}});
    registry.add("unlike", "Unlike", Hist);
    registry.add("likep", "Likep", Hist);
    registry.add("liken", "Liken", Hist);
  }

  template <typename T>
  bool selectedTrack(const T& track)
  {
    if (!track.isPrimaryTrack())
      return false;
    return true;
  }

  template <typename T>
  bool selectedPair(TLorentzVector& mother, const T& track1, const T& track2)
  {
    d1.SetXYZM(track1.px(), track1.py(), track1.pz(), mass1);
    d2.SetXYZM(track2.px(), track2.py(), track2.pz(), mass2);
    mother = d1 + d2;

    return true;
  }

  void process(EventCandidate const& collision, TrackCandidates const& tracks)
  {
    if (test && collision.globalIndex() != 1)
      return;

    auto posDauthers = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthers = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (collision.globalIndex() == 0) {
      if (test) {
        LOGF(info, "BAD pos=%lld neg=%lld, Z vertex position: %f [cm], %d, mult:%f.0", posDauthers.size(), negDauthers.size(), collision.posZ(),
             collision.globalIndex(), multiplicity);
      }
      return;
    }

    multiplicity = collision.multFT0A() + collision.multFT0C();

    if (refresh > 0 && collision.globalIndex() % refresh == refresh_index)
      LOGF(info, "pos=%lld neg=%lld, Z vertex position: %f [cm], %d, mult:%f.0", posDauthers.size(), negDauthers.size(), collision.posZ(),
           collision.globalIndex(), multiplicity);

    for (const auto& trk : posDauthers) {
      registry.fill(HIST("hNsigmaPos"), trk.tpcNSigmaKa());
    }

    for (const auto& trk : negDauthers) {
      registry.fill(HIST("hNsigmaNeg"), trk.tpcNSigmaKa());
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(posDauthers, negDauthers))) {

      if (!selectedTrack(track1))
        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;
      registry.fill(HIST("unlike"), mother.Mag(), mother.Pt(), multiplicity);
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDauthers, posDauthers))) {

      if (!selectedTrack(track1))
        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;
      registry.fill(HIST("likep"), mother.Mag(), mother.Pt(), multiplicity);
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negDauthers, negDauthers))) {

      if (!selectedTrack(track1))
        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;
      registry.fill(HIST("liken"), mother.Mag(), mother.Pt(), multiplicity);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<rsnanalysisTHnSparse>(cfgc)};
}
