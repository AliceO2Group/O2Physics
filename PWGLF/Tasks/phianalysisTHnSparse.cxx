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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/PID.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phianalysisTHnSparse {

  SliceCache cache;

  Configurable<int> verboselevel{"verbose-level", 0, "Verbose level"};
  Configurable<int> refresh{"print-refresh", 0, "Freqency of print event information."};
  Configurable<int> refresh_index{"print-refresh-index", 0, "Freqency of print event information index."};
  Configurable<float> tpcnSigma1{"tpcnSigma1", 3.0f, "TPC NSigma cut of the first particle."};
  Configurable<bool> ignorezeroevent{"ignore-zero-event", true, "Flag if zero event is skipped"};
  Configurable<float> tpcnSigma2{"tpcnSigma2", 3.0f, "TPC NSigma cut of the second particle."};
  Configurable<int> dauther1{"dauther1", 3, "Particle type of the first dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> dauther2{"dauther2", 3, "Particle type of the second dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<float> zVertex{"zvertex", 10.0f, "Z vertex range."};
  Configurable<float> rapidityCut{"rapidity-max", 0.5, "Rapidity cut."};

  ConfigurableAxis invaxis{"invAxis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptAxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis mcposZ{"mcposZ", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicityAxis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidityAxis", {10., -1.0 * rapidityCut, rapidityCut}, "Rapidity axis binning."};
  ConfigurableAxis nsigmatrackaxis{"nsigmatrackaxis", {300, -15., 15.}, "NSigma axis binning."};
  ConfigurableAxis nsigmaaxis1{"nsigmaAxis1", {1, 0., tpcnSigma1}, "NSigma axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxis2{"nsigmaAxis2", {1, 0., tpcnSigma2}, "NSigma axis binning in THnSparse."};

  HistogramRegistry registry{"registry",
                             {{"hNsigmaPos", "hNsigmaPos", {HistType::kTH1F, {nsigmatrackaxis}}},
                              {"hNsigmaNeg", "hNsigmaNeg", {HistType::kTH1F, {nsigmatrackaxis}}},
                              {"motherGen", "motherGen", {HistType::kTH1F, {ptaxis}}},
                              {"motherTrue", "motherTrue", {HistType::kTH1F, {ptaxis}}},
                              {"motherBgr", "motherBgr", {HistType::kTH1F, {ptaxis}}},
                              {"mcTrueposZ", "mcTrueposZ", {HistType::kTH1F, {mcposZ}}},
                              {"mcGenposZ", "mcGenposZ", {HistType::kTH1F, {mcposZ}}}}};

  // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
  float mass1 = o2::track::PID::getMass(dauther1);
  float mass2 = o2::track::PID::getMass(dauther2);

  Service<o2::framework::O2DatabasePDG> pdg;

  float multiplicity;
  float multiplicityMC;
  int n = 0;

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < zVertex);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;

  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullKa>;

  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullKa, aod::McTrackLabels>;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma1);
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma2);

  Partition<TrackCandidatesMC> positiveMC = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma1);
  Partition<TrackCandidatesMC> negativeMC = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigma2);

  TLorentzVector d1, d2, mother;

  void init(o2::framework::InitContext&)
  {
    AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
    AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
    AxisSpec mAxis = {multiplicityaxis, "N", "m"};
    AxisSpec yAxis = {rapidityaxis, "y", "y"};
    AxisSpec nsigmatrackaxis1 = {nsigmaaxis1, fmt::format("nSigma particle 1({})", mass1), "ns1"};
    AxisSpec nsigmatrackaxis2 = {nsigmaaxis1, fmt::format("nSigma particle 2({})", mass2), "ns2"};
    HistogramConfigSpec pairHisto({HistType::kTHnSparseF, {invAxis, ptAxis, mAxis, nsigmatrackaxis1, nsigmatrackaxis2, yAxis}});
    registry.add("unlike", "Unlike", pairHisto);
    registry.add("likep", "Likep", pairHisto);
    registry.add("liken", "Liken", pairHisto);

    registry.add("unlikeTrue", "UnlikeTrue", pairHisto);
    registry.add("unlikeGen", "UnlikeGen", pairHisto);
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

    if (std::abs(mother.Rapidity()) > 0.5)
      return false;

    return true;
  }

  void processData(EventCandidate const& collision, TrackCandidates const& tracks)
  {
    auto posDauthers = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthers = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (ignorezeroevent && collision.globalIndex() == 0) {
      if (verboselevel > 0)
        LOGF(info, "BAD pos=%lld neg=%lld, Z vertex position: %f [cm], %d, mult:%f.0", posDauthers.size(), negDauthers.size(), collision.posZ(),
             collision.globalIndex(), multiplicity);

      return;
    }

    multiplicity = collision.multFT0A() + collision.multFT0C();

    if (verboselevel > 0 && refresh > 0 && collision.globalIndex() % refresh == refresh_index)
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

      if (verboselevel > 1)
        LOGF(info, "Unlike-sign: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

      registry.fill(HIST("unlike"), mother.Mag(), mother.Pt(), multiplicity, std::abs(track1.tpcNSigmaKa()), std::abs(track2.tpcNSigmaKa()), mother.Rapidity());
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDauthers, posDauthers))) {

      if (!selectedTrack(track1))
        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;

      if (verboselevel > 1)
        LOGF(info, "Like-sign positive: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

      registry.fill(HIST("likep"), mother.Mag(), mother.Pt(), multiplicity, std::abs(track1.tpcNSigmaKa()), std::abs(track2.tpcNSigmaKa()), mother.Rapidity());
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negDauthers, negDauthers))) {

      if (!selectedTrack(track1))
        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;

      if (verboselevel > 1)
        LOGF(info, "Like-sign negative: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

      registry.fill(HIST("liken"), mother.Mag(), mother.Pt(), multiplicity, std::abs(track1.tpcNSigmaKa()), std::abs(track2.tpcNSigmaKa()), mother.Rapidity());
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processData, "Process Event for Data", true);

  void processTrue(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    auto posDauthersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }

    if (std::abs(collision.posZ()) < zVertex) {
      registry.fill(HIST("mcTrueposZ"), collision.posZ());

      multiplicityMC = collision.multFT0A() + collision.multFT0C();

      for (auto& [track1, track2] : combinations(o2::soa::CombinationsUpperIndexPolicy(posDauthersMC, negDauthersMC))) {

        if (!track1.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        if (!track2.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }

        if (!selectedTrack(track1))
          continue;

        if (!selectedTrack(track2))
          continue;

        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());

        if (!(track1PDG == 321 && track2PDG == 321)) {
          continue;
        }
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
              continue;

            if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
              continue;

            if (std::abs(mothertrack1.y()) > rapidityCut)
              continue;

            if (std::abs(mothertrack2.y()) > rapidityCut)
              continue;

            if (std::abs(mothertrack1.pdgCode()) != 333)
              continue;

            registry.fill(HIST("motherTrue"), mothertrack1.pt());
            n++;
            if (verboselevel > 1)
              LOGF(info, "True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());

            if (!selectedPair(mother, mctrack1, mctrack2))
              continue;
            registry.fill(HIST("unlikeTrue"), mother.Mag(), mother.Pt(), multiplicityMC, std::abs(track1.tpcNSigmaKa()), std::abs(track2.tpcNSigmaKa()), mother.Rapidity());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  int numberofEntries = 0;

  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) < zVertex) {
      registry.fill(HIST("mcGenposZ"), mcCollision.posZ());

      int nuberofPhi = 0;

      for (auto& particle : mcParticles) {
        if (std::abs(particle.y()) > rapidityCut)
          continue;

        if (particle.pdgCode() == 333) {
          auto daughters = particle.daughters_as<aod::McParticles>();
          if (daughters.size() != 2)
            continue;

          auto daup = false;
          auto daun = false;

          for (auto& dau : daughters) {
            if (!dau.isPhysicalPrimary())
              continue;

            if (dau.pdgCode() == +321) {
              daup = true;
              d1.SetXYZM(dau.px(), dau.py(), dau.pz(), mass1);
            } else if (dau.pdgCode() == -321) {
              daun = true;
              d2.SetXYZM(dau.px(), dau.py(), dau.pz(), mass2);
            }
          }
          if (!daup && !daun)
            continue;

          mother = d1 + d2;

          registry.fill(HIST("unlikeGen"), mother.Mag(), mother.Pt(), multiplicityMC, tpcnSigma1 / 2.0, tpcnSigma2 / 2.0, mother.Rapidity());
          registry.fill(HIST("motherGen"), particle.pt());

          nuberofPhi++;
          numberofEntries++;

          if (verboselevel > 1)
            LOGF(info, "Gen:  %d, #Phi =%d, mother=%d (%ld), Inv.mass:%f, Pt= %f", numberofEntries, nuberofPhi, particle.pdgCode(), particle.globalIndex(), mother.Mag(), mother.Pt());
        } else {
          registry.fill(HIST("motherBgr"), particle.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processGen, "Process generated.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisTHnSparse>(cfgc)};
}
