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

#include <TLorentzVector.h>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "PWGLF/Utils/rsnOutput.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phianalysisTHnSparse {

  SliceCache cache;

  Configurable<bool> produceTrue{"produce-true", false, "Produce True and Gen histograms."};
  Configurable<int> verboselevel{"verbose-level", 0, "Verbose level"};
  Configurable<int> refresh{"print-refresh", 0, "Freqency of print event information."};
  Configurable<int> refresh_index{"print-refresh-index", 0, "Freqency of print event information index."};
  Configurable<bool> ignorezeroevent{"ignore-zero-event", true, "Flag if zero event is skipped"};
  Configurable<float> tpcnSigmaPos{"tpc-ns-pos", 3.0f, "TPC NSigma cut of the positive particle."};
  Configurable<float> tpcnSigmaNeg{"tpc-ns-neg", 3.0f, "TPC NSigma cut of the negative particle."};
  Configurable<int> dautherPos{"dauther-type-pos", 3, "Particle type of the positive dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> dautherNeg{"dauther-type-neg", 3, "Particle type of the negative dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<float> zVertexCut{"zvertex-cut", 10.0f, "Z vertex range."};
  Configurable<float> rapidityCut{"rapidity-max", 0.5, "Rapidity cut."};
  Configurable<int> motherPDG{"mother-pdg", 333, "PDG code of mother particle."};
  Configurable<int> dautherPosPDG{"dauther-pdg-pos", 321, "PDG code of positive dauther particle."};
  Configurable<int> dautherNegPDG{"dauther-pdg-neg", 321, "PDG code of negative dauther particle."};

  Configurable<std::vector<std::string>> sparseAxes{"sparse-axes", std::vector<std::string>{o2::analysis::rsn::PariAxis::names}, "Axes."};

  ConfigurableAxis invaxis{"invAxis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"ptAxis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis posZ{"posZ", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicityAxis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidityAxis", {10., -1.0 * rapidityCut, rapidityCut}, "Rapidity axis binning."};
  ConfigurableAxis nsigmaaxisPos{"nsigmaAxisPos", {1, 0., tpcnSigmaPos}, "NSigma of positive particle axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxisNeg{"nsigmaAxisNeg", {1, 0., tpcnSigmaNeg}, "NSigma of negative particle axis binning in THnSparse."};

  // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
  float massPos = o2::track::PID::getMass(dautherPos);
  float massNeg = o2::track::PID::getMass(dautherNeg);

  AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
  AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
  AxisSpec mAxis = {multiplicityaxis, "N", "mu"};
  AxisSpec yAxis = {rapidityaxis, "y", "y"};
  AxisSpec nsigmatrackaxisPos = {nsigmaaxisPos, fmt::format("nSigma of positive particle ({})", massPos), "ns1"};
  AxisSpec nsigmatrackaxisNeg = {nsigmaaxisNeg, fmt::format("nSigma of negative particle ({})", massNeg), "ns2"};

  // All axes has to have same order as defined enum o2::analysis::rsn::PairAxisType (name from AxisSpec is taken to compare in o2::analysis::rsn::Output::init())
  std::vector<AxisSpec> allAxes = {invAxis, ptAxis, mAxis, nsigmatrackaxisPos, nsigmatrackaxisNeg, yAxis};
  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg;
  float multiplicity;
  float multiplicityMC;
  int n = 0;
  double* pointPair = nullptr;
  double* pointEvent = nullptr;
  TLorentzVector d1, d2, mother;

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < zVertexCut);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullKa>;

  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::McCollisionLabels>>;
  using TrackCandidatesMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullKa, aod::McTrackLabels>;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigmaPos);
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigmaNeg);

  Partition<TrackCandidatesMC> positiveMC = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigmaPos);
  Partition<TrackCandidatesMC> negativeMC = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < tpcnSigmaNeg);

  void init(o2::framework::InitContext&)
  {
    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    pointEvent = new double[static_cast<int>(o2::analysis::rsn::EventType::all)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, produceTrue, &registry);
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
    d1.SetXYZM(track1.px(), track1.py(), track1.pz(), massPos);
    d2.SetXYZM(track2.px(), track2.py(), track2.pz(), massNeg);
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

    if (std::abs(collision.posZ()) > zVertexCut)
      return;

    pointEvent[0] = collision.posZ();
    rsnOutput->fill(o2::analysis::rsn::EventType::zvertex, pointEvent);

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthers, negDauthers))) {

      if (!selectedTrack(track1))

        continue;
      if (!selectedTrack(track2))
        continue;

      if (!selectedPair(mother, track1, track2))
        continue;

      if (verboselevel > 1)
        LOGF(info, "Unlike-sign: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = mother.Mag();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = mother.Pt();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = multiplicity;
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = std::abs(track1.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = std::abs(track2.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = mother.Rapidity();

      rsnOutput->fillUnlike(pointPair);
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

      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = mother.Mag();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = mother.Pt();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = multiplicity;
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = std::abs(track1.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = std::abs(track2.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = mother.Rapidity();

      rsnOutput->fillLikepp(pointPair);
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

      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = mother.Mag();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = mother.Pt();
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = multiplicity;
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = std::abs(track1.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = std::abs(track2.tpcNSigmaKa());
      pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = mother.Rapidity();

      rsnOutput->fillLikemm(pointPair);
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processData, "Process Event for Data", true);

  void processTrue(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    if (!produceTrue)
      return;

    auto posDauthersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }

    if (std::abs(collision.posZ()) > zVertexCut)
      return;

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

      if (!(track1PDG == dautherPosPDG && track2PDG == dautherNegPDG)) {
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

          if (std::abs(mothertrack1.pdgCode()) != motherPDG)
            continue;

          n++;
          if (verboselevel > 1)
            LOGF(info, "True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());

          if (!selectedPair(mother, mctrack1, mctrack2))
            continue;

          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = mother.Mag();
          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = mother.Pt();
          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = multiplicity;
          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = std::abs(track1.tpcNSigmaKa());
          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = std::abs(track2.tpcNSigmaKa());
          pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = mother.Rapidity();

          rsnOutput->fillUnliketrue(pointPair);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  int numberofEntries = 0;

  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    if (!produceTrue)
      return;

    if (std::abs(mcCollision.posZ()) > zVertexCut)
      return;

    int nuberofPhi = 0;

    for (auto& particle : mcParticles) {
      if (std::abs(particle.y()) > rapidityCut)
        continue;

      if (particle.pdgCode() == motherPDG) {
        auto daughters = particle.daughters_as<aod::McParticles>();
        if (daughters.size() != 2)
          continue;

        auto daup = false;
        auto daun = false;

        for (auto& dau : daughters) {
          if (!dau.isPhysicalPrimary())
            continue;

          if (dau.pdgCode() == dautherPosPDG) {
            daup = true;
            d1.SetXYZM(dau.px(), dau.py(), dau.pz(), massPos);
          } else if (dau.pdgCode() == -dautherNegPDG) {
            daun = true;
            d2.SetXYZM(dau.px(), dau.py(), dau.pz(), massNeg);
          }
        }
        if (!daup && !daun)
          continue;

        mother = d1 + d2;

        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = mother.Mag();
        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = mother.Pt();
        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = multiplicityMC;
        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = std::abs(tpcnSigmaPos / 2.0);
        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = std::abs(tpcnSigmaNeg / 2.0);
        pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = mother.Rapidity();

        rsnOutput->fillUnlikegen(pointPair);

        nuberofPhi++;
        numberofEntries++;

        if (verboselevel > 1)
          LOGF(info, "Gen:  %d, #Phi =%d, mother=%d (%ld), Inv.mass:%f, Pt= %f", numberofEntries, nuberofPhi, particle.pdgCode(), particle.globalIndex(), mother.Mag(), mother.Pt());
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
