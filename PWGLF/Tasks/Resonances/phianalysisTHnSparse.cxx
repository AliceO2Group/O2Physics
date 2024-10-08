// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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
/// \since April 3, 2024

#include <TLorentzVector.h>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
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

  struct : ConfigurableGroup {
    Configurable<bool> QA{"produce-qa", false, "Produce QA histograms."};
    Configurable<bool> True{"produce-true", false, "Produce True and Gen histograms."};
    Configurable<bool> Likesign{"produce-likesign", false, "Produce Like sign histograms."};
    Configurable<bool> eventMixing{"produce-event-mixing", false, "Produce Event Mixing histograms."};
  } produce;

  struct : ConfigurableGroup {
    Configurable<int> verboselevel{"verbose-level", 0, "Verbose level"};
    Configurable<int> refresh{"print-refresh", 0, "Freqency of print event information."};
    Configurable<int> refresh_index{"print-refresh-index", 0, "Freqency of print event information index."};
  } verbose;

  Configurable<int> dautherPos{"dauther-type-pos", 3, "Particle type of the positive dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> dautherNeg{"dauther-type-neg", 3, "Particle type of the negative dauther according to ReconstructionDataFormats/PID.h (Default = Kaon)"};
  Configurable<int> motherPDG{"mother-pdg", 333, "PDG code of mother particle."};
  Configurable<int> dautherPosPDG{"dauther-pdg-pos", 321, "PDG code of positive dauther particle."};
  Configurable<int> dautherNegPDG{"dauther-pdg-neg", 321, "PDG code of negative dauther particle."};

  struct : ConfigurableGroup {
    Configurable<float> tpcnSigmaPos{"tpc-ns-pos", 3.0f, "TPC NSigma cut of the positive particle."};
    Configurable<float> tpcnSigmaNeg{"tpc-ns-neg", 3.0f, "TPC NSigma cut of the negative particle."};
    Configurable<float> vZ{"zvertex-cut", 10.0f, "Z vertex range."};
    Configurable<float> y{"rapidity-cut", 0.5, "Rapidity cut (maximum)."};
    Configurable<float> etatrack{"eta-track-cut", 0.8, "Eta cut for track."};
    Configurable<float> etapair{"eta-pair-cut", 0.8, "Eta cut for pair."};
    Configurable<float> pt{"pt-cut", 0.15f, "Cut: Minimal value of tracks pt."};
    Configurable<float> dcaXY{"dcaXY-cut", 1.0f, "Cut: Maximal value of tracks DCA XY."};
    Configurable<float> dcaZ{"dcaZ-cut", 1.0f, "Cut: Maximal value of tracks DCA Z."};
    Configurable<int> tpcNClsFound{"tpc-ncl-found-cut", 70, "Cut: Minimal value of found TPC clasters"};
  } cut;

  Configurable<std::vector<std::string>> sparseAxes{"sparse-axes", std::vector<std::string>{o2::analysis::rsn::PairAxis::names}, "Axes."};
  Configurable<std::vector<std::string>> sysAxes{"systematics-axes", std::vector<std::string>{o2::analysis::rsn::SystematicsAxis::names}, "Axes."};

  ConfigurableAxis invaxis{"inv-axis", {130, 0.97, 1.1}, "Invariant mass axis binning."};
  ConfigurableAxis ptaxis{"pt-axis", {20, 0., 20.}, "Pt axis binning."};
  ConfigurableAxis vzaxis{"vz-axis", {40, -20., 20.}, "Z vertex position axis binning."};
  ConfigurableAxis multiplicityaxis{"multiplicity-axis", {50, 0., 5000.}, "Multiplicity axis binning."};
  ConfigurableAxis centralityaxis{"centrality-axis", {20, 0., 100.}, "Centrality axis binning."};
  ConfigurableAxis etaaxis{"eta-axis", {16., -1.0 * static_cast<float>(cut.etatrack), static_cast<float>(cut.etatrack)}, "Pseudorapidity axis binning."};
  ConfigurableAxis rapidityaxis{"rapidity-axis", {10., -1.0 * static_cast<float>(cut.y), static_cast<float>(cut.y)}, "Rapidity axis binning."};
  ConfigurableAxis nsigmaaxisPos{"nsigma-pos-axis", {1, 0., static_cast<float>(cut.tpcnSigmaPos)}, "NSigma of positive particle axis binning in THnSparse."};
  ConfigurableAxis nsigmaaxisNeg{"nsigma-neg-axis", {1, 0., static_cast<float>(cut.tpcnSigmaNeg)}, "NSigma of negative particle axis binning in THnSparse."};

  // mixing
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  Configurable<int> numberofMixedEvents{"number-of-mixed-events", 5, "Number of events that should be mixed."};
  ConfigurableAxis axisVertexMixing{"vertex-axis-mixing", {5, -10, 10}, "Z vertex axis for bin"};
  ConfigurableAxis axisMultiplicityMixing{"multiplicity-axis-mixing", {5, 0, 5000}, "TPC multiplicity for bin"};

  // defined in DataFormats/Reconstruction/include/ReconstructionDataFormats/PID.h
  float massPos = o2::track::PID::getMass(dautherPos);
  float massNeg = o2::track::PID::getMass(dautherNeg);

  // Axes specifications
  AxisSpec posZaxis = {400, -20., 20., "V_z (cm)"};
  AxisSpec dcaXYaxis = {120, -3.0, 3.0, "DCA_{xy} (cm)"};
  AxisSpec dcaZaxis = {120, -3.0, 3.0, "DCA_{z} (cm)"};

  HistogramRegistry registry{"registry"};
  o2::analysis::rsn::Output* rsnOutput = nullptr;

  Service<o2::framework::O2DatabasePDG> pdg;

  int n = 0;
  double* pointPair = nullptr;
  double* pointSys = nullptr;
  TLorentzVector d1, d2, mother;
  bool QA = false;
  float etapair = 0.8;
  float tpcnSigmaPos = 3.0f;
  float tpcnSigmaNeg = 3.0f;
  int tpcNClsFound = 70;

  Filter triggerFilter = (o2::aod::evsel::sel8 == true);
  Filter vtxFilter = (nabs(o2::aod::collision::posZ) < static_cast<float>(cut.vZ));

  Filter ptFilter = nabs(aod::track::pt) > static_cast<float>(cut.pt);
  Filter etaFilter = nabs(aod::track::eta) < static_cast<float>(cut.etatrack);
  Filter dcaFilter = (nabs(o2::aod::track::dcaXY) < static_cast<float>(cut.dcaXY)) && (nabs(o2::aod::track::dcaZ) < static_cast<float>(cut.dcaZ));

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>>;
  using EventCandidate = EventCandidates::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa>>;

  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::McCollisionLabels>>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::McTrackLabels>>;

  Partition<TrackCandidates> positive = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaPos)));
  Partition<TrackCandidates> negative = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaNeg)));

  Partition<TrackCandidatesMC> positiveMC = (aod::track::signed1Pt > 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaPos)));
  Partition<TrackCandidatesMC> negativeMC = (aod::track::signed1Pt < 0.0f) && (nabs(o2::aod::pidtpc::tpcNSigmaKa) < std::abs(static_cast<float>(cut.tpcnSigmaNeg)));

  void init(o2::framework::InitContext&)
  {
    // Sparse axes
    AxisSpec invAxis = {invaxis, "Inv. mass (GeV/c^{2})", "im"};
    AxisSpec ptAxis = {ptaxis, "p_{T} (GeV/c)", "pt"};
    AxisSpec muAxis = {multiplicityaxis, "N", "mu"};
    AxisSpec mumAxis = {multiplicityaxis, "N", "mum"};
    AxisSpec ceAxis = {centralityaxis, "N", "ce"};
    AxisSpec cemAxis = {centralityaxis, "N", "cem"};
    AxisSpec etaAxis = {etaaxis, "#eta", "eta"};
    AxisSpec yAxis = {rapidityaxis, "y", "y"};
    AxisSpec nsAxisPos = {nsigmaaxisPos, fmt::format("nSigma of positive particle ({})", massPos), "ns1"};
    AxisSpec nsAxisNeg = {nsigmaaxisNeg, fmt::format("nSigma of negative particle ({})", massNeg), "ns2"};
    AxisSpec vzAxis = {vzaxis, "V_{z} (cm)", "vz"};
    AxisSpec vzmAxis = {vzaxis, "V_{z} (cm)", "vzm"};

    // Systematics axes
    std::vector<double> tpcNClsFound_bins = {65., 66., 67., 68., 69., 70., 71., 72., 73.};
    AxisSpec tpcNClsFoundAxis = {tpcNClsFound_bins, "TPC NCl", "ncl"};

    // All axes has to have same order as defined enum o2::analysis::rsn::PairAxisType (name from AxisSpec is taken to compare in o2::analysis::rsn::Output::init())
    std::vector<AxisSpec> allAxes = {invAxis, ptAxis, muAxis, ceAxis, nsAxisPos, nsAxisNeg, etaAxis, yAxis, vzAxis, mumAxis, cemAxis, vzmAxis};
    std::vector<AxisSpec> allAxes_sys = {tpcNClsFoundAxis};

    QA = static_cast<bool>(produce.QA);
    etapair = static_cast<float>(cut.etapair);
    tpcnSigmaPos = static_cast<float>(cut.tpcnSigmaPos);
    tpcnSigmaNeg = static_cast<float>(cut.tpcnSigmaNeg);
    tpcNClsFound = static_cast<int>(cut.tpcNClsFound);

    pointPair = new double[static_cast<int>(o2::analysis::rsn::PairAxisType::unknown)];
    pointSys = new double[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::unknown)];
    rsnOutput = new o2::analysis::rsn::OutputSparse();
    rsnOutput->init(sparseAxes, allAxes, sysAxes, allAxes_sys, static_cast<bool>(produce.True), static_cast<bool>(produce.eventMixing), static_cast<bool>(produce.Likesign), &registry);

    if (QA) {
      // Event QA
      registry.add("QAEvent/hVtxZ", "", kTH1F, {posZaxis});
      registry.add("QAEvent/hCent", "", kTH1F, {{20, 0., 100.}});
      registry.add("QAEvent/hMult", "", kTH1F, {{50, 0., 10000.}});
      registry.add("QAEvent/s4Size", "", kTHnSparseF, {{30, 0., 30.}, {30, 0., 30.}, muAxis, vzAxis});
      registry.add("QAEvent/s2Mult_Vz", "", kTHnSparseF, {axisMultiplicityMixing, axisVertexMixing});
      // Track QA
      registry.add("QATrack/unlikepm/beforeSelection/hTrack1pt", "", kTH1F, {ptAxis});
      registry.add("QATrack/unlikepm/beforeSelection/hTrackDCAxy", "", kTH1F, {dcaXYaxis});
      registry.add("QATrack/unlikepm/beforeSelection/hTrackDCAz", "", kTH1F, {dcaZaxis});
      registry.add("QATrack/unlikepm/beforeSelection/hTrack1eta", "", kTH1F, {{100, -1.0, 1.0}});
      registry.add("QATrack/unlikepm/beforeSelection/hTrack1tpcNClsFound", "", kTH1F, {{110, 50., 160.}});

      registry.add("QATrack/unlikepm/afterSelection/hTrack1pt", "", kTH1F, {ptAxis});
      registry.add("QATrack/unlikepm/afterSelection/hTrackDCAxy", "", kTH1F, {dcaXYaxis});
      registry.add("QATrack/unlikepm/afterSelection/hTrackDCAz", "", kTH1F, {dcaZaxis});
      registry.add("QATrack/unlikepm/afterSelection/hTrack1eta", "", kTH1F, {{100, -1.0, 1.0}});
      registry.add("QATrack/unlikepm/afterSelection/hTrack1tpcNClsFound", "", kTH1F, {{110, 50., 160.}});

      registry.add("QATrack/unlikepm/TPCPID/h2TracknSigma", "", kTH2F, {{120, -6, 6}, {120, -6, 6}});

      registry.add("QAMC/hInvMassTrueFalse", "", kTH1F, {invAxis});

      // Mixing QA
      registry.add("QAMixing/s4Mult_Vz", "", kTHnSparseF, {axisMultiplicityMixing, axisMultiplicityMixing, axisVertexMixing, axisVertexMixing});
      registry.add("QAMixing/s2Mult_Vz", "", kTHnSparseF, {axisMultiplicityMixing, axisVertexMixing});
    }

    pointSys[static_cast<int>(o2::analysis::rsn::SystematicsAxisType::ncl)] = tpcNClsFound;
    rsnOutput->fillSystematics(pointSys);
  }

  template <typename T>
  bool selectedTrack(const T& track)
  {
    if (track.tpcNClsFound() < tpcNClsFound)
      return false;
    if (!track.isPrimaryTrack())
      return false;
    if (!track.isPVContributor())
      return false;
    return true;
  }

  template <typename T>
  bool selectedPair(TLorentzVector& mother, const T& track1, const T& track2)
  {
    d1.SetXYZM(track1.px(), track1.py(), track1.pz(), massPos);
    d2.SetXYZM(track2.px(), track2.py(), track2.pz(), massNeg);
    mother = d1 + d2;

    if (std::abs(mother.Eta()) > etapair)
      return false;
    return true;
  }

  template <typename T>
  float GetMultiplicity(const T& collision)
  {
    float multiplicity = collision.multFT0C() + collision.multFT0A();
    return multiplicity;
  }
  template <typename T>
  float GetCentrality(const T& collision)
  {
    float centrality = collision.centFT0M();
    return centrality;
  }

  double* FillPointPair(double im, double pt, double mu, double ce, double ns1, double ns2, double eta, double y, double vz, double mum, double cem, double vzm)
  {
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::im)] = im;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::pt)] = pt;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mu)] = mu;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ce)] = ce;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns1)] = ns1;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::ns2)] = ns2;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::eta)] = eta;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::y)] = y;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::vz)] = vz;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::mum)] = mum;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::cem)] = cem;
    pointPair[static_cast<int>(o2::analysis::rsn::PairAxisType::vzm)] = vzm;

    return pointPair;
  }

  void processData(EventCandidate const& collision, TrackCandidates const& /*tracks*/)
  {
    auto posDauthers = positive->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthers = negative->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (QA) {
      registry.fill(HIST("QAEvent/s4Size"), posDauthers.size(), negDauthers.size(), GetMultiplicity(collision), collision.posZ());
      registry.fill(HIST("QAEvent/s2Mult_Vz"), GetMultiplicity(collision), collision.posZ());
    }

    if (static_cast<int>(verbose.verboselevel) > 0 && static_cast<int>(verbose.refresh) > 0 && collision.globalIndex() % static_cast<int>(verbose.refresh) == static_cast<int>(verbose.refresh_index))
      LOGF(info, "pos=%lld neg=%lld, Z vertex position: %f [cm], %d, mult:%f.0", posDauthers.size(), negDauthers.size(), collision.posZ(),
           collision.globalIndex(), GetMultiplicity(collision));

    if (QA) {
      registry.fill(HIST("QAEvent/hVtxZ"), collision.posZ());
      registry.fill(HIST("QAEvent/hMult"), GetMultiplicity(collision));
      registry.fill(HIST("QAEvent/hCent"), GetCentrality(collision));
    }

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthers, negDauthers))) {
      if (QA) {
        registry.fill(HIST("QATrack/unlikepm/beforeSelection/hTrack1pt"), track1.pt());
        registry.fill(HIST("QATrack/unlikepm/beforeSelection/hTrackDCAxy"), track1.dcaXY());
        registry.fill(HIST("QATrack/unlikepm/beforeSelection/hTrackDCAz"), track1.dcaZ());
        registry.fill(HIST("QATrack/unlikepm/beforeSelection/hTrack1eta"), track1.eta());
        registry.fill(HIST("QATrack/unlikepm/beforeSelection/hTrack1tpcNClsFound"), track1.tpcNClsFound());
      }

      if (!selectedTrack(track1))
        continue;

      if (!selectedTrack(track2))
        continue;

      if (QA) {
        registry.fill(HIST("QATrack/unlikepm/afterSelection/hTrack1pt"), track1.pt());
        registry.fill(HIST("QATrack/unlikepm/afterSelection/hTrackDCAxy"), track1.dcaXY());
        registry.fill(HIST("QATrack/unlikepm/afterSelection/hTrackDCAz"), track1.dcaZ());
        registry.fill(HIST("QATrack/unlikepm/afterSelection/hTrack1eta"), track1.eta());
        registry.fill(HIST("QATrack/unlikepm/afterSelection/hTrack1tpcNClsFound"), track1.tpcNClsFound());
      }

      if (!selectedPair(mother, track1, track2))
        continue;

      if (QA) {
        registry.fill(HIST("QATrack/unlikepm/TPCPID/h2TracknSigma"), track1.tpcNSigmaKa(), track2.tpcNSigmaKa());
      }

      pointPair = FillPointPair(mother.Mag(),
                                mother.Pt(),
                                GetMultiplicity(collision),
                                GetCentrality(collision),
                                (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                mother.Eta(),
                                mother.Rapidity(),
                                collision.posZ(),
                                0,
                                0,
                                0);
      rsnOutput->fillUnlikepm(pointPair);
    }

    if (static_cast<bool>(produce.Likesign)) {

      for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(posDauthers, posDauthers))) {
        if (!selectedTrack(track1))
          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign positive: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

        pointPair = FillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  GetMultiplicity(collision),
                                  GetCentrality(collision),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  collision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillLikepp(pointPair);
      }

      for (auto& [track1, track2] : combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(negDauthers, negDauthers))) {
        if (!selectedTrack(track1))
          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Like-sign negative: d1=%ld , d2=%ld , mother=%f", track1.globalIndex(), track2.globalIndex(), mother.Mag());

        pointPair = FillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  GetMultiplicity(collision),
                                  GetCentrality(collision),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  collision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillLikemm(pointPair);
      }
    }
  }
  PROCESS_SWITCH(phianalysisTHnSparse, processData, "Process Event for Data", true);

  void processTrue(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& /*tracks*/, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!static_cast<bool>(produce.True))
      return;

    auto posDauthersMC = positiveMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negDauthersMC = negativeMC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!collision.has_mcCollision()) {
      LOGF(warning, "No MC collision for this collision, skip...");
      return;
    }

    if (std::abs(collision.posZ()) > static_cast<float>(cut.vZ))
      return;

    for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersMC, negDauthersMC))) {

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
      n = 0;
      for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
        for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {

          if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
            continue;

          if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
            continue;

          if (std::abs(mothertrack1.y()) > static_cast<float>(cut.y))
            continue;

          if (std::abs(mothertrack2.y()) > static_cast<float>(cut.y))
            continue;

          if (std::abs(mothertrack1.pdgCode()) != motherPDG)
            continue;

          if (static_cast<int>(verbose.verboselevel) > 1) {
            LOGF(info, "True: %d, d1=%d (%ld), d2=%d (%ld), mother=%d (%ld)", n, mctrack1.pdgCode(), mctrack1.globalIndex(), mctrack2.pdgCode(), mctrack2.globalIndex(), mothertrack1.pdgCode(), mothertrack1.globalIndex());
            LOGF(info, "%d px: %f, py=%f, pz=%f, px: %f, py=%f, pz=%f", n, mctrack1.px(), mctrack1.py(), mctrack1.pz(), mctrack2.px(), mctrack2.py(), mctrack2.pz());
          }

          if (!selectedPair(mother, mctrack1, mctrack2))
            continue;

          if (n > 0) {
            if (QA)
              registry.fill(HIST("QAMC/hInvMassTrueFalse"), mother.Mag());
            continue;
          }

          pointPair = FillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    GetMultiplicity(collision),
                                    GetCentrality(collision),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    collision.posZ(),
                                    0,
                                    0,
                                    0);

          rsnOutput->fillUnliketrue(pointPair);
          n++;
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processTrue, "Process Event for MC reconstruction.", false);

  int numberofEntries = 0;

  void processGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {

    if (!static_cast<bool>(produce.True))
      return;

    if (std::abs(mcCollision.posZ()) > static_cast<float>(cut.vZ))
      return;

    int nuberofPhi = 0;

    for (auto& particle : mcParticles) {
      if (std::abs(particle.y()) > static_cast<float>(cut.y))
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

        pointPair = FillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  0,
                                  0,
                                  std::abs(static_cast<float>(cut.tpcnSigmaPos) / 2.0),
                                  std::abs(static_cast<float>(cut.tpcnSigmaNeg) / 2.0),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  mcCollision.posZ(),
                                  0,
                                  0,
                                  0);

        rsnOutput->fillUnlikegen(pointPair);

        nuberofPhi++;
        numberofEntries++;

        if (static_cast<int>(verbose.verboselevel) > 1)
          LOGF(info, "Gen:  %d, #Phi =%d, mother=%d (%ld), Inv.mass:%f, Pt= %f", numberofEntries, nuberofPhi, particle.pdgCode(), particle.globalIndex(), mother.Mag(), mother.Pt());
      }
    }
  }

  PROCESS_SWITCH(phianalysisTHnSparse, processGen, "Process generated.", false);

  int id;

  void processMixed(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    if (!static_cast<bool>(produce.eventMixing))
      return;
    BinningType binning{{axisVertexMixing, axisMultiplicityMixing}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binning, numberofMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [c1, tracks1, c2, tracks2] : pair) {

      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (QA) {
        if (!(id == c1.globalIndex())) {
          registry.fill(HIST("QAMixing/s2Mult_Vz"), GetMultiplicity(c1), c1.posZ());
          id = c1.globalIndex();
        }
      }

      auto posDauthersc1 = positive->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto posDauthersc2 = positive->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);
      auto negDauthersc1 = negative->sliceByCached(aod::track::collisionId, c1.globalIndex(), cache);
      auto negDauthersc2 = negative->sliceByCached(aod::track::collisionId, c2.globalIndex(), cache);

      if (QA)
        registry.fill(HIST("QAMixing/s4Mult_Vz"), GetMultiplicity(c1), GetMultiplicity(c2), c1.posZ(), c2.posZ());

      for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, negDauthersc2))) {

        if (!selectedTrack(track1))

          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        pointPair = FillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  GetMultiplicity(c1),
                                  GetCentrality(c1),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  c1.posZ(),
                                  GetMultiplicity(c2),
                                  GetCentrality(c2),
                                  c2.posZ());

        rsnOutput->fillMixingpm(pointPair);
      }

      if (static_cast<bool>(produce.Likesign)) {

        for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc1, posDauthersc2))) {

          if (!selectedTrack(track1))

            continue;
          if (!selectedTrack(track2))
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;

          pointPair = FillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    GetMultiplicity(c1),
                                    GetCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    GetMultiplicity(c2),
                                    GetCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingpp(pointPair);
        }

        for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(negDauthersc1, negDauthersc2))) {

          if (!selectedTrack(track1))

            continue;
          if (!selectedTrack(track2))
            continue;

          if (!selectedPair(mother, track1, track2))
            continue;
          pointPair = FillPointPair(mother.Mag(),
                                    mother.Pt(),
                                    GetMultiplicity(c1),
                                    GetCentrality(c1),
                                    (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                    (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                    mother.Eta(),
                                    mother.Rapidity(),
                                    c1.posZ(),
                                    GetMultiplicity(c2),
                                    GetCentrality(c2),
                                    c2.posZ());

          rsnOutput->fillMixingmm(pointPair);
        }
      }

      for (auto& [track1, track2] : combinations(o2::soa::CombinationsFullIndexPolicy(posDauthersc2, negDauthersc1))) {

        if (!selectedTrack(track1))

          continue;
        if (!selectedTrack(track2))
          continue;

        if (!selectedPair(mother, track1, track2))
          continue;

        pointPair = FillPointPair(mother.Mag(),
                                  mother.Pt(),
                                  GetMultiplicity(c1),
                                  GetCentrality(c1),
                                  (tpcnSigmaPos > 0) ? std::abs(track1.tpcNSigmaKa()) : track1.tpcNSigmaKa(),
                                  (tpcnSigmaNeg > 0) ? std::abs(track2.tpcNSigmaKa()) : track2.tpcNSigmaKa(),
                                  mother.Eta(),
                                  mother.Rapidity(),
                                  c1.posZ(),
                                  GetMultiplicity(c2),
                                  GetCentrality(c2),
                                  c2.posZ());

        rsnOutput->fillMixingmp(pointPair);
      }
    }
  }
  PROCESS_SWITCH(phianalysisTHnSparse, processMixed, "Process Mixing Event.", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisTHnSparse>(cfgc)};
}
