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

/// \file rhoanalysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Nasir Mehdi Malik

#include "RecoDecay.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <array>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct rhoanalysis {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<int> cMultiplicityType{"cMultiplicityType", 2, " Set 1 for MultFT0M, 2 for centFT0M"};
  Configurable<bool> isUnlike{"isUnlike", true, "Use unlike charged tracks"};
  Configurable<bool> isLike{"isLike", true, "use same charge tracks"};
  Configurable<bool> isMixed{"isMixed", false, "Use mixed events"};
  Configurable<bool> isMC{"isMC", false, "Run MC"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec ptAxis = {140, 0., 7., "#it{p}_{T} (GeV/c)"};
    AxisSpec dcaxyAxis = {100, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {100, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};

    AxisSpec multAxis = {2000, 0., 2000, "mult"};
    if (cMultiplicityType == 2)
      multAxis = {100, 0., 100., "multiplicity percentile"};
    AxisSpec ptAxiss = {100, 0.0f, 10.0f, "#it{p}_{T}  (GeV/#it{c})"};
    AxisSpec invmassAxis = {180, 0.2, 2.0, "{M}_{#pi #pi}} (GeV/#it{c}^2)"};
    AxisSpec rapidityAxis = {100, -1.0, 1.0, "y"};
    if (isUnlike) {
      histos.add("hmultFT0M", "Multipicity FT0M", kTH1F, {multAxis});

      histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("ptpi", "pi pT distribution ", kTH1F, {ptAxis});
      histos.add("dcaXY", "dcaXY ", kTH1F, {dcaxyAxis});
      histos.add("dcaZ", "dcaZ ", kTH1F, {dcazAxis});
      histos.add("hNsigmaPionTPC", "NsigmaPion TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});

      histos.add("hNsigmaPionTPCvspT_beforecut", "NsigmaPion TPC distribution before cut Vs pT", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTPCvspT_aftercut", "NsigmaPion TPC distribution after cuts Vs pT", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});
      histos.add("hNsigmaPionTOF", "NsigmaPion TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});

      histos.add("hrhoUnlike", "Invariant mass of rho", kTHnSparseF, {invmassAxis, ptAxiss, multAxis});
    }

    if (isLike) {
      histos.add("hrhoLikepp", "Invariant mass of rho like sign pp", kTHnSparseF, {invmassAxis, ptAxiss, multAxis});
      histos.add("hrhoLikemm", "Invariant mass of rho like sign mm", kTHnSparseF, {invmassAxis, ptAxiss, multAxis});
    }
    if (isMixed) {
      histos.add("hrhoMixed", "Invariant mass of rho mixed event", kTHnSparseF, {invmassAxis, ptAxiss, multAxis});
    }

    if (isMC) {

      histos.add("hNsigmaPionTPCvspT", "NsigmaPion TPC distribution Vs pT", kTH2F, {ptAxis, {100, -10.0f, 10.0f}});
      histos.add("hrhoRec", "Invariant mass of rho mc", kTHnSparseF, {invmassAxis, ptAxiss, multAxis});

      histos.add("hrhoGen", "Gen mass of rho mc", kTHnSparseF, {ptAxiss, {30, 100, 130, "pdg code"}});
    }
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC)
      return true;
    else
      return false;
  }

  template <typename T>
  float Multiplicity(const T& event)
  {
    float multiplicity = 0.0;
    if (cMultiplicityType == 2)
      multiplicity = event.centFT0M();
    else
      multiplicity = event.multFT0M();

    return multiplicity;
  }

  enum class FillMode {
    Unlike,
    Mixed,
    LikeSign
  };

  std::array<float, 3> pvec0;
  std::array<float, 3> pvec1;

  double mass{0.};
  double recMass{0.};
  double pT{0.};
  double rapidity;

  Service<o2::framework::O2DatabasePDG> Ipdg;

  template <typename T1, typename T2>
  void FillHistogam(const T1& candidate1, const T2& candidate2, FillMode fillMode, float multiplicity = 0.0, float massd1 = 0., float massd2 = 0.)
  {

    auto _px1 = candidate1.px();
    auto _py1 = candidate1.py();
    auto _pz1 = candidate1.pz();
    auto _px2 = candidate2.px();
    auto _py2 = candidate2.py();
    auto _pz2 = candidate2.pz();

    pvec0 = std::array{_px1, _py1, _pz1};
    pvec1 = std::array{_px2, _py2, _pz2};

    auto arrMom = std::array{pvec0, pvec1};
    mass = RecoDecay::m(arrMom, std::array{massd1, massd2});

    rapidity = RecoDecay::y(std::array{_px1 + _px2, _py1 + _py2, _pz1 + _pz2}, mass);
    if (std::abs(rapidity) >= 0.5)
      return;
    pT = RecoDecay::pt(std::array{_px1 + _px2, _py1 + _py2});

    switch (fillMode) {
      case FillMode::LikeSign:
        if (candidate1.sign() > 0 && candidate2.sign() > 0)
          histos.fill(HIST("hrhoLikepp"), mass, pT, multiplicity);
        else
          histos.fill(HIST("hrhoLikemm"), mass, pT, multiplicity);
        break;
      case FillMode::Mixed:
        histos.fill(HIST("hrhoMixed"), mass, pT, multiplicity);
        break;
      case FillMode::Unlike:
        histos.fill(HIST("hrhoUnlike"), mass, pT, multiplicity);
        break;
    }
  }

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta) && (nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using TrackPi = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi>>;
  using Event = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::FT0Mults, aod::CentFT0Ms>>;

  void processUnlike(Event::iterator const& events, TrackPi const& tracks)
  {
    if (!events.sel8())
      return;
    float massPiplus = Ipdg->Mass(kPiPlus);
    float massPiminus = Ipdg->Mass(kPiMinus);
    float multiplicity = Multiplicity(events);

    histos.fill(HIST("hmultFT0M"), multiplicity);
    histos.fill(HIST("hNcontributor"), events.numContrib());
    histos.fill(HIST("hVtxZ"), events.posZ());

    for (auto& track1 : tracks) {

      auto _pt1 = track1.pt();

      histos.fill(HIST("hNsigmaPionTPCvspT_beforecut"), _pt1, track1.tpcNSigmaPi());

      for (auto& track2 : tracks) {

        if (track1.sign() * track2.sign() > 0)
          continue;

        if (!(selectionPID(track1)))
          continue;
        if (!(selectionPID(track2)))
          continue;

        histos.fill(HIST("hNsigmaPionTPCvspT_aftercut"), _pt1, track1.tpcNSigmaPi());
        histos.fill(HIST("ptpi"), _pt1);
        histos.fill(HIST("dcaXY"), track1.dcaXY());
        histos.fill(HIST("dcaZ"), track1.dcaZ());
        histos.fill(HIST("hNsigmaPionTPC"), track1.tpcNSigmaPi());
        histos.fill(HIST("hNsigmaPionTOF"), track1.tofNSigmaPi());

        FillHistogam(track1, track2, FillMode::Unlike, multiplicity, massPiplus, massPiminus);
      }
    }
  }
  PROCESS_SWITCH(rhoanalysis, processUnlike, "Process Same event", isUnlike);

  void processlike(Event::iterator const& events, TrackPi const& tracks)
  {
    if (!events.sel8())
      return;
    float massPiplus = Ipdg->Mass(kPiPlus);
    float massPiminus = Ipdg->Mass(kPiMinus);
    float multiplicity = Multiplicity(events);

    for (auto& track1 : tracks) {

      for (auto& track2 : tracks) {

        if (track1.sign() * track2.sign() < 0)
          continue;
        if (!(selectionPID(track1)))
          continue;
        if (!(selectionPID(track2)))
          continue;
        FillHistogam(track1, track2, FillMode::LikeSign, multiplicity, massPiplus, massPiminus);
      }
    }
  }
  PROCESS_SWITCH(rhoanalysis, processlike, "Process Same event with like sign track", isLike);

  ConfigurableAxis axisVertex{"axisVertex", {VARIABLE_WIDTH, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12}, "vertex axis for bin"};

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 2.75, 5.25, 7.75, 12.75, 17.75, 22.75, 27.75, 32.75, 37.75, 42.75, 47.75, 52.75, 57.75, 62.75, 67.75, 72.75, 77.75, 82.75, 87.75, 92.75, 97.75, 102.75, 107.75, 112.75, 117.75, 122.75, 127.75, 132.75, 137.75, 142.75, 147.75, 152.75, 157.75, 162.75, 167.75, 172.75, 177.75, 182.75, 187.75, 192.75, 197.75, 202.75, 207.75, 212.75, 217.75, 222.75, 227.75, 232.75, 237.75, 242.75, 247.75, 252.75, 257.75, 262.75, 267.75, 272.75, 277.75, 282.75, 287.75, 292.75, 297.75, 302.75, 307.75, 312.75, 317.75, 322.75, 327.75, 332.75, 337.75, 342.75, 347.75, 352.75, 357.75, 362.75, 367.75, 372.75, 377.75, 382.75, 387.75, 392.75, 397.75, 402.75, 407.75, 412.75, 417.75, 422.75, 427.75, 432.75, 437.75, 442.75, 447.75, 452.75, 457.75, 462.75, 467.75, 472.75, 477.75, 482.75, 487.75, 492.75, 497.75, 502.75, 507.75, 512.75, 517.75, 522.75, 527.75, 532.75, 537.75, 542.75, 547.75, 552.75, 557.75, 562.75, 567.75, 572.75, 577.75, 582.75, 587.75, 592.75, 597.75, 602.75, 607.75, 612.75, 617.75, 622.75, 627.75, 632.75, 637.75, 642.75, 647.75, 652.75, 657.75, 662.75, 667.75, 672.75, 677.75, 682.75, 687.75, 692.75, 697.75, 702.75, 707.75, 712.75, 717.75, 722.75, 727.75, 732.75, 737.75, 742.75, 747.75, 752.75, 757.75, 762.75, 767.75, 772.75, 777.75, 782.75, 787.75, 792.75, 797.75, 802.75, 807.75, 812.75, 817.75, 822.75, 827.75, 832.75, 837.75, 842.75, 847.75, 852.75, 857.75, 862.75, 867.75, 872.75, 877.75, 882.75, 887.75, 892.75, 897.75, 902.75, 907.75, 912.75, 917.75, 922.75, 927.75, 932.75, 937.75, 942.75, 947.75, 952.75, 957.75, 962.75, 967.75, 972.75, 977.75, 982.75, 987.75, 992.75, 997.75, 1002.75, 1007.75, 1012.75, 1017.75, 1022.75, 1027.75, 1032.75, 1037.75, 1042.75, 1047.75, 1052.75}, "multiplicity axis for histograms"};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;

  void processMixedEvent(Event const& events, TrackPi const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default)
    SameKindPair<Event, TrackPi, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1, events, tracksTuple, &cache};
    float massPiplus = Ipdg->Mass(kPiPlus);
    float massPiminus = Ipdg->Mass(kPiMinus);
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8())
        continue;
      if (!c2.sel8())
        continue;
      float multiplicity = Multiplicity(c1);

      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0)
          continue;
        if (!(selectionPID(track1)))
          continue;
        if (!(selectionPID(track2)))
          continue;

        FillHistogam(track1, track2, FillMode::Mixed, multiplicity, massPiplus, massPiminus);
      }
    }
  }
  PROCESS_SWITCH(rhoanalysis, processMixedEvent, "Process Mixed event", isMixed);

  using EventMC = soa::Join<aod::Collisions, aod::TPCMults, aod::FT0Mults, aod::CentFT0Ms, aod::McCollisionLabels>;
  using TrackMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksDCA, aod::pidTPCFullPi, aod::McTrackLabels>>;
  void processMC(EventMC::iterator const& events, TrackMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!events.has_mcCollision())
      return;

    float massPi = Ipdg->Mass(kPiPlus);
    float multiplicity = Multiplicity(events);

    if (std::abs(events.mcCollision().posZ()) > cfgCutVertex) {
      return;
    }
    for (auto track1 : tracks) {
      histos.fill(HIST("hNsigmaPionTPCvspT"), track1.pt(), track1.tpcNSigmaPi());
      if (abs(track1.tpcNSigmaPi()) > nsigmaCutCombined)
        continue;
      for (auto track2 : tracks) {

        if (abs(track2.tpcNSigmaPi()) > nsigmaCutCombined)
          continue;

        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());

        if (!(track1PDG == 211 && track2PDG == 211)) {
          continue;
        }
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
              continue;

            if (std::abs(mothertrack1.y()) > 0.5) {
              continue;
            }

            if (std::abs(mothertrack1.pdgCode()) != 113) {
              continue;
            } //

            pvec0 = std::array{track1.px(), track1.py(), track1.pz()};
            pvec1 = std::array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = std::array{pvec0, pvec1};
            recMass = RecoDecay::m(arrMomrec, std::array{massPi, massPi});
            histos.fill(HIST("hrhoRec"), recMass, mothertrack1.pt(), multiplicity);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(rhoanalysis, processMC, "Process mc reconstructed", isMC);
  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {

    if (std::abs(mcCollision.posZ()) > cfgCutVertex)
      return;
    for (auto& mcParticle : mcParticles) {

      if (std::abs(mcParticle.y()) > 0.5)
        continue;

      if (std::abs(mcParticle.pdgCode()) != 113)
        continue;

      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();

      if (kDaughters.size() != 2)
        continue;

      for (auto kCurrentDaughter : kDaughters) {

        if (!kCurrentDaughter.isPhysicalPrimary())
          continue;

        if ((kCurrentDaughter.pdgCode() == 211 || kCurrentDaughter.pdgCode() == -211)) {

          histos.fill(HIST("hrhoGen"), mcParticle.pt(), mcParticle.pdgCode());
        }
      }
    }
  }
  PROCESS_SWITCH(rhoanalysis, processGen, "Process Generated", isMC);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<rhoanalysis>(cfgc, TaskName{"lf-rhoanalysis"})};
}
