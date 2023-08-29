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
// Analysis task for anti-lithium4 analysis

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

namespace
{
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"He3"};

constexpr float lithiumMass = 3.74976;
constexpr float he3Mass = o2::constants::physics::MassHelium3;
constexpr float protonMass = o2::constants::physics::MassProton;
constexpr float piMass = o2::constants::physics::MassPionCharged;
constexpr int lithium4PDG = 1000030040;
constexpr int protonPDG = 2212;
constexpr int he3PDG = 1000020030;

} // namespace

struct lithium4analysis {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  // bethe bloch parameters
  std::array<float, 6> mBBparamsHe;
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for He3"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  void init(o2::framework::InitContext&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{2001, -0.5, 2000.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 2000.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("h2NsigmaHe3TPC", "NsigmaHe3 TPC distribution", kTH2F, {{20, -5.0f, 5.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaHe3TOF", "NsigmaHe3 TOF distribution", kTH2F, {{20, -5.0f, 5.0f}, {200, -10.0f, 10.0f}});
    histos.add("h2NsigmaProtonTPC", "NsigmaProton TPC distribution", kTH2F, {{20, -3.0f, 3.0f}, {200, -5.0f, 5.0f}});
    histos.add("h2NsigmaProtonTOF", "NsigmaProton TOF distribution", kTH2F, {{20, -3.0f, 3.0f}, {200, -10.0f, 10.0f}});

    if (!isMC) {
      histos.add("h3LithiumInvMassUnlikeSign", "Invariant mass of Lithium4 Unlike Sign", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3LithiumInvMassLikeSignPP", "Invariant mass of Lithium4 Like Sign positive", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3LithiumInvMassLikeSignMM", "Invariant mass of Lithium4 Like Sign negative", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3LithiumInvMassMixed", "Invariant mass of Lithium4 Mixed", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{2, 0.0f, 2.0f}});
      histos.add("h1LitGen", "Lithium4 Gen", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h2LitRec", "Lithium4 Rec", kTH2F, {{100, 0.0f, 10.0f}, {200, -0.1, 0.1}});
    }

    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = cfgBetheBlochParams->get("He3", Form("p%i", i));
    }
    mBBparamsHe[5] = cfgBetheBlochParams->get("He3", "resolution");
  }

  double rapidity;
  double genMass, recMass, resolution;
  double mass{0.};
  double pT{0.};
  array<float, 3> momHe3;
  array<float, 3> momPr;
  array<float, 3> momLith;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut && !(candidate.isGlobalTrack() || candidate.isPVContributor())) {
      return false;
    } else {
      if (!candidate.isGlobalTrackWoDCA() || std::abs(candidate.dcaXY()) > cfgCutDCAxy || std::abs(candidate.dcaZ()) > cfgCutDCAz || !candidate.isPVContributor()) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selectionPIDProton(const T& candidate)
  {
    if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      histos.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
      histos.fill(HIST("h2NsigmaProtonTOF"), candidate.p(), candidate.tofNSigmaPr());
      return true;
    } else if (std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      histos.fill(HIST("h2NsigmaProtonTPC"), candidate.tpcInnerParam(), candidate.tpcNSigmaPr());
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDHe3(const T& candidate)
  {
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() * 2 / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);
    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    float nSigmaHe3 = static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
    if (std::abs(nSigmaHe3) < nsigmaCutTPC) {
      histos.fill(HIST("h2NsigmaHe3TPC"), candidate.tpcInnerParam() * candidate.sign(), nSigmaHe3);
      return true;
    }
    return false;
  }

  template <typename T1, typename T2>
  void FillinvMass(const T1& candidateHe3, const T2& candidatePr, float multiplicity, bool unlike, bool mix, bool likesign, float massd1, float massd2)
  {
    momHe3 = array{2 * candidateHe3.px(), 2 * candidateHe3.py(), 2 * candidateHe3.pz()};
    momPr = array{candidatePr.px(), candidatePr.py(), candidatePr.pz()};
    array momLith = array{momHe3[0] + momPr[0], momHe3[1] + momPr[1], momHe3[2] + momPr[2]};
    int he3Sign = candidateHe3.sign();
    int protonSign = candidatePr.sign();

    mass = RecoDecay::m(array{momHe3, momPr}, array{massd1, massd2});
    pT = RecoDecay::pt(momLith);
    rapidity = RecoDecay::y(momLith, lithiumMass);

    if (std::abs(rapidity) < 1 && he3Sign * protonSign < 0 && unlike) {
      histos.fill(HIST("h3LithiumInvMassUnlikeSign"), multiplicity, pT, mass);
    } else if (std::abs(rapidity) < 1 && he3Sign * protonSign < 0 && mix) {
      histos.fill(HIST("h3LithiumInvMassMixed"), multiplicity, he3Sign * pT, mass); // anti-lithium4 with negative pT
    } else if (std::abs(rapidity) < 1 && he3Sign * protonSign > 0 && likesign) {
      if (he3Sign > 0 && protonSign > 0) {
        histos.fill(HIST("h3LithiumInvMassLikeSignPP"), multiplicity, pT, mass);
      } else {
        histos.fill(HIST("h3LithiumInvMassLikeSignMM"), multiplicity, pT, mass);
      }
    }
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullPr, aod::pidTOFFullPr,
                                                    aod::McTrackLabels>>;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 2.750, 5.250, 7.750, 12.750, 17.750, 22.750, 27.750, 32.750, 37.750, 42.750, 47.750, 52.750, 57.750, 62.750, 67.750, 72.750, 77.750, 82.750, 87.750, 92.750, 97.750, 250.1}, "multiplicity axis for histograms"};

  // using BinningType = BinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicityClass}, true};

  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningType binningOnPositions{{axisVertex, axisMultiplicityClass}, true};

  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    float multiplicity;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (!cfgMultFT0)
      multiplicity = collision.centFT0M() - 0.5;
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hNcontributor"), collision.numContrib());
    histos.fill(HIST("hVtxZ"), collision.posZ());

    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionPIDHe3(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        bool unlike = true;
        bool mix = false;
        bool likesign = true;
        if (!selectionPIDProton(track2)) {
          continue;
        }
        FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, he3Mass, protonMass);
      }
    }
  }
  PROCESS_SWITCH(lithium4analysis, processSameEvent, "Process Same event", false);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }

      float multiplicity;
      if (cfgMultFT0)
        multiplicity = c1.multZeqFT0A() + c1.multZeqFT0C();
      if (!cfgMultFT0)
        multiplicity = c1.centFT0M() - 0.5;

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        bool likesign = false;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (selectionPIDHe3(t1) && selectionPIDProton(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, he3Mass, protonMass);
        } else if (selectionPIDHe3(t2) && selectionPIDProton(t1)) {
          FillinvMass(t2, t1, multiplicity, unlike, mix, likesign, he3Mass, protonMass);
        }
      }
    }
  }
  PROCESS_SWITCH(lithium4analysis, processMixedEvent, "Process Mixed event", false);

  void processMCGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 0.5);
      for (auto& mcParticle : mcParticles) {

        if (std::abs(mcParticle.y()) > 1) {
          continue;
        }
        if (mcParticle.pdgCode() != lithium4PDG) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }
        auto daughtHe3 = false;
        auto daughtPr = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (std::abs(kCurrentDaughter.pdgCode()) == he3PDG) {
            daughtHe3 = true;
          } else if (std::abs(kCurrentDaughter.pdgCode()) == protonPDG) {
            daughtPr = true;
          }
        }
        if (daughtHe3 && daughtPr) {
          histos.fill(HIST("h1LitGen"), mcParticle.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(lithium4analysis, processMCGen, "Process Generated", false);
  void processMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles)
  {
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
      return;
    }
    float multiplicity;
    if (cfgMultFT0)
      multiplicity = collision.multZeqFT0A() + collision.multZeqFT0C();
    if (!cfgMultFT0)
      multiplicity = collision.centFT0M() - 0.5;

    histos.fill(HIST("hMC"), 1.5);
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hNcontributor"), collision.numContrib());
    histos.fill(HIST("hVtxZ"), collision.posZ());

    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }

      if (!selectionPIDHe3(track1)) {
        continue;
      }

      for (auto track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }

        if (!selectionPIDProton(track2)) {
          continue;
        }

        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrackHe3 = track1.mcParticle();
        const auto mctrackPr = track2.mcParticle();

        if (!mctrackHe3.isPhysicalPrimary() || !mctrackPr.isPhysicalPrimary()) {
          continue;
        }

        if (std::abs(mctrackHe3.pdgCode()) != he3PDG || std::abs(mctrackPr.pdgCode()) != protonPDG) {
          continue;
        }

        for (auto& mothertrack : mctrackHe3.mothers_as<aod::McParticles>()) {
          for (auto& mothertrackPr : mctrackPr.mothers_as<aod::McParticles>()) {

            if (mothertrack != mothertrackPr || mothertrack.pdgCode() != lithium4PDG) {
              continue;
            }

            if (std::abs(mothertrack.y()) > 1) {
              continue;
            }

            FillinvMass(track1, track2, multiplicity, true, false, false, he3Mass, protonMass);

            momHe3 = array{2 * track1.px(), 2 * track1.py(), 2 * track1.pz()};
            momPr = array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = array{momHe3, momPr};

            auto motherP = mothertrack.p();
            auto motherE = mothertrack.e();
            genMass = std::sqrt(motherE * motherE - motherP * motherP);

            recMass = RecoDecay::m(arrMomrec, array{he3Mass, protonMass});
            histos.fill(HIST("h2LitRec"), mothertrack.pt(), recMass - genMass);
          }
        }
      }
    }
  }

  //   PROCESS_SWITCH(phianalysisrun3, processRec, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lithium4analysis>(cfgc, TaskName{"lithium4analysis"})};
}
