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
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include <TH1F.h>
#include <TRandom.h>
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
#include <Math/Vector4D.h>
#include <vector>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phik0shortanalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry eventHist{"eventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry K0SHist{"K0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiHist{"PhiHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry PhiK0SHist{"PhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurables for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1, "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.98, "V0 CosPA"};
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4.0, "NSigmaTPCPion"};

  // Configurables for drawing  TPCPion vs pT
  Configurable<float> lowmK0S{"lowmK0S", 0.48, "Lower limit on mK0Short to draw TPCPion vs pT"};
  Configurable<float> upmK0S{"upmK0S", 0.52, "Upper limit on mK0Short to draw TPCPion vs pT"};

  // Configurables for phi selection
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minimum pt cut"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "Cut on charge"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};

  // Configurables for invariant mass histograms filling
  Configurable<double> cfgFirstCutonDeltay{"cgfFirstCutonDeltay", 0.5, "First upper bound on Deltay selection"};
  Configurable<double> cfgSecondCutonDeltay{"cgfSecondCutonDeltay", 0.2, "Second upper bound on Deltay selection"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0SmassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec PhimassAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec deltayAxis = {16, 0.0f, 0.8f, "#it{#Deltay}"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};

    // Histograms
    // Number of events per selection
    eventHist.add("hEventSelection", "hEVentSelection", kTH1F, {{5, -0.5f, 4.5f}});
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    eventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a #phi");

    // Event information
    eventHist.add("hVertexZRec", "hVertexZRec", kTH1F, {vertexZAxis});
    eventHist.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});

    // K0s topological/PID cuts
    K0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    K0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    K0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {{100, -5.f, 5.f}, ptAxis});
    K0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {{100, -5.f, 5.f}, ptAxis});

    // Phis tpological/PID cuts
    PhiHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhiHist.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhiHist.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    PhiHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    PhiHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});

    // 2D mass for Phi and K0S
    PhiK0SHist.add("h3PhiK0SInvMassSameEventInclusive", "2D Invariant mass of Phi and K0Short for Same Event Inclusive", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassSameEventFirstCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < FirstCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassSameEventSecondCut", "2D Invariant mass of Phi and K0Short for Same Event Deltay < SecondCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventInclusive", "2D Invariant mass of Phi and K0Short for Mixed Event Inclusive", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventFirstCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < FirstCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
    PhiK0SHist.add("h3PhiK0SInvMassMixedEventSecondCut", "2D Invariant mass of Phi and K0Short for Mixed Event Deltay < SecondCut", kTH3F, {multAxis, K0SmassAxis, PhimassAxis});
  }

  // Constants
  double massKa = o2::constants::physics::MassKPlus;

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv && nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv && aod::v0data::dcaV0daughters < v0setting_dcav0dau);

  // Defining the type of the event
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;

  // Defining the type of the Phi daughter tracks
  using PhiDaughterCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;

  // Defining the type of the V0s
  using V0Candidates = soa::Filtered<aod::V0Datas>;

  // Defining the type of the V0 daughter tracks
  using V0DaughterCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;

  // Defining the binning policy for mixed event
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  SliceCache cache;
  Partition<PhiDaughterCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<PhiDaughterCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  // Event selection and QA filling
  template <typename T>
  bool acceptEvent(const T& collision, bool isFillEventSelectionQA)
  {
    if (isFillEventSelectionQA) {
      eventHist.fill(HIST("hEventSelection"), 0); // all collisions
    }
    if (!collision.sel8()) {
      return false;
    }
    if (isFillEventSelectionQA) {
      eventHist.fill(HIST("hEventSelection"), 1); // sel8 collisions
    }
    if (std::abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    if (isFillEventSelectionQA) {
      eventHist.fill(HIST("hEventSelection"), 2); // vertex-Z selected
      eventHist.fill(HIST("hVertexZRec"), collision.posZ());
    }
    if (!collision.isInelGt0()) {
      return false;
    }
    if (isFillEventSelectionQA) {
      eventHist.fill(HIST("hEventSelection"), 3); // INEL>0 collisions
    }
    return true;
  }

  // Track selection for Phi reconstruction
  // Topological
  template <typename T>
  bool selectionTrack(const T& track)
  {
    if (std::abs(track.pt()) < cMinPtcut) {
      return false;
    }
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut) {
      return false;
    }
    if (cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (cfgPVContributor && !track.isPVContributor()) {
      return false;
    }
    return true;
  }

  // PID
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!isNoTOF && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    if (!isNoTOF && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (isNoTOF && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    return false;
  }

  // Reconstruct the Phi
  template <typename T1, typename T2>
  TLorentzVector recMother(const T1& candidate1, const T2& candidate2, float masscand1, float masscand2)
  {
    TLorentzVector daughter1, daughter2, mother;

    daughter1.SetXYZM(candidate1.px(), candidate1.py(), candidate1.pz(), masscand1); // set the daughter1 4-momentum
    daughter2.SetXYZM(candidate2.px(), candidate2.py(), candidate2.pz(), masscand2); // set the daughter2 4-momentum
    mother = daughter1 + daughter2;                                                  // calculate the mother 4-momentum

    return mother;
  }

  // V0 selection
  template <typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& daughter1, const T2& daughter2)
  {
    if (v0.v0cosPA() < v0setting_cospa) {
      return false;
    }
    if (v0.v0radius() < v0setting_radius) {
      return false;
    }
    if (TMath::Abs(daughter1.tpcNSigmaPi()) > NSigmaTPCPion) {
      return false;
    }
    if (TMath::Abs(daughter2.tpcNSigmaPi()) > NSigmaTPCPion) {
      return false;
    }
    return true;
  }

  // Fill 2D invariant mass histogram for strange hadrons and Phi
  void fillInvMass2D(TLorentzVector V0, std::vector<TLorentzVector> listPhi, float multiplicity, bool same, bool mix, double weightInclusive, double weightLtFirstCut, double weightLtSecondCut)
  {
    double massV0 = V0.M();
    double rapidityV0 = V0.Rapidity();

    for (unsigned int phitag = 0; phitag < listPhi.size(); phitag++) {
      double massPhi = listPhi[phitag].M();
      double rapidityPhi = listPhi[phitag].Rapidity();
      double deltay = std::abs(rapidityV0 - rapidityPhi);

      if (same) { // same event
        PhiK0SHist.fill(HIST("h3PhiK0SInvMassSameEventInclusive"), multiplicity, massV0, massPhi, weightInclusive);
        if (deltay <= cfgFirstCutonDeltay) {
          PhiK0SHist.fill(HIST("h3PhiK0SInvMassSameEventFirstCut"), multiplicity, massV0, massPhi, weightLtFirstCut);
          if (deltay <= cfgSecondCutonDeltay) {
            PhiK0SHist.fill(HIST("h3PhiK0SInvMassSameEventSecondCut"), multiplicity, massV0, massPhi, weightLtSecondCut);
          }
        }
      }
      if (mix) { // mixed event
        PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventInclusive"), multiplicity, massV0, massPhi, weightInclusive);
        if (deltay <= cfgFirstCutonDeltay) {
          PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventFirstCut"), multiplicity, massV0, massPhi, weightLtFirstCut);
          if (deltay <= cfgSecondCutonDeltay) {
            PhiK0SHist.fill(HIST("h3PhiK0SInvMassMixedEventSecondCut"), multiplicity, massV0, massPhi, weightLtSecondCut);
          }
        }
      }
    }
  }

  void processSE(EventCandidates::iterator const& collision, PhiDaughterCandidates const&, V0Candidates const& V0s, V0DaughterCandidates const&)
  {
    // Check if the event selection is passed
    bool isFillEventSelectionQA = true;
    if (!acceptEvent(collision, isFillEventSelectionQA)) {
      return;
    }

    float multiplicity = collision.centFT0M();
    eventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;
    bool isFilledhPhiCand = false;

    bool same = true;
    bool mix = false;

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<V0DaughterCandidates>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DaughterCandidates>();

      // Cut on V0 dynamic columns
      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack)) {
        continue;
      }

      K0SHist.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      K0SHist.fill(HIST("hV0CosPA"), v0.v0cosPA());
      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (lowmK0S < v0.mK0Short() && v0.mK0Short() < upmK0S) {
        K0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
        K0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
      }

      TLorentzVector recK0S;
      recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
      if (recK0S.Rapidity() > 0.8) {
        continue;
      }
      std::vector<TLorentzVector> listrecPhi;
      double countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

      // Phi reconstruction
      // Loop over positive candidates
      for (auto track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrack(track1) || !selectionPID(track1)) {
          continue; // topological and PID selection
        }

        if (!isFilledhPhiCand) {
          PhiHist.fill(HIST("hEta"), track1.eta());
          PhiHist.fill(HIST("hDcaxy"), track1.dcaXY());
          PhiHist.fill(HIST("hDcaz"), track1.dcaZ());
          PhiHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
          PhiHist.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
        }

        auto track1ID = track1.globalIndex();

        // Loop over all negative candidates
        for (auto track2 : negThisColl) {
          if (!selectionTrack(track2) || !selectionPID(track2)) {
            continue; // topological and PID selection
          }
          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID) {
            continue; // condition to avoid double counting of pair
          }

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8) {
            continue;
          }
          listrecPhi.push_back(recPhi);
          countInclusive++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) < cfgFirstCutonDeltay) {
            countLtFirstCut++;
            if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) < cfgSecondCutonDeltay) {
              countLtSecondCut++;
            }
          }

          if (!isCountedPhi && isFillEventSelectionQA) {
            eventHist.fill(HIST("hEventSelection"), 4); // at least a Phi in the event
            isCountedPhi = true;
          }
        }
      }

      isFilledhPhiCand = true;

      double weightInclusive = 1 / countInclusive;
      double weightLtFirstCut = 1 / countLtFirstCut;
      double weightLtSecondCut = 1 / countLtSecondCut;
      fillInvMass2D(recK0S, listrecPhi, multiplicity, same, mix, weightInclusive, weightLtFirstCut, weightLtSecondCut);
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processSE, "Process Same Event", true);

  void processME(EventCandidates const& collisions, PhiDaughterCandidates const&, V0Candidates const& V0s, V0DaughterCandidates const&)
  {
    // Mixing the events with similar vertex z and multiplicity
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      // Check if the event selection is passed
      bool isFillEventSelectionQA = false;
      if (!acceptEvent(collision1, isFillEventSelectionQA) || !acceptEvent(collision2, isFillEventSelectionQA)) {
        continue;
      }

      float multiplicity = collision1.centFT0M();

      // Defining V0s from collision1
      auto V0ThisColl = V0s.sliceByCached(aod::v0::collisionId, collision1.globalIndex(), cache);

      // Defining positive and negative tracks for phi reconstruction from collision1 and collision2, respectively
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      bool same = false;
      bool mix = true;

      for (const auto& v0 : V0ThisColl) {
        const auto& posDaughterTrack = v0.posTrack_as<V0DaughterCandidates>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DaughterCandidates>();

        // Cut on V0 dynamic columns
        if (!selectionV0(v0, posDaughterTrack, negDaughterTrack)) {
          continue;
        }

        TLorentzVector recK0S;
        recK0S.SetXYZM(v0.px(), v0.py(), v0.pz(), v0.mK0Short());
        if (recK0S.Rapidity() > 0.8) {
          continue;
        }
        std::vector<TLorentzVector> listrecPhi;
        double countInclusive = 0, countLtFirstCut = 0, countLtSecondCut = 0;

        // Combinatorial background simulation
        for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
          if (!selectionTrack(track1) || !selectionPID(track1) || !selectionTrack(track2) || !selectionPID(track2)) {
            continue; // topological and PID selection
          }

          TLorentzVector recPhi;
          recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Rapidity() > 0.8) {
            continue;
          }
          listrecPhi.push_back(recPhi);
          countInclusive++;
          if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) < cfgFirstCutonDeltay) {
            countLtFirstCut++;
            if (std::abs(recK0S.Rapidity() - recPhi.Rapidity()) < cfgSecondCutonDeltay) {
              countLtSecondCut++;
            }
          }
        }

        double weightInclusive = 1 / countInclusive;
        double weightLtFirstCut = 1 / countLtFirstCut;
        double weightLtSecondCut = 1 / countLtSecondCut;
        fillInvMass2D(recK0S, listrecPhi, multiplicity, same, mix, weightInclusive, weightLtFirstCut, weightLtSecondCut);
      }
    }
  }

  PROCESS_SWITCH(phik0shortanalysis, processME, "Process Mixed Event", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phik0shortanalysis>(cfgc, TaskName{"lf-phik0shortanalysis"})};
}
