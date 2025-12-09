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
// Preliminary QA analysis task for resonances
// (1) For Run3
// (2) Event and track selection need to be optimized
// (3) particle = 0 --> phi
// (4) particle = 1 --> kstar
// (5) particle = 2 --> lambdastar
// (6) 4 process function (a) Data same event (b) Data mixed event (c) MC generated (d) MC reconstructed

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct resonanceqa {
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
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> isEtaAssym{"isEtaAssym", false, "isEtaAssym"};
  // particle
  Configurable<int> cfgparticletype{"cfgparticletype", 0, "Resonance particle type: 0(kstar), 1(phi), 2(Lambdastar)"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  void init(o2::framework::InitContext&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{500, 0.0, 500.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{1000, 0.0f, 1000.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaPionTPC", "NsigmaPion TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaPionTOF", "NsigmaPion TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaProtonTPC", "NsigmaProton TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaProtonTOF", "NsigmaProton TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    if (cfgparticletype == 1 && !isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3PhiInvMassRotational", "Invariant mass of Phi meson Rotational", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      if (isEtaAssym) {
        histos.add("h3PhiInvMassUnlikeSignAside", "Invariant mass of Phi meson Unlike Sign A side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
        histos.add("h3PhiInvMassLikeSignAside", "Invariant mass of Phi meson Like Sign A side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
        histos.add("h3PhiInvMassMixedAside", "Invariant mass of Phi meson Mixed A side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
        histos.add("h3PhiInvMassUnlikeSignCside", "Invariant mass of Phi meson Unlike Sign C side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
        histos.add("h3PhiInvMassLikeSignCside", "Invariant mass of Phi meson Like Sign C side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
        histos.add("h3PhiInvMassMixedCside", "Invariant mass of Phi meson Mixed C side", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
      }
    } else if (cfgparticletype == 0 && !isMC) {
      histos.add("h3KstarInvMassUnlikeSign", "Invariant mass of Kstar meson Unlike Sign", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {180, 0.6, 1.5}});
      histos.add("h3KstarInvMassLikeSignPP", "Invariant mass of Kstar meson Like Sign positive", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {180, 0.6, 1.5}});
      histos.add("h3KstarInvMassLikeSignMM", "Invariant mass of Kstar meson Like Sign negative", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {180, 0.6, 1.5}});
      histos.add("h3KstarInvMassRotational", "Invariant mass of Kstar meson Rotational", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {180, 0.6, 1.5}});
      histos.add("h3KstarInvMassMixed", "Invariant mass of Kstar meson Mixed", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {180, 0.6, 1.5}});
    } else if (cfgparticletype == 2 && !isMC) {
      histos.add("h3LambdastarInvMassUnlikeSign", "Invariant mass of Lambdastar baryon Unlike Sign", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {160, 1.2, 2.0}});
      histos.add("h3LambdastarInvMassLikeSignPP", "Invariant mass of Lambdastar baryon Like Sign positive", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {160, 1.2, 2.0}});
      histos.add("h3LambdastarInvMassLikeSignMM", "Invariant mass of Lambdastar baryon Like Sign negative", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {160, 1.2, 2.0}});
      histos.add("h3LambdastarInvMassRotational", "Invariant mass of Lambdastar baryon Rotational", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {160, 1.2, 2.0}});
      histos.add("h3LambdastarInvMassMixed", "Invariant mass of Lambdastar baryon Mixed", kTH3F, {{500, 0.0f, 500.0f}, {100, 0.0f, 10.0f}, {160, 1.2, 2.0}});
    } else if (cfgparticletype == 1 && isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{2, 0.0f, 2.0f}});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h2PhiRec", "Phi meson Rec", kTH2F, {{100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    } else if (cfgparticletype == 0 && isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{2, 0.0f, 2.0f}});
      histos.add("h1KstarGen", "Kstar meson Gen", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h1KstarRec", "Kstar meson Rec", kTH1F, {{100, 0.0f, 10.0f}});
    } else if (cfgparticletype == 2 && isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{2, 0.0f, 2.0f}});
      histos.add("h1LambdastarGen", "Lambdastar meson Gen", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h1LambdastarRec", "Lambdastar meson Rec", kTH1F, {{100, 0.0f, 10.0f}});
    }
  }
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;
  double massPr = 0.938272088f;
  double rapidity;
  double genMass, recMass, resolution;
  double mass{0.};
  double mass_roat{0.};
  double pT{0.};
  double pT_roat{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  array<float, 3> pvec1_roat;
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!candidate.isGlobalTrack()) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (candidate.hasTOF()) {
      if (PID == 0 && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      } else if (PID == 1 && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      } else if (PID == 2 && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
        return true;
      }
    } else {
      if (PID == 0 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      } else if (PID == 1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      } else if (PID == 2 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }
  template <typename T1, typename T2>
  void FillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, bool rotational, bool likesign, float massd1, float massd2)
  {
    pvec0 = array{candidate1.px(), candidate1.py(), candidate1.pz()};
    pvec1 = array{candidate2.px(), candidate2.py(), candidate2.pz()};
    pvec1_roat = array{-candidate2.px(), -candidate2.py(), candidate2.pz()};
    auto arrMom = array{pvec0, pvec1};
    auto arrMom_roat = array{pvec0, pvec1_roat};
    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    mass = RecoDecay::m(arrMom, array{massd1, massd2});
    mass_roat = RecoDecay::m(arrMom_roat, array{massd1, massd2});
    pT = RecoDecay::pt(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py()});
    pT_roat = RecoDecay::pt(array{candidate1.px() - candidate2.px(), candidate1.py() - candidate2.py()});
    rapidity = RecoDecay::y(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py(), candidate1.pz() + candidate2.pz()}, mass);
    if (isEtaAssym && unlike && track1Sign * track2Sign < 0) {
      if (candidate1.eta() > 0.2 && candidate1.eta() < 0.8 && candidate2.eta() > 0.2 && candidate2.eta() < 0.8) {
        histos.fill(HIST("h3PhiInvMassUnlikeSignAside"), multiplicity, pT, mass);
      } else if (candidate1.eta() > -0.6 && candidate1.eta() < 0.0 && candidate2.eta() > -0.6 && candidate2.eta() < 0.0) {
        histos.fill(HIST("h3PhiInvMassUnlikeSignCside"), multiplicity, pT, mass);
      }
    }
    if (isEtaAssym && mix && track1Sign * track2Sign < 0) {
      if (candidate1.eta() > 0.2 && candidate1.eta() < 0.8 && candidate2.eta() > 0.2 && candidate2.eta() < 0.8) {
        histos.fill(HIST("h3PhiInvMassMixedAside"), multiplicity, pT, mass);
      } else if (candidate1.eta() > -0.6 && candidate1.eta() < 0.0 && candidate2.eta() > -0.6 && candidate2.eta() < 0.0) {
        histos.fill(HIST("h3PhiInvMassMixedCside"), multiplicity, pT, mass);
      }
    }
    if (isEtaAssym && likesign && track1Sign * track2Sign > 0) {
      if (candidate1.eta() > 0.2 && candidate1.eta() < 0.8 && candidate2.eta() > 0.2 && candidate2.eta() < 0.8) {
        histos.fill(HIST("h3PhiInvMassLikeSignAside"), multiplicity, pT, mass);
      } else if (candidate1.eta() > -0.6 && candidate1.eta() < 0.0 && candidate2.eta() > -0.6 && candidate2.eta() < 0.0) {
        histos.fill(HIST("h3PhiInvMassLikeSignCside"), multiplicity, pT, mass);
      }
    }
    if (std::abs(rapidity) < 0.5 && !isEtaAssym) {
      if (track1Sign * track2Sign < 0 && unlike) /// unlike sign
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, pT, mass);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
        } else if (cfgparticletype == 2) {
          histos.fill(HIST("h3LambdastarInvMassUnlikeSign"), multiplicity, pT, mass);
        }
      }
      if (track1Sign * track2Sign < 0 && mix) /// mix
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, pT, mass);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
        } else if (cfgparticletype == 2) {
          histos.fill(HIST("h3LambdastarInvMassMixed"), multiplicity, pT, mass);
        }
      }
      if (track1Sign * track2Sign < 0 && rotational) // rotational
      {
        if (cfgparticletype == 0) {
          histos.fill(HIST("h3KstarInvMassRotational"), multiplicity, pT, mass_roat);
        } else if (cfgparticletype == 1) {
          histos.fill(HIST("h3PhiInvMassRotational"), multiplicity, pT, mass_roat);
        } else if (cfgparticletype == 2) {
          histos.fill(HIST("h3LambdastarInvMassRotational"), multiplicity, pT, mass_roat);
        }
      }
      if (track1Sign * track2Sign > 0 && likesign) /// like sign
      {
        if (track1Sign > 0 && track2Sign > 0) {
          if (cfgparticletype == 0) {
            histos.fill(HIST("h3KstarInvMassLikeSignPP"), multiplicity, pT, mass);
          } else if (cfgparticletype == 1) {
            histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
          } else if (cfgparticletype == 2) {
            histos.fill(HIST("h3LambdastarInvMassLikeSignPP"), multiplicity, pT, mass);
          }
        } else {
          if (cfgparticletype == 0) {
            histos.fill(HIST("h3KstarInvMassLikeSignMM"), multiplicity, pT, mass);
          } else if (cfgparticletype == 1) {
            histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
          } else if (cfgparticletype == 2) {
            histos.fill(HIST("h3LambdastarInvMassLikeSignMM"), multiplicity, pT, mass);
          }
        }
      }
    }
  }
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                  aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::TPCMults, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                    aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                    aod::pidTPCFullPr, aod::pidTOFFullPr,
                                                    aod::McTrackLabels>>;
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0.0, 2.750, 5.250, 7.750, 12.750, 17.750, 22.750, 27.750, 32.750, 37.750, 42.750, 47.750, 52.750, 57.750, 62.750, 67.750, 72.750, 77.750, 82.750, 87.750, 92.750, 97.750, 250.1}, "multiplicity axis for histograms"};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default)
  SameKindPair<EventCandidates, TrackCandidates, BinningType> pair{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    float multiplicity = collision.multTPC();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hNcontributor"), collision.numContrib());
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      histos.fill(HIST("hNsigmaPionTPC"), track1.tpcNSigmaPi());
      histos.fill(HIST("hNsigmaPionTOF"), track1.tofNSigmaPi());
      histos.fill(HIST("hNsigmaProtonTPC"), track1.tpcNSigmaPr());
      histos.fill(HIST("hNsigmaProtonTOF"), track1.tofNSigmaPr());
      auto track1ID = track1.globalIndex();
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        bool unlike = true;
        bool mix = false;
        bool rotational = true;
        bool likesign = true;
        if (cfgparticletype == 0) {
          if (selectionPID(track1, 0) && selectionPID(track2, 1)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massPi, massKa);
          } else if (selectionPID(track1, 1) && selectionPID(track2, 0)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massKa, massPi);
          }
        } else if (cfgparticletype == 1) {
          if (selectionPID(track1, 1) && selectionPID(track2, 1)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massKa, massKa);
          }
        } else if (cfgparticletype == 2) {
          if (selectionPID(track1, 1) && selectionPID(track2, 2)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massKa, massPr);
          } else if (selectionPID(track1, 2) && selectionPID(track2, 1)) {
            FillinvMass(track1, track2, multiplicity, unlike, mix, rotational, likesign, massPr, massKa);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(resonanceqa, processSameEvent, "Process Same event", false);
  void processMixedEvent(EventCandidates const&, TrackCandidates const&)
  {
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      float multiplicity = c1.multTPC();
      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        bool rotational = false;
        bool likesign = false;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (cfgparticletype == 0) {
          if (selectionPID(t1, 0) && selectionPID(t2, 1)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massPi, massKa);
          } else if (selectionPID(t1, 1) && selectionPID(t2, 0)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massKa, massPi);
          }
        } else if (cfgparticletype == 1) {
          if (selectionPID(t1, 1) && selectionPID(t2, 1)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massKa, massKa);
          }
        } else if (cfgparticletype == 2) {
          if (selectionPID(t1, 1) && selectionPID(t2, 2)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massKa, massPr);
          } else if (selectionPID(t1, 2) && selectionPID(t2, 1)) {
            FillinvMass(t1, t2, multiplicity, unlike, mix, rotational, likesign, massPr, massKa);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(resonanceqa, processMixedEvent, "Process Mixed event", false);
  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 0.5);
      for (auto& mcParticle : mcParticles) {

        if (std::abs(mcParticle.y()) > 0.5) {
          continue;
        }
        if (cfgparticletype == 0 && std::abs(mcParticle.pdgCode()) != 313) {
          continue;
        }
        if (cfgparticletype == 1 && mcParticle.pdgCode() != 333) {
          continue;
        }
        if (cfgparticletype == 2 && std::abs(mcParticle.pdgCode()) != 3124) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +321 && cfgparticletype == 1) {
            daughtp = true;
          } else if (kCurrentDaughter.pdgCode() == -321 && cfgparticletype == 1) {
            daughtm = true;
          } else if ((kCurrentDaughter.pdgCode() == 321 || kCurrentDaughter.pdgCode() == 211) && cfgparticletype == 0) {
            daughtp = true;
          } else if ((kCurrentDaughter.pdgCode() == -321 || kCurrentDaughter.pdgCode() == -211) && cfgparticletype == 0) {
            daughtm = true;
          } else if ((kCurrentDaughter.pdgCode() == 321 || kCurrentDaughter.pdgCode() == 2212) && cfgparticletype == 2) {
            daughtp = true;
          } else if ((kCurrentDaughter.pdgCode() == -321 || kCurrentDaughter.pdgCode() == -2212) && cfgparticletype == 2) {
            daughtm = true;
          }
        }
        if (daughtp && daughtm) {
          if (cfgparticletype == 1) {
            histos.fill(HIST("h1PhiGen"), mcParticle.pt());
          } else if (cfgparticletype == 0) {
            histos.fill(HIST("h1KstarGen"), mcParticle.pt());
          } else if (cfgparticletype == 2) {
            histos.fill(HIST("h1LambdastarGen"), mcParticle.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(resonanceqa, processGen, "Process Generated", false);
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("hMC"), 1.5);
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }
      auto track1ID = track1.globalIndex();
      for (auto track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID <= track1ID) {
          continue;
        }
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        if (cfgparticletype == 0 && !((track1PDG == 321 && track2PDG == 211) || (track1PDG == 211 && track2PDG == 321))) {
          continue;
        }
        if (cfgparticletype == 1 && !(track1PDG == 321 && track2PDG == 321)) {
          continue;
        }
        if (cfgparticletype == 2 && !((track1PDG == 321 && track2PDG == 2212) || (track1PDG == 2212 && track2PDG == 321))) {
          continue;
        }
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            if (mothertrack1 != mothertrack2) {
              continue;
            }
            if (std::abs(mothertrack1.y()) > 0.5) {
              continue;
            }
            if (cfgparticletype == 0) {
              if ((selectionPID(track1, 0) && selectionPID(track2, 1)) || (selectionPID(track1, 1) && selectionPID(track2, 0))) {
                if (std::abs(mothertrack1.pdgCode()) != 313) {
                  continue;
                }
                histos.fill(HIST("h1KstarRec"), mothertrack1.pt());
              }
            } else if (cfgparticletype == 2) {
              if ((selectionPID(track1, 1) && selectionPID(track2, 2)) || (selectionPID(track1, 2) && selectionPID(track2, 1))) {
                if (std::abs(mothertrack1.pdgCode()) != 3124) {
                  continue;
                }
                histos.fill(HIST("h1LambdastarRec"), mothertrack1.pt());
              }
            } else if (cfgparticletype == 1) {
              if (selectionPID(track1, 1) && selectionPID(track2, 1)) {
                if (std::abs(mothertrack1.pdgCode()) != 333) {
                  continue;
                }
                pvec0 = array{track1.px(), track1.py(), track1.pz()};
                pvec1 = array{track2.px(), track2.py(), track2.pz()};
                auto arrMomrec = array{pvec0, pvec1};
                recMass = RecoDecay::m(arrMomrec, array{massKa, massKa});
                histos.fill(HIST("h2PhiRec"), mothertrack1.pt(), recMass);
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(resonanceqa, processRec, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<resonanceqa>(cfgc, TaskName{"resonanceqa"})};
}
