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

#include <Math/Vector4D.h>
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
struct phianalysisrun3 {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<bool> applyEvsel{"applyEvsel", false, "applyEvsel"};
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
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> isITSOnlycut{"isITSOnlycut", true, "isITSOnlycut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", true, "TPC Time frame boundary cut"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  void init(o2::framework::InitContext&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{201, -0.5, 200.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    if (!isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRotation", "Invariant mass of Phi meson Rotation", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      if (isEtaAssym) {
        histos.add("h3PhiInvMassUnlikeSignAside", "Invariant mass of Phi meson Unlike Sign A side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassLikeSignAside", "Invariant mass of Phi meson Like Sign A side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassMixedAside", "Invariant mass of Phi meson Mixed A side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassUnlikeSignCside", "Invariant mass of Phi meson Unlike Sign C side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassLikeSignCside", "Invariant mass of Phi meson Like Sign C side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassMixedCside", "Invariant mass of Phi meson Mixed C side", kTH3F, {{201, -0.5, 200.5}, {100, 0.0f, 10.0f}, {200, 0.9, 1.1}});
      }
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{6, 0.0f, 6.0f}});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("h3PhiRec", "Phi meson Rec", kTH3F, {{100, 0.0f, 10.0f}, {100, 0.0f, 10.0f}, {200, -0.1, 0.1}});
    }
  }

  double massKa = o2::constants::physics::MassKPlus;
  double rapidity;
  double genMass, recMass, resolution;
  double mass{0.};
  double massrotation{0.};
  double pT{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  array<float, 3> pvec1rotation;
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (iscustomDCAcut && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (isITSOnlycut && !(candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }
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
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }
  template <typename T1, typename T2>
  void FillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, bool likesign, bool rotation, float massd1, float massd2)
  {
    pvec0 = array{candidate1.px(), candidate1.py(), candidate1.pz()};
    pvec1 = array{candidate2.px(), candidate2.py(), candidate2.pz()};
    pvec1rotation = array{-candidate2.px(), -candidate2.py(), candidate2.pz()};
    auto arrMom = array{pvec0, pvec1};
    auto arrMomrotation = array{pvec0, pvec1rotation};
    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    mass = RecoDecay::m(arrMom, array{massd1, massd2});
    massrotation = RecoDecay::m(arrMomrotation, array{massd1, massd2});
    pT = RecoDecay::pt(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py()});
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

    // default filling
    if (std::abs(rapidity) < 0.5 && !isEtaAssym && track1Sign * track2Sign < 0) {
      if (unlike) {
        histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
      }
      if (mix) {
        histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
      }
      if (rotation) {
        histos.fill(HIST("h3PhiInvMassRotation"), multiplicity, pT, massrotation);
      }
    }
    if (std::abs(rapidity) < 0.5 && !isEtaAssym && track1Sign * track2Sign > 0 && likesign) {
      if (track1Sign > 0 && track2Sign > 0) {
        histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
      } else {
        histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
      }
    }
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  // using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::McCollisionLabels>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                    aod::McTrackLabels>>;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

  // using BinningType = BinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicityClass}, true};

  // using BinningTypeTPCMultiplicity =  ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  // using BinningTypeCentrality = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true};

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (applyEvsel && !collision.sel8()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    float multiplicity;
    if (cfgMultFT0)
      multiplicity = collision.centFT0C();
    if (!cfgMultFT0)
      multiplicity = collision.numContrib();
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
      auto track1ID = track1.index();
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        bool unlike = true;
        bool mix = false;
        bool likesign = true;
        bool rotation = true;
        if (isITSOnlycut) {
          FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
        if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
          FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3, processSameEvent, "Process Same event", false);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (applyEvsel && !c1.sel8()) {
        continue;
      }
      if (applyEvsel && !c2.sel8()) {
        continue;
      }
      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      float multiplicity;
      if (cfgMultFT0)
        multiplicity = c1.centFT0C();
      if (!cfgMultFT0)
        multiplicity = c1.numContrib();

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        bool likesign = false;
        bool rotation = false;
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (isITSOnlycut) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
        if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3, processMixedEvent, "Process Mixed event", false);
  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    histos.fill(HIST("hMC"), 0.5);
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 1.5);
    }
    int Nchinel = 0;
    for (auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          Nchinel = Nchinel + 1;
        }
      }
    }
    if (Nchinel > 0 && std::abs(mcCollision.posZ()) < cfgCutVertex)
      histos.fill(HIST("hMC"), 2.5);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
    histos.fill(HIST("hMC"), 3.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("hMC"), 4.5);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) > 0.5) {
        continue;
      }
      if (mcParticle.pdgCode() != 333) {
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
        if (kCurrentDaughter.pdgCode() == +321) {
          daughtp = true;
        } else if (kCurrentDaughter.pdgCode() == -321) {
          daughtm = true;
        }
      }
      if (daughtp && daughtm) {
        histos.fill(HIST("h1PhiGen"), mcParticle.pt());
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3, processGen, "Process Generated", false);
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex || !collision.sel8()) {
      return;
    }
    histos.fill(HIST("hMC"), 5.5);
    auto oldindex = -999;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (!track1.has_mcParticle()) {
        continue;
      }
      auto track1ID = track1.index();
      for (auto track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
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
        if (!(track1PDG == 321 && track2PDG == 321)) {
          continue;
        }
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            if (std::abs(mothertrack1.y()) > 0.5) {
              continue;
            }
            if (std::abs(mothertrack1.pdgCode()) != 333) {
              continue;
            }
            if (!isITSOnlycut && !(selectionPID(track1) && selectionPID(track2))) {
              continue;
            }
            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            pvec0 = array{track1.px(), track1.py(), track1.pz()};
            pvec1 = array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = array{pvec0, pvec1};
            auto motherP = mothertrack1.p();
            auto motherE = mothertrack1.e();
            genMass = std::sqrt(motherE * motherE - motherP * motherP);
            recMass = RecoDecay::m(arrMomrec, array{massKa, massKa});
            auto recpt = TMath::Sqrt((track1.px() + track2.px()) * (track1.px() + track2.px()) + (track1.py() + track2.py()) * (track1.py() + track2.py()));
            histos.fill(HIST("h3PhiRec"), mothertrack1.pt(), recpt, recMass - genMass);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3, processRec, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisrun3>(cfgc, TaskName{"phianalysisrun3"})};
}
