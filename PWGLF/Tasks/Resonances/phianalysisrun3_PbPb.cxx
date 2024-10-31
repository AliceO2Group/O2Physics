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

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phianalysisrun3_PbPb {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 2.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> fillOccupancy{"fillOccupancy", true, "fill Occupancy"};
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 2500, "Occupancy cut"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> isEtaAssym{"isEtaAssym", false, "isEtaAssym"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> isITSOnlycut{"isITSOnlycut", true, "isITSOnlycut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};

  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  void init(o2::framework::InitContext&)
  {
    std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
    AxisSpec occupancyAxis = {occupancyBinning, "Occupancy"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy distribution", kTH1F, {occupancyAxis});
    if (!isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRotation", "Invariant mass of Phi meson Rotation", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      if (isEtaAssym) {
        histos.add("h3PhiInvMassUnlikeSignAside", "Invariant mass of Phi meson Unlike Sign A side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassLikeSignAside", "Invariant mass of Phi meson Like Sign A side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassMixedAside", "Invariant mass of Phi meson Mixed A side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassUnlikeSignCside", "Invariant mass of Phi meson Unlike Sign C side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassLikeSignCside", "Invariant mass of Phi meson Like Sign C side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
        histos.add("h3PhiInvMassMixedCside", "Invariant mass of Phi meson Mixed C side", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      }
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiGen1", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("Centrec", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h3PhiRec3", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h1PhiRec1", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1Phimassgen", "Phi meson gen", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phipt", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
    }

    // DCA QA
    histos.add("QAbefore/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAbefore/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAafter/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAafter/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAbefore/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
    histos.add("QAbefore/TPC_Nsigma_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAafter/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
    histos.add("QAafter/TPC_Nsigma_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {{200, 0.0, 20.0}, {100, -6, 6}}});
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
    if (iscustomDCAcut && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    if (isITSOnlycut && !(candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster)) {
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
  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
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
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                    aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                    aod::McTrackLabels>>;
  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;
  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  SliceCache cache;

  // using BinningType = BinningPolicy<aod::collision::PosZ, aod::mult::MultFT0M<aod::mult::MultFT0A, aod::mult::MultFT0C>>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicityClass}, true};

  // using BinningTypeTPCMultiplicity =  ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  // using BinningTypeCentrality = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  // BinningType binningOnPositions{{axisVertex, axisMultiplicity}, true};
  ROOT::Math::PxPyPzMVector PhiMesonMother, KaonPlus, KaonMinus;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && occupancy < cfgOccupancyCut) // occupancy info is available for this collision (*)
    {
      return;
    }
    float multiplicity;
    if (cfgMultFT0)
      multiplicity = collision.centFT0C();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("QAbefore/TPC_Nsigma_all"), track1.pt(), track1.tpcNSigmaKa());
      histos.fill(HIST("QAbefore/TOF_Nsigma_all"), track1.pt(), track1.tofNSigmaKa());
      histos.fill(HIST("QAbefore/trkDCAxy"), track1.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAz"), track1.dcaZ());
      histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());

      auto track1ID = track1.globalIndex();
      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
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
          histos.fill(HIST("QAafter/TPC_Nsigma_all"), track1.pt(), track1.tpcNSigmaKa());
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), track1.pt(), track1.tofNSigmaKa());
          histos.fill(HIST("QAafter/trkDCAxy"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAz"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
        if (!isITSOnlycut && !ispTdepPID && selectionPID(track1) && selectionPID(track2)) {
          histos.fill(HIST("QAafter/TPC_Nsigma_all"), track1.pt(), track1.tpcNSigmaKa());
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), track1.pt(), track1.tofNSigmaKa());
          histos.fill(HIST("QAafter/trkDCAxy"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAz"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
        if (!isITSOnlycut && ispTdepPID && selectionPIDpTdependent(track1) && selectionPIDpTdependent(track2)) {
          histos.fill(HIST("QAafter/TPC_Nsigma_all"), track1.pt(), track1.tpcNSigmaKa());
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), track1.pt(), track1.tofNSigmaKa());
          histos.fill(HIST("QAafter/trkDCAxy"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAz"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processSameEvent, "Process Same event", false);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (!c1.sel8()) {
        continue;
      }
      if (!c2.sel8()) {
        continue;
      }
      if (additionalEvSel2 && (!c1.selection_bit(aod::evsel::kNoSameBunchPileup) || !c1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      if (additionalEvSel2 && (!c2.selection_bit(aod::evsel::kNoSameBunchPileup) || !c2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      if (additionalEvSel3 && (!c1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      if (additionalEvSel3 && (!c2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      int occupancy1 = c1.trackOccupancyInTimeRange();
      int occupancy2 = c2.trackOccupancyInTimeRange();
      if (fillOccupancy && occupancy1 < cfgOccupancyCut && occupancy2 < cfgOccupancyCut) // occupancy info is available for this collision (*)
      {
        return;
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
        if (!isITSOnlycut && !ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
        if (!isITSOnlycut && ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent, "Process Mixed event", false);

  void processMC(CollisionMCTrueTable::iterator const& /*TrueCollision*/, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("hMC"), 2);
      return;
    }
    for (auto& RecCollision : RecCollisions) {
      histos.fill(HIST("hMC"), 3);
      if (!RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }
      if (timFrameEvsel && (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        histos.fill(HIST("hMC"), 5);
        continue;
      }
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("Centrec"), centrality);
      auto oldindex = -999;
      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      for (auto track1 : Rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track1)) {
          continue;
        }
        if (ispTdepPID && !selectionPIDpTdependent(track1)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        for (auto track2 : Rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track2)) {
            continue;
          }
          if (ispTdepPID && !selectionPIDpTdependent(track2)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
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
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());
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
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (TMath::Abs(mothertrack1.y()) > confRapidity) {
                continue;
              }
              if (TMath::Abs(mothertrack1.pdgCode()) != 333) {
                continue;
              }
              if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
                continue;
              }
              if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              if (track1.sign() > 0 && track2.sign() < 0) {
                KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              if (track1.sign() < 0 && track2.sign() > 0) {
                KaonMinus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonPlus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              PhiMesonMother = KaonPlus + KaonMinus;

              if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
                continue;
              }
              histos.fill(HIST("h1PhiRec1"), PhiMesonMother.pt());
              histos.fill(HIST("h2PhiRec2"), PhiMesonMother.pt(), centrality);
              histos.fill(HIST("h1Phimassrec"), PhiMesonMother.M());
              histos.fill(HIST("h3PhiRec3"), PhiMesonMother.pt(), centrality, PhiMesonMother.M());
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
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
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            if (!genacceptancecut) {
              daughtp = true;
            }
            KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -321) {
            if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            if (!genacceptancecut) {
              daughtm = true;
            }
            KaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          PhiMesonMother = KaonPlus + KaonMinus;
          histos.fill(HIST("h1PhiGen"), PhiMesonMother.pt());
          histos.fill(HIST("h2PhiGen2"), PhiMesonMother.pt(), centrality);
          histos.fill(HIST("h1Phimassgen"), PhiMesonMother.M());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMC, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisrun3_PbPb>(cfgc, TaskName{"phianalysisrun3_PbPb"})};
}
