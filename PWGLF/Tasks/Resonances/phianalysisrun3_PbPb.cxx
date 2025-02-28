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
#include <vector>

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
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 2.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> fillOccupancy{"fillOccupancy", true, "fill Occupancy"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<int> nBkgRotations{"nBkgRotations", 3, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
  Configurable<bool> PDGcheck{"PDGcheck", true, "PDGcheck"};
  Configurable<bool> Reco{"Reco", true, "Reco"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0, 3.5, 5.67, 7.45, 8.85, 10.0, 11.21, 12.26, 13.28, 14.23, 15.27}, "Binning of the impact parameter axis"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutMaxOccupancy{"cfgCutMaxOccupancy", 2000.0f, "Accepted maximum Occupancy"};
  Configurable<bool> cfgApplyOccupancyCut{"cfgApplyOccupancyCut", false, "Apply maximum Occupancy"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 3000, "Occupancy cut"};

  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  void init(o2::framework::InitContext&)
  {
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hOccupancy", "Occupancy distribution", kTH1F, {{500, 0, 50000}});
    if (!isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRot", "Invariant mass of Phi meson Rotation", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassSame", "Invariant mass of Phi meson same", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("EL1", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("EL2", "MC Event statistics", kTH1F, {centAxis});
      histos.add("ES1", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("ES3", "MC Event statistics", kTH1F, {impactParAxis});
      histos.add("ES2", "MC Event statistics", kTH1F, {centAxis});
      histos.add("ES4", "MC Event statistics", kTH1F, {centAxis});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiGen1", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("Centrec", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("Centgen", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h3PhiRec3", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3Phi1Rec3", "Phi meson Rec", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3PhiGen3", "Phi meson Gen", kTH3F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassMixedMC", "Invariant mass of Phi meson Mixed", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassSameMC", "Invariant mass of Phi meson same", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h3PhiInvMassRotMC", "Invariant mass of Phi meson Rotation", kTH3F, {{200, 0.0, 200.0}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
      histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h2PhiGen1", "Phi meson gen", kTH2F, {ptAxis, impactParAxis});
      histos.add("h1PhiRec1", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1Phimassgen", "Phi meson gen", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimasssame", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassmix", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrot", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phi1massrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phipt", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hOccupancy1", "Occupancy distribution", kTH1F, {{500, 0, 50000}});
      histos.add("hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterRec", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterGenCen", "Impact parameter of generated MC events", kTH2F, {impactParAxis, centAxis});
      histos.add("hImpactParameterRecCen", "Impact parameter of generated MC events", kTH2F, {impactParAxis, centAxis});
      histos.add("TOF_Nsigma_MC", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      histos.add("TPC_Nsigma_MC", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
      if (doprocessEvtLossSigLossMC) {
        histos.add("QAevent/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
        histos.add("QAevent/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
        histos.add("QAevent/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{120, 0.0f, 120.0f}, impactParAxis});
        histos.add("QAevent/phigenBeforeEvtSel", "phi before event selections", kTH2F, {ptAxis, impactParAxis});
        histos.add("QAevent/phigenAfterEvtSel", "phi after event selections", kTH2F, {ptAxis, impactParAxis});
      }
    }

    // DCA QA
    histos.add("QAbefore/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAbefore/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAafter/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    histos.add("QAafter/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {{150, 0.0f, 1.0f}});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAbefore/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAbefore/TPC_Nsigma_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {{100, -6, 6}, {100, -6, 6}}});
    histos.add("QAafter/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
    histos.add("QAafter/TPC_Nsigma_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3D, {{200, -12, 12}, {200, 0.0, 200.0}, {200, 0.0f, 20.0f}}});
  }

  double massKa = o2::constants::physics::MassKPlus;
  double rapidity;
  double genMass, recMass, resolution;
  ROOT::Math::PxPyPzMVector phiMother, daughter1, daughter2;
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
    if (!candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.hasTOF() && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }
  template <typename CollType>
  bool myEventSelections(const CollType& collision)
  {
    if (std::abs(collision.posZ()) > cfgCutVertex)
      return false;
    if (!collision.sel8())
      return false;
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)))
      return false;
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
      return false;
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy))
      return false;
    return true;
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
  void FillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, float massd1, float massd2)
  {
    pvec0 = array{candidate1.px(), candidate1.py(), candidate1.pz()};
    pvec1 = array{candidate2.px(), candidate2.py(), candidate2.pz()};
    auto arrMom = array{pvec0, pvec1};
    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    mass = RecoDecay::m(arrMom, array{massd1, massd2});
    pT = RecoDecay::pt(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py()});
    rapidity = RecoDecay::y(array{candidate1.px() + candidate2.px(), candidate1.py() + candidate2.py(), candidate1.pz() + candidate2.pz()}, mass);

    // default filling
    if (std::abs(rapidity) < 0.5 && track1Sign * track2Sign < 0) {
      if (unlike) {
        histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
      }
      if (mix) {
        histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
      }
    }
  }
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

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

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
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
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    float multiplicity{-1};
    if (cfgMultFT0)
      multiplicity = collision.centFT0C();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("QAbefore/TPC_Nsigma_all"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
      histos.fill(HIST("QAbefore/TOF_Nsigma_all"), track1.tofNSigmaKa(), multiplicity, track1.pt());
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
        if (!ispTdepPID && selectionPID(track1) && selectionPID(track2)) {
          histos.fill(HIST("QAafter/TPC_Nsigma_all"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), track1.tofNSigmaKa(), multiplicity, track1.pt());
          histos.fill(HIST("QAafter/trkDCAxy"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAz"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          FillinvMass(track1, track2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(track1) && selectionPIDpTdependent(track2)) {
          histos.fill(HIST("QAafter/TPC_Nsigma_all"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), track1.tofNSigmaKa(), multiplicity, track1.pt());
          histos.fill(HIST("QAafter/trkDCAxy"), track1.dcaXY());
          histos.fill(HIST("QAafter/trkDCAz"), track1.dcaZ());
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
          FillinvMass(track1, track2, multiplicity, unlike, mix, massKa, massKa);
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
      if (fillOccupancy && (occupancy1 > cfgCutOccupancy)) {
        continue;
      }
      if (fillOccupancy && (occupancy2 > cfgCutOccupancy)) {
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
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && selectionPID(t1) && selectionPID(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
        if (ispTdepPID && selectionPIDpTdependent(t1) && selectionPIDpTdependent(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEvent, "Process Mixed event", false);
  void processRotEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
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
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    float multiplicity{-1};
    if (cfgMultFT0)
      multiplicity = collision.centFT0C();
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hOccupancy"), occupancy);
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("QAbefore/TPC_Nsigma_all"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
      histos.fill(HIST("QAbefore/TOF_Nsigma_all"), track1.tofNSigmaKa(), multiplicity, track1.pt());
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
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (track1.sign() * track2.sign() < 0) {
          KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassSame"), multiplicity, PhiMesonMother.pt(), PhiMesonMother.M());
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (track1.sign() * track2.sign() < 0) {
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              KaonPlus = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            PhiMesonMother = KaonPlus + KaonMinus;
            if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h3PhiInvMassRot"), multiplicity, PhiMesonMother.pt(), PhiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRotEvent, "Process Rot event", false);
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
      if (additionalEvSel2 && (!RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup) || !RecCollision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      int occupancy = RecCollision.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
        continue;
      }
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }
      histos.fill(HIST("hMC"), 7);
      auto centrality = RecCollision.centFT0C();
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
              if (PDGcheck && TMath::Abs(mothertrack1.pdgCode()) != 333) {
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
              if (track1.sign() * track2.sign() < 0) {
                KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
                KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              }
              PhiMesonMother = KaonPlus + KaonMinus;

              if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
                continue;
              }
              histos.fill(HIST("h1PhiRec1"), PhiMesonMother.pt());
              histos.fill(HIST("h2PhiRec2"), PhiMesonMother.pt(), centrality);
              histos.fill(HIST("h1Phimassrec"), PhiMesonMother.M());
              histos.fill(HIST("h3PhiRec3"), PhiMesonMother.pt(), centrality, PhiMesonMother.M());
              histos.fill(HIST("Centrec"), centrality);
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.y()) > confRapidity) {
          continue;
        }
        if (PDGcheck && mcParticle.pdgCode() != 333) {
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
          histos.fill(HIST("Centgen"), centrality);
          histos.fill(HIST("h1Phimassgen"), PhiMesonMother.M());
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMC, "Process Reconstructed", false);
  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {

    histos.fill(HIST("hMC"), 0.5);
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 1.5);
    }
    float imp = mcCollision.impactParameter();
    histos.fill(HIST("hImpactParameterGen"), imp);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    auto multiplicity = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      int occupancy = collision.trackOccupancyInTimeRange();
      if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
        continue;
      }
      histos.fill(HIST("hOccupancy1"), occupancy);
      multiplicity = collision.centFT0C();
      histos.fill(HIST("Centgen"), multiplicity);
      histos.fill(HIST("hImpactParameterGenCen"), imp, multiplicity);

      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
      histos.fill(HIST("hMC"), 2.5);
    }
    SelectedEvents.resize(nevts);

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
    histos.fill(HIST("EL1"), imp);
    histos.fill(HIST("EL2"), multiplicity);
    if (Reco && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("ES1"), imp);
    histos.fill(HIST("ES2"), multiplicity);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      if (PDGcheck && mcParticle.pdgCode() != 333) {
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
          KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        } else if (kCurrentDaughter.pdgCode() == -321) {
          daughtm = true;
          KaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
        }
      }
      if (daughtp && daughtm) {
        PhiMesonMother = KaonPlus + KaonMinus;
        histos.fill(HIST("h1PhiGen"), PhiMesonMother.pt());
        histos.fill(HIST("h2PhiGen2"), PhiMesonMother.pt(), multiplicity);
        histos.fill(HIST("h2PhiGen1"), PhiMesonMother.pt(), imp);
        histos.fill(HIST("h1Phimassgen"), PhiMesonMother.M());
        histos.fill(HIST("h3PhiGen3"), PhiMesonMother.pt(), multiplicity, PhiMesonMother.M());
      }
    }
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processGen, "Process Generated", false);
  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (std::abs(collision.mcCollision().posZ()) > cfgCutVertex || !collision.sel8()) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    int occupancy = collision.trackOccupancyInTimeRange();
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    auto multiplicity = collision.centFT0C();
    histos.fill(HIST("Centrec"), multiplicity);
    float imp = collision.mcCollision().impactParameter();
    histos.fill(HIST("hImpactParameterRec"), imp);
    histos.fill(HIST("hImpactParameterRecCen"), imp, multiplicity);
    histos.fill(HIST("ES3"), imp);
    histos.fill(HIST("ES4"), multiplicity);
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
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);

        phiMother = daughter1 + daughter2;
        histos.fill(HIST("h1Phi1massrec"), phiMother.M());
        histos.fill(HIST("h3Phi1Rec3"), phiMother.pt(), multiplicity, phiMother.M());
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
            if (std::abs(mothertrack1.y()) >= 0.5) {
              continue;
            }
            if (PDGcheck && std::abs(mothertrack1.pdgCode()) != 333) {
              continue;
            }
            if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
              continue;
            }
            if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
              continue;
            }

            histos.fill(HIST("TPC_Nsigma_MC"), track1.tpcNSigmaKa(), multiplicity, track1.pt());
            histos.fill(HIST("TOF_Nsigma_MC"), track1.tofNSigmaKa(), multiplicity, track1.pt());
            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            if (track1.sign() * track2.sign() < 0) {
              KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            PhiMesonMother = KaonPlus + KaonMinus;

            if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h1PhiRec1"), PhiMesonMother.pt());
            histos.fill(HIST("h2PhiRec2"), PhiMesonMother.pt(), multiplicity);
            histos.fill(HIST("h1Phimassrec"), PhiMesonMother.M());
            histos.fill(HIST("h3PhiRec3"), PhiMesonMother.pt(), multiplicity, PhiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRec, "Process Reconstructed", false);
  void processSameEventMC(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const& /*mcParticles*/, aod::McCollisions const& /*mcCollisions*/)
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
    if (fillOccupancy && (occupancy > cfgCutOccupancy)) {
      return;
    }
    float multiplicity{-1};
    multiplicity = collision.centFT0C();
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
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
        if (!ispTdepPID && (!selectionPID(track1) || !selectionPID(track2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(track1) || !selectionPIDpTdependent(track2))) {
          continue;
        }
        if (track1.sign() * track2.sign() < 0) {
          KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        }
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassSameMC"), multiplicity, PhiMesonMother.pt(), PhiMesonMother.M());
        histos.fill(HIST("h1Phimasssame"), PhiMesonMother.M());
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            if (track1.sign() * track2.sign() < 0) {
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              KaonPlus = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            PhiMesonMother = KaonPlus + KaonMinus;
            if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
              continue;
            }
            histos.fill(HIST("h3PhiInvMassRotMC"), multiplicity, PhiMesonMother.pt(), PhiMesonMother.M());
            histos.fill(HIST("h1Phimassrot"), PhiMesonMother.M());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processSameEventMC, "Process Same event", false);
  void processMixedEventMC(EventCandidatesMC const& recCollisions, TrackCandidatesMC const& RecTracks, aod::McParticles const&)
  {

    auto tracksTuple = std::make_tuple(RecTracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidatesMC, TrackCandidatesMC, BinningTypeVertexContributor> pairs{binningOnPositions, cfgNoMixedEvents, -1, recCollisions, tracksTuple, &cache};

    for (auto& [c1, tracks1, c2, tracks2] : pairs) {
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
      if (fillOccupancy && (occupancy1 > cfgCutOccupancy)) {
        continue;
      }
      if (fillOccupancy && (occupancy2 > cfgCutOccupancy)) {
        continue;
      }
      auto multiplicity = c1.centFT0C();
      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        histos.fill(HIST("hMC"), 6.5);
        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }
        if (!ispTdepPID && (!selectionPID(t1) || !selectionPID(t2))) {
          continue;
        }
        if (ispTdepPID && (!selectionPIDpTdependent(t1) || !selectionPIDpTdependent(t2))) {
          continue;
        }
        if (t1.sign() * t2.sign() < 0) {
          KaonPlus = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa);
        }
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) > confRapidity) {
          continue;
        }
        histos.fill(HIST("h3PhiInvMassMixedMC"), multiplicity, PhiMesonMother.pt(), PhiMesonMother.M());
        histos.fill(HIST("h1Phimassmix"), PhiMesonMother.M());
      }
    }
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processMixedEventMC, "Process Mixed event MC", true);
  void processEvtLossSigLossMC(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& recCollisions)
  {

    // Event loss estimation
    auto impactPar = mcCollision.impactParameter();
    histos.fill(HIST("QAevent/hImpactParameterGen"), impactPar);

    bool isSel = false;
    auto centrality = -999.;
    for (const auto& RecCollision : recCollisions) {
      if (!myEventSelections(RecCollision))
        continue;
      centrality = RecCollision.centFT0C();
      isSel = true;
    }

    if (isSel) {
      histos.fill(HIST("QAevent/hImpactParameterRec"), impactPar);
      histos.fill(HIST("QAevent/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if (std::abs(mcPart.y()) >= 0.5 || std::abs(mcPart.pdgCode()) != 333)
        continue;

      // signal loss estimation
      histos.fill(HIST("QAevent/phigenBeforeEvtSel"), mcPart.pt(), impactPar);
      if (isSel) {
        // signal loss estimation
        histos.fill(HIST("QAevent/phigenAfterEvtSel"), mcPart.pt(), impactPar);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(phianalysisrun3_PbPb, processEvtLossSigLossMC, "Process Signal Loss, Event Loss", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisrun3_PbPb>(cfgc, TaskName{"phianalysisrun3_PbPb"})};
}
