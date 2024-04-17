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
#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "TF1.h"

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
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phianalysisrun3_PbPb {
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
  Configurable<bool> cfgMultFT0{"cfgMultFT0", true, "cfgMultFT0"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> isITSOnlycut{"isITSOnlycut", true, "isITSOnlycut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<double> confRapidity{"confRapidity", 0.5, "Rapidity cut"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", true, "TPC Time frame boundary cut"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> genacceptancecut{"genacceptancecut", true, "use acceptance cut for generated"};
  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  void init(o2::framework::InitContext&)
  {
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{90, 0.0f, 90.0f}, {600, -0.5f, 5999.5f}});
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
      histos.add("hMC", "MC Event statistics", kTH1F, {{6, 0.0f, 6.0f}});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("Centrec", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("Centgen", "MC Centrality", kTH1F, {{200, 0.0, 200.0}});
      histos.add("h2PhiRec2", "Phi meson Rec", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h2PhiGen2", "Phi meson gen", kTH2F, {{200, 0.0f, 20.0f}, {200, 0.0, 200.0}});
      histos.add("h1PhiRec1", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("h1Phimassgen", "Phi meson gen", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phimassrec", "Phi meson Rec", kTH1F, {{200, 0.9, 1.1}});
      histos.add("h1Phipt", "Phi meson Rec", kTH1F, {{200, 0.0f, 20.0f}});
    }
    if (additionalQAplots) {
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
    // Event selection cut additional - Alex
    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
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
  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& multiplicity)
  {
    if (collision.alias_bit(kTVXinTRD)) {
      // TRD triggered
      // return 0;
    }
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(multiplicity))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(multiplicity))
      return 0;
    // if (multTrk < fMultCutLow->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(centrality))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;

    return 1;
  }

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
  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.p() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.p() >= 0.5 && candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
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
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  // using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  // using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::McCollisionLabels>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
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
    if (!collision.sel8()) {
      return;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    auto multiplicity = collision.centFT0C();
    auto multTPC = collision.multNTracksPV();
    if (additionalEvsel && !eventSelected(collision, multiplicity)) {
      return;
    }
    histos.fill(HIST("hFTOCvsTPC"), multiplicity, multTPC);
    histos.fill(HIST("hCentrality"), multiplicity);
    histos.fill(HIST("hNcontributor"), collision.numContrib());
    histos.fill(HIST("hVtxZ"), collision.posZ());
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }

      histos.fill(HIST("QAbefore/TPC_Nsigma_all"), track1.pt(), track1.tpcNSigmaKa());
      histos.fill(HIST("QAbefore/TOF_Nsigma_all"), track1.pt(), track1.tofNSigmaKa());
      histos.fill(HIST("QAbefore/trkDCAxy"), track1.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAz"), track1.dcaZ());
      histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());

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
      if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      auto multiplicity = c1.centFT0C();
      auto multiplicity2 = c2.centFT0C();
      if (additionalEvsel && !eventSelected(c1, multiplicity)) {
        // printf("Mix = %d\n", 4);
        continue;
      }
      if (additionalEvsel && !eventSelected(c2, multiplicity2)) {
        // printf("Mix = %d\n", 5);
        continue;
      }

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
    auto multiplicity = 0;
    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      multiplicity = collision.centFT0C();
      histos.fill(HIST("Centgen"), multiplicity);
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
      if (std::abs(mcParticle.y()) > confRapidity) {
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
        } else if (kCurrentDaughter.pdgCode() == -321) {

          if (genacceptancecut && kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
            daughtm = true;
          }
          if (!genacceptancecut) {
            daughtm = true;
          }
        }
      }
      if (daughtp && daughtm) {
        histos.fill(HIST("h1PhiGen"), mcParticle.pt());
        histos.fill(HIST("h2PhiGen2"), mcParticle.pt(), multiplicity);
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
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    auto multiplicity = collision.centFT0C();
    histos.fill(HIST("Centrec"), multiplicity);
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
            if (std::abs(mothertrack1.y()) >= confRapidity) {
              continue;
            }
            if (std::abs(mothertrack1.pdgCode()) != 333) {
              continue;
            }
            if (!isITSOnlycut && !ispTdepPID && !(selectionPID(track1) && selectionPID(track2))) {
              continue;
            }
            if (!isITSOnlycut && ispTdepPID && !(selectionPIDpTdependent(track1) && selectionPIDpTdependent(track2))) {
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
            histos.fill(HIST("h1PhiRec1"), mothertrack1.pt());
            histos.fill(HIST("h2PhiRec2"), mothertrack1.pt(), multiplicity);
            histos.fill(HIST("h1Phimassgen"), genMass);
            histos.fill(HIST("h1Phimassrec"), recMass);
            histos.fill(HIST("h1Phipt"), recpt);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(phianalysisrun3_PbPb, processRec, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phianalysisrun3_PbPb>(cfgc, TaskName{"phianalysisrun3_PbPb"})};
}
