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
/// \brief kaon kaon analysis for higher mass resonances (code taken from phianalysisrun3)
/// \author Sawan (sawan.sawan@cern.ch)

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
#include "TF1.h"
#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

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
struct kaonkaonAnalysisRun3 {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
  Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", true, "TPC Time frame boundary cut"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selcection"};
  Configurable<bool> otherQAplots{"otherQAplots", true, "Other QA plots"};
  Configurable<bool> QAPID{"QAPID", true, "QA PID plots"};
  Configurable<bool> QAevents{"QAevents", true, "QA events"};
  Configurable<bool> cfgMultFT0M{"cfgMultFT0M", true, "true for pp (FT0M estimator) and false for PbPb (FT0C estimator)"};

  // Event selection cuts - Alex (Temporary, need to fix!)
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  // track
  Configurable<int> rotational_cut{"rotational_cut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> isEtaAssym{"isEtaAssym", false, "isEtaAssym"};
  Configurable<bool> iscustomDCAcut{"iscustomDCAcut", false, "iscustomDCAcut"};
  Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<bool> isITSOnlycut{"isITSOnlycut", true, "isITSOnlycut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<float> cmultLow{"cmultLow", 0.0f, "Low centrality percentile"};
  Configurable<float> cmultHigh{"cmultHigh", 150.0f, "High centrality percentile"};
  Configurable<int> cmultBins{"cmultBins", 150, "Number of centrality bins"};
  Configurable<float> cpTlow{"cpTlow", 0.0f, "Low pT"};
  Configurable<float> cpThigh{"cpThigh", 10.0f, "High pT"};
  Configurable<int> cpTbins{"cpTbins", 100, "Number of pT bins"};
  Configurable<float> cMasslow{"cMasslow", 0.9f, "Low mass"};
  Configurable<float> cMasshigh{"cMasshigh", 2.5f, "High mass"};
  Configurable<int> cMassbins{"cMassbins", 320, "Number of mass bins"};
  Configurable<int> c_nof_rotations{"c_nof_rotations", 3, "Number of random rotations in the rotational background"};
  ConfigurableAxis axisdEdx{"axisdEdx", {20000, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {2000, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis axisMultdist{"axisMultdist", {3500, 0, 70000}, "Multiplicity distribution"};

  // different frames
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};

  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec axisMult{cmultBins, cmultLow, cmultHigh, "Multiplicity"};
    AxisSpec axisPt{cpTbins, cpTlow, cpThigh, "pT (GeV/c)"};
    AxisSpec axisMass{cMassbins, cMasslow, cMasshigh, "Invariant mass (GeV/c^2)"};
    const AxisSpec thnAxisPOL{configThnAxisPOL, "Frame axis"};

    //  THnSparses
    std::array<int, 4> sparses = {activateTHnSparseCosThStarHelicity, activateTHnSparseCosThStarProduction, activateTHnSparseCosThStarBeam, activateTHnSparseCosThStarRandom};

    // std::array<int, 1> sparses = {activateTHnSparseCosThStarHelicity};

    if (std::accumulate(sparses.begin(), sparses.end(), 0) == 0) {
      LOGP(fatal, "No output THnSparses enabled");
    } else {
      if (activateTHnSparseCosThStarHelicity) {
        LOGP(info, "THnSparse with cosThStar w.r.t. helicity axis active.");
      }
      if (activateTHnSparseCosThStarProduction) {
        LOGP(info, "THnSparse with cosThStar w.r.t. production axis active.");
      }
      if (activateTHnSparseCosThStarBeam) {
        LOGP(info, "THnSparse with cosThStar w.r.t. beam axis active. (Gottified jackson frame)");
      }
      if (activateTHnSparseCosThStarRandom) {
        LOGP(info, "THnSparse with cosThStar w.r.t. random axis active.");
      }
    }

    if (QAevents) {
      histos.add("hmutiplicity", "Multiplicity percentile distribution", kTH1F, {axisMult});
      histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
      histos.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {axisMultdist});
      histos.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {axisMultdist});
      histos.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {axisMultdist});
    }

    if (QAPID) {
      histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hNsigmaKaonTPC_before", "NsigmaKaon TPC distribution", kTH2F, {{axisPt}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTOF_before", "NsigmaKaon TOF distribution", kTH2F, {{axisPt}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaKaonTPC_after", "NsigmaKaon TPC distribution", kTH2F, {{axisPt}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaKaonTOF_after", "NsigmaKaon TOF distribution", kTH2F, {{axisPt}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaKaonTOF_TPC_before", "NsigmaKaon TOF-TPC distribution", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaKaonTOF_TPC_after", "NsigmaKaon TOF-TPC distribution", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
    }

    if (otherQAplots) {
      histos.add("hpTvsRapidity", "pT vs Rapidity", kTH2F, {{100, 0.0f, 10.0f}, {300, -1.5f, 1.5f}});
      histos.add("hFTOCvsTPCNoCut", "Mult correlation FT0C vs. TPC without any cut", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});

      // chi2 plots
      histos.add("Chi2perclusterITS", "Chi2 / cluster for the ITS track segment", kTH1F, {{50, 0.0f, 50.0f}});
      histos.add("Chi2perclusterTPC", "Chi2 / cluster for the TPC track segment", kTH1F, {{50, 0.0f, 50.0f}});
      histos.add("Chi2perclusterTRD", "Chi2 / cluster for the TRD track segment", kTH1F, {{50, 0.0f, 50.0f}});
      histos.add("Chi2perclusterTOF", "Chi2 / cluster for the TOF track segment", kTH1F, {{50, 0.0f, 50.0f}});
    }
    if (!isMC) {
      histos.add("h3PhiInvMassUnlikeSign", "KK Unlike Sign", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      histos.add("h3PhiInvMassLikeSignPP", "KK Like Sign +", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      histos.add("h3PhiInvMassLikeSignMM", "KK Like Sign -", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      histos.add("h3PhiInvMassMixed", "KK Mixed", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      histos.add("h3PhiInvMassRotation", "KK Rotation", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
    } else if (isMC) {
      histos.add("hMC", "MC Event statistics", kTH1F, {{6, 0.0f, 6.0f}});
      histos.add("h1PhiGen", "Phi meson Gen", kTH1F, {axisPt});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {axisPt});
      histos.add("h3PhiRec", "Phi meson Rec", kTHnSparseF, {axisPt, axisPt, {200, -0.1, 0.1}}, true);
    }
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
  double massrotation1{0.};
  double massrotation2{0.};
  double pT{0.};
  array<float, 3> pvec0;
  array<float, 3> pvec1;
  array<float, 3> pvec1rotation;
  array<float, 3> pvec2rotation;
  ROOT::Math::PxPyPzMVector daughter1, daughter2;

  template <typename Collision>
  bool eventselection(Collision const& collision, const float& multiplicity)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    // if (collision.alias_bit(kTVXinTRD)) {
    //   // TRD triggered
    //   // return 0;
    // }
    auto multNTracksPV = collision.multNTracksPV();
    if (additionalEvsel && multNTracksPV < fMultPVCutLow->Eval(multiplicity)) {
      return false;
    }
    if (additionalEvsel && multNTracksPV > fMultPVCutHigh->Eval(multiplicity)) {
      return false;
    }
    // if (multTrk < fMultCutLow->Eval(multiplicity))
    //  return 0;
    // if (multTrk > fMultCutHigh->Eval(multiplicity))
    //  return 0;
    // if (multTrk > fMultMultPVCut->Eval(multNTracksPV))
    //  return 0;
    return true;
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
  template <typename T1, typename T2, typename T3>
  void FillinvMass(const T1& candidate1, const T2& candidate2, const T3& framecalculation, float multiplicity, bool unlike, bool mix, bool likesign, bool rotation, float massd1, float massd2)
  {
    TRandom* rn = new TRandom();
    int track1Sign = candidate1.sign();
    int track2Sign = candidate2.sign();
    TLorentzVector vec1, vec2, vec3, vec4, vec5;
    vec1.SetPtEtaPhiM(candidate1.pt(), candidate1.eta(), candidate1.phi(), massd1);
    vec2.SetPtEtaPhiM(candidate2.pt(), candidate2.eta(), candidate2.phi(), massd2);
    vec3 = vec1 + vec2;
    // daughter1 = ROOT::Math::PxPyPzMVector(candidate1.px(), candidate1.py(), candidate1.pz(), massd1); // Kplus
    // daughter2 = ROOT::Math::PxPyPzMVector(candidate2.px(), candidate2.py(), candidate2.pz(), massd2); // Kminus
    double rapidity = vec3.Rapidity();

    if (otherQAplots) {
      histos.fill(HIST("Chi2perclusterITS"), candidate1.itsChi2NCl());
      histos.fill(HIST("Chi2perclusterITS"), candidate2.itsChi2NCl());
      histos.fill(HIST("Chi2perclusterTPC"), candidate1.tpcChi2NCl());
      histos.fill(HIST("Chi2perclusterTPC"), candidate2.tpcChi2NCl());
      histos.fill(HIST("Chi2perclusterTRD"), candidate1.trdChi2());
      histos.fill(HIST("Chi2perclusterTRD"), candidate2.trdChi2());
      histos.fill(HIST("Chi2perclusterTOF"), candidate1.tofChi2());
      histos.fill(HIST("Chi2perclusterTOF"), candidate2.tofChi2());
    }
    if (QAPID) {
      histos.fill(HIST("dE_by_dx_TPC"), candidate1.p(), candidate1.tpcSignal());
      histos.fill(HIST("dE_by_dx_TPC"), candidate2.p(), candidate2.tpcSignal());
    }

    // polarization calculations
    // auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
    // auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
    // ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massd1); // Kaon

    // ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(vec3.Px(), vec3.Py(), vec3.Pz(), vec3.M()); // mass of KaonKaon pair
    // ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                                             // boost mother to center of mass frame
    // ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);                                                     // boost the frame of daughter same as mother
    // ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();                                                      // get the 3 vector of daughter in the frame of mother

    // default filling
    // if (activateTHnSparseCosThStarHelicity) {
    // ROOT::Math::XYZVector helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
    // auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));
    if (std::abs(rapidity) < 0.5 && track1Sign * track2Sign < 0) {
      if (unlike) {
        histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, vec3.Pt(), vec3.M(), framecalculation);
      }
      if (mix) {
        histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, vec3.Pt(), vec3.M(), framecalculation);
      }

      if (rotation) {
        for (int i = 0; i < c_nof_rotations; i++) {
          float theta2 = rn->Uniform(TMath::Pi() - TMath::Pi() / rotational_cut, TMath::Pi() + TMath::Pi() / rotational_cut);
          vec4.SetPtEtaPhiM(candidate1.pt(), candidate1.eta(), candidate1.phi() + theta2, massd1); // for rotated background
          vec5 = vec4 + vec2;
          histos.fill(HIST("h3PhiInvMassRotation"), multiplicity, vec5.Pt(), vec5.M(), framecalculation);
        }
      }
    }
    if (std::abs(rapidity) < 0.5 && track1Sign * track2Sign > 0 && likesign) {
      if (track1Sign > 0 && track2Sign > 0) {
        histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, vec3.Pt(), vec3.M(), framecalculation);
      } else {
        histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, vec3.Pt(), vec3.M(), framecalculation);
      }
    }
    // }
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  // using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::McCollisionLabels>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!eventselection(collision, collision.centFT0M())) {
      return;
    }
    float multiplicity;
    if (cfgMultFT0M == true)
      multiplicity = collision.centFT0M();
    else
      multiplicity = collision.centFT0C();
    if (QAevents) {
      histos.fill(HIST("hmutiplicity"), multiplicity);
      histos.fill(HIST("hVtxZ"), collision.posZ());
      histos.fill(HIST("hNcontributor"), collision.numContrib());
      histos.fill(HIST("multdist_FT0M"), collision.multFT0M());
      histos.fill(HIST("multdist_FT0A"), collision.multFT0A());
      histos.fill(HIST("multdist_FT0C"), collision.multFT0C());
    }
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      if (QAPID) {
        histos.fill(HIST("hEta"), track1.eta());
        histos.fill(HIST("hDcaxy"), track1.dcaXY());
        histos.fill(HIST("hDcaz"), track1.dcaZ());
        histos.fill(HIST("hNsigmaKaonTPC_before"), track1.pt(), track1.tpcNSigmaKa());
        histos.fill(HIST("hNsigmaKaonTOF_before"), track1.pt(), track1.tofNSigmaKa());
        histos.fill(HIST("hNsigmaKaonTOF_TPC_before"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
      }
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

        // calculation of event planes
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa); // Kplus
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa); // Kminus

        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
        ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massKa); // Kaon
        TLorentzVector lv1, lv2, lv3;
        lv1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), massKa);
        lv2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), massKa);
        lv3 = lv1 + lv2;

        ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KaonKaon pair
        ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                                         // boost mother to center of mass frame
        ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);                                                 // boost the frame of daughter same as mother
        ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();                                                  // get the 3 vector of daughter in the frame of mother

        bool unlike = true;
        bool mix = false;
        bool likesign = true;
        bool rotation = true;
        if (activateTHnSparseCosThStarHelicity) {
          ROOT::Math::XYZVector helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
          auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));

          if (isITSOnlycut) {
            FillinvMass(track1, track2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
          if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
            FillinvMass(track1, track2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
        } else if (activateTHnSparseCosThStarProduction) {
          ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
          auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));

          if (isITSOnlycut) {
            FillinvMass(track1, track2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
          if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
            FillinvMass(track1, track2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
        } else if (activateTHnSparseCosThStarBeam) {
          ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
          auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

          if (isITSOnlycut) {
            FillinvMass(track1, track2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
          if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
            FillinvMass(track1, track2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
        } else if (activateTHnSparseCosThStarRandom) {
          ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
          auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

          if (isITSOnlycut) {
            FillinvMass(track1, track2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
          if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
            FillinvMass(track1, track2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          }
        }

        // if (!isITSOnlycut && selectionPID(track1) && selectionPID(track2)) {
        //   // histos.fill(HIST("hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        //   // histos.fill(HIST("hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
        //   // histos.fill(HIST("hNsigmaKaonTOF_TPC_after"), track1.tofNSigmaKa(), track1.tpcNSigmaKa());
        // }
      }
    }
  }

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};
  using BinningTypeVertexContributor1 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeVertexContributor2 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  PROCESS_SWITCH(kaonkaonAnalysisRun3, processSameEvent, "Process Same event", false);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    //////// currently mixing the event with similar TPC multiplicity ////////
    BinningTypeVertexContributor1 binningOnPositions1{{axisVertex, axisMultiplicity}, true};
    BinningTypeVertexContributor2 binningOnPositions2{{axisVertex, axisMultiplicity}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor1> pair1{binningOnPositions1, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor2> pair2{binningOnPositions2, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    if (cfgMultFT0M == true) {
      for (auto& [c1, tracks1, c2, tracks2] : pair1) {
        float multiplicity = c1.centFT0M();

        if (!eventselection(c1, multiplicity)) {
          continue;
        }
        if (!eventselection(c2, multiplicity)) {
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

          // calculation of event planes
          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa); // Kplus
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa); // Kminus

          auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
          auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
          ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massKa); // Kaon
          TLorentzVector lv1, lv2, lv3;
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massKa);
          lv3 = lv1 + lv2;

          ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KaonKaon pair
          ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                                         // boost mother to center of mass frame
          ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);                                                 // boost the frame of daughter same as mother
          ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();                                                  // get the 3 vector of daughter in the frame of mother

          if (activateTHnSparseCosThStarHelicity) {
            ROOT::Math::XYZVector helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
            auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarProduction) {
            ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
            auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarBeam) {
            ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
            auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarRandom) {
            ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
            auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          }

          // if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
          //   histos.fill(HIST("hNsigmaKaonTPC_after"), t1.pt(), t1.tpcNSigmaKa());
          //   histos.fill(HIST("hNsigmaKaonTOF_after"), t1.pt(), t1.tofNSigmaKa());
          //   histos.fill(HIST("hNsigmaKaonTOF_TPC_after"), t1.tofNSigmaKa(), t1.tpcNSigmaKa());
          // }
          // if (isITSOnlycut) {
          //   FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          // }
          // if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
          //   FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          // }
        }
      }
    } else {
      for (auto& [c1, tracks1, c2, tracks2] : pair2) {
        float multiplicity = c1.centFT0C();

        if (!eventselection(c1, multiplicity)) {
          continue;
        }
        if (!eventselection(c2, multiplicity)) {
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

          // calculation of event planes
          daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa); // Kplus
          daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa); // Kminus

          auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
          auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
          ROOT::Math::PxPyPzMVector fourVecDau = ROOT::Math::PxPyPzMVector(daughter1.Px(), daughter1.Py(), daughter1.Pz(), massKa); // Kaon
          TLorentzVector lv1, lv2, lv3;
          lv1.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
          lv2.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massKa);
          lv3 = lv1 + lv2;

          ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KaonKaon pair
          ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                                         // boost mother to center of mass frame
          ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau);                                                 // boost the frame of daughter same as mother
          ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();                                                  // get the 3 vector of daughter in the frame of mother

          if (activateTHnSparseCosThStarHelicity) {
            ROOT::Math::XYZVector helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
            auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarHelicity, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarProduction) {
            ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
            auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarProduction, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarBeam) {
            ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
            auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarBeam, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          } else if (activateTHnSparseCosThStarRandom) {
            ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
            auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

            if (isITSOnlycut) {
              FillinvMass(t1, t2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
            if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
              FillinvMass(t1, t2, cosThetaStarRandom, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
            }
          }

          // if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
          //   histos.fill(HIST("hNsigmaKaonTPC_after"), t1.pt(), t1.tpcNSigmaKa());
          //   histos.fill(HIST("hNsigmaKaonTOF_after"), t1.pt(), t1.tofNSigmaKa());
          //   histos.fill(HIST("hNsigmaKaonTOF_TPC_after"), t1.tofNSigmaKa(), t1.tpcNSigmaKa());
          // }

          // if (isITSOnlycut) {
          //   FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          // }
          // if (!isITSOnlycut && selectionPID(t1) && selectionPID(t2)) {
          //   FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, rotation, massKa, massKa);
          // }
        }
      }
    }
  }

  PROCESS_SWITCH(kaonkaonAnalysisRun3, processMixedEvent, "Process Mixed event", false);
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

  PROCESS_SWITCH(kaonkaonAnalysisRun3, processGen, "Process Generated", false);
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

  PROCESS_SWITCH(kaonkaonAnalysisRun3, processRec, "Process Reconstructed", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kaonkaonAnalysisRun3>(cfgc, TaskName{"kaonkaonAnalysisRun3"})};
}
