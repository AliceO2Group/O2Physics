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
// (4) particle = 1 --> Phi
// (5) particle = 2 --> lambdastar
// (6) 4 process function (a) Data same event (b) Data mixed event (c) MC generated (d) MC reconstructed
/// \brief kaon kaon analysis for higher mass resonances (code taken from phianalysisrun3)
/// \author Sawan (sawan.sawan@cern.ch)

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

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;

struct kaonkaonAnalysisRun3 {
  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  RCTFlagsChecker rctChecker;

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry hInvMass{"hInvMass", {}, OutputObjHandlingPolicy::AnalysisObject};

  // For histograms
  Configurable<bool> calcLikeSign{"calcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> calcRotational{"calcRotational", false, "Calculate Rotational"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // Configurable<bool> piluprejection{"piluprejection", false, "Pileup rejection"};
  // Configurable<bool> goodzvertex{"goodzvertex", false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference."};
  // Configurable<bool> itstpctracks{"itstpctracks", false, "selects collisions with at least one ITS-TPC track,"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", true, "TPC Time frame boundary cut"};
  Configurable<bool> otherQAplots{"otherQAplots", true, "Other QA plots"};
  Configurable<bool> QAPID{"QAPID", true, "QA PID plots"};
  Configurable<bool> QAevents{"QAevents", true, "QA events"};
  Configurable<bool> cfgMultFT0M{"cfgMultFT0M", true, "true for pp (FT0M estimator) and false for PbPb (FT0C estimator)"};

  // track
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
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
  Configurable<bool> isITSOnlycut{"isITSOnlycut", false, "isITSOnlycut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};

  Configurable<int> cRotations{"cRotations", 3, "Number of random rotations in the rotational background"};
  ConfigurableAxis axisdEdx{"axisdEdx", {1, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {1, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis axisMultdist{"axisMultdist", {1, 0, 70000}, "Multiplicity distribution"};

  // different frames
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};
  ConfigurableAxis invMassKKAxis{"invMassKKAxis", {200, 1.0f, 3.0f}, "KK pair invariant mass axis"};
  ConfigurableAxis ptAxisKK{"ptAxisKK", {200, 0.0f, 20.0f}, "KK pair pT axis"};
  ConfigurableAxis multAxis{"multAxis", {110, 0.0f, 110.0f}, "THnSparse multiplicity axis"};

  // MC
  Configurable<bool> isMC{"isMC", false, "Run MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};
  TRandom* rn = new TRandom();

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    AxisSpec axisMult{multAxis, "Multiplicity"};
    AxisSpec axisPt{ptAxisKK, "pT (GeV/c)"};
    AxisSpec axisMass{invMassKKAxis, "Invariant mass (GeV/c^2)"};
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
      hInvMass.add("h3PhiInvMassUnlikeSign", "KK Unlike Sign", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      hInvMass.add("h3PhiInvMassLikeSignPP", "KK Like Sign +", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      hInvMass.add("h3PhiInvMassLikeSignMM", "KK Like Sign -", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      hInvMass.add("h3PhiInvMassMixed", "KK Mixed", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
      hInvMass.add("h3PhiInvMassRotated", "KK Rotation", kTHnSparseF, {axisMult, axisPt, axisMass, thnAxisPOL}, true);
    } else if (isMC) {
      hInvMass.add("h1PhiGen", "Phi meson Gen", kTH1F, {axisMult, axisPt});
      hInvMass.add("h3PhiRec", "Phi meson Rec", kTHnSparseF, {axisMult, axisPt, axisMass});
      histos.add("hMC", "MC Event statistics", kTH1F, {{6, 0.0f, 6.0f}});
      histos.add("h1PhiRecsplit", "Phi meson Rec split", kTH1F, {axisPt});
      histos.add("Recmutiplicity", "Reconstructed multiplicity distribution", kTH1F, {axisMult});
      histos.add("Genmutiplicity", "Generated multiplicity distribution", kTH1F, {axisMult});
    }
  }

  double massKa = o2::constants::physics::MassKPlus;
  double rapidity{0.0}, mass{0.}, massrotation1{0.}, massrotation2{0.}, pT{0.};
  float theta2;
  ROOT::Math::PxPyPzMVector daughter1, daughter2, daughterRot, mother, motherRot, daughterSelected, fourVecDauCM, daughterRotCM;
  ROOT::Math::XYZVector randomVec, beamVec, normalVec;
  bool isMix = false;

  template <typename Collision>
  bool eventselection(Collision const& collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return false;
    }
    // if (piluprejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    //   return false;
    // }
    // if (goodzvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    //   return false;
    // }
    // if (itstpctracks && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    //   return false;
    // }
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

  template <typename T1, typename T2>
  void fillInvMass(const T1& daughter1, const T1& daughter2, const T1& mother, float multiplicity, bool isMix, const T2& track1, const T2& track2)
  {
    ROOT::Math::Boost boost{mother.BoostToCM()};
    fourVecDauCM = boost(daughter1);

    if (std::abs(mother.Rapidity()) < 0.5) {
      if (activateTHnSparseCosThStarHelicity) {
        auto cosThetaStarHelicity = mother.Vect().Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(mother.Vect().Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);

            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);

              daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

              motherRot = daughterRot + daughter2;

              ROOT::Math::Boost boost2{motherRot.BoostToCM()};
              daughterRotCM = boost2(daughterRot);

              auto cosThetaStarHelicityRot = motherRot.Vect().Dot(daughterRotCM.Vect()) / (std::sqrt(daughterRotCM.Vect().Mag2()) * std::sqrt(motherRot.Vect().Mag2()));

              if (calcRotational)
                hInvMass.fill(HIST("h3PhiInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarHelicityRot);
            }
          } else {
            hInvMass.fill(HIST("h3PhiInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign) {
              if (track1.sign() * track2.sign() > 0) {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
              } else {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarHelicity);
              }
            }
          }
        }

      } else if (activateTHnSparseCosThStarProduction) {
        normalVec = ROOT::Math::XYZVector(mother.Py(), -mother.Px(), 0.f);
        auto cosThetaStarProduction = normalVec.Dot(fourVecDauCM.Vect()) / (std::sqrt(fourVecDauCM.Vect().Mag2()) * std::sqrt(normalVec.Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

              motherRot = daughterRot + daughter2;
              if (calcRotational)
                hInvMass.fill(HIST("h3PhiInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarProduction);
            }
          } else {
            hInvMass.fill(HIST("h3PhiInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign) {
              if (track1.sign() * track2.sign() > 0) {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
              } else {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarProduction);
              }
            }
          }
        }
      } else if (activateTHnSparseCosThStarBeam) {
        beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        auto cosThetaStarBeam = beamVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

              motherRot = daughterRot + daughter2;
              if (calcRotational)
                hInvMass.fill(HIST("h3PhiInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarBeam);
            }
          } else {
            hInvMass.fill(HIST("h3PhiInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
          }
        } else {
          if (calcLikeSign) {
            if (track1.sign() * track2.sign() > 0) {
              hInvMass.fill(HIST("h3PhiInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
            } else {
              hInvMass.fill(HIST("h3PhiInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarBeam);
            }
          }
        }
      } else if (activateTHnSparseCosThStarRandom) {
        auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
        auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);

        randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        auto cosThetaStarRandom = randomVec.Dot(fourVecDauCM.Vect()) / std::sqrt(fourVecDauCM.Vect().Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            hInvMass.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
            for (int i = 0; i < cRotations; i++) {
              theta2 = rn->Uniform(0, o2::constants::math::PI);
              daughterRot = ROOT::Math::PxPyPzMVector(daughter1.Px() * std::cos(theta2) - daughter1.Py() * std::sin(theta2), daughter1.Px() * std::sin(theta2) + daughter1.Py() * std::cos(theta2), daughter1.Pz(), daughter1.M());

              motherRot = daughterRot + daughter2;
              if (calcRotational)
                hInvMass.fill(HIST("h3PhiInvMassRotated"), multiplicity, motherRot.Pt(), motherRot.M(), cosThetaStarRandom);
            }
          } else {
            hInvMass.fill(HIST("h3PhiInvMassMixed"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
          }
        } else {
          if (!isMix) {
            if (calcLikeSign) {
              if (track1.sign() * track2.sign() > 0) {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignPP"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
              } else {
                hInvMass.fill(HIST("h3PhiInvMasslikeSignMM"), multiplicity, mother.Pt(), mother.M(), cosThetaStarRandom);
              }
            }
          }
        }
      }
    }
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  // using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::McCollisionLabels>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (rctCut.requireRCTFlagChecker && !rctChecker(collision)) {
      return;
    }
    if (!eventselection(collision)) {
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

        if (!selectionPID(track1)) // Track 1 is checked with Kaon
          continue;
        if (!selectionPID(track2)) // Track 2 is checked with Pion
          continue;

        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        mother = daughter1 + daughter2; // Kstar meson
        isMix = false;
        fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, track1, track2);
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
    BinningTypeVertexContributor1 binningOnPositions1{{axisVertex, axisMultiplicity}, true}; // for pp
    BinningTypeVertexContributor2 binningOnPositions2{{axisVertex, axisMultiplicity}, true}; // for PbPb
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor1> pair1{binningOnPositions1, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor2> pair2{binningOnPositions2, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto& [c1, tracks1, c2, tracks2] : pair1) {
      float multiplicity = c1.centFT0M();

      if (rctCut.requireRCTFlagChecker && !rctChecker(c1)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker && !rctChecker(c2)) {
        continue;
      }

      if (!eventselection(c1)) {
        continue;
      }
      if (!eventselection(c2)) {
        continue;
      }

      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

        if (!selectionTrack(t1)) {
          continue;
        }
        if (!selectionTrack(t2)) {
          continue;
        }
        if (!selectionPair(t1, t2)) {
          continue;
        }

        if (!selectionPID(t1))
          continue;
        if (!selectionPID(t2))
          continue;

        daughter1 = ROOT::Math::PxPyPzMVector(t1.px(), t1.py(), t1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(t2.px(), t2.py(), t2.pz(), massKa);
        mother = daughter1 + daughter2; // Kstar meson

        isMix = true;

        if (std::abs(mother.Rapidity()) < 0.5) {
          fillInvMass(daughter1, daughter2, mother, multiplicity, isMix, t1, t2);
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
    auto multiplicity = -1;
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
      multiplicity = collision.centFT0M();
    }
    SelectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
    histos.fill(HIST("hMC"), 3.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("Genmutiplicity"), multiplicity);
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
        hInvMass.fill(HIST("h1PhiGen"), multiplicity, mcParticle.pt());
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
    auto multiplicity = collision.centFT0M();
    histos.fill(HIST("Recmutiplicity"), multiplicity);
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
            if (!(selectionPID(track1))) {
              continue;
            }
            if (!(selectionPID(track2))) {
              continue;
            }
            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
              continue;
            }

            if (track1.sign() * track2.sign() < 0) {
              daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
            }
            mother = daughter1 + daughter2;

            if (TMath::Abs(mother.Rapidity()) >= 0.5) {
              continue;
            }
            hInvMass.fill(HIST("h3PhiRec"), multiplicity, mother.Pt(), mother.M());
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
