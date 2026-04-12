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
/// \file f0980analysis.cxx
/// \brief f0(980) analysis in pp 13.6 TeV
/// \author Yunseul Bae (ybae@cern.ch), Junlee Kim (jikim1290@gmail.com)
/// \since 01/07/2024

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "TVector2.h"

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;

struct f0980analysis {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoMCCollisions>;
  using ResoMCColsEP = o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoMCCollisions, o2::aod::ResoEvtPlCollisions>;

  //  Event selections
  Configurable<float> cfgMinpT{"cfgMinpT", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 0.2, "Maximum longitudinal DCA"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};

  //  Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", true, "Global track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};

  //  PID
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"};
  Configurable<double> cMaxTPCnSigmaPionWoTOF{"cMaxTPCnSigmaPionWoTOF", 2.0, "TPC nSigma cut without TOF for Pion"};
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};
  Configurable<int> selectType{"selectType", 3, "PID selection type"};

  //  Axis
  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis pTAxis{"pTAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 95.0, 100.0, 105.0, 110.0}, "Centrality Binning"};

  Configurable<bool> cfgFindRT{"cfgFindRT", false, "boolean for RT analysis"};
  Configurable<bool> cfgFindEP{"cfgFindEP", false, "boolean for Event plane analysis"};

  Configurable<bool> cRecoINELgt0{"cRecoINELgt0", true, "check if INEL>0 for reco events"};

  Configurable<bool> cfgQAEPLT{"cfgQAEPLT", false, "Fill QA histograms for Event Plane and Leading Track"};
  Configurable<bool> cfgQACutflow{"cfgQACutflow", true, "Fill cutflow QA histograms"};
  Configurable<bool> cfgQAPairQA{"cfgQAPairQA", true, "Fill pair-level QA histograms"};
  Configurable<bool> cfgQATrackStages{"cfgQATrackStages", true, "Fill track QA at each stage"};
  Configurable<bool> cfgQASelection{"cfgQASelection", true, "Fill QA histograms for Selection"};
  Configurable<bool> cfgQAMCTrue{"cfgQAMCTrue", false, "Fill QA histograms for MC True Selection"};

  //  MC event selection
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> lptBinning = {0, 5.0, 13.0, 20.0, 50.0, 1000.0};

    AxisSpec rtAxis = {3, 0, 3};
    AxisSpec lptAxis = {lptBinning};                    //  Minimum leading hadron pT selection
    AxisSpec epAxis = {10, 0, o2::constants::math::PI}; //  Event Plane
    AxisSpec epqaAxis = {200, -o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec epResAxis = {200, -2, 2};

    AxisSpec pidqaAxis = {60, -6, 6, "#sigma"};
    AxisSpec pTqaAxis = {200, 0, 20, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec phiqaAxis = {72, 0, o2::constants::math::TwoPI, "#Phi"}; //  Azimuthal angle axis
    AxisSpec etaqaAxis = {150, -2, 2, "#eta"};                        //  Pseudorapidity axis
    AxisSpec rapqaAxis = {60, -1.5, 1.5, "#it{y}"};                   //  Rapidity axis
    AxisSpec dcaxyAxis = {200, -5.0, 5.0, "DCA_{xy} (cm)"};           //  DCAxy axis
    AxisSpec dcazAxis = {200, -5.0, 5.0, "DCA_{z} (cm)"};             //  DCAz axis
    AxisSpec invMass1DAxis = {400, 0.2, 2.2, "M_{#pi#pi} (GeV/#it{c}^{2})"};

    AxisSpec evtCutflowAxis = {5, -1.5, 3.5, "Event cutflow step"};
    AxisSpec trkCutflowAxis = {8, -0.5, 7.5, "Track cutflow step"};
    AxisSpec pairRapAxis = {80, -2.0, 2.0, "y_{#pi#pi}"};
    AxisSpec mcSelAxis = {5, -1.5, 3.5, "MC event selection step"};

    if (cfgFindRT) {
      histos.add("hInvMass_f0980_US_RT", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSpp_RT", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSmm_RT", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
    } else if (cfgFindEP) {
      histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
    } else {
      histos.add("hInvMass_f0980_US", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
      histos.add("hInvMass_f0980_LSpp", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
      histos.add("hInvMass_f0980_LSmm", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
    }

    if (cfgQAEPLT) {
      // Event Plane QA
      histos.add("QA/EPhist", "", {HistType::kTH2F, {centAxis, epqaAxis}});
      histos.add("QA/hEPResAB", "", {HistType::kTH2F, {centAxis, epResAxis}});
      histos.add("QA/hEPResBC", "", {HistType::kTH2F, {centAxis, epResAxis}});
      histos.add("QA/hEPResAC", "", {HistType::kTH2F, {centAxis, epResAxis}});
      //  Leading track pT QA
      histos.add("QA/LTpt", "", {HistType::kTH3F, {pTqaAxis, centAxis, phiqaAxis}});
    }
    if (cfgQASelection) {
      //  Event QA
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});
      //  General track QA
      if (cfgQACutflow) {
        histos.add("QAtrack/hTrackCutflow", "Track cutflow; step; counts", HistType::kTH1F, {trkCutflowAxis});
      }
      if (cfgQATrackStages) {
        // stage: 0=input(initializer tracks), 1=base(kin/quality), 2=type(global/primary), 3=DCA, 4=PID(final pion candidates)
        histos.add("QAtrack/stage/hPt_stage", "track pT vs stage; pT; stage", {HistType::kTH2F, {pTqaAxis, trkCutflowAxis}});
        histos.add("QAtrack/stage/hEta_stage", "track eta vs stage; eta; stage", {HistType::kTH2F, {etaqaAxis, trkCutflowAxis}});
        histos.add("QAtrack/stage/hPhi_stage", "track phi vs stage; phi; stage", {HistType::kTH2F, {phiqaAxis, trkCutflowAxis}});
      }
      histos.add("QAtrack/TrackPt", "", {HistType::kTH1F, {pTqaAxis}});
      histos.add("QAtrack/TrackEta", "", {HistType::kTH1F, {etaqaAxis}});
      histos.add("QAtrack/TrackPhi", "", {HistType::kTH1F, {phiqaAxis}});
      //  Track selection QA
      histos.add("QAtrack/trkDCAxy_BC", "DCA_{xy} for pion tracks (before cuts)", HistType::kTH2F, {pTqaAxis, dcaxyAxis});
      histos.add("QAtrack/trkDCAz_BC", "DCA_{z} for pion tracks (before cuts)", HistType::kTH2F, {pTqaAxis, dcazAxis});
      histos.add("QAtrack/trkDCAxy", "DCA_{xy} for pion tracks (after cuts)", HistType::kTH2F, {pTqaAxis, dcaxyAxis});
      histos.add("QAtrack/trkDCAz", "DCA_{z} for pion tracks (after cuts)", HistType::kTH2F, {pTqaAxis, dcazAxis});
      //  PID QA
      histos.add("QApid/Nsigma_TPC_BC", "TPC n#sigma^{#pi} (before PID cuts); p_{T} (GeV/c); n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QApid/Nsigma_TOF_BC", "TOF n#sigma^{#pi} (before PID cuts); p_{T} (GeV/c); n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QApid/Nsigma_TPC_TOF_BC", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});
      histos.add("QApid/Nsigma_TPC", "TPC n#sigma^{#pi} (after PID cuts); p_{T} (GeV/c); n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QApid/Nsigma_TOF", "TOF n#sigma^{#pi} (after PID cuts); p_{T} (GeV/c); n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QApid/Nsigma_TPC_TOF", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});
      if (cfgQAPairQA) {
        histos.add("QApair/hYpipi_BC", "y_{#pi#pi} (before y cut); y; counts", HistType::kTH1F, {pairRapAxis});
        histos.add("QApair/hYpipi_AC", "y_{#pi#pi} (after y cut); y; counts", HistType::kTH1F, {pairRapAxis});
        // low-dimensional mass QA (for quick visual sanity)
        histos.add("QApair/hMass_US", "M_{#pi#pi} US; M; counts", HistType::kTH1F, {invMass1DAxis});
        histos.add("QApair/hMass_LSpp", "M_{#pi#pi} ++; M; counts", HistType::kTH1F, {invMass1DAxis});
        histos.add("QApair/hMass_LSmm", "M_{#pi#pi} --; M; counts", HistType::kTH1F, {invMass1DAxis});
      }
    }

    if (doprocessMCRec) {
      histos.add("MCL/hMass_f0980_REC", "M_{REC} f0", HistType::kTH3F, {massAxis, pTqaAxis, centAxis});
      histos.add("MCL/hpT_f0980_REC", "REC truth-matched f0; p_{T}; Centrality (%)", HistType::kTH2F, {pTqaAxis, centAxis});
      histos.add("MCL/hpT_f0980_REC_truePt", "REC truth-matched f0; p_{T}^{true mother}; Centrality (%)", HistType::kTH2F, {pTqaAxis, centAxis});
    }
    if (doprocessMCTrue) {
      // histos.add("MCL/hpT_f0980_GEN", "GEN f0; p_{T} (GeV/c); Centrality (%)", HistType::kTH2F, {pTqaAxis, centAxis});
      histos.add("MCL/hpT_f0980_GEN", "GEN f0; step; p_{T} (GeV/c); centFT0M (%)", HistType::kTH3F, {mcSelAxis, pTqaAxis, centAxis});
      if (cfgQAMCTrue) {
        histos.add("QAMCTrue/f0_pt_y", "Generated f0 ; #it{p}_{T} (GeV/#it{c}) ; #it{y}", HistType::kTH2F, {pTqaAxis, rapqaAxis});
        histos.add("QAMCTrue/f0_pt_cent", "Generated f0 ; #it{p}_{T} (GeV/#it{c}); Centrality (%)", HistType::kTH2F, {pTqaAxis, centAxis});
        histos.add("QAMCTrue/f0YieldSteps", "MC true f0 counts per event-selection; step; counts", HistType::kTH1F, {evtCutflowAxis});
        // histos.add("QAMCTrue/trueINELgt0", "MC true INEL>0 flag; flag; counts", HistType::kTH1F, {inelAxis});
      }
    }
    histos.print();
  }

  double massPi = MassPionCharged;

  static constexpr float OneThird = 1.0f / 3.0f;
  static constexpr float PIthird = o2::constants::math::PI * OneThird;
  static constexpr float TWOPIthird = o2::constants::math::TwoPI * OneThird;

  int rtIndex(double pairphi, double lhphi)
  {
    double dphi = std::fabs(TVector2::Phi_mpi_pi(lhphi - pairphi));

    if (dphi < PIthird)
      return 0;
    if (dphi < TWOPIthird && dphi > PIthird)
      return 1;
    if (dphi > TWOPIthird)
      return 2;

    return -1;
  }

  template <typename TrackType>
  bool selTrackBase(const TrackType& track)
  {
    if (std::abs(track.pt()) < cfgMinpT)
      return false;
    if (std::fabs(track.eta()) > cfgMaxEta)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;

    return true;
  }

  template <typename TrackType>
  bool selTrackType(const TrackType& track)
  {
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;

    return true;
  }

  template <typename TrackType>
  bool selTrackDCA(const TrackType& track)
  {
    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    return true;
  }

  template <typename TrackType>
  bool selTrack(const TrackType& track)
  {
    if (!selTrackBase(track))
      return false;
    if (!selTrackType(track))
      return false;
    if (!selTrackDCA(track))
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;

    return true;
  }

  template <typename TrackType>
  bool selPion(const TrackType track)
  {
    switch (selectType) {
      case 0:
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
          return false;
        break;
      case 1:
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion)
          return false;
        break;
      case 2:
        if (track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi() >= nsigmaCutCombinedPion * nsigmaCutCombinedPion)
          return false;
        break;
      case 3:
        if (track.hasTOF()) {
          if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
            return false;
        } else {
          if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPionWoTOF)
            return false;
        }
        break;
    }
    return true;
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks,
                      aod::McParticles const* mcParticles = nullptr,
                      float centOverride = NAN)
  {
    const float cent = std::isfinite(centOverride) ? centOverride : collision.cent();
    if (cfgQASelection) {
      histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), cent);
    }
    double lhpT = 0.;
    double lhphi = 0.;
    double relphi = 0.;
    if (cfgFindRT) {
      for (const auto& trk : dTracks) {
        if (trk.pt() > lhpT) {
          lhpT = trk.pt();
          lhphi = trk.phi();
        }
      }
      if (cfgQAEPLT)
        histos.fill(HIST("QA/LTpt"), lhpT, cent, lhphi);
    } else if (cfgFindEP) {
      if (cfgQAEPLT) {
        histos.fill(HIST("QA/EPhist"), cent, collision.evtPl());
        histos.fill(HIST("QA/hEPResAB"), cent, collision.evtPlResAB());
        histos.fill(HIST("QA/hEPResBC"), cent, collision.evtPlResBC());
        histos.fill(HIST("QA/hEPResAC"), cent, collision.evtPlResAC());
      }
    }

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> pion1, pion2, reco;
    // Single-track QA
    if (cfgQASelection) {
      for (const auto& trk : dTracks) {
        if (cfgQACutflow) {
          histos.fill(HIST("QAtrack/hTrackCutflow"), 0); // all tracks (initializer track)
        }
        if (cfgQATrackStages) {
          histos.fill(HIST("QAtrack/stage/hPt_stage"), trk.pt(), 0);
          histos.fill(HIST("QAtrack/stage/hEta_stage"), trk.eta(), 0);
          histos.fill(HIST("QAtrack/stage/hPhi_stage"), trk.phi(), 0);
        }
        histos.fill(HIST("QAtrack/trkDCAxy_BC"), trk.pt(), trk.dcaXY());
        histos.fill(HIST("QAtrack/trkDCAz_BC"), trk.pt(), trk.dcaZ());

        if (!selTrackBase(trk)) {
          continue;
        }
        if (cfgQACutflow) {
          histos.fill(HIST("QAtrack/hTrackCutflow"), 1); // after kinematics cut
        }
        if (cfgQATrackStages) {
          histos.fill(HIST("QAtrack/stage/hPt_stage"), trk.pt(), 1);
          histos.fill(HIST("QAtrack/stage/hEta_stage"), trk.eta(), 1);
          histos.fill(HIST("QAtrack/stage/hPhi_stage"), trk.phi(), 1);
        }
        histos.fill(HIST("QApid/Nsigma_TPC_BC"), trk.pt(), trk.tpcNSigmaPi());
        if (trk.hasTOF()) {
          histos.fill(HIST("QApid/Nsigma_TOF_BC"), trk.pt(), trk.tofNSigmaPi());
          histos.fill(HIST("QApid/Nsigma_TPC_TOF_BC"), trk.tpcNSigmaPi(), trk.tofNSigmaPi());
        }

        if (!selTrackType(trk)) {
          continue;
        }
        if (cfgQACutflow) {
          histos.fill(HIST("QAtrack/hTrackCutflow"), 2); // after track type
        }
        if (cfgQATrackStages) {
          histos.fill(HIST("QAtrack/stage/hPt_stage"), trk.pt(), 2);
          histos.fill(HIST("QAtrack/stage/hEta_stage"), trk.eta(), 2);
          histos.fill(HIST("QAtrack/stage/hPhi_stage"), trk.phi(), 2);
        }

        if (!selTrackDCA(trk)) {
          continue;
        }
        histos.fill(HIST("QAtrack/trkDCAxy"), trk.pt(), trk.dcaXY());
        histos.fill(HIST("QAtrack/trkDCAz"), trk.pt(), trk.dcaZ());
        if (cfgQACutflow) {
          histos.fill(HIST("QAtrack/hTrackCutflow"), 3); // after DCA
        }
        if (cfgQATrackStages) {
          histos.fill(HIST("QAtrack/stage/hPt_stage"), trk.pt(), 3);
          histos.fill(HIST("QAtrack/stage/hEta_stage"), trk.eta(), 3);
          histos.fill(HIST("QAtrack/stage/hPhi_stage"), trk.phi(), 3);
        }

        if (!selPion(trk)) {
          continue;
        }
        if (cfgQACutflow) {
          histos.fill(HIST("QAtrack/hTrackCutflow"), 4); //  after PID for pion selection
        }
        if (cfgQATrackStages) {
          histos.fill(HIST("QAtrack/stage/hPt_stage"), trk.pt(), 4);
          histos.fill(HIST("QAtrack/stage/hEta_stage"), trk.eta(), 4);
          histos.fill(HIST("QAtrack/stage/hPhi_stage"), trk.phi(), 4);
        }
        histos.fill(HIST("QApid/Nsigma_TPC"), trk.pt(), trk.tpcNSigmaPi());
        if (trk.hasTOF()) {
          histos.fill(HIST("QApid/Nsigma_TOF"), trk.pt(), trk.tofNSigmaPi());
          histos.fill(HIST("QApid/Nsigma_TPC_TOF"), trk.tpcNSigmaPi(), trk.tofNSigmaPi());
        }

        histos.fill(HIST("QAtrack/TrackPt"), trk.pt());
        histos.fill(HIST("QAtrack/TrackEta"), trk.eta());
        histos.fill(HIST("QAtrack/TrackPhi"), trk.phi());
      }
    }

    for (const auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(dTracks, dTracks))) {
      if (!selTrack(trk1) || !selTrack(trk2))
        continue;
      if (!selPion(trk1) || !selPion(trk2))
        continue;

      pion1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      pion2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      reco = pion1 + pion2;
      if (cfgQASelection && cfgQAPairQA) {
        histos.fill(HIST("QApair/hYpipi_BC"), reco.Rapidity());
      }

      if (reco.Rapidity() > cfgMaxRap || reco.Rapidity() < cfgMinRap)
        continue;
      if (cfgQASelection && cfgQAPairQA) {
        histos.fill(HIST("QApair/hYpipi_AC"), reco.Rapidity());
      }

      if (cfgFindEP) {
        relphi = TVector2::Phi_0_2pi(reco.Phi() - collision.evtPl());
        if (relphi > o2::constants::math::PI) {
          relphi -= o2::constants::math::PI;
        }
      }

      if (trk1.sign() * trk2.sign() < 0) {
        if (cfgQASelection && cfgQAPairQA) {
          histos.fill(HIST("QApair/hMass_US"), reco.M());
        }
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_US_RT"), reco.M(), reco.Pt(), cent, rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_US_EPA"), reco.M(), reco.Pt(), cent, relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_US"), reco.M(), reco.Pt(), cent);
        }

        if constexpr (IsMC) {
          // truth-matched f0 -> pi+pi- (signal-only REC)
          if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kPiPlus)
            continue;
          if (trk1.motherId() != trk2.motherId())
            continue;
          if (std::abs(trk1.motherPDG()) != 9010221)
            continue;

          histos.fill(HIST("MCL/hMass_f0980_REC"), reco.M(), reco.Pt(), cent);
          histos.fill(HIST("MCL/hpT_f0980_REC"), reco.Pt(), cent);
          if (mcParticles && trk1.motherId() >= 0) {
            auto mother = mcParticles->iteratorAt(trk1.motherId());
            histos.fill(HIST("MCL/hpT_f0980_REC_truePt"), mother.pt(), cent);
          }
        }
      } else if (trk1.sign() > 0 && trk2.sign() > 0) {
        if (cfgQASelection && cfgQAPairQA) {
          histos.fill(HIST("QApair/hMass_LSpp"), reco.M());
        }
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_LSpp_RT"), reco.M(), reco.Pt(), cent, rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), reco.M(), reco.Pt(), cent, relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_LSpp"), reco.M(), reco.Pt(), cent);
        }
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        if (cfgQASelection && cfgQAPairQA) {
          histos.fill(HIST("QApair/hMass_LSmm"), reco.M());
        }
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_LSmm_RT"), reco.M(), reco.Pt(), cent, rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), reco.M(), reco.Pt(), cent, relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_LSmm"), reco.M(), reco.Pt(), cent);
        }
      }
    }
  }

  void processData(o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoEvtPlCollisions>::iterator const& resoCollision,
                   o2::aod::ResoCollisionColls const& collisionIndex,
                   o2::soa::Join<o2::aod::Collisions, o2::aod::PVMults> const& collisions,
                   o2::aod::ResoTracks const& resotracks)
  {
    if (cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId();
      if (collId < 0)
        return;
      auto coll = collisions.iteratorAt(collId);
      if (!coll.isInelGt0())
        return;
    }
    fillHistograms<false>(resoCollision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processData, "Process Event for data", true);

  void processMCRec(ResoMCColsEP::iterator const& collision,
                    o2::aod::ResoCollisionColls const& collisionIndex,
                    o2::soa::Join<o2::aod::ResoCollisionCandidatesMC, o2::aod::PVMults> const& collisionsMC,
                    o2::soa::Join<o2::aod::ResoTracks, o2::aod::ResoMCTracks> const& resotracks,
                    o2::aod::McParticles const& mcParticles,
                    o2::soa::Join<o2::aod::McCollisions, o2::aod::McCentFT0Ms> const&)
  {
    auto linkRow = collisionIndex.iteratorAt(collision.globalIndex());
    auto collId = linkRow.collisionId();
    if (collId < 0) {
      return;
    }
    auto coll = collisionsMC.iteratorAt(collId);
    if (cRecoINELgt0 && !coll.isInelGt0())
      return;
    auto mcColl = coll.mcCollision_as<o2::soa::Join<o2::aod::McCollisions, o2::aod::McCentFT0Ms>>();
    const float centFT0M = mcColl.centFT0M();
    if (!collision.isInAfterAllCuts() || (std::abs(collision.posZ()) > cZvertCutMC)) {
      return;
    }

    fillHistograms<true>(collision, resotracks, &mcParticles, centFT0M);
  }
  PROCESS_SWITCH(f0980analysis, processMCRec, "Process Event for MC (REC)", false);

  void processMCTrue(ResoMCCols::iterator const& resoCollision,
                     o2::aod::ResoCollisionColls const& collisionIndex,
                     o2::aod::ResoMCParents const& resoParents,
                     o2::aod::ResoCollisionCandidatesMC const& collisionsMC,
                     o2::soa::Join<o2::aod::McCollisions, o2::aod::McCentFT0Ms> const&)
  {
    auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
    auto collId = linkRow.collisionId();
    if (collId < 0) {
      return;
    }
    auto coll = collisionsMC.iteratorAt(collId);
    if (cRecoINELgt0 && !coll.isInelGt0())
      return;
    auto mcColl = coll.mcCollision_as<o2::soa::Join<o2::aod::McCollisions, o2::aod::McCentFT0Ms>>();
    const float centFT0M = mcColl.centFT0M();

    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != 9010221)
        continue;
      if (!part.producedByGenerator())
        continue;
      if (part.y() < cfgMinRap || part.y() > cfgMaxRap) {
        continue;
      }
      const bool decayPions = (std::abs(part.daughterPDG1()) == kPiPlus && std::abs(part.daughterPDG2()) == kPiPlus);
      if (!decayPions)
        continue;

      // no event selection: baseline
      histos.fill(HIST("MCL/hpT_f0980_GEN"), -1, part.pt(), centFT0M);
      histos.fill(HIST("QAMCTrue/f0YieldSteps"), -1);
      // |zvtx|<10
      if (resoCollision.isVtxIn10()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 0, part.pt(), centFT0M);
        histos.fill(HIST("QAMCTrue/f0YieldSteps"), 0);
      }
      // |zvtx|<10 && sel8
      if (resoCollision.isVtxIn10() && resoCollision.isInSel8()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 1, part.pt(), centFT0M);
        histos.fill(HIST("QAMCTrue/f0YieldSteps"), 1);
      }
      // |zvtx|<10 && TVX
      if (resoCollision.isVtxIn10() && resoCollision.isTriggerTVX()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 2, part.pt(), centFT0M);
        histos.fill(HIST("QAMCTrue/f0YieldSteps"), 2);
      }
      // after all event cuts
      if (resoCollision.isInAfterAllCuts()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 3, part.pt(), centFT0M);
        histos.fill(HIST("QAMCTrue/f0YieldSteps"), 3);
      }
      if (cfgQAMCTrue) {
        histos.fill(HIST("QAMCTrue/f0_pt_y"), part.pt(), part.y());
        histos.fill(HIST("QAMCTrue/f0_pt_cent"), part.pt(), centFT0M);
      }
    }
  };
  PROCESS_SWITCH(f0980analysis, processMCTrue, "Process Event for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<f0980analysis>(cfgc, TaskName{"lf-f0980analysis"})};
}
