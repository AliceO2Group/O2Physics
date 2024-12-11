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

/// \file   flowSP.cxx
/// \author Noor Koster
/// \since  01/12/2024
/// \brief  task to evaluate flow with respect to spectator plane.

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/RecoDecay.h"

#include "PWGCF/DataModel/SPTableZDC.h"
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowSP {

  O2_DEFINE_CONFIGURABLE(cfgDCAxy, float, 0.2, "Cut on DCA in the transverse direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgDCAz, float, 2, "Cut on DCA in the longitudinal direction (cm)");
  O2_DEFINE_CONFIGURABLE(cfgNcls, float, 70, "Cut on number of TPC clusters found");
  O2_DEFINE_CONFIGURABLE(cfgPtmin, float, 0.2, "minimum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgPtmax, float, 10, "maximum pt (GeV/c)");
  O2_DEFINE_CONFIGURABLE(cfgEta, float, 0.8, "eta cut");
  O2_DEFINE_CONFIGURABLE(cfgVtxZ, float, 10, "vertex cut (cm)");
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, true, "Bool to enable Additional Event Cut");
  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalTrackCut, bool, true, "Bool to enable Additional Track Cut");
  O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 60, "Maximum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 10, "Minimum cenrality for selected events");

  O2_DEFINE_CONFIGURABLE(cfgDoubleTrackFunction, bool, true, "Include track cut at low pt");
  O2_DEFINE_CONFIGURABLE(cfgTrackCutSize, float, 0.06, "Spread of track cut");
  O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 500, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, true, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy&& nabs(aod::track::dcaZ) < cfgDCAz;
  using UsedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::SPTableZDC, aod::Qvectors>>;
  using UsedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

  //  Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry"};

  // Event selection cuts - Alex
  TF1* fPhiCutLow = nullptr;
  TF1* fPhiCutHigh = nullptr;
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};

    AxisSpec phiModAxis = {100, 0, constants::math::PI / 9, "fmod(#varphi,#pi/9)"};
    AxisSpec ptAxis = {ptbinning, "#it{p}_{T} GeV/#it{c}"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    registry.add<TH1>("hSPplaneA", "hSPplaneA", kTH1D, {{100, -3, 3}});
    registry.add<TH1>("hSPplaneC", "hSPplaneC", kTH1D, {{100, -3, 3}});
    registry.add<TH1>("hSPplaneA-C", "hSPplaneA-C", kTH1D, {{100, -3, 3}});
    registry.add<TH1>("hCent", "hCent", kTH1D, {{80, 0, 80}});

    registry.add<TH1>("hqIm", "hqIm", kTH1D, {{100, -2, 2}});
    registry.add<TH1>("hqRe", "hqRe", kTH1D, {{100, -2, 2}});

    registry.add<TProfile>("hCosdPhi", "hCosdPhi; Centrality(%); #LT Cos( #Psi^{A}-#Psi^{C})#GT", kTProfile, {{80, 0, 80}});
    registry.add<TProfile>("hSindPhi", "hSindPhi; Centrality(%); #LT Sin( #Psi^{A}-#Psi^{C})#GT", kTProfile, {{80, 0, 80}});
    registry.add<TProfile>("hSPlaneRes", "hSPlaneRes; Centrality(%); ", kTProfile, {{80, 0, 80}});

    registry.add("pt_phi_bef", "", {HistType::kTH2D, {ptAxis, phiModAxis}});
    registry.add("pt_phi_aft", "", {HistType::kTH2D, {ptAxis, phiModAxis}});

    registry.add<TProfile2D>("v1_eta", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1A_eta", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1C_eta", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1AC_eta", "", kTProfile2D, {{10, -.8, .8}, ptAxis});

    registry.add<TProfile2D>("v1_eta_odd", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even", "", kTProfile2D, {{10, -.8, .8}, ptAxis});

    registry.add<TProfile2D>("v1_eta_odd_dev", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even_dev", "", kTProfile2D, {{10, -.8, .8}, ptAxis});

    registry.add<TProfile2D>("v1_eta_odd_dev_pos", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even_dev_pos", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_odd_pos", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even_pos", "", kTProfile2D, {{10, -.8, .8}, ptAxis});

    registry.add<TProfile2D>("v1_eta_odd_dev_neg", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even_dev_neg", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_odd_neg", "", kTProfile2D, {{10, -.8, .8}, ptAxis});
    registry.add<TProfile2D>("v1_eta_even_neg", "", kTProfile2D, {{10, -.8, .8}, ptAxis});

    registry.add<TProfile2D>("v2_cent", "", kTProfile2D, {{80, 0, 80}, ptAxis});
    registry.add<TProfile2D>("v2A_cent", "", kTProfile2D, {{80, 0, 80}, ptAxis});
    registry.add<TProfile2D>("v2C_cent", "", kTProfile2D, {{80, 0, 80}, ptAxis});
    registry.add<TProfile2D>("v2AC_cent", "", kTProfile2D, {{80, 0, 80}, ptAxis});

    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{11, 0, 11}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "occupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(9, "after Mult cuts");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(10, "after Cent cuts");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(11, "isSelected");

    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
    }

    if (cfgUseAdditionalTrackCut) {
      fPhiCutLow = new TF1("fPhiCutLow", "0.06/x+pi/18.0-0.06", 0, 100);
      fPhiCutHigh = new TF1("fPhiCutHigh", "0.1/x+pi/18.0+0.06", 0, 100);
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // static o2::parameters::GRPObject* grpo = nullptr;
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int& multTrk, const float& centrality)
  {
    if (cfgTVXinTRD) {
      if (collision.alias_bit(kTVXinTRD)) {
        // TRD triggered
        // "CMTVX-B-NOPF-TRD,minbias_TVX"
        return 0;
      }
      registry.fill(HIST("hEventCount"), 3.5);
    }

    if (cfgNoSameBunchPileupCut) {
      if (!collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        // rejects collisions which are associated with the same "found-by-T0" bunch crossing
        // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
        return 0;
      }
      registry.fill(HIST("hEventCount"), 4.5);
    }
    if (cfgIsGoodZvtxFT0vsPV) {
      if (!collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
        // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
        // use this cut at low multiplicities with caution
        return 0;
      }
      registry.fill(HIST("hEventCount"), 5.5);
    }
    if (cfgNoCollInTimeRangeStandard) {
      if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        //  Rejection of the collisions which have other events nearby
        return 0;
      }
      registry.fill(HIST("hEventCount"), 6.5);
    }

    if (cfgIsVertexITSTPC) {
      if (!collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
        // selects collisions with at least one ITS-TPC track, and thus rejects vertices built from ITS-only tracks
        return 0;
      }
      registry.fill(HIST("hEventCount"), 7.5);
    }

    float vtxz = -999;
    if (collision.numContrib() > 1) {
      vtxz = collision.posZ();
      float zRes = std::sqrt(collision.covZZ());
      if (zRes > 0.25 && collision.numContrib() < 20)
        vtxz = -999;
    }
    // auto multV0A = collision.multFV0A();
    // auto multT0A = collision.multFT0A();
    // auto multT0C = collision.multFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    if (vtxz > 10 || vtxz < -10)
      return 0;
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;
    if (multTrk < fMultCutLow->Eval(centrality))
      return 0;
    if (multTrk > fMultCutHigh->Eval(centrality))
      return 0;

    registry.fill(HIST("hEventCount"), 8.5);

    if (centrality > cfgCentMax || centrality < cfgCentMin)
      return 0;

    registry.fill(HIST("hEventCount"), 9.5);

    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int& field)
  {
    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, o2::constants::math::PI / 9.0);
    registry.fill(HIST("pt_phi_bef"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    registry.fill(HIST("pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  void process(UsedCollisions::iterator const& collision, aod::BCsWithTimestamps const&, UsedTracks const& tracks)
  {
    // Hier sum over collisions and get ZDC data.
    registry.fill(HIST("hEventCount"), .5);

    if (!collision.sel8())
      return;
    registry.fill(HIST("hEventCount"), 1.5);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    auto field = (cfgMagField == 99999) ? getMagneticField(bc.timestamp()) : cfgMagField;

    auto centrality = collision.centFT0C();
    // auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (!eventSelected(collision, tracks.size(), centrality))
      return;

    if (collision.isSelected()) {
      registry.fill(HIST("hEventCount"), 10.5);

      registry.fill(HIST("hCent"), centrality);

      double qxA = collision.qxA();
      double qyA = collision.qyA();
      double qxC = collision.qxC();
      double qyC = collision.qyC();

      double psiA = 1.0 * std::atan2(qyA, qxA);
      registry.fill(HIST("hSPplaneA"), psiA, 1);

      double psiC = 1.0 * std::atan2(qyC, qxC);
      registry.fill(HIST("hSPplaneC"), psiC, 1);

      registry.fill(HIST("hSPplaneA-C"), psiA - psiC, 1);

      registry.fill(HIST("hCosdPhi"), centrality, std::cos(psiA - psiC));
      if (std::cos(psiA - psiC) < 0)
        registry.fill(HIST("hSPlaneRes"), centrality, std::sqrt(-1. * std::cos(psiA - psiC)));

      registry.fill(HIST("hSindPhi"), centrality, std::sin(psiA - psiC));

      auto qxAqxC = qxA * qxC;
      auto qyAqyC = qyA * qyC;

      for (const auto& track : tracks) {
        if (!trackSelected(track, field))
          continue;

        bool pos;
        if (track.sign() == 0.0)
          continue;
        if (track.sign() > 0) {
          pos = true;
        } else {
          pos = false;
        }

        // constrain angle to 0 -> [0,0+2pi]
        auto phi = RecoDecay::constrainAngle(track.phi(), 0);

        auto ux = std::cos(phi);
        auto uy = std::sin(phi);

        // auto uxQxA = ux * qxA;
        // auto uyQyA = uy * qyA;
        // auto uxyQxyA = uxQxA + uyQyA;
        // auto uxQxC = ux * qxC;
        // auto uyQyC = uy * qyC;
        // auto uxyQxyC = uxQxC + uyQyC;

        auto oddv1 = ux * (qxA - qxC) + uy * (qyA - qyC);
        auto evenv1 = ux * (qxA + qxC) + uy * (qyA + qyC);

        auto oddv1Dev = ux * (qxA - qxC) / std::sqrt(std::abs(qxAqxC)) + uy * (qyA - qyC) / std::sqrt(std::abs(qyAqyC));
        auto evenv1Dev = ux * (qxA + qxC) / std::sqrt(std::abs(qxAqxC)) + uy * (qyA + qyC) / std::sqrt(std::abs(qyAqyC));

        double v1A = std::cos(phi - psiA);
        double v1C = std::cos(phi - psiC);

        double v1AC = std::cos(phi - (psiA - psiC));

        registry.fill(HIST("v1_eta"), track.eta(), track.pt(), (1. / std::sqrt(2)) * (v1A - v1C));
        registry.fill(HIST("v1A_eta"), track.eta(), track.pt(), (v1A));
        registry.fill(HIST("v1C_eta"), track.eta(), track.pt(), (v1C));
        registry.fill(HIST("v1AC_eta"), track.eta(), track.pt(), (v1AC));

        registry.fill(HIST("v1_eta_odd"), track.eta(), track.pt(), oddv1);
        registry.fill(HIST("v1_eta_even"), track.eta(), track.pt(), evenv1);

        registry.fill(HIST("v1_eta_odd_dev"), track.eta(), track.pt(), oddv1Dev);
        registry.fill(HIST("v1_eta_even_dev"), track.eta(), track.pt(), evenv1Dev);

        if (pos) {
          registry.fill(HIST("v1_eta_odd_pos"), track.eta(), track.pt(), oddv1);
          registry.fill(HIST("v1_eta_even_pos"), track.eta(), track.pt(), evenv1);

          registry.fill(HIST("v1_eta_odd_dev_pos"), track.eta(), track.pt(), oddv1Dev);
          registry.fill(HIST("v1_eta_even_dev_pos"), track.eta(), track.pt(), evenv1Dev);
        } else {
          registry.fill(HIST("v1_eta_odd_neg"), track.eta(), track.pt(), oddv1);
          registry.fill(HIST("v1_eta_even_neg"), track.eta(), track.pt(), evenv1);

          registry.fill(HIST("v1_eta_odd_dev_neg"), track.eta(), track.pt(), oddv1Dev);
          registry.fill(HIST("v1_eta_even_dev_neg"), track.eta(), track.pt(), evenv1Dev);
        }

        double v2A = std::cos(2 * (phi - psiA));
        double v2C = std::cos(2 * (phi - psiC));
        double v2AC = std::cos(2 * (phi - (psiA - psiC)));

        registry.fill(HIST("v2_cent"), centrality, track.pt(), (1. / std::sqrt(2)) * (v2A - v2C));
        registry.fill(HIST("v2A_cent"), centrality, track.pt(), (v2A));
        registry.fill(HIST("v2C_cent"), centrality, track.pt(), (v2C));
        registry.fill(HIST("v2AC_cent"), centrality, track.pt(), (v2AC));
      }

      float qIm = collision.qvecIm()[0];
      float qRe = collision.qvecRe()[0];

      registry.fill(HIST("hqIm"), qIm);
      registry.fill(HIST("hqRe"), qRe);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSP>(cfgc),
  };
}
