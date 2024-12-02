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

// author: Noor Koster noor.koster@cern.ch

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

#include "PWGCF/DataModel/SPTableZDC.h"
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct flowAnalysisSP {

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
  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::SPTableZDC, aod::Qvectors>>;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;

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

    registry.add<TProfile>("v1_eta", "", kTProfile, {{10, -.8, .8}});
    registry.add<TProfile>("v1A_eta", "", kTProfile, {{10, -.8, .8}});
    registry.add<TProfile>("v1C_eta", "", kTProfile, {{10, -.8, .8}});
    registry.add<TProfile>("v1AC_eta", "", kTProfile, {{10, -.8, .8}});

    registry.add<TProfile>("v2_cent", "", kTProfile, {{80, 0, 80}});
    registry.add<TProfile>("v2A_cent", "", kTProfile, {{80, 0, 80}});
    registry.add<TProfile>("v2C_cent", "", kTProfile, {{80, 0, 80}});
    registry.add<TProfile>("v2AC_cent", "", kTProfile, {{80, 0, 80}});

    registry.add("hEventCount", "Number of Event;; Count", {HistType::kTH1D, {{10, 0, 10}}});
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(3, "occupancy");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(4, "kTVXinTRD");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(5, "kNoSameBunchPileup");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(6, "kIsGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeStandard");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(8, "kIsVertexITSTPC");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(9, "after Mult cuts");
    registry.get<TH1>(HIST("hEventCount"))->GetXaxis()->SetBinLabel(10, "isSelected");

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
      float zRes = TMath::Sqrt(collision.covZZ());
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

    return 1;
  }

  template <typename TTrack>
  bool trackSelected(TTrack track, const int& field)
  {
    double phimodn = track.phi();
    if (field < 0) // for negative polarity field
      phimodn = TMath::TwoPi() - phimodn;
    if (track.sign() < 0) // for negative charge
      phimodn = TMath::TwoPi() - phimodn;
    if (phimodn < 0)
      LOGF(warning, "phi < 0: %g", phimodn);

    phimodn += TMath::Pi() / 18.0; // to center gap in the middle
    phimodn = fmod(phimodn, TMath::Pi() / 9.0);
    registry.fill(HIST("pt_phi_bef"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    registry.fill(HIST("pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  void process(myCollisions::iterator const& collision, aod::BCsWithTimestamps const&, myTracks const& tracks)
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
      registry.fill(HIST("hEventCount"), 9.5);

      registry.fill(HIST("hCent"), centrality);

      double qxA = collision.qxA();
      double qyA = collision.qyA();
      double qxC = collision.qxC();
      double qyC = collision.qyC();

      double Psi_A = 1.0 * TMath::ATan2(qyA, qxA);
      registry.fill(HIST("hSPplaneA"), Psi_A, 1);

      double Psi_C = 1.0 * TMath::ATan2(qyC, qxC);
      registry.fill(HIST("hSPplaneC"), Psi_C, 1);

      registry.fill(HIST("hSPplaneA-C"), Psi_A - Psi_C, 1);

      registry.fill(HIST("hCosdPhi"), centrality, TMath::Cos(Psi_A - Psi_C));
      if (TMath::Cos(Psi_A - Psi_C) < 0)
        registry.fill(HIST("hSPlaneRes"), centrality, TMath::Sqrt(-1. * TMath::Cos(Psi_A - Psi_C)));
      registry.fill(HIST("hSindPhi"), centrality, TMath::Sin(Psi_A - Psi_C));

      for (auto& track : tracks) {
        if (!trackSelected(track, field))
          continue;

        double v1A = TMath::Cos(track.phi() - Psi_A);
        double v1C = TMath::Cos(track.phi() - Psi_C);

        double v1AC = TMath::Cos(track.phi() - (Psi_A - Psi_C));

        registry.fill(HIST("v1_eta"), track.eta(), (1. / TMath::Sqrt(2)) * (v1A - v1C));
        registry.fill(HIST("v1A_eta"), track.eta(), (v1A));
        registry.fill(HIST("v1C_eta"), track.eta(), (v1C));
        registry.fill(HIST("v1AC_eta"), track.eta(), (v1AC));

        double v2A = TMath::Cos(2 * (track.phi() - Psi_A));
        double v2C = TMath::Cos(2 * (track.phi() - Psi_C));
        double v2AC = TMath::Cos(2 * (track.phi() - (Psi_A - Psi_C)));

        registry.fill(HIST("v2_cent"), centrality, (1. / TMath::Sqrt(2)) * (v2A - v2C));
        registry.fill(HIST("v2A_cent"), centrality, (v2A));
        registry.fill(HIST("v2C_cent"), centrality, (v2C));
        registry.fill(HIST("v2AC_cent"), centrality, (v2AC));
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
    adaptAnalysisTask<flowAnalysisSP>(cfgc),
  };
}
