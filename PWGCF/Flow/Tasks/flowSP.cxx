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
  O2_DEFINE_CONFIGURABLE(cfgCentMax, float, 90, "Maximum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgCentMin, float, 0, "Minimum cenrality for selected events");
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, false, "Fill NUA weights")

  O2_DEFINE_CONFIGURABLE(cfgDoubleTrackFunction, bool, true, "Include track cut at low pt");
  O2_DEFINE_CONFIGURABLE(cfgTrackCutSize, float, 0.06, "Spread of track cut");
  O2_DEFINE_CONFIGURABLE(cfgMaxOccupancy, int, 500, "Maximum occupancy of selected events");
  O2_DEFINE_CONFIGURABLE(cfgNoSameBunchPileupCut, bool, true, "kNoSameBunchPileupCut");
  O2_DEFINE_CONFIGURABLE(cfgIsGoodZvtxFT0vsPV, bool, true, "kIsGoodZvtxFT0vsPV");
  O2_DEFINE_CONFIGURABLE(cfgNoCollInTimeRangeStandard, bool, true, "kNoCollInTimeRangeStandard");
  O2_DEFINE_CONFIGURABLE(cfgDoOccupancySel, bool, true, "Bool for event selection on detector occupancy");
  O2_DEFINE_CONFIGURABLE(cfgMultCut, bool, true, "Use additional evenr cut on mult correlations");
  O2_DEFINE_CONFIGURABLE(cfgTVXinTRD, bool, false, "Use kTVXinTRD (reject TRD triggered events)");
  O2_DEFINE_CONFIGURABLE(cfgIsVertexITSTPC, bool, true, "Selects collisions with at least one ITS-TPC track");


  ConfigurableAxis axisDCAz{"axisDCAz", {200, -.5, .5}, "DCA_{z} (cm)"};
  ConfigurableAxis axisDCAxy{"axisDCAxy", {200, -.5, .5}, "DCA_{xy} (cm)"};
  ConfigurableAxis axisPhiMod = {"axisPhiMod", {100, 0, constants::math::PI / 9}, "fmod(#varphi,#pi/9)"};
  ConfigurableAxis axisPhi = {"axisPhi", {60, 0, constants::math::TwoPI}, "#varphi"};
  ConfigurableAxis axisEta = {"axisEta", {64, -1, 1}, "#eta"}; 
  ConfigurableAxis axisEtaV1 = {"axisEtaV1", {8, -.8, .8}, "#eta"}; 
  ConfigurableAxis axisVz = {"axisVz", {40, -10, 10}, "v_{z}"}; 
  ConfigurableAxis axisCent = {"axisCent", {90,0,90}, "Centrality(%)"}; 
  ConfigurableAxis axisPhiPlane = {"axisPhiPlane", {100, -constants::math::PI,constants::math::PI}, "#Psi"}; 


  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZ;
  Filter trackFilter = nabs(aod::track::eta) < cfgEta && aod::track::pt > cfgPtmin&& aod::track::pt < cfgPtmax && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && nabs(aod::track::dcaXY) < cfgDCAxy&& nabs(aod::track::dcaZ) < cfgDCAz;
  using UsedCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::SPTableZDC>>;
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
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<double> ptbinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
    AxisSpec axisPt = {ptbinning, "#it{p}_{T} GeV/#it{c}"};
    AxisSpec nchAxis = {4000, 0, 4000, "N_{ch}"};
    AxisSpec t0cAxis = {70, 0, 70000, "N_{ch} (T0C)"};
    AxisSpec t0aAxis = {200, 0, 200, "N_{ch}"};
    AxisSpec multpvAxis = {4000, 0, 4000, "N_{ch} (PV)"};

    registry.add<TH1>("hSPplaneA", "hSPplaneA", kTH1D, {axisPhiPlane});
    registry.add<TH1>("hSPplaneC", "hSPplaneC", kTH1D, {axisPhiPlane});
    registry.add<TH1>("hSPplaneFull", "hSPplaneFull", kTH1D, {axisPhiPlane});

    registry.add<TProfile>("hCosPhiACosPhiC", "hCosPhiACosPhiC; Centrality(%); #LT Cos(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hSinPhiASinPhiC", "hSinPhiASinPhiC; Centrality(%); #LT Sin(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hSinPhiACosPhiC", "hSinPhiACosPhiC; Centrality(%); #LT Sin(#Psi^{A})Cos(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hCosPhiASinsPhiC", "hCosPhiASinsPhiC; Centrality(%); #LT Cos(#Psi^{A})Sin(#Psi^{C})#GT", kTProfile, {axisCent});
    registry.add<TProfile>("hFullEvPlaneRes", "hFullEvPlaneRes; Centrality(%); #sqrt{2#LT Cos(#Psi^{A} - #Psi^{C})#GT} ", kTProfile, {axisCent});

    registry.add("QA/after/hCent", "", {HistType::kTH1D, {axisCent}});
    registry.add("QA/after/globalTracks_centT0C", "", {HistType::kTH2D, {axisCent, nchAxis}});
    registry.add("QA/after/PVTracks_centT0C", "", {HistType::kTH2D, {axisCent, multpvAxis}});
    registry.add("QA/after/globalTracks_PVTracks", "", {HistType::kTH2D, {multpvAxis, nchAxis}});
    registry.add("QA/after/globalTracks_multT0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
    registry.add("QA/after/globalTracks_multV0A", "", {HistType::kTH2D, {t0aAxis, nchAxis}});
    registry.add("QA/after/multV0A_multT0A", "", {HistType::kTH2D, {t0aAxis, t0aAxis}});
    registry.add("QA/after/multT0C_centT0C", "", {HistType::kTH2D, {axisCent, t0cAxis}});
    registry.addClone("QA/after/", "QA/before/");
    
    registry.add("QA/after/pt_phi_bef", "", {HistType::kTH2D, {axisPt, axisPhiMod}});
    registry.add("QA/after/pt_phi_aft", "", {HistType::kTH2D, {axisPt, axisPhiMod}});
    registry.add("QA/after/hPhi_Eta_vz", "", {HistType::kTH3D, {axisPhi, axisEta, axisVz}});
    registry.add("QA/after/hDCAxy", "", {HistType::kTH1D, {axisDCAxy}});
    registry.add("QA/after/hDCAz", "", {HistType::kTH1D, {axisDCAz}});

    // track properties per centrality and per eta, pt bin 
    registry.add<TProfile3D>("uxqxA_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyA_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uxqxC_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyC_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1A_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1C_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1Full_eta_pt", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});

    registry.add<TProfile3D>("uxqxA_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyA_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uxqxC_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyC_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1A_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1C_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1Full_eta_pt_pos", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});

    registry.add<TProfile3D>("uxqxA_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyA_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uxqxC_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("uyqyC_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1A_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1C_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});
    registry.add<TProfile3D>("v1Full_eta_pt_neg", "", kTProfile3D, {axisEtaV1, axisPt, axisCent});

    registry.add<TProfile>("qAqCX", "", kTProfile, {axisCent});
    registry.add<TProfile>("qAqCY", "", kTProfile, {axisCent});
    registry.add<TProfile>("qAqCXY", "", kTProfile, {axisCent});

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
    registry.fill(HIST("QA/after/pt_phi_bef"), track.pt(), phimodn);
    if (phimodn < fPhiCutHigh->Eval(track.pt()) && phimodn > fPhiCutLow->Eval(track.pt()))
      return false; // reject track
    registry.fill(HIST("QA/after/pt_phi_aft"), track.pt(), phimodn);
    return true;
  }

  template <typename CollisionObject, typename TracksObject>
  inline void FillEventQA(CollisionObject collision, TracksObject tracks, bool before)
  {
    if(before){
      registry.fill(HIST("QA/before/hCent"), collision.centFT0C());
      registry.fill(HIST("QA/before/globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("QA/before/PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("QA/before/globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("QA/before/globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("QA/before/globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("QA/before/multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("QA/before/multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    } else {
      registry.fill(HIST("QA/after/hCent"), collision.centFT0C());
      registry.fill(HIST("QA/after/globalTracks_centT0C"), collision.centFT0C(), tracks.size());
      registry.fill(HIST("QA/after/PVTracks_centT0C"), collision.centFT0C(), collision.multNTracksPV());
      registry.fill(HIST("QA/after/globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
      registry.fill(HIST("QA/after/globalTracks_multT0A"), collision.multFT0A(), tracks.size());
      registry.fill(HIST("QA/after/globalTracks_multV0A"), collision.multFV0A(), tracks.size());
      registry.fill(HIST("QA/after/multV0A_multT0A"), collision.multFT0A(), collision.multFV0A());
      registry.fill(HIST("QA/after/multT0C_centT0C"), collision.centFT0C(), collision.multFT0C());
    }
    return;
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

    FillEventQA(collision, tracks, true);

    if (collision.isSelected()) {

      registry.fill(HIST("hEventCount"), 10.5);

      FillEventQA(collision, tracks, false);

      double qxA = collision.qxA();
      double qyA = collision.qyA();
      double qxC = collision.qxC();
      double qyC = collision.qyC();

      double psiA = 1.0 * std::atan2(qyA, qxA);
      registry.fill(HIST("hSPplaneA"), psiA, 1);

      double psiC = 1.0 * std::atan2(qyC, qxC);
      registry.fill(HIST("hSPplaneC"), psiC, 1);

      // https://twiki.cern.ch/twiki/pub/ALICE/DirectedFlowAnalysisNote/v1_ZDC_ALICE_INT_NOTE_version02.pdf
      double psiFull = 1.0 * std::atan2(qyA+qyC, qxA + qxC);
      registry.fill(HIST("hSPplaneFull"), psiFull, 1);

      registry.fill(HIST("hCosPhiACosPhiC"), centrality, std::cos(psiA)*std::cos(psiC));
      registry.fill(HIST("hSinPhiASinPhiC"), centrality, std::sin(psiA)*std::sin(psiC));
      registry.fill(HIST("hSinPhiACosPhiC"), centrality, std::sin(psiA)*std::cos(psiC));
      registry.fill(HIST("hCosPhiASinsPhiC"), centrality, std::cos(psiA)*std::sin(psiC));

      if(std::cos(psiA) - std::cos(psiC) < 0)
        registry.fill(HIST("hFullEvPlaneRes"), centrality, std::sqrt(-2. * (std::cos(psiA) - std::cos(psiC))));

      registry.fill(HIST("qAqCXY"), centrality, qxA * qxC + qyA * qyC);
      registry.fill(HIST("qAqCX"), centrality, qxA * qxC);
      registry.fill(HIST("qAqCY"), centrality, qyA * qyC);

      for (const auto& track : tracks) {
        if (!trackSelected(track, field))
          continue;
        
        if (track.sign() == 0.0)
          continue;
        bool pos = (track.sign() > 0) ? true : false;

        // constrain angle to 0 -> [0,0+2pi]
        auto phi = RecoDecay::constrainAngle(track.phi(), 0);

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        auto ux = std::cos(phi);
        auto uy = std::sin(phi);
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        double uxqxA = ux*qxA;
        double uyqyA = uy*qyA;

        registry.fill(HIST("uxqxA_eta_pt"), track.eta(), track.pt(), centrality, uxqxA);
        registry.fill(HIST("uyqyA_eta_pt"), track.eta(), track.pt(), centrality, uyqyA);

        double uxqxC = ux*qxC;
        double uyqyC = uy*qyC;

        registry.fill(HIST("uxqxC_eta_pt"), track.eta(), track.pt(), centrality, uxqxC);
        registry.fill(HIST("uyqyC_eta_pt"), track.eta(), track.pt(), centrality, uyqyC);
       
        double v1A = std::cos(phi - psiA);
        double v1C = std::cos(phi - psiC);
        double v1Full = std::cos(phi - psiFull);

        registry.fill(HIST("v1A_eta_pt"), track.eta(), track.pt(), centrality, v1A);
        registry.fill(HIST("v1C_eta_pt"), track.eta(), track.pt(), centrality, v1C);
        registry.fill(HIST("v1Full_eta_pt"), track.eta(), track.pt(), centrality, v1Full);

        if (pos) {
          registry.fill(HIST("uxqxA_eta_pt_pos"), track.eta(), track.pt(), centrality, uxqxA);
          registry.fill(HIST("uyqyA_eta_pt_pos"), track.eta(), track.pt(), centrality, uyqyA);
          registry.fill(HIST("uxqxC_eta_pt_pos"), track.eta(), track.pt(), centrality, uxqxC);
          registry.fill(HIST("uyqyC_eta_pt_pos"), track.eta(), track.pt(), centrality, uyqyC);
          registry.fill(HIST("v1A_eta_pt_pos"), track.eta(), track.pt(), centrality, v1A);
          registry.fill(HIST("v1C_eta_pt_pos"), track.eta(), track.pt(), centrality, v1C);
          registry.fill(HIST("v1Full_eta_pt_pos"), track.eta(), track.pt(), centrality, v1Full);
        } else {
          registry.fill(HIST("uxqxA_eta_pt_neg"), track.eta(), track.pt(), centrality, uxqxA);
          registry.fill(HIST("uyqyA_eta_pt_neg"), track.eta(), track.pt(), centrality, uyqyA);
          registry.fill(HIST("uxqxC_eta_pt_neg"), track.eta(), track.pt(), centrality, uxqxC);
          registry.fill(HIST("uyqyC_eta_pt_neg"), track.eta(), track.pt(), centrality, uyqyC);
          registry.fill(HIST("v1A_eta_pt_neg"), track.eta(), track.pt(), centrality, v1A);
          registry.fill(HIST("v1C_eta_pt_neg"), track.eta(), track.pt(), centrality, v1C);
          registry.fill(HIST("v1Full_eta_pt_neg"), track.eta(), track.pt(), centrality, v1Full);
        }
        //QA plots
        registry.fill(HIST("QA/after/hPhi_Eta_vz"), track.phi(), track.eta(), collision.posZ());
        registry.fill(HIST("QA/after/hDCAxy"), track.dcaXY());
        registry.fill(HIST("QA/after/hDCAz"), track.dcaZ());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowSP>(cfgc),
  };
}

