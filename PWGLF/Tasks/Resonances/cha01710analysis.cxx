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

/// \file cha01710analysis.cxx
/// \brief charged a01710 resonance analysis
/// \author Junlee Kim (jikim1290@gmail.com)

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TRandom.h>
#include <TVector2.h>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Cha01710analysis {
  enum CentSel {
    kFT0C = 0,
    kFT0M
  };

  enum V0MassRegion {
    kReject = 0,
    kSignal,
    kSideband
  };

  TRandom* rn = new TRandom();

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> cfgUrl{"cfgUrl", "http://alice-ccdb.cern.ch", "CCDB URL"};

  struct : ConfigurableGroup {
    Configurable<int> cfgCentEst{"cfgCentEst", 1, "0: FT0C, 1: FT0M"};
    Configurable<float> cfgZVertexMax{"cfgZVertexMax", 10.f, "maximum |zPV| (cm)"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8"};
    Configurable<bool> cfgRequireGoodZvtx{"cfgRequireGoodZvtx", true, "require kIsGoodZvtxFT0vsPV"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", true, "require kNoSameBunchPileup"};
  } eventCuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0.f, "minimum K0S pT"};
    Configurable<float> cfgV0PtMax{"cfgV0PtMax", 100.f, "maximum K0S pT"};
    Configurable<float> cfgV0RapidityMax{"cfgV0RapidityMax", 0.5f, "maximum |y(K0S)|"};
    Configurable<float> cfgV0DcaDaughtersMax{"cfgV0DcaDaughtersMax", 1.f, "maximum daughter DCA"};
    Configurable<float> cfgV0CosPAMin{"cfgV0CosPAMin", 0.97f, "minimum V0 cosine pointing angle"};
    Configurable<float> cfgV0RadiusMin{"cfgV0RadiusMin", 0.5f, "minimum V0 radius (cm)"};
    Configurable<float> cfgV0CtauMax{"cfgV0CtauMax", 15.f, "maximum K0S c tau (cm)"};
    Configurable<float> cfgV0DaughterDcaPVMin{"cfgV0DaughterDcaPVMin", 0.06f, "minimum daughter |DCA to PV|"};
    Configurable<float> cfgV0DaughterEtaMax{"cfgV0DaughterEtaMax", 0.8f, "maximum daughter |eta|"};
    Configurable<float> cfgV0DaughterPtMin{"cfgV0DaughterPtMin", 0.1f, "minimum daughter pT"};
    Configurable<int> cfgV0DaughterTPCNClsMin{"cfgV0DaughterTPCNClsMin", 70, "minimum daughter TPC clusters"};
    Configurable<float> cfgV0DaughterTPCNSigmaPiMax{"cfgV0DaughterTPCNSigmaPiMax", 5.f, "maximum daughter |TPC nSigma(pi)|"};
    Configurable<float> cfgV0DaughterTOFNSigmaPiMax{"cfgV0DaughterTOFNSigmaPiMax", 5.f, "maximum daughter |TOF nSigma(pi)|"};
    Configurable<float> cfgKsMassWindow{"cfgKsMassWindow", 0.01, "K short mass window"};
    Configurable<bool> cfgRejectLambda{"cfgRejectLambda", true, "reject Lambda/anti-Lambda competitors"};
    Configurable<float> cfgLambdaMassWindow{"cfgLambdaMassWindow", 0.005f, "Lambda rejection window (GeV/c2)"};
  } v0Cuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "minimum charged-kaon pT"};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "maximum charged-kaon |eta|"};
    Configurable<int> cfgTPCCrossedRowsMin{"cfgTPCCrossedRowsMin", 70, "minimum crossed TPC rows"};
    Configurable<float> cfgCrossedRowsRatioMin{"cfgCrossedRowsRatioMin", 0.8f, "minimum crossed/findable ratio"};
    Configurable<float> cfgTPCChi2Max{"cfgTPCChi2Max", 4.f, "maximum TPC chi2/cluster"};
    Configurable<float> cfgITSChi2Max{"cfgITSChi2Max", 36.f, "maximum ITS chi2/cluster"};
    Configurable<float> cfgDcaXYMax{"cfgDcaXYMax", 0.5f, "maximum |DCAxy| (cm)"};
    Configurable<float> cfgDcaZMax{"cfgDcaZMax", 2.f, "maximum |DCAz| (cm)"};
    Configurable<float> cfgTPCNSigmaMax{"cfgTPCNSigmaMax", 3.f, "maximum |TPC nSigma(K)| without TOF"};
    Configurable<float> cfgCombinedNSigmaMax{"cfgCombinedNSigmaMax", 3.f, "maximum combined TPC-TOF nSigma(K)"};
  } kaonCuts;

  struct : ConfigurableGroup {
    Configurable<bool> cfgRequireQVector{"cfgRequireQVector", true, "reject events with an invalid Q-vector"};
    Configurable<int> cfgQvecHarmonic{"cfgQvecHarmonic", 2, "event-plane harmonic"};
    Configurable<int> cfgQvecNumDetectors{"cfgQvecNumDetectors", 7, "number of detectors in the Q-vector table"};
    Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0M", "detector used for the event plane"};
    Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "event-plane reference A"};
    Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "event-plane reference B"};
    Configurable<float> cfgQvecAmplitudeMin{"cfgQvecAmplitudeMin", 1.e-4f, "minimum accepted Q-vector amplitude"};
  } epConfig;

  Configurable<float> cfgMotherRapidityMax{"cfgMotherRapidityMax", 0.5f, "maximum |y(K0S K)|"};
  Configurable<int> cfgNRotations{"cfgNRotations", 3, "number of deterministic rotations per candidate"};

  ConfigurableAxis cfgAxisMass{"cfgAxisMass", {300, 1.0, 2.5}, "M(K^{0}_{S}K^{#pm}) (GeV/c^{2})"};
  ConfigurableAxis cfgAxisPt{"cfgAxisPt", {200, 0., 20.}, "pT(K^{0}_{S}K^{#pm}) (GeV/c)"};
  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100.}, "FT0M centrality (%)"};
  ConfigurableAxis cfgAxisEP{"cfgAxisEP", {6, 0., constants::math::TwoPI}, "n#Delta#varphi = n(#varphi-#Psi_{n})"};

  Filter collisionFilter = nabs(aod::collision::posZ) < eventCuts.cfgZVertexMax;
  Filter trackFilter = nabs(aod::track::eta) < kaonCuts.cfgEtaMax && aod::track::pt > kaonCuts.cfgPtMin;

  using Collisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::Qvectors>>;
  using Tracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                         aod::pidTPCFullPi, aod::pidTOFFullPi,
                                         aod::pidTPCFullKa, aod::pidTOFFullKa>>;
  using V0s = soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas>;

  int detId = 0;
  int refAId = 4;
  int refBId = 5;
  int minHarm = 2;

  float centrality = -1.;
  float rotmin = 5.0 / 6.0;
  float rotmax = 7.0 / 6.0;

  template <typename T>
  int getDetId(T const& name)
  {
    if (name.value == "FT0C") {
      return 0;
    }
    if (name.value == "FT0A") {
      return 1;
    }
    if (name.value == "FT0M") {
      return 2;
    }
    if (name.value == "FV0A") {
      return 3;
    }
    if (name.value == "TPCpos") {
      return 4;
    }
    if (name.value == "TPCneg") {
      return 5;
    }
    LOGF(warning, "Unknown Q-vector detector %s; using FT0C", name.value.c_str());
    return 0;
  }

  void init(InitContext const&)
  {
    histos.add("Event/hCutFlow", "event cut flow", HistType::kTH1F, {{6, -0.5, 5.5}});
    histos.add("Event/hCentDist", "", HistType::kTH1F, {{120, 0, 120}});
    histos.add("Event/hPVzDist", "", HistType::kTH1F, {{200, -12, 12}});

    histos.add("EventPlane/hPsiDet", "event plane, detector", HistType::kTH2F, {cfgAxisCent, {180, -constants::math::PI, constants::math::PI}});
    histos.add("EventPlane/hPsiRefA", "event plane, reference A", HistType::kTH2F, {cfgAxisCent, {180, -constants::math::PI, constants::math::PI}});
    histos.add("EventPlane/hPsiRefB", "event plane, reference B", HistType::kTH2F, {cfgAxisCent, {180, -constants::math::PI, constants::math::PI}});
    histos.add("EventPlane/hResolutionDetRefA", "cos[n(PsiDet-PsiRefA)]", HistType::kTH2F, {cfgAxisCent, {102, -1.02, 1.02}});
    histos.add("EventPlane/hResolutionDetRefB", "cos[n(PsiDet-PsiRefB)]", HistType::kTH2F, {cfgAxisCent, {102, -1.02, 1.02}});
    histos.add("EventPlane/hResolutionRefARefB", "cos[n(PsiRefA-PsiRefB)]", HistType::kTH2F, {cfgAxisCent, {102, -1.02, 1.02}});

    histos.add("V0/hMassBefore", "K0S mass before selection", HistType::kTH2F, {{200, 0.45, 0.55}, cfgAxisPt});
    histos.add("V0/hMassSelected", "selected K0S signal", HistType::kTH2F, {{200, 0.45, 0.55}, cfgAxisPt});

    histos.add("Kaon/hTPCNSigma", "charged kaon TPC PID", HistType::kTH2F, {cfgAxisPt, {120, -6., 6.}});
    histos.add("Kaon/hTOFNSigma", "charged kaon TOF PID", HistType::kTH2F, {cfgAxisPt, {120, -6., 6.}});

    histos.add("Pair/hSignalPlus", "K0S K+ same-event versus EP", HistType::kTHnSparseF, {cfgAxisMass, cfgAxisPt, cfgAxisCent, cfgAxisEP});
    histos.add("Pair/hSignalMinus", "K0S K- same-event versus EP", HistType::kTHnSparseF, {cfgAxisMass, cfgAxisPt, cfgAxisCent, cfgAxisEP});
    histos.add("Pair/hSidebandPlus", "K0S-sideband K+ versus EP", HistType::kTHnSparseF, {cfgAxisMass, cfgAxisPt, cfgAxisCent, cfgAxisEP});
    histos.add("Pair/hSidebandMinus", "K0S-sideband K- versus EP", HistType::kTHnSparseF, {cfgAxisMass, cfgAxisPt, cfgAxisCent, cfgAxisEP});
    histos.add("Pair/hRotatedPlus", "rotated K+ background versus EP", HistType::kTH3F, {cfgAxisMass, cfgAxisPt, cfgAxisCent});
    histos.add("Pair/hRotatedMinus", "rotated K- background versus EP", HistType::kTH3F, {cfgAxisMass, cfgAxisPt, cfgAxisCent});
    histos.add("Pair/hMassVsK0SMass", "pair mass versus reconstructed K0S mass", HistType::kTH2F, {cfgAxisMass, {200, 0.43, 0.57}});

    detId = getDetId(epConfig.cfgQvecDetName);
    refAId = getDetId(epConfig.cfgQvecRefAName);
    refBId = getDetId(epConfig.cfgQvecRefBName);
    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(warning, "Q-vector detector names are not distinct; falling back to FT0M, TPCpos, TPCneg");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }
    if (epConfig.cfgQvecHarmonic < minHarm) {
      LOGF(fatal, "cfgQvecHarmonic must be >= 2");
    }

    ccdb->setURL(cfgUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  template <typename C>
  bool selectEvent(C const& collision)
  {
    histos.fill(HIST("Event/hCutFlow"), 0.);
    if (eventCuts.cfgRequireSel8 && !collision.sel8()) {
      return false;
    }
    histos.fill(HIST("Event/hCutFlow"), 1.);
    if (eventCuts.cfgRequireGoodZvtx && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("Event/hCutFlow"), 2.);
    if (eventCuts.cfgRequireNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("Event/hCutFlow"), 3.);

    if (epConfig.cfgRequireQVector &&
        (collision.qvecAmp()[detId] < epConfig.cfgQvecAmplitudeMin ||
         collision.qvecAmp()[refAId] < epConfig.cfgQvecAmplitudeMin ||
         collision.qvecAmp()[refBId] < epConfig.cfgQvecAmplitudeMin)) {
      return false;
    }
    histos.fill(HIST("Event/hCutFlow"), 4.);
    return true;
  }

  template <typename T>
  bool selectKaon(T const& track)
  {
    if (!track.isPVContributor() || !track.isGlobalTrackWoDCA() || !track.isPrimaryTrack()) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < kaonCuts.cfgTPCCrossedRowsMin ||
        track.tpcCrossedRowsOverFindableCls() < kaonCuts.cfgCrossedRowsRatioMin ||
        track.tpcChi2NCl() > kaonCuts.cfgTPCChi2Max || track.itsChi2NCl() > kaonCuts.cfgITSChi2Max ||
        std::abs(track.dcaXY()) > kaonCuts.cfgDcaXYMax || std::abs(track.dcaZ()) > kaonCuts.cfgDcaZMax) {
      return false;
    }
    histos.fill(HIST("Kaon/hTPCNSigma"), track.pt(), track.tpcNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST("Kaon/hTOFNSigma"), track.pt(), track.tofNSigmaKa());
      return std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()) < kaonCuts.cfgCombinedNSigmaMax;
    }
    return std::abs(track.tpcNSigmaKa()) < kaonCuts.cfgTPCNSigmaMax;
  }

  template <typename T>
  bool selectPionDaughter(T const& track)
  {
    if (!track.hasTPC() || track.tpcNClsFound() < v0Cuts.cfgV0DaughterTPCNClsMin ||
        track.pt() < v0Cuts.cfgV0DaughterPtMin || std::abs(track.eta()) > v0Cuts.cfgV0DaughterEtaMax ||
        std::abs(track.tpcNSigmaPi()) > v0Cuts.cfgV0DaughterTPCNSigmaPiMax) {
      return false;
    }
    return true;
  }

  template <typename C, typename V>
  V0MassRegion selectV0(C const& collision, V const& v0)
  {
    histos.fill(HIST("V0/hMassBefore"), v0.mK0Short(), v0.pt());
    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * constants::physics::MassK0Short;
    if (v0.pt() < v0Cuts.cfgV0PtMin || v0.pt() > v0Cuts.cfgV0PtMax ||
        v0.dcaV0daughters() > v0Cuts.cfgV0DcaDaughtersMax || v0.v0cosPA() < v0Cuts.cfgV0CosPAMin ||
        v0.v0radius() < v0Cuts.cfgV0RadiusMin || std::abs(ctau) > v0Cuts.cfgV0CtauMax) {
      return V0MassRegion::kReject;
    }
    if (v0Cuts.cfgRejectLambda &&
        (std::abs(v0.mLambda() - constants::physics::MassLambda0) < v0Cuts.cfgLambdaMassWindow ||
         std::abs(v0.mAntiLambda() - constants::physics::MassLambda0) < v0Cuts.cfgLambdaMassWindow)) {
      return V0MassRegion::kReject;
    }
    auto pos = v0.template posTrack_as<Tracks>();
    auto neg = v0.template negTrack_as<Tracks>();
    if (!selectPionDaughter(pos) || !selectPionDaughter(neg)) {
      return V0MassRegion::kReject;
    }
    ROOT::Math::PxPyPzMVector k0(v0.px(), v0.py(), v0.pz(), constants::physics::MassK0Short);
    if (std::abs(k0.Rapidity()) > v0Cuts.cfgV0RapidityMax) {
      return V0MassRegion::kReject;
    }
    float dm = std::abs(v0.mK0Short() - constants::physics::MassK0Short);
    if (dm < v0Cuts.cfgKsMassWindow) {
      histos.fill(HIST("V0/hMassSelected"), v0.mK0Short(), v0.pt());
      return V0MassRegion::kSignal;
    }
    if (dm > v0Cuts.cfgKsMassWindow) {
      histos.fill(HIST("V0/hMassSelected"), v0.mK0Short(), v0.pt());
      return V0MassRegion::kSideband;
    }
    return V0MassRegion::kReject;
  }

  template <typename CollisionType, typename TracksType, typename V0Type>
  void fillHistograms(CollisionType const& collision, TracksType const& dTracks, V0Type const& v0s)
  {
    int harmonic = epConfig.cfgQvecHarmonic;
    int qVecDetInd = detId * 4 + 3 + (harmonic - 2) * epConfig.cfgQvecNumDetectors * 4;
    int qVecRefAInd = refAId * 4 + 3 + (harmonic - 2) * epConfig.cfgQvecNumDetectors * 4;
    int qVecRefBInd = refBId * 4 + 3 + (harmonic - 2) * epConfig.cfgQvecNumDetectors * 4;

    float eventPlaneDet = std::atan2(collision.qvecIm()[qVecDetInd], collision.qvecRe()[qVecDetInd]) / static_cast<double>(harmonic);
    float eventPlaneRefA = std::atan2(collision.qvecIm()[qVecRefAInd], collision.qvecRe()[qVecRefAInd]) / static_cast<double>(harmonic);
    float eventPlaneRefB = std::atan2(collision.qvecIm()[qVecRefBInd], collision.qvecRe()[qVecRefBInd]) / static_cast<double>(harmonic);

    histos.fill(HIST("EventPlane/hPsiDet"), centrality, eventPlaneDet);
    histos.fill(HIST("EventPlane/hPsiRefA"), centrality, eventPlaneRefA);
    histos.fill(HIST("EventPlane/hPsiRefB"), centrality, eventPlaneRefB);
    histos.fill(HIST("EventPlane/hResolutionDetRefA"), centrality, std::cos(harmonic * (eventPlaneDet - eventPlaneRefA)));
    histos.fill(HIST("EventPlane/hResolutionDetRefB"), centrality, std::cos(harmonic * (eventPlaneDet - eventPlaneRefB)));
    histos.fill(HIST("EventPlane/hResolutionRefARefB"), centrality, std::cos(harmonic * (eventPlaneRefA - eventPlaneRefB)));

    for (const auto& v0 : v0s) {
      auto region = selectV0(collision, v0);
      if (region == kReject) {
        continue;
      }
      auto pos = v0.template posTrack_as<Tracks>();
      auto neg = v0.template negTrack_as<Tracks>();

      ROOT::Math::PxPyPzMVector k0(v0.px(), v0.py(), v0.pz(), constants::physics::MassK0Short);

      for (const auto& track : dTracks) {
        if (track.globalIndex() == pos.globalIndex() || track.globalIndex() == neg.globalIndex() || !selectKaon(track)) {
          continue;
        }
        ROOT::Math::PxPyPzMVector kaon(track.px(), track.py(), track.pz(), constants::physics::MassKaonCharged);
        auto mother = k0 + kaon;
        if (std::abs(mother.Rapidity()) > cfgMotherRapidityMax) {
          continue;
        }

        float relPhi = TVector2::Phi_0_2pi((mother.Phi() - eventPlaneDet) * harmonic);
        histos.fill(HIST("Pair/hMassVsK0SMass"), mother.M(), v0.mK0Short());
        if (region != kSignal)
          continue;
        if (track.sign() > 0) {
          histos.fill(HIST("Pair/hSignalPlus"), mother.M(), mother.Pt(), centrality, relPhi);
        } else if (track.sign() < 0) {
          histos.fill(HIST("Pair/hSignalMinus"), mother.M(), mother.Pt(), centrality, relPhi);
        }
        for (int i = 0; i < cfgNRotations; ++i) {
          auto randomPhi = rn->Uniform(o2::constants::math::PI * rotmin, o2::constants::math::PI * rotmax);
          randomPhi += kaon.Phi();
          auto kaonRot = ROOT::Math::PxPyPzMVector(kaon.Pt() * std::cos(randomPhi), kaon.Pt() * std::sin(randomPhi), track.pz(), o2::constants::physics::MassKaonCharged);
          auto motherRot = k0 + kaonRot;
          if (std::abs(motherRot.Rapidity()) > cfgMotherRapidityMax)
            continue;

          if (track.sign() > 0) {
            histos.fill(HIST("Pair/hRotatedPlus"), motherRot.M(), motherRot.Pt(), centrality);
          } else if (track.sign() < 0) {
            histos.fill(HIST("Pair/hRotatedMinus"), motherRot.M(), motherRot.Pt(), centrality);
          }
        }
      }
    }
  }

  void processDataSame(Collisions::iterator const& collision, Tracks const& tracks, V0s const& v0s)
  {
    if (!selectEvent(collision)) {
      return;
    }
    if (eventCuts.cfgCentEst == kFT0M) {
      centrality = collision.centFT0M();
    } else if (eventCuts.cfgCentEst == kFT0C) {
      centrality = collision.centFT0C();
    } else {
      centrality = collision.centFT0M();
    }

    histos.fill(HIST("Event/hCentDist"), centrality);
    histos.fill(HIST("Event/hPVzDist"), collision.posZ());

    fillHistograms(collision, tracks, v0s);
  };
  PROCESS_SWITCH(Cha01710analysis, processDataSame, "Process Event for data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Cha01710analysis>(cfgc)};
}
