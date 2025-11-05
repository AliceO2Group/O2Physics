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

/// \author Junlee Kim (jikim1290@gmail.com)

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <TMath.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambdaTwoPartPolarization {
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using V0TrackCandidate = aod::V0Datas;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<float> cfgCentSel{"cfgCentSel", 100., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<bool> cfgEvtSel{"cfgEvtSel", true, "event selection flag"};
  Configurable<bool> cfgPVSel{"cfgPVSel", true, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 10.0, "Additional PV selection range for syst"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", true, "flag for additional pileup selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPrToPVMin{"cfgDCAPrToPVMin", 0.05, "minimum DCA to PV for proton track"};
  Configurable<float> cfgDCAPiToPVMin{"cfgDCAPiToPVMin", 0.1, "minimum DCA to PV for pion track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.99, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 50, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};

  Configurable<float> cfgHypMassWindow{"cfgHypMassWindow", 0.005, "single lambda mass selection"};

  Configurable<bool> cfgEffCor{"cfgEffCor", false, "flag to apply efficiency correction"};
  Configurable<std::string> cfgEffCorPath{"cfgEffCorPath", "", "path for pseudo efficiency correction"};

  Configurable<bool> cfgAccCor{"cfgAccCor", false, "flag to apply acceptance correction"};
  Configurable<std::string> cfgAccCorPath{"cfgAccCorPath", "", "path for pseudo acceptance correction"};

  Configurable<bool> cfgRotBkg{"cfgRotBkg", true, "flag to construct rotational backgrounds"};
  Configurable<int> cfgNRotBkg{"cfgNRotBkg", 10, "the number of rotational backgrounds"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 10, "Number of mixed events per event"};

  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 100}, "Centrality interval"};
  ConfigurableAxis RapAxis{"RapAxis", {10, -0.5, 0.5}, "Rapidity axis"};
  ConfigurableAxis detaAxis{"dyAxis", {20, -1, 1}, "relative rapidity axis"};
  ConfigurableAxis dphiAxis{"dphiAxis", {20, -constants::math::PI * 0.5, constants::math::PI * 1.5}, "relative azimuth axis"};

  ConfigurableAxis cosSigAxis{"cosSigAxis", {110, -1.05, 1.05}, "Signal cosine axis"};
  ConfigurableAxis cosAccAxis{"cosAccAxis", {110, -7.05, 7.05}, "Accepatance cosine axis"};

  ConfigurableAxis vertexAxis{"vertexAxis", {5, -10, 10}, "vertex axis for mixing"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  float centrality;
  float dphi;
  float weight;

  TProfile2D* EffMap = nullptr;
  TProfile2D* AccMap = nullptr;

  void init(o2::framework::InitContext&)
  {
    AxisSpec centQaAxis = {100, 0.0, 100.0};
    AxisSpec PVzQaAxis = {300, -15.0, 15.0};
    AxisSpec epAxis = {6, 0.0, 2.0 * constants::math::PI};

    AxisSpec pidAxis = {100, -10, 10};

    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {20, 0, 20, "basis"};

    if (cfgQAv0) {
      histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});
      histos.add("QA/PVzDist", "", {HistType::kTH1F, {PVzQaAxis}});

      histos.add("QA/nsigma_tpc_pt_ppr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_ppi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
    }

    histos.add("Ana/Signal", "", {HistType::kTHnSparseF, {ptAxis, ptAxis, detaAxis, dphiAxis, centAxis, cosSigAxis}});
    histos.add("Ana/Acceptance", "", {HistType::kTHnSparseF, {ptAxis, centAxis, RapAxis, cosAccAxis}});

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  ROOT::Math::PxPyPzMVector ProtonVec1, PionVec1, LambdaVec1, ProtonBoostedVec1, PionBoostedVec1;
  ROOT::Math::PxPyPzMVector ProtonVec2, PionVec2, LambdaVec2, ProtonBoostedVec2, PionBoostedVec2;
  int V01Tag;
  int V02Tag;
  double costhetastar1;
  double costhetastar2;

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }
    if (cfgCentSel < centrality) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return 0;
    }
    if (cfgAddEvtSelPileup && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TCollision, typename V0>
  bool SelectionV0(TCollision const& collision, V0 const& candidate, int LambdaTag)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (LambdaTag) {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPrToPVMin)
        return false;
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPiToPVMin)
        return false;
    } else if (!LambdaTag) {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPiToPVMin)
        return false;
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPrToPVMin)
        return false;
    }
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.yLambda() < cfgV0EtaMin)
      return false;
    if (candidate.yLambda() > cfgV0EtaMax)
      return false;
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime)
      return false;

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  template <typename C1, typename C2, typename V01, typename V02>
  void FillHistograms(C1 const& c1, C2 const& c2, V01 const& V01s, V02 const& V02s)
  {
    for (auto& v01 : V01s) {
      auto postrack_v01 = v01.template posTrack_as<TrackCandidates>();
      auto negtrack_v01 = v01.template negTrack_as<TrackCandidates>();

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack_v01, 0) && isSelectedV0Daughter(negtrack_v01, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack_v01, 0) && isSelectedV0Daughter(postrack_v01, 1)) {
        aLambdaTag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!SelectionV0(c1, v01, LambdaTag))
        continue;

      if (LambdaTag) {
        ProtonVec1 = ROOT::Math::PxPyPzMVector(v01.pxpos(), v01.pypos(), v01.pzpos(), massPr);
        PionVec1 = ROOT::Math::PxPyPzMVector(v01.pxneg(), v01.pyneg(), v01.pzneg(), massPi);
        V01Tag = 0;
      }
      if (aLambdaTag) {
        ProtonVec1 = ROOT::Math::PxPyPzMVector(v01.pxneg(), v01.pyneg(), v01.pzneg(), massPr);
        PionVec1 = ROOT::Math::PxPyPzMVector(v01.pxpos(), v01.pypos(), v01.pzpos(), massPi);
        V01Tag = 1;
      }
      LambdaVec1 = ProtonVec1 + PionVec1;
      LambdaVec1.SetM(massLambda);

      ROOT::Math::Boost boost1{LambdaVec1.BoostToCM()};
      ProtonBoostedVec1 = boost1(ProtonVec1);

      costhetastar1 = ProtonBoostedVec1.Pz() / ProtonBoostedVec1.P();

      histos.fill(HIST("Ana/Acceptance"), v01.pt(), centrality, v01.yLambda(), costhetastar1 * costhetastar1);

      for (auto& v02 : V02s) {
        if (v01.v0Id() <= v02.v0Id() && doprocessDataSame)
          continue;
        auto postrack_v02 = v02.template posTrack_as<TrackCandidates>();
        auto negtrack_v02 = v02.template negTrack_as<TrackCandidates>();

        LambdaTag = 0;
        aLambdaTag = 0;

        if (isSelectedV0Daughter(postrack_v02, 0) && isSelectedV0Daughter(negtrack_v02, 1)) {
          LambdaTag = 1;
        }
        if (isSelectedV0Daughter(negtrack_v02, 0) && isSelectedV0Daughter(postrack_v02, 1)) {
          aLambdaTag = 1;
        }

        if (LambdaTag == aLambdaTag)
          continue;

        if (!SelectionV0(c2, v02, LambdaTag))
          continue;

        if (doprocessDataSame) {
          if (postrack_v01.globalIndex() == postrack_v02.globalIndex() || postrack_v01.globalIndex() == negtrack_v02.globalIndex() || negtrack_v01.globalIndex() == postrack_v02.globalIndex() || negtrack_v01.globalIndex() == negtrack_v02.globalIndex())
            continue; // no shared decay products
        }

        if (LambdaTag) {
          ProtonVec2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), massPr);
          PionVec2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), massPi);
          V02Tag = 0;
        }
        if (aLambdaTag) {
          ProtonVec2 = ROOT::Math::PxPyPzMVector(v02.pxneg(), v02.pyneg(), v02.pzneg(), massPr);
          PionVec2 = ROOT::Math::PxPyPzMVector(v02.pxpos(), v02.pypos(), v02.pzpos(), massPi);
          V02Tag = 1;
        }
        LambdaVec2 = ProtonVec2 + PionVec2;
        LambdaVec2.SetM(massLambda);

        ROOT::Math::Boost boost2{LambdaVec2.BoostToCM()};
        ProtonBoostedVec2 = boost2(ProtonVec2);

        costhetastar2 = ProtonBoostedVec2.Pz() / ProtonBoostedVec2.P();

        weight = 1.0;
        weight *= cfgEffCor ? 1.0 / EffMap->GetBinContent(EffMap->GetXaxis()->FindBin(v01.pt()), EffMap->GetYaxis()->FindBin(centrality)) : 1.;
        weight *= cfgAccCor ? 1.0 / AccMap->GetBinContent(AccMap->GetXaxis()->FindBin(v01.pt()), AccMap->GetYaxis()->FindBin(v01.yLambda())) : 1.;
        weight *= cfgEffCor ? 1.0 / EffMap->GetBinContent(EffMap->GetXaxis()->FindBin(v02.pt()), EffMap->GetYaxis()->FindBin(centrality)) : 1.;
        weight *= cfgAccCor ? 1.0 / AccMap->GetBinContent(AccMap->GetXaxis()->FindBin(v02.pt()), AccMap->GetYaxis()->FindBin(v02.yLambda())) : 1.;

        if (V01Tag != V02Tag) {
          weight *= -1.0;
        }

        dphi = TVector2::Phi_0_2pi(v01.phi() - v02.phi());
        if (dphi > constants::math::PI * 1.5) {
          dphi -= constants::math::PI * 2.0;
        }

        histos.fill(HIST("Ana/Signal"), v01.pt(), v02.pt(), v01.yLambda() - v02.yLambda(), dphi, centrality, costhetastar1 * costhetastar2 * weight);
      }
    }
  }

  void processDataSame(EventCandidates::iterator const& collision,
                       TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s,
                       aod::BCsWithTimestamps const&)
  {
    if (cfgCentEst == 1) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == 2) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision) && cfgEvtSel) {
      return;
    }

    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/PVzDist"), collision.posZ(), 1.0);

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgEffCor) {
      EffMap = ccdb->getForTimeStamp<TProfile2D>(cfgEffCorPath.value, bc.timestamp());
    }
    if (cfgAccCor) {
      AccMap = ccdb->getForTimeStamp<TProfile2D>(cfgAccCorPath.value, bc.timestamp());
    }

    FillHistograms(collision, collision, V0s, V0s);
  }
  PROCESS_SWITCH(lambdaTwoPartPolarization, processDataSame, "Process event for same data", true);

  SliceCache cache;
  Preslice<aod::V0Datas> tracksPerCollisionV0 = aod::v0data::collisionId;

  using BinningTypeT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningTypeT0C colBinningT0C{{vertexAxis, centAxis}, true};

  void processDataMixedT0C(EventCandidates const& collisions,
                           TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    for (auto& [c1, c2] : selfCombinations(colBinningT0C, cfgNoMixedEvents, -1, collisions, collisions)) {

      if (c1.index() == c2.index())
        continue;

      centrality = c1.centFT0C();
      if (cfgAccCor) {
        auto bc = c1.bc_as<aod::BCsWithTimestamps>();
        AccMap = ccdb->getForTimeStamp<TProfile2D>(cfgAccCorPath.value, bc.timestamp());
      }
      if (!eventSelected(c1))
        continue;
      if (!eventSelected(c2))
        continue;

      auto tracks1 = V0s.sliceBy(tracksPerCollisionV0, c1.globalIndex());
      auto tracks2 = V0s.sliceBy(tracksPerCollisionV0, c2.globalIndex());

      FillHistograms(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(lambdaTwoPartPolarization, processDataMixedT0C, "Process event for mixed data in PbPb", false);

  using BinningTypeT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeT0M colBinningT0M{{vertexAxis, centAxis}, true};

  void processDataMixedT0M(EventCandidates const& collisions,
                           TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s, aod::BCsWithTimestamps const&)
  {
    for (auto& [c1, c2] : selfCombinations(colBinningT0M, cfgNoMixedEvents, -1, collisions, collisions)) {

      if (c1.index() == c2.index())
        continue;

      centrality = c1.centFT0M();
      if (cfgAccCor) {
        auto bc = c1.bc_as<aod::BCsWithTimestamps>();
        AccMap = ccdb->getForTimeStamp<TProfile2D>(cfgAccCorPath.value, bc.timestamp());
      }
      if (!eventSelected(c1))
        continue;
      if (!eventSelected(c2))
        continue;

      auto tracks1 = V0s.sliceBy(tracksPerCollisionV0, c1.globalIndex());
      auto tracks2 = V0s.sliceBy(tracksPerCollisionV0, c2.globalIndex());

      FillHistograms(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(lambdaTwoPartPolarization, processDataMixedT0M, "Process event for mixed data in pp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdaTwoPartPolarization>(cfgc, TaskName{"lf-lambdaTwoPartPolarization"})};
}
