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

#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambdalambda {
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
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", false, "flag for additional pileup selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
  Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};

  Configurable<float> cfgHypMassWindow{"cfgHypMassWindow", 0.02, "single lambda mass selection"};
  Configurable<float> cfgV0V0RapMax{"cfgV0V0RapMax", 0.5, "rapidity selection for V0V0"};

  Configurable<bool> cfgV0V0Sel{"cfgV0V0Sel", false, "application of V0V0 selections"};
  Configurable<float> cfgV0V0Radius{"cfgV0V0Radius", 1.0, "maximum radius of v0v0"};
  Configurable<float> cfgV0V0CPA{"cfgV0V0CPA", 0.6, "minimum CPA of v0v0"};
  Configurable<float> cfgV0V0Distance{"cfgV0V0Distance", 1, "minimum distance of v0v0"};
  Configurable<float> cfgV0V0DCA{"cfgV0V0DCA", 1.0, "maximum DCA of v0v0"};

  Configurable<bool> cfgEffCor{"cfgEffCor", false, "flag to apply efficiency correction"};
  Configurable<std::string> cfgEffCorPath{"cfgEffCorPath", "", "path for pseudo efficiency correction"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 10, "Number of mixed events per event"};

  ConfigurableAxis massAxis{"massAxis", {110, 2.22, 2.33}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 10, 20, 50, 100}, "Centrality interval"};
  ConfigurableAxis vertexAxis{"vertexAxis", {10, -10, 10}, "vertex axis for mixing"};

  ConfigurableAxis RadiusAxis{"RadiusAxis", {100, 0, 5}, "radius of v0v0"};
  ConfigurableAxis CPAAxis{"CPAAxis", {102, -1.02, 1.02}, "CPA of v0v0"};
  ConfigurableAxis DistanceAxis{"DistanceAxis", {100, 0, 10}, "distance of v0v0"};
  ConfigurableAxis DCAAxis{"DCAAxis", {100, 0, 5}, "DCA of v0v0"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  float centrality;
  TProfile2D* EffMap = nullptr;

  void init(o2::framework::InitContext&)
  {
    AxisSpec centQaAxis = {80, 0.0, 80.0};
    AxisSpec PVzQaAxis = {300, -15.0, 15.0};
    AxisSpec combAxis = {3, -0.5, 2.5};

    histos.add("Radius_V0V0_full", "", {HistType::kTHnSparseF, {massAxis, ptAxis, RadiusAxis, combAxis}});
    histos.add("CPA_V0V0_full", "", {HistType::kTHnSparseF, {massAxis, ptAxis, CPAAxis, combAxis}});
    histos.add("Distance_V0V0_full", "", {HistType::kTHnSparseF, {massAxis, ptAxis, DistanceAxis, combAxis}});
    histos.add("DCA_V0V0_full", "", {HistType::kTHnSparseF, {massAxis, ptAxis, DCAAxis, combAxis}});

    histos.add("Radius_V0V0_sel", "", {HistType::kTHnSparseF, {massAxis, ptAxis, RadiusAxis, combAxis}});
    histos.add("CPA_V0V0_sel", "", {HistType::kTHnSparseF, {massAxis, ptAxis, CPAAxis, combAxis}});
    histos.add("Distance_V0V0_sel", "", {HistType::kTHnSparseF, {massAxis, ptAxis, DistanceAxis, combAxis}});
    histos.add("DCA_V0V0_sel", "", {HistType::kTHnSparseF, {massAxis, ptAxis, DCAAxis, combAxis}});

    histos.add("h_InvMass_same", "", {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, combAxis}});
    histos.add("h_InvMass_mixed", "", {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, combAxis}});
    if (cfgQAv0) {
      histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});
      histos.add("QA/PVzDist", "", {HistType::kTH1F, {PVzQaAxis}});
    }

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
  ROOT::Math::PxPyPzMVector RecoV01, RecoV02, RecoV0V0;

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
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
  bool SelectionV0(TCollision const& collision, V0 const& candidate)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (std::abs(candidate.dcapostopv()) < cfgDCAPosToPVMin)
      return false;
    if (std::abs(candidate.dcanegtopv()) < cfgDCANegToPVMin)
      return false;
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

  template <typename V01, typename V02>
  bool isSelectedV0V0(V01 const& v01, V02 const& v02)
  {
    if (getDCAofV0V0(v01, v02) > cfgV0V0DCA)
      return false;
    if (getCPA(v01, v02) < cfgV0V0CPA)
      return false;
    if (getDistance(v01, v02) < cfgV0V0Distance)
      return false;
    if (getRadius(v01, v02) > cfgV0V0Radius)
      return false;

    return true;
  }

  template <typename V01, typename V02>
  float getDCAofV0V0(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos, v01mom, v02mom;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    v01mom.SetXYZ(v01.px(), v01.py(), v01.pz());
    v02mom.SetXYZ(v02.px(), v02.py(), v02.pz());

    ROOT::Math::XYZVector posdiff = v02pos - v01pos;
    ROOT::Math::XYZVector cross = v01mom.Cross(v02mom);
    if (std::sqrt(cross.Mag2()) < 1e-6)
      return 999.;
    return std::abs(posdiff.Dot(cross)) / std::sqrt(cross.Mag2());
  }

  template <typename V01, typename V02>
  float getCPA(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01mom, v02mom;
    v01mom.SetXYZ(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    v02mom.SetXYZ(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());
    return v01mom.Dot(v02mom);
  }

  template <typename V01, typename V02>
  float getDistance(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    ROOT::Math::XYZVector posdiff = v02pos - v01pos;
    return std::sqrt(posdiff.Mag2());
  }

  template <typename V01, typename V02>
  float getRadius(V01 const& v01, V02 const& v02)
  {
    ROOT::Math::XYZVector v01pos, v02pos, v01mom, v02mom;
    v01pos.SetXYZ(v01.x(), v01.y(), v01.z());
    v02pos.SetXYZ(v02.x(), v02.y(), v02.z());
    v01mom.SetXYZ(v01.px() / v01.p(), v01.py() / v01.p(), v01.pz() / v01.p());
    v02mom.SetXYZ(v02.px() / v02.p(), v02.py() / v02.p(), v02.pz() / v02.p());
    ROOT::Math::XYZVector posdiff = v02pos - v01pos;

    float d = 1. - TMath::Power(v01mom.Dot(v02mom), 2);
    if (d < 1e-5)
      return 999;
    float t = posdiff.Dot(v01mom - v01mom.Dot(v02mom) * v02mom) / d;
    float s = -posdiff.Dot(v02mom - v01mom.Dot(v02mom) * v01mom) / d;
    ROOT::Math::XYZVector dca = v01pos + v02pos + t * v01mom + s * v02mom;
    dca /= 2.;
    return std::sqrt(dca.Mag2());
  }

  template <typename C1, typename C2, typename V01, typename V02>
  void FillHistograms(C1 const& c1, C2 const& c2, V01 const& V01s, V02 const& V02s)
  {
    for (auto& v01 : V01s) {
      auto postrack_v01 = v01.template posTrack_as<TrackCandidates>();
      auto negtrack_v01 = v01.template negTrack_as<TrackCandidates>();

      int LambdaTag = 0;
      int aLambdaTag = 0;
      int V01Tag = -2;

      if (isSelectedV0Daughter(postrack_v01, 0) && isSelectedV0Daughter(negtrack_v01, 1)) {
        LambdaTag = 1;
        V01Tag = 0;
      }
      if (isSelectedV0Daughter(negtrack_v01, 0) && isSelectedV0Daughter(postrack_v01, 1)) {
        aLambdaTag = 1;
        V01Tag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!SelectionV0(c1, v01))
        continue;

      if (LambdaTag) {
        if (std::abs(massLambda - v01.mLambda()) > cfgHypMassWindow)
          continue;
        RecoV01 = ROOT::Math::PxPyPzMVector(v01.px(), v01.py(), v01.pz(), v01.mLambda());
      } else if (aLambdaTag) {
        if (std::abs(massLambda - v01.mAntiLambda()) > cfgHypMassWindow)
          continue;
        RecoV01 = ROOT::Math::PxPyPzMVector(v01.px(), v01.py(), v01.pz(), v01.mAntiLambda());
      }

      for (auto& v02 : V02s) {
        if (v01.v0Id() <= v02.v0Id() && doprocessDataSame)
          continue;
        auto postrack_v02 = v02.template posTrack_as<TrackCandidates>();
        auto negtrack_v02 = v02.template negTrack_as<TrackCandidates>();

        LambdaTag = 0;
        aLambdaTag = 0;
        int V02Tag = -2;

        if (isSelectedV0Daughter(postrack_v02, 0) && isSelectedV0Daughter(negtrack_v02, 1)) {
          LambdaTag = 1;
          V02Tag = 0;
        }
        if (isSelectedV0Daughter(negtrack_v02, 0) && isSelectedV0Daughter(postrack_v02, 1)) {
          aLambdaTag = 1;
          V02Tag = 1;
        }

        if (LambdaTag == aLambdaTag)
          continue;

        if (!SelectionV0(c2, v02))
          continue;

        if (doprocessDataSame) {
          if (postrack_v01.globalIndex() == postrack_v02.globalIndex() || postrack_v01.globalIndex() == negtrack_v02.globalIndex() || negtrack_v01.globalIndex() == postrack_v02.globalIndex() || negtrack_v01.globalIndex() == negtrack_v02.globalIndex())
            continue; // no shared decay products
        }

        if (LambdaTag) {
          if (std::abs(massLambda - v02.mLambda()) > cfgHypMassWindow)
            continue;
          RecoV02 = ROOT::Math::PxPyPzMVector(v02.px(), v02.py(), v02.pz(), v02.mLambda());
        } else if (aLambdaTag) {
          if (std::abs(massLambda - v02.mAntiLambda()) > cfgHypMassWindow)
            continue;
          RecoV02 = ROOT::Math::PxPyPzMVector(v02.px(), v02.py(), v02.pz(), v02.mAntiLambda());
        }

        RecoV0V0 = RecoV01 + RecoV02;
        if (std::abs(RecoV0V0.Rapidity()) > cfgV0V0RapMax)
          continue;

        histos.fill(HIST("Radius_V0V0_full"), RecoV0V0.M(), RecoV0V0.Pt(), getRadius(v01, v02), V01Tag + V02Tag);
        histos.fill(HIST("CPA_V0V0_full"), RecoV0V0.M(), RecoV0V0.Pt(), getCPA(v01, v02), V01Tag + V02Tag);
        histos.fill(HIST("Distance_V0V0_full"), RecoV0V0.M(), RecoV0V0.Pt(), getDistance(v01, v02), V01Tag + V02Tag);
        histos.fill(HIST("DCA_V0V0_full"), RecoV0V0.M(), RecoV0V0.Pt(), getDCAofV0V0(v01, v02), V01Tag + V02Tag);

        if (isSelectedV0V0(v01, v02)) {
          histos.fill(HIST("Radius_V0V0_sel"), RecoV0V0.M(), RecoV0V0.Pt(), getRadius(v01, v02), V01Tag + V02Tag);
          histos.fill(HIST("CPA_V0V0_sel"), RecoV0V0.M(), RecoV0V0.Pt(), getCPA(v01, v02), V01Tag + V02Tag);
          histos.fill(HIST("Distance_V0V0_sel"), RecoV0V0.M(), RecoV0V0.Pt(), getDistance(v01, v02), V01Tag + V02Tag);
          histos.fill(HIST("DCA_V0V0_sel"), RecoV0V0.M(), RecoV0V0.Pt(), getDCAofV0V0(v01, v02), V01Tag + V02Tag);
        }

        if (cfgV0V0Sel && !isSelectedV0V0(v01, v02))
          continue;

        if (doprocessDataSame) {
          histos.fill(HIST("h_InvMass_same"), RecoV0V0.M(), RecoV0V0.Pt(), centrality, V01Tag + V02Tag);
        }
        if (doprocessDataMixed) {
          histos.fill(HIST("h_InvMass_mixed"), RecoV0V0.M(), RecoV0V0.Pt(), centrality, V01Tag + V02Tag);
        }
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
    if (!eventSelected(collision)) {
      return;
    }

    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/PVzDist"), collision.posZ(), 1.0);

    if (cfgEffCor) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      EffMap = ccdb->getForTimeStamp<TProfile2D>(cfgEffCorPath.value, bc.timestamp());
    }
    FillHistograms(collision, collision, V0s, V0s);
  }
  PROCESS_SWITCH(lambdalambda, processDataSame, "Process Event for same data", true);

  SliceCache cache;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  void processDataMixed(EventCandidates const& collisions,
                        TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s)
  {
    auto tracksTuple = std::make_tuple(V0s);
    BinningTypeVertexContributor binningOnPositions{{vertexAxis, centAxis}, true};
    SameKindPair<EventCandidates, V0TrackCandidate, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      if (cfgCentEst == 1) {
        centrality = c1.centFT0C();
      } else if (cfgCentEst == 2) {
        centrality = c1.centFT0M();
      }
      if (!eventSelected(c1))
        continue;
      if (!eventSelected(c2))
        continue;
      if (c1.bcId() == c2.bcId())
        continue;

      FillHistograms(c1, c2, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(lambdalambda, processDataMixed, "Process Event for mixed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdalambda>(cfgc, TaskName{"lf-lambdalambda"})};
}
