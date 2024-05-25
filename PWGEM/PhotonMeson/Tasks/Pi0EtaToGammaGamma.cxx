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
//
// ========================
//
// This code loops over photons and makes pairs for neutral mesons analyses.
//    Please write to: daiki.sekihata@cern.ch

#include <cstring>
#include <iterator>

#include "TString.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzRotation.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "CCDB/BasicCCDBManager.h"
#include "Tools/ML/MlResponse.h"
#include "Tools/ML/model.h"

#include "Common/Core/RecoDecay.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/PHOSPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Utils/EMTrack.h"
#include "PWGEM/PhotonMeson/Utils/EventMixingHandler.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::photonpair;
using namespace o2::aod::pwgem::photon;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec, aod::EMEventsNgPCM, aod::EMEventsNgPHOS, aod::EMEventsNgEMC, aod::EMEventsNee>;
using MyCollision = MyCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyDalitzEEs = soa::Join<aod::DalitzEEs, aod::DalitzEEEMEventIds>;
using MyDalitzEE = MyDalitzEEs::iterator;

using MyDalitzMuMus = soa::Join<aod::DalitzMuMus, aod::DalitzMuMuEMEventIds>;
using MyDalitzMuMu = MyDalitzMuMus::iterator;

using MyPrimaryElectrons = soa::Join<aod::EMPrimaryElectrons, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBit>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

using MyPrimaryMuons = soa::Join<aod::EMPrimaryMuons, aod::EMPrimaryMuonEMEventIds, aod::EMPrimaryMuonsPrefilterBit>;
using MyPrimaryMuon = MyPrimaryMuons::iterator;

using MyEMCClusters = soa::Join<aod::SkimEMCClusters, aod::EMCEMEventIds>;
using MyEMCCluster = MyEMCClusters::iterator;

using MyPHOSClusters = soa::Join<aod::PHOSClusters, aod::PHOSEMEventIds>;
using MyPHOSCluster = MyEMCClusters::iterator;

namespace o2::aod::pwgem::photonmeson::nmgghistogram
{
void addhistograms(HistogramRegistry* fRegistry, bool doFlow, const char* pairname = "#gamma#gamma")
{
  // event info
  auto hCollisionCounter = fRegistry->add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
  hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
  hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
  hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
  hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
  hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
  hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
  hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
  hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
  hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

  fRegistry->add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
  fRegistry->add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
  fRegistry->add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
  fRegistry->add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
  fRegistry->add("Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
  fRegistry->add("Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);

  if (doFlow) {
    fRegistry->add("Event/before/hQ2xFT0M_CentFT0C", "hQ2xFT0M_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0M_CentFT0C", "hQ2yFT0M_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0M}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0A_CentFT0C", "hQ2xFT0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0A_CentFT0C", "hQ2yFT0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFT0C_CentFT0C", "hQ2xFT0C_CentFT0C;centrality FT0C (%);Q_{2,x}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFT0C_CentFT0C", "hQ2yFT0C_CentFT0C;centrality FT0C (%);Q_{2,y}^{FT0C}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xFV0A_CentFT0C", "hQ2xFV0A_CentFT0C;centrality FT0C (%);Q_{2,x}^{FV0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yFV0A_CentFT0C", "hQ2yFV0A_CentFT0C;centrality FT0C (%);Q_{2,y}^{FV0A}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBPos_CentFT0C", "hQ2xBPos_CentFT0C;centrality FT0C (%);Q_{2,x}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBPos_CentFT0C", "hQ2yBPos_CentFT0C;centrality FT0C (%);Q_{2,y}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2xBNeg_CentFT0C", "hQ2xBNeg_CentFT0C;centrality FT0C (%);Q_{2,x}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2yBNeg_CentFT0C", "hQ2yBNeg_CentFT0C;centrality FT0C (%);Q_{2,y}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hEP2FT0M_CentFT0C", "2nd harmonics event plane FT0M;centrality FT0C (%);#Psi_{2}^{FT0M} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0A_CentFT0C", "2nd harmonics event plane FT0A;centrality FT0C (%);#Psi_{2}^{FT0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FT0C_CentFT0C", "2nd harmonics event plane FT0C;centrality FT0C (%);#Psi_{2}^{FT0C} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2FV0A_CentFT0C", "2nd harmonics event plane FV0A;centrality FT0C (%);#Psi_{2}^{FV0A} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BPos_CentFT0C", "2nd harmonics event plane BPos;centrality FT0C (%);#Psi_{2}^{BPos} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);
    fRegistry->add("Event/before/hEP2BNeg_CentFT0C", "2nd harmonics event plane BNeg;centrality FT0C (%);#Psi_{2}^{BNeg} (rad.)", kTH2F, {{110, 0, 110}, {180, -M_PI_2, +M_PI_2}}, false);

    fRegistry->add("Event/before/hQ2FT0MQ2BPos_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0MQ2BNeg_CentFT0C", "Q_{2}^{FT0M} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0M} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2BPosQ2BNeg_CentFT0C", "Q_{2}^{BPos} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{BPos} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false); // this is common for FT0M, FT0A, FT0C, FV0A resolution.

    fRegistry->add("Event/before/hQ2FT0CQ2BPos_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0CQ2BNeg_CentFT0C", "Q_{2}^{FT0C} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0C} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hQ2FT0AQ2BPos_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FT0AQ2BNeg_CentFT0C", "Q_{2}^{FT0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FT0A} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);

    fRegistry->add("Event/before/hQ2FV0AQ2BPos_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BPos};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BPos}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
    fRegistry->add("Event/before/hQ2FV0AQ2BNeg_CentFT0C", "Q_{2}^{FV0A} #upoint Q_{2}^{BNeg};centrality FT0C (%);Q_{2}^{FV0A} #upoint Q_{2}^{BNeg}", kTH2F, {{110, 0, 110}, {200, -10, +10}}, false);
  }

  fRegistry->addClone("Event/before/", "Event/after/");

  // pair info

  std::vector<double> ptbins;
  for (int i = 0; i < 50; i++) {
    ptbins.emplace_back(0.1 * (i - 0) + 0.0); // from 0 to 5 GeV/c, every 0.1 GeV/c
  }
  for (int i = 50; i < 60; i++) {
    ptbins.emplace_back(0.5 * (i - 50) + 5.0); // from 5 to 10 GeV/c, evety 0.5 GeV/c
  }
  for (int i = 60; i < 71; i++) {
    ptbins.emplace_back(1.0 * (i - 60) + 10.0); // from 10 to 20 GeV/c, evety 1 GeV/c
  }
  const AxisSpec axis_pt{ptbins, Form("p_{T,%s} (GeV/c)", pairname)};

  const AxisSpec axis_mass{400, 0, 0.8, Form("m_{%s} (GeV/c^{2})", pairname)};
  const AxisSpec sp2_mass{50, -5, 5, Form("u_{2}^{%s} #upoint Q_{2}", pairname)};

  if (doFlow) {
    std::string_view sp_names[4] = {"FT0M", "FT0A", "FT0C", "FV0A"};
    for (int i = 0; i < 4; i++) {
      fRegistry->add(Form("Pair/same/hs_same_SPQ2%s", sp_names[i].data()), "2photon", kTHnSparseF, {axis_mass, axis_pt, sp2_mass}, true);
    }
    fRegistry->add("Pair/mix/hMggPt", "2photon", kTH2F, {axis_mass, axis_pt}, true);
  } else {
    fRegistry->add("Pair/same/hMggPt", "2photon", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry->addClone("Pair/same/", "Pair/mix/");
  }
}

template <const int ev_id, typename TCollision>
void fillEventInfo(HistogramRegistry* fRegistry, TCollision const& collision, const bool doFlow, const float weight = 1.f)
{
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 1.0);
  if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 2.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 3.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 4.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 5.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 6.0);
  }
  if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 7.0);
  }
  if (collision.sel8()) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 8.0);
  }
  if (abs(collision.posZ()) < 10.0) {
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCollisionCounter"), 9.0);
  }
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hZvtx"), collision.posZ());

  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0A"), collision.centFT0A());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0C"), collision.centFT0C());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0M"), collision.centFT0M());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
  fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());

  if (doFlow) {
    std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
    std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
    std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
    std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};
    std::array<float, 2> q2bpos = {collision.q2xbpos(), collision.q2ybpos()};
    std::array<float, 2> q2bneg = {collision.q2xbneg(), collision.q2ybneg()};

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0M_CentFT0C"), collision.centFT0C(), collision.q2xft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0M_CentFT0C"), collision.centFT0C(), collision.q2yft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0A_CentFT0C"), collision.centFT0C(), collision.q2xft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0A_CentFT0C"), collision.centFT0C(), collision.q2yft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFT0C_CentFT0C"), collision.centFT0C(), collision.q2xft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFT0C_CentFT0C"), collision.centFT0C(), collision.q2yft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xFV0A_CentFT0C"), collision.centFT0C(), collision.q2xfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yFV0A_CentFT0C"), collision.centFT0C(), collision.q2yfv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBPos_CentFT0C"), collision.centFT0C(), collision.q2xbpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBPos_CentFT0C"), collision.centFT0C(), collision.q2ybpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2xBNeg_CentFT0C"), collision.centFT0C(), collision.q2xbneg());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2yBNeg_CentFT0C"), collision.centFT0C(), collision.q2ybneg());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0M_CentFT0C"), collision.centFT0C(), collision.ep2ft0m());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0A_CentFT0C"), collision.centFT0C(), collision.ep2ft0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FT0C_CentFT0C"), collision.centFT0C(), collision.ep2ft0c());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2FV0A_CentFT0C"), collision.centFT0C(), collision.ep2fv0a());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BPos_CentFT0C"), collision.centFT0C(), collision.ep2bpos());
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hEP2BNeg_CentFT0C"), collision.centFT0C(), collision.ep2bneg());

    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0MQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0MQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0CQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FT0CQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FV0AQ2BPos_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bpos));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2FV0AQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bneg));
    fRegistry->fill(HIST("Event/") + HIST(event_types[ev_id]) + HIST("hQ2BPosQ2BNeg_CentFT0C"), collision.centFT0C(), RecoDecay::dotProd(q2bpos, q2bneg));
  }
}

template <int ev_id, PairType pairtype, typename TCollision, typename TDiphoton>
bool fillPairInfo(HistogramRegistry* fRegistry, TCollision const& collision, TDiphoton const& diphoton, const bool doFlow, const float weight = 1.f)
{

  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};
  if (doFlow) {
    if constexpr (ev_id == 0) { // same
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};
      std::array<float, 2> u2_gg = {static_cast<float>(std::cos(2 * diphoton.Phi())), static_cast<float>(std::sin(2 * diphoton.Phi()))};

      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0M"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0m));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0A"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0a));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0C"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0c));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FV0A"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2fv0a));
    } else { // mix
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMggPt"), diphoton.M(), diphoton.Pt());
    }
  } else {
    fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMggPt"), diphoton.M(), diphoton.Pt());
  }
  return true;
}

} // namespace o2::aod::pwgem::photonmeson::nmgghistogram

template <PairType pairtype, typename... Types>
struct Pi0EtaToGammaGamma_v2 {
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};
  Configurable<bool> cfgDoFlow{"cfgDoFlow", false, "flag to analyze vn"};
  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f, 999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfEPBins{"ConfEPBins", {VARIABLE_WIDTH, 0.0f, M_PI / 4, M_PI / 2, M_PI}, "Mixing bins - event plane angle"};

  EMEventCut fEMEventCut;
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
  Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
  Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
  Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
  Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
  Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    Configurable<bool> cfg_require_v0_on_wwire_ib{"cfg_require_v0_on_wwire_ib", false, "flag to select V0s on W wires ITSib"};
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.9, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.99, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 3.0, "max distance btween 2 legs"};
    Configurable<bool> cfg_require_v0_with_correct_xz{"cfg_require_v0_with_correct_xz", true, "flag to select V0s with correct xz"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};

    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 10, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 20, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.5, "max mass"};
    Configurable<float> cfg_min_pair_dca3d{"cfg_min_pair_dca3d", 0.0, "min pair dca3d in sigma"};
    Configurable<float> cfg_max_pair_dca3d{"cfg_max_pair_dca3d", 1e+10, "max pair dca3d in sigma"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_apply_pf{"cfg_apply_pf", false, "flag to apply phiv prefilter"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", true, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", false, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.9, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq), "pid scheme [kTOFreq : 0, kTPChadrej : 1, kTPChadrejORTOFreq : 2, kTPConly : 3]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaMu{"cfg_min_TPCNsigmaMu", -0.0, "min. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaMu{"cfg_max_TPCNsigmaMu", +0.0, "max. TPC n sigma for muon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3.0, "min. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3.0, "max. TPC n sigma for kaon exclusion"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3.0, "min. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3.0, "max. TPC n sigma for proton exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};

    // CCDB configuration for PID ML
    Configurable<std::string> BDTLocalPathGamma{"BDTLocalPathGamma", "pid_ml_xgboost.onnx", "Path to the local .onnx file"};

    Configurable<std::string> BDTPathCCDB{"BDTPathCCDB", "Users/d/dsekihat/pwgem/pidml/", "Path on CCDB"};
    Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
    Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
    Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
  } dileptoncuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<bool> requireCaloReadout{"requireCaloReadout", true, "Require calorimeters readout when analyzing EMCal/PHOS"};
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
    Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
    Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
    Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
  } emccuts;

  PHOSPhotonCut fPHOSCut;
  struct : ConfigurableGroup {
    std::string prefix = "phoscut_group";
    Configurable<float> cfg_min_Ecluster{"cfg_min_Ecluster", 0.3, "Minimum cluster energy for PHOS in GeV"};
  } phoscuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  static constexpr std::string_view event_types[2] = {"before/", "after/"};
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};

  std::vector<float> zvtx_bin_edges;
  std::vector<float> cent_bin_edges;

  void init(InitContext&)
  {
    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());

    cent_bin_edges = std::vector<float>(ConfCentBins.value.begin(), ConfCentBins.value.end());
    cent_bin_edges.erase(cent_bin_edges.begin());
    emh1 = new o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>(ndepth);
    emh2 = new o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>(ndepth);

    if constexpr (pairtype == PairType::kPCMDalitzEE) {
      o2::aod::pwgem::photonmeson::nmgghistogram::addhistograms(&fRegistry, cfgDoFlow, "ee#gamma");
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      o2::aod::pwgem::photonmeson::nmgghistogram::addhistograms(&fRegistry, cfgDoFlow, "#mu#mu#gamma");
    } else {
      o2::aod::pwgem::photonmeson::nmgghistogram::addhistograms(&fRegistry, cfgDoFlow, "#gamma#gamma");
    }
    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();
    DefineEMCCut();
    DefinePHOSCut();

    if constexpr (pairtype == kEMCEMC) {
      fRegistry.addClone("Pair/same/", "Pair/rotation/");
    }
  }

  ~Pi0EtaToGammaGamma_v2()
  {
    delete emh1;
    emh1 = 0x0;
    delete emh2;
    emh2 = 0x0;

    used_photonIds.clear();
    used_photonIds.shrink_to_fit();
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-cfgZvtxMax, +cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(cfgRequireGoodZvtxFT0vsPV);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, 1e10f);
    fV0PhotonCut.SetV0EtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetTrackPtRange(pcmcuts.cfg_min_pt_v0 * 0.4, 1e+10f);
    fV0PhotonCut.SetTrackEtaRange(-pcmcuts.cfg_max_eta_v0, +pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    if (pcmcuts.cfg_reject_v0_on_itsib) {
      fV0PhotonCut.SetNClustersITS(2, 4);
    } else {
      fV0PhotonCut.SetNClustersITS(0, 7);
    }
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetIsWithinBeamPipe(pcmcuts.cfg_require_v0_with_correct_xz);

    if (pcmcuts.cfg_require_v0_with_itstpc) {
      fV0PhotonCut.SetRequireITSTPC(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 40);
    }
    if (pcmcuts.cfg_require_v0_with_itsonly) {
      fV0PhotonCut.SetRequireITSonly(true);
      fV0PhotonCut.SetMaxPCA(1.0);
      fV0PhotonCut.SetRxyRange(4, 24);
    }
    if (pcmcuts.cfg_require_v0_with_tpconly) {
      fV0PhotonCut.SetRequireTPConly(true);
      fV0PhotonCut.SetMaxPCA(3.0);
      fV0PhotonCut.SetRxyRange(36, 90);
    }
    if (pcmcuts.cfg_require_v0_on_wwire_ib) {
      fV0PhotonCut.SetMaxPCA(0.3);
      fV0PhotonCut.SetOnWwireIB(true);
      fV0PhotonCut.SetOnWwireOB(false);
      fV0PhotonCut.SetRxyRange(7, 14);
    }
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.SetPairDCARange(dileptoncuts.cfg_min_pair_dca3d, dileptoncuts.cfg_max_pair_dca3d); // in sigma
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.ApplyPrefilter(dileptoncuts.cfg_apply_pf);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    // for track
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMeanClusterSizeITSob(0, 16);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaMuRange(dileptoncuts.cfg_min_TPCNsigmaMu, dileptoncuts.cfg_max_TPCNsigmaMu);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTPCNsigmaKaRange(dileptoncuts.cfg_min_TPCNsigmaKa, dileptoncuts.cfg_max_TPCNsigmaKa);
    fDileptonCut.SetTPCNsigmaPrRange(dileptoncuts.cfg_min_TPCNsigmaPr, dileptoncuts.cfg_max_TPCNsigmaPr);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);

    if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) { // please call this at the end of DefineDileptonCut
      o2::ml::OnnxModel* eid_bdt = new o2::ml::OnnxModel();
      if (dileptoncuts.loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        std::map<std::string, std::string> metadata;
        bool retrieveSuccessGamma = ccdbApi.retrieveBlob(dileptoncuts.BDTPathCCDB.value, ".", metadata, dileptoncuts.timestampCCDB.value, false, dileptoncuts.BDTLocalPathGamma.value);
        if (retrieveSuccessGamma) {
          eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
        } else {
          LOG(fatal) << "Error encountered while fetching/loading the Gamma model from CCDB! Maybe the model doesn't exist yet for this runnumber/timestamp?";
        }
      } else {
        eid_bdt->initModel(dileptoncuts.BDTLocalPathGamma.value, dileptoncuts.enableOptimizations.value);
      }

      fDileptonCut.SetPIDModel(eid_bdt);
    } // end of PID ML
  }

  void DefineEMCCut()
  {
    const float a = emccuts.EMC_TM_Eta->at(0);
    const float b = emccuts.EMC_TM_Eta->at(1);
    const float c = emccuts.EMC_TM_Eta->at(2);

    const float d = emccuts.EMC_TM_Phi->at(0);
    const float e = emccuts.EMC_TM_Phi->at(1);
    const float f = emccuts.EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);

    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);

    fEMCCut.SetTrackMatchingEta([&a, &b, &c](float pT) { return a + pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([&d, &e, &f](float pT) { return d + pow(pT + e, f); });

    fEMCCut.SetMinEoverP(emccuts.EMC_Eoverp);
    fEMCCut.SetUseExoticCut(emccuts.EMC_UseExoticCut);
  }

  void DefinePHOSCut()
  {
    fPHOSCut.SetEnergyRange(phoscuts.cfg_min_Ecluster, 1e+10);
  }

  template <typename TCollision, typename TTracks, typename TCut>
  std::vector<EMTrack> createULSPairs(TCollision const& collision, TTracks const& posTracks, TTracks const& negTracks, TCut const& cut, const float mlepton)
  {
    std::vector<EMTrack> pairs;
    pairs.reserve(posTracks.size() * negTracks.size());

    for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks, negTracks))) { // ULS

      if (pos.trackId() == ele.trackId()) { // this is protection against pairing identical 2 tracks.
        continue;
      }

      if (dileptoncuts.cfg_pid_scheme == static_cast<int>(DalitzEECut::PIDSchemes::kPIDML)) {
        if (!cut.template IsSelectedTrack<true>(pos, collision) || !cut.template IsSelectedTrack<true>(ele, collision)) {
          continue;
        }
      } else { // cut-based
        if (!cut.template IsSelectedTrack<false>(pos, collision) || !cut.template IsSelectedTrack<false>(ele, collision)) {
          continue;
        }
      }

      if (!cut.IsSelectedPair(pos, ele, collision.bz())) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), mlepton);
      ROOT::Math::PtEtaPhiMVector v2(ele.pt(), ele.eta(), ele.phi(), mlepton);
      ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
      float dca_pos_3d = pos.dca3DinSigma();
      float dca_ele_3d = ele.dca3DinSigma();
      float dca_ee_3d = std::sqrt((dca_pos_3d * dca_pos_3d + dca_ele_3d * dca_ele_3d) / 2.);
      pairs.emplace_back(EMTrack(collision.globalIndex(), ndilepton, v12.Pt(), v12.Eta(), v12.Phi(), v12.M(), 0, dca_ee_3d));
      ndilepton++;
    }
    return pairs;
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, EMCPhotonCut const& cut, aod::SkimEMCMTs const&)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const float rotationAngle = M_PI / 2.0; // rotaion angle 90 degree
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    for (auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!cut.template IsSelected<aod::SkimEMCMTs>(photon)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector photon3(photon.pt(), photon.eta(), photon.phi(), 0.);
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));

      if (openingAngle1 > emccuts.minOpenAngle && abs(mother1.Rapidity()) < maxY) {
        fRegistry.fill(HIST("Pair/rotation/hMggPt"), mother1.M(), mother1.Pt());
      }
      if (openingAngle2 > emccuts.minOpenAngle && abs(mother2.Rapidity()) < maxY) {
        fRegistry.fill(HIST("Pair/rotation/hMggPt"), mother2.M(), mother2.Pt());
      }
    }
  }

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::emeventId;
  Preslice<MyEMCClusters> perCollision_emc = aod::emccluster::emeventId;
  Preslice<MyPHOSClusters> perCollision_phos = aod::phoscluster::emeventId;

  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectron::emeventId;
  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaEl) < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaEl);

  Preslice<MyPrimaryMuons> perCollision_muon = aod::emprimarymuon::emeventId;
  Partition<MyPrimaryMuons> muons_pos = o2::aod::emprimarymuon::sign > int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt&& nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu&& o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);
  Partition<MyPrimaryMuons> muons_neg = o2::aod::emprimarymuon::sign < int8_t(0) && static_cast<float>(dileptoncuts.cfg_min_pt_track) < o2::aod::track::pt && nabs(o2::aod::track::eta) < static_cast<float>(dileptoncuts.cfg_max_eta_track) && static_cast<float>(dileptoncuts.cfg_min_TPCNsigmaMu) < o2::aod::pidtpc::tpcNSigmaMu && o2::aod::pidtpc::tpcNSigmaMu < static_cast<float>(dileptoncuts.cfg_max_TPCNsigmaMu);

  o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>* emh1 = nullptr;
  o2::aod::pwgem::photonmeson::utils::EventMixingHandler<std::tuple<int, int, int>, std::pair<int, int64_t>, EMTrack>* emh2 = nullptr;
  std::map<std::tuple<int, int, int>, std::vector<std::pair<int, int64_t>>> map_mix_bins; // zvtx, centrality, event plane bins -> pair<df index, vector of collisionId>
  std::map<std::pair<int, int64_t>, std::vector<EMTrack>> map_photons_to_collision;       // pair<df index, collisionId> -> vector of track
  std::vector<std::pair<int, int>> used_photonIds;                                        // <ndf, trackId>

  template <typename TCollisions, typename TPhotons1, typename TPhotons2, typename TSubInfos1, typename TSubInfos2, typename TPreslice1, typename TPreslice2, typename TCut1, typename TCut2, typename TTracksMatchedWithEMC, typename TTracksMatchedWithPHOS>
  void runPairing(TCollisions const& collisions,
                  TPhotons1 const& photons1, TPhotons2 const& photons2,
                  TSubInfos1 const& subinfos1, TSubInfos2 const& subinfos2,
                  TPreslice1 const& perCollision1, TPreslice2 const& perCollision2,
                  TCut1 const& cut1, TCut2 const& cut2,
                  TTracksMatchedWithEMC const& tracks_emc, TTracksMatchedWithPHOS const& tracks_phos)
  {
    for (auto& collision : collisions) {
      int ndiphoton = 0;
      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == PairType::kEMCEMC || pairtype == PairType::kPCMEMC) && ((!collision.alias_bit(triggerAliases::kTVXinEMC) && emccuts.requireCaloReadout) || collision.ncollsPerBC() != 1)) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      o2::aod::pwgem::photonmeson::nmgghistogram::fillEventInfo<0>(&fRegistry, collision, cfgDoFlow);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      o2::aod::pwgem::photonmeson::nmgghistogram::fillEventInfo<1>(&fRegistry, collision, cfgDoFlow);
      fRegistry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      fRegistry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      int zbin = lower_bound(zvtx_bin_edges.begin(), zvtx_bin_edges.end(), collision.posZ()) - zvtx_bin_edges.begin() - 1;
      if (zbin < 0) {
        zbin = 0;
      } else if (static_cast<int>(zvtx_bin_edges.size()) - 2 < zbin) {
        zbin = static_cast<int>(zvtx_bin_edges.size()) - 2;
      }

      float centrality = centralities[cfgCentEstimator];
      int centbin = lower_bound(cent_bin_edges.begin(), cent_bin_edges.end(), centrality) - cent_bin_edges.begin() - 1;
      if (centbin < 0) {
        centbin = 0;
      } else if (static_cast<int>(cent_bin_edges.size()) - 2 < centbin) {
        centbin = static_cast<int>(cent_bin_edges.size()) - 2;
      }
      int epbin = 0;

      std::tuple<int, int, int64_t> key_bin = std::make_tuple(zbin, centbin, epbin);
      std::pair<int, int64_t> key_df_collision = std::make_pair(ndf, collision.globalIndex());

      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) { // same kinds pairing
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (abs(v12.Rapidity()) > maxY) {
            continue;
          }
          o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);

          if constexpr (pairtype == PairType::kEMCEMC) {
            RotationBackground<MyEMCClusters>(v12, v1, v2, photons2_per_collision, g1.globalIndex(), g2.globalIndex(), cut1, tracks_emc);
          }

          std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
          std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.globalIndex());

          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id1);
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
            emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.globalIndex(), g2.pt(), g2.eta(), g2.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id2);
          }
          ndiphoton++;
        } // end of pairing loop
      } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto positrons_per_collision = positrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto electrons_per_collision = electrons->sliceByCached(o2::aod::emprimaryelectron::emeventId, collision.globalIndex(), cache);
        auto photons2_per_collision = createULSPairs(collision, positrons_per_collision, electrons_per_collision, cut2, o2::constants::physics::MassElectron);

        for (auto& g1 : photons1_per_collision) {
          for (auto& g2 : photons2_per_collision) {
            if (!cut1.template IsSelected<TSubInfos1>(g1)) {
              continue;
            }
            ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
            ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
            ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
            if (abs(v12.Rapidity()) > maxY) {
              continue;
            }
            o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);

            std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
            std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.trackId());
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
              used_photonIds.emplace_back(pair_tmp_id1);
            }
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
              emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.trackId(), g2.pt(), g2.eta(), g2.phi(), g2.mass(), g2.sign(), g2.dca3DinSigma()));
              used_photonIds.emplace_back(pair_tmp_id2);
            }
            ndiphoton++;
          } // end of g2 loop
        }   // end of g1 loop
      } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto muons_pos_per_collision = muons_pos->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
        auto muons_neg_per_collision = muons_neg->sliceByCached(o2::aod::emprimarymuon::emeventId, collision.globalIndex(), cache);
        auto photons2_per_collision = createULSPairs(collision, muons_pos_per_collision, muons_neg_per_collision, cut2, o2::constants::physics::MassMuon);

        for (auto& g1 : photons1_per_collision) {
          for (auto& g2 : photons2_per_collision) {
            if (!cut1.template IsSelected<TSubInfos1>(g1)) {
              continue;
            }
            ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
            ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), g2.mass());
            ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
            if (abs(v12.Rapidity()) > maxY) {
              continue;
            }
            o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);

            std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
            std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.trackId());
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
              emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
              used_photonIds.emplace_back(pair_tmp_id1);
            }
            if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
              emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.trackId(), g2.pt(), g2.eta(), g2.phi(), g2.mass(), g2.sign(), g2.dca3DinSigma()));
              used_photonIds.emplace_back(pair_tmp_id2);
            }
            ndiphoton++;
          }    // end of g2 loop
        }      // end of g1 loop
      } else { // PCM-EMC, PCM-PHOS. Nightmare. don't run these pairs.
        auto photons1_per_collision = photons1.sliceBy(perCollision1, collision.globalIndex());
        auto photons2_per_collision = photons2.sliceBy(perCollision2, collision.globalIndex());

        for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_per_collision, photons2_per_collision))) {
          if (!cut1.template IsSelected<TSubInfos1>(g1) || !cut2.template IsSelected<TSubInfos2>(g2)) {
            continue;
          }
          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          if (abs(v12.Rapidity()) > maxY) {
            continue;
          }
          o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<0, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
          std::pair<int, int> pair_tmp_id1 = std::make_pair(ndf, g1.globalIndex());
          std::pair<int, int> pair_tmp_id2 = std::make_pair(ndf, g2.globalIndex());

          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id1) == used_photonIds.end()) {
            emh1->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g1.globalIndex(), g1.pt(), g1.eta(), g1.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id1);
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), pair_tmp_id2) == used_photonIds.end()) {
            emh2->AddTrackToEventPool(key_df_collision, EMTrack(collision.globalIndex(), g2.globalIndex(), g2.pt(), g2.eta(), g2.phi(), 0, 0, 0));
            used_photonIds.emplace_back(pair_tmp_id2);
          }
          ndiphoton++;
        } // end of pairing loop
      }   // end of pairing in same event

      // event mixing
      if (!cfgDoMix || !(ndiphoton > 0)) {
        continue;
      }

      // make a vector of selected photons in this collision.
      auto selected_photons1_in_this_event = emh1->GetTracksPerCollision(key_df_collision);
      auto selected_photons2_in_this_event = emh2->GetTracksPerCollision(key_df_collision);

      auto collisionIds1_in_mixing_pool = emh1->GetCollisionIdsFromEventPool(key_bin);
      auto collisionIds2_in_mixing_pool = emh2->GetCollisionIdsFromEventPool(key_bin);

      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) { // same kinds pairing
        selected_photons1_in_this_event.insert(selected_photons1_in_this_event.end(), selected_photons2_in_this_event.begin(), selected_photons2_in_this_event.end());

        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          photons1_from_event_pool.insert(photons1_from_event_pool.end(), photons2_from_event_pool.begin(), photons2_from_event_pool.end());
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        }      // end of loop over mixed event pool
      } else { //[photon1 from event1, photon2 from event2] and [photon1 from event2, photon2 from event1]

        for (auto& mix_dfId_collisionId : collisionIds2_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons2_from_event_pool = emh2->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), ngamma = %d | event pool (%d, %d), nll = %d", ndf, collision.globalIndex(), selected_photons1_in_this_event.size(), mix_dfId, mix_collisionId, photons2_from_event_pool.size());

          for (auto& g1 : selected_photons1_in_this_event) {
            for (auto& g2 : photons2_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v2.SetM(g2.mass());
              }
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        } // end of loop over mixed event pool

        for (auto& mix_dfId_collisionId : collisionIds1_in_mixing_pool) {
          int mix_dfId = mix_dfId_collisionId.first;
          int64_t mix_collisionId = mix_dfId_collisionId.second;

          if (collision.globalIndex() == mix_collisionId && ndf == mix_dfId) { // this never happens. only protection.
            continue;
          }

          auto photons1_from_event_pool = emh1->GetTracksPerCollision(mix_dfId_collisionId);
          // LOGF(info, "Do event mixing: current event (%d, %d), nll = %d | event pool (%d, %d), ngamma = %d", ndf, collision.globalIndex(), selected_photons2_in_this_event.size(), mix_dfId, mix_collisionId, photons1_from_event_pool.size());

          for (auto& g1 : selected_photons2_in_this_event) {
            for (auto& g2 : photons1_from_event_pool) {
              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) { //[photon from event1, dilepton from event2] and [photon from event2, dilepton from event1]
                v1.SetM(g1.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              o2::aod::pwgem::photonmeson::nmgghistogram::fillPairInfo<1, pairtype>(&fRegistry, collision, v12, cfgDoFlow);
            }
          }
        } // end of loop over mixed event pool
      }

      if (ndiphoton > 0) {
        emh1->AddCollisionIdAtLast(key_bin, key_df_collision);
        emh2->AddCollisionIdAtLast(key_bin, key_df_collision);
      }

    } // end of collision loop
  }

  int ndf = 0;
  int ndilepton = 0;
  void processAnalysis(MyCollisions const& collisions, Types const&... args)
  {
    ndilepton = 0;
    // LOGF(info, "ndf = %d", ndf);
    if constexpr (pairtype == PairType::kPCMPCM) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      runPairing(collisions, v0photons, v0photons, v0legs, v0legs, perCollision_pcm, perCollision_pcm, fV0PhotonCut, fV0PhotonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimaryelectrons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimaryelectrons, v0legs, emprimaryelectrons, perCollision_pcm, perCollision_electron, fV0PhotonCut, fDileptonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      auto v0photons = std::get<0>(std::tie(args...));
      auto v0legs = std::get<1>(std::tie(args...));
      auto emprimarymuons = std::get<2>(std::tie(args...));
      // LOGF(info, "electrons.size() = %d, positrons.size() = %d", electrons.size(), positrons.size());
      runPairing(collisions, v0photons, emprimarymuons, v0legs, emprimarymuons, perCollision_pcm, perCollision_muon, fV0PhotonCut, fDileptonCut, nullptr, nullptr);
    } else if constexpr (pairtype == PairType::kEMCEMC) {
      auto emcclusters = std::get<0>(std::tie(args...));
      auto emcmatchedtracks = std::get<1>(std::tie(args...));
      runPairing(collisions, emcclusters, emcclusters, nullptr, nullptr, perCollision_emc, perCollision_emc, fEMCCut, fEMCCut, emcmatchedtracks, nullptr);
    } else if constexpr (pairtype == PairType::kPHOSPHOS) {
      auto phosclusters = std::get<0>(std::tie(args...));
      runPairing(collisions, phosclusters, phosclusters, nullptr, nullptr, perCollision_phos, perCollision_phos, fPHOSCut, fPHOSCut, nullptr, nullptr);
    }
    // else if constexpr (pairtype == PairType::kPCMEMC) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto emcclusters = std::get<2>(std::tie(args...));
    //   auto emcmatchedtracks = std::get<3>(std::tie(args...));
    //   runPairing(collisions, v0photons, emcclusters, v0legs, nullptr, perCollision_pcm, perCollision_emc, fV0PhotonCut, fEMCCut, emcmatchedtracks, nullptr);
    // } else if constexpr (pairtype == PairType::kPCMPHOS) {
    //   auto v0photons = std::get<0>(std::tie(args...));
    //   auto v0legs = std::get<1>(std::tie(args...));
    //   auto phosclusters = std::get<2>(std::tie(args...));
    //   runPairing(collisions, v0photons, phosclusters, v0legs, nullptr, perCollision_pcm, perCollision_phos, fV0PhotonCut, fPHOSCut, nullptr, nullptr);
    // }
    ndf++;
  }
  PROCESS_SWITCH(Pi0EtaToGammaGamma_v2, processAnalysis, "process pair analysis", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(Pi0EtaToGammaGamma_v2, processDummy, "Dummy function", true);
};

struct Pi0EtaToGammaGamma {

  Configurable<bool> cfgDoFlow{"cfgDoFlow", false, "flag to analyze vn"};

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999, "max. centrality"};

  Configurable<float> maxY{"maxY", 0.9, "maximum rapidity for reconstructed particles"};
  Configurable<std::string> fConfigPCMCuts{"cfgPCMCuts", "qc,nocut", "Comma separated list of V0 photon cuts"};
  Configurable<std::string> fConfigDalitzEECuts{"cfgDalitzEECuts", "mee_0_120_tpchadrejortofreq_lowB,mee_120_500_tpchadrejortofreq_lowB,mee_0_500_tpchadrejortofreq_lowB", "Comma separated list of Dalitz ee cuts"};
  Configurable<std::string> fConfigDalitzMuMuCuts{"cfgDalitzMuMuCuts", "mmumu_0_500_tpctof_lowB", "Comma separated list of Dalitz mumu cuts"};
  Configurable<std::string> fConfigPHOSCuts{"cfgPHOSCuts", "test02,test03", "Comma separated list of PHOS photon cuts"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "nocut", "Comma separated list of pair cuts"};

  Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
  Configurable<std::string> fConfigEMCCuts{"cfgEMCCuts", "custom,standard,nocut", "Comma separated list of EMCal photon cuts"};

  // Configurable for EMCal cuts
  Configurable<bool> requireCaloReadout{"requireCaloReadout", true, "Require calorimeters readout when analyzing EMCal/PHOS"};
  Configurable<float> EMC_minTime{"EMC_minTime", -20., "Minimum cluster time for EMCal time cut"};
  Configurable<float> EMC_maxTime{"EMC_maxTime", +25., "Maximum cluster time for EMCal time cut"};
  Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
  Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
  Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
  Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
  Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
  Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
  Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};

  Configurable<std::string> fConfigEMEventCut{"cfgEMEventCut", "minbias", "em event cut"}; // only 1 event cut per wagon
  EMEventCut fEMEventCut;
  static constexpr std::string_view event_types[2] = {"before", "after"};

  OutputObj<THashList> fOutputEvent{"Event"};
  OutputObj<THashList> fOutputPair{"Pair"}; // 2-photon pair
  THashList* fMainList = new THashList();

  std::vector<V0PhotonCut> fPCMCuts;
  std::vector<DalitzEECut> fDalitzEECuts;
  std::vector<DalitzEECut> fDalitzMuMuCuts;
  std::vector<PHOSPhotonCut> fPHOSCuts;
  std::vector<EMCPhotonCut> fEMCCuts;
  std::vector<PairCut> fPairCuts;
  std::vector<std::string> fPairNames;

  void init(InitContext& context)
  {
    if (context.mOptions.get<bool>("processPCMPCM")) {
      fPairNames.push_back("PCMPCM");
    }
    if (context.mOptions.get<bool>("processPHOSPHOS")) {
      fPairNames.push_back("PHOSPHOS");
    }
    if (context.mOptions.get<bool>("processEMCEMC")) {
      fPairNames.push_back("EMCEMC");
    }
    if (context.mOptions.get<bool>("processPCMPHOS")) {
      fPairNames.push_back("PCMPHOS");
    }
    if (context.mOptions.get<bool>("processPCMEMC")) {
      fPairNames.push_back("PCMEMC");
    }
    if (context.mOptions.get<bool>("processPCMDalitzEE")) {
      fPairNames.push_back("PCMDalitzEE");
    }
    if (context.mOptions.get<bool>("processPCMDalitzMuMu")) {
      fPairNames.push_back("PCMDalitzMuMu");
    }
    if (context.mOptions.get<bool>("processPHOSEMC")) {
      fPairNames.push_back("PHOSEMC");
    }

    DefinePCMCuts();
    DefineDalitzEECuts();
    DefineDalitzMuMuCuts();
    DefinePHOSCuts();
    DefineEMCCuts();
    DefinePairCuts();
    addhistograms();

    TString ev_cut_name = fConfigEMEventCut.value;
    fEMEventCut = *eventcuts::GetCut(ev_cut_name.Data());

    fOutputEvent.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Event")));
    fOutputPair.setObject(reinterpret_cast<THashList*>(fMainList->FindObject("Pair")));
  }

  template <typename TCuts1, typename TCuts2, typename TCuts3>
  void add_pair_histograms(THashList* list_pair, const std::string pairname, TCuts1 const& cuts1, TCuts2 const& cuts2, TCuts3 const& cuts3)
  {
    for (auto& cut1 : cuts1) {
      for (auto& cut2 : cuts2) {
        std::string cutname1 = cut1.GetName();
        std::string cutname2 = cut2.GetName();

        if ((pairname == "PCMPCM" || pairname == "PHOSPHOS" || pairname == "EMCEMC") && (cutname1 != cutname2))
          continue;

        THashList* list_pair_subsys = reinterpret_cast<THashList*>(list_pair->FindObject(pairname.data()));
        std::string photon_cut_name = cutname1 + "_" + cutname2;
        o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys, photon_cut_name.data());
        THashList* list_pair_subsys_photoncut = reinterpret_cast<THashList*>(list_pair_subsys->FindObject(photon_cut_name.data()));

        for (auto& cut3 : cuts3) {
          std::string pair_cut_name = cut3.GetName();
          o2::aod::pwgem::photon::histogram::AddHistClass(list_pair_subsys_photoncut, pair_cut_name.data());
          THashList* list_pair_subsys_paircut = reinterpret_cast<THashList*>(list_pair_subsys_photoncut->FindObject(pair_cut_name.data()));

          if (cfgDoFlow) {
            o2::aod::pwgem::photon::histogram::DefineHistograms(list_pair_subsys_paircut, "gammagamma_mass_pt", Form("%s,%s", pairname.data(), ",qvector"));
          } else {
            o2::aod::pwgem::photon::histogram::DefineHistograms(list_pair_subsys_paircut, "gammagamma_mass_pt", pairname.data());
          }
        } // end of cut3 loop pair cut
      }   // end of cut2 loop
    }     // end of cut1 loop
  }

  static constexpr std::string_view pairnames[9] = {"PCMPCM", "PHOSPHOS", "EMCEMC", "PCMPHOS", "PCMEMC", "PCMDalitzEE", "PCMDalitzMuMu", "PHOSEMC", "DalitzEEDalitzEE"};
  void addhistograms()
  {
    fMainList->SetOwner(true);
    fMainList->SetName("fMainList");

    // create sub lists first.
    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Event");
    THashList* list_ev = reinterpret_cast<THashList*>(fMainList->FindObject("Event"));

    o2::aod::pwgem::photon::histogram::AddHistClass(fMainList, "Pair");
    THashList* list_pair = reinterpret_cast<THashList*>(fMainList->FindObject("Pair"));

    for (auto& pairname : fPairNames) {
      LOGF(info, "Enabled pairs = %s", pairname.data());

      THashList* list_ev_pair = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev, pairname.data()));
      for (const auto& evtype : event_types) {
        THashList* list_ev_type = reinterpret_cast<THashList*>(o2::aod::pwgem::photon::histogram::AddHistClass(list_ev_pair, evtype.data()));

        if (cfgDoFlow) {
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", "qvector");
        } else {
          o2::aod::pwgem::photon::histogram::DefineHistograms(list_ev_type, "Event", "");
        }
      }

      o2::aod::pwgem::photon::histogram::AddHistClass(list_pair, pairname.data());

      if (pairname == "PCMPCM") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPCMCuts, fPairCuts);
      }
      if (pairname == "PHOSPHOS") {
        add_pair_histograms(list_pair, pairname, fPHOSCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "EMCEMC") {
        add_pair_histograms(list_pair, pairname, fEMCCuts, fEMCCuts, fPairCuts);
      }
      if (pairname == "PCMPHOS") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fPHOSCuts, fPairCuts);
      }
      if (pairname == "PCMEMC") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fEMCCuts, fPairCuts);
      }
      if (pairname == "PCMDalitzEE") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fDalitzEECuts, fPairCuts);
      }
      if (pairname == "PCMDalitzMuMu") {
        add_pair_histograms(list_pair, pairname, fPCMCuts, fDalitzMuMuCuts, fPairCuts);
      }
      if (pairname == "PHOSEMC") {
        add_pair_histograms(list_pair, pairname, fPHOSCuts, fEMCCuts, fPairCuts);
      }

    } // end of pair name loop
  }

  void DefinePCMCuts()
  {
    TString cutNamesStr = fConfigPCMCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPCMCuts.push_back(*pcmcuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PCM cuts = %d", fPCMCuts.size());
  }

  void DefineDalitzEECuts()
  {
    TString cutNamesStr = fConfigDalitzEECuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fDalitzEECuts.push_back(*dalitzeecuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of DalitzEE cuts = %d", fDalitzEECuts.size());
  }
  void DefineDalitzMuMuCuts()
  {
    TString cutNamesStr = fConfigDalitzMuMuCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fDalitzMuMuCuts.push_back(*dalitzeecuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of DalitzMuMu cuts = %d", fDalitzMuMuCuts.size());
  }

  void DefinePHOSCuts()
  {
    TString cutNamesStr = fConfigPHOSCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPHOSCuts.push_back(*phoscuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of PHOS cuts = %d", fPHOSCuts.size());
  }

  void DefineEMCCuts()
  {
    const float a = EMC_TM_Eta->at(0);
    const float b = EMC_TM_Eta->at(1);
    const float c = EMC_TM_Eta->at(2);

    const float d = EMC_TM_Phi->at(0);
    const float e = EMC_TM_Phi->at(1);
    const float f = EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);

    TString cutNamesStr = fConfigEMCCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        if (std::strcmp(cutname, "custom") == 0) {
          EMCPhotonCut* custom_cut = new EMCPhotonCut(cutname, cutname);
          custom_cut->SetMinE(EMC_minE);
          custom_cut->SetMinNCell(EMC_minNCell);
          custom_cut->SetM02Range(EMC_minM02, EMC_maxM02);
          custom_cut->SetTimeRange(EMC_minTime, EMC_maxTime);

          custom_cut->SetTrackMatchingEta([&a, &b, &c](float pT) {
            return a + pow(pT + b, c);
          });
          custom_cut->SetTrackMatchingPhi([&d, &e, &f](float pT) {
            return d + pow(pT + e, f);
          });

          custom_cut->SetMinEoverP(EMC_Eoverp);
          custom_cut->SetUseExoticCut(EMC_UseExoticCut);
          fEMCCuts.push_back(*custom_cut);
        } else {
          fEMCCuts.push_back(*emccuts::GetCut(cutname));
        }
      }
    }
    LOGF(info, "Number of EMCal cuts = %d", fEMCCuts.size());
  }

  void DefinePairCuts()
  {
    TString cutNamesStr = fConfigPairCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        const char* cutname = objArray->At(icut)->GetName();
        LOGF(info, "add cut : %s", cutname);
        fPairCuts.push_back(*paircuts::GetCut(cutname));
      }
    }
    LOGF(info, "Number of Pair cuts = %d", fPairCuts.size());
  }

  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  Preslice<MyDalitzEEs> perCollision_dalitzee = aod::dalitzee::emeventId;
  Preslice<MyDalitzMuMus> perCollision_dalitzmumu = aod::dalitzmumu::emeventId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::skimmedcluster::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  template <PairType pairtype, typename TG1, typename TG2, typename TCut1, typename TCut2>
  bool IsSelectedPair(TG1 const& g1, TG2 const& g2, TCut1 const& cut1, TCut2 const& cut2)
  {
    bool is_selected_pair = false;
    if constexpr (pairtype == PairType::kPCMPCM) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, aod::V0Legs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPHOSPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<int, int>(g1, g2, cut1, cut2); // dummy, because track matching is not ready.
    } else if constexpr (pairtype == PairType::kEMCEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::SkimEMCMTs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMPHOS) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, int>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMDalitzEE) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, MyPrimaryElectrons>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<aod::V0Legs, MyPrimaryMuons>(g1, g2, cut1, cut2);
    } else if constexpr (pairtype == PairType::kPHOSEMC) {
      is_selected_pair = o2::aod::photonpair::IsSelectedPair<int, aod::SkimEMCMTs>(g1, g2, cut1, cut2);
    } else {
      is_selected_pair = true;
    }
    return is_selected_pair;
  }

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons, typename TEMCMTs>
  void SameEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& /*legs*/, TEMPrimaryElectrons const& /*emprimaryelectrons*/, TEMPrimaryMuons const& /*emprimarymuons*/, TEMCMTs const& emcmatchedtracks)
  {
    THashList* list_ev_pair_before = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[0].data()));
    THashList* list_ev_pair_after = static_cast<THashList*>(fMainList->FindObject("Event")->FindObject(pairnames[pairtype].data())->FindObject(event_types[1].data()));
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));
    double values[3] = {0, 0, 0};

    for (auto& collision : collisions) {
      if ((pairtype == PairType::kPHOSPHOS || pairtype == PairType::kPCMPHOS) && !collision.alias_bit(triggerAliases::kTVXinPHOS)) {
        continue;
      }
      if ((pairtype == PairType::kEMCEMC || pairtype == PairType::kPCMEMC) && ((!collision.alias_bit(triggerAliases::kTVXinEMC) && requireCaloReadout) || collision.ncollsPerBC() != 1)) {
        continue;
      }

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      if (cfgDoFlow) {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent_Cent_Qvec>(list_ev_pair_before, "", collision);
      } else {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_before, "", collision);
      }

      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      if (cfgDoFlow) {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent_Cent_Qvec>(list_ev_pair_after, "", collision);
      } else {
        o2::aod::pwgem::photon::histogram::FillHistClass<EMHistType::kEvent>(list_ev_pair_after, "", collision);
      }

      reinterpret_cast<TH1F*>(list_ev_pair_before->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      reinterpret_cast<TH1F*>(list_ev_pair_after->FindObject("hCollisionCounter"))->Fill("accepted", 1.f);
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};

      auto photons1_coll = photons1.sliceBy(perCollision1, collision.globalIndex());
      auto photons2_coll = photons2.sliceBy(perCollision2, collision.globalIndex());

      if constexpr (pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) {
        for (auto& cut : cuts1) {
          for (auto& paircut : paircuts) {

            for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons1_coll, photons2_coll))) {
              if (!IsSelectedPair<pairtype>(g1, g2, cut, cut)) {
                continue;
              }

              if (!paircut.IsSelected(g1, g2)) {
                continue;
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }

              if (cfgDoFlow) {
                values[0] = v12.M(), values[1] = v12.Pt();
                std::array<float, 2> u2_gg = {static_cast<float>(std::cos(2 * v12.Phi())), static_cast<float>(std::sin(2 * v12.Phi()))};
                values[2] = RecoDecay::dotProd(u2_gg, q2ft0m);
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0M"))->Fill(values);
                values[2] = RecoDecay::dotProd(u2_gg, q2ft0a);
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0A"))->Fill(values);
                values[2] = RecoDecay::dotProd(u2_gg, q2ft0c);
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0C"))->Fill(values);
                values[2] = RecoDecay::dotProd(u2_gg, q2fv0a);
                reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FV0A"))->Fill(values);
              } else {
                reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same"))->Fill(v12.M(), v12.Pt());
                if constexpr (pairtype == PairType::kEMCEMC) {
                  RotationBackground<aod::SkimEMCClusters>(v12, v1, v2, photons2_coll, g1.globalIndex(), g2.globalIndex(), cut, paircut, emcmatchedtracks);
                }
              }

            } // end of combination
          }   // end of pair cut loop
        }     // end of cut loop

      } else { // different subsystem pairs
        for (auto& cut1 : cuts1) {
          for (auto& cut2 : cuts2) {
            for (auto& paircut : paircuts) {
              for (auto& [g1, g2] : combinations(CombinationsFullIndexPolicy(photons1_coll, photons2_coll))) {

                if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) {
                  auto pos_pv = g2.template posTrack_as<MyPrimaryElectrons>();
                  auto ele_pv = g2.template negTrack_as<MyPrimaryElectrons>();
                  std::tuple<MyPrimaryElectron, MyPrimaryElectron, float> pair2 = std::make_tuple(pos_pv, ele_pv, collision.bz());
                  if (!IsSelectedPair<pairtype>(g1, pair2, cut1, cut2)) {
                    continue;
                  }
                } else {
                  if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
                    continue;
                  }
                }

                if (!paircut.IsSelected(g1, g2)) {
                  continue;
                }

                if constexpr (pairtype == PairType::kPCMPHOS || pairtype == PairType::kPCMEMC) {
                  auto pos = g1.template posTrack_as<aod::V0Legs>();
                  auto ele = g1.template negTrack_as<aod::V0Legs>();

                  for (auto& v0leg : {pos, ele}) {
                    float deta = v0leg.eta() - g2.eta();
                    float dphi = TVector2::Phi_mpi_pi(TVector2::Phi_0_2pi(v0leg.phi()) - TVector2::Phi_0_2pi(g2.phi()));
                    float Ep = g2.e() / v0leg.p();
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hdEtadPhi"))->Fill(dphi, deta);
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hdEtaPt"))->Fill(v0leg.pt(), deta);
                    reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hdPhiPt"))->Fill(v0leg.pt(), dphi);
                    if (pow(deta / 0.02, 2) + pow(dphi / 0.4, 2) < 1) {
                      reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hEp_E"))->Fill(g2.e(), Ep);
                    }
                  }

                  if constexpr (pairtype == PairType::kPCMPHOS) {
                    if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos, g2, 0.02, 0.4, 0.2) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele, g2, 0.02, 0.4, 0.2)) {
                      continue;
                    }
                  } else if constexpr (pairtype == PairType::kPCMEMC) {
                    if (o2::aod::photonpair::DoesV0LegMatchWithCluster(pos, g2, 0.02, 0.4, 0.5) || o2::aod::photonpair::DoesV0LegMatchWithCluster(ele, g2, 0.02, 0.4, 0.5)) {
                      continue;
                    }
                  }
                }

                ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
                ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
                if constexpr (pairtype == PairType::kPCMDalitzEE) {
                  v2.SetM(g2.mass());
                  auto pos_sv = g1.template posTrack_as<aod::V0Legs>();
                  auto ele_sv = g1.template negTrack_as<aod::V0Legs>();
                  auto pos_pv = g2.template posTrack_as<MyPrimaryElectrons>();
                  auto ele_pv = g2.template negTrack_as<MyPrimaryElectrons>();
                  if (pos_sv.trackId() == pos_pv.trackId() || ele_sv.trackId() == ele_pv.trackId()) {
                    continue;
                  }
                } else if constexpr (pairtype == PairType::kPCMDalitzMuMu) {
                  v2.SetM(g2.mass());
                  auto pos_sv = g1.template posTrack_as<aod::V0Legs>();
                  auto ele_sv = g1.template negTrack_as<aod::V0Legs>();
                  auto pos_pv = g2.template posTrack_as<MyPrimaryMuons>();
                  auto ele_pv = g2.template negTrack_as<MyPrimaryMuons>();
                  if (pos_sv.trackId() == pos_pv.trackId() || ele_sv.trackId() == ele_pv.trackId()) {
                    continue;
                  }
                }

                ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
                if (abs(v12.Rapidity()) > maxY) {
                  continue;
                }
                if (cfgDoFlow) {
                  values[0] = v12.M(), values[1] = v12.Pt();
                  std::array<float, 2> u_gg = {static_cast<float>(v12.Px() / v12.Pt()), static_cast<float>(v12.Py() / v12.Pt())};
                  values[2] = RecoDecay::dotProd(u_gg, q2ft0m);
                  reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0M"))->Fill(values);
                  values[2] = RecoDecay::dotProd(u_gg, q2ft0a);
                  reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0A"))->Fill(values);
                  values[2] = RecoDecay::dotProd(u_gg, q2ft0c);
                  reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FT0C"))->Fill(values);
                  values[2] = RecoDecay::dotProd(u_gg, q2fv0a);
                  reinterpret_cast<THnSparseF*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hs_Same_SPQ2FV0A"))->Fill(values);
                } else {
                  reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same"))->Fill(v12.M(), v12.Pt());
                }

              } // end of combination
            }   // end of pair cut loop
          }     // end of cut2 loop
        }       // end of cut1 loop
      }
    } // end of collision loop
  }

  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 999.f}, "Mixing bins - centrality"};
  using BinningType_M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningType_A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningType_C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningType_M colBinning_M{{ConfVtxBins, ConfCentBins}, true};
  BinningType_A colBinning_A{{ConfVtxBins, ConfCentBins}, true};
  BinningType_C colBinning_C{{ConfVtxBins, ConfCentBins}, true};

  template <PairType pairtype, typename TEvents, typename TPhotons1, typename TPhotons2, typename TPreslice1, typename TPreslice2, typename TCuts1, typename TCuts2, typename TPairCuts, typename TLegs, typename TEMPrimaryElectrons, typename TEMPrimaryMuons, typename TEMCMTs, typename TMixedBinning>
  void MixedEventPairing(TEvents const& collisions, TPhotons1 const& photons1, TPhotons2 const& photons2, TPreslice1 const& perCollision1, TPreslice2 const& perCollision2, TCuts1 const& cuts1, TCuts2 const& cuts2, TPairCuts const& paircuts, TLegs const& /*legs*/, TEMPrimaryElectrons const& /*emprimaryelectrons*/, TEMPrimaryMuons const& /*emprimarymuons*/, TEMCMTs const& /*emcmatchedtracks*/, TMixedBinning const& colBinning)
  {
    THashList* list_pair_ss = static_cast<THashList*>(fMainList->FindObject("Pair")->FindObject(pairnames[pairtype].data()));
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ndepth, -1, collisions, collisions)) { // internally, CombinationsStrictlyUpperIndexPolicy(collisions, collisions) is called.

      // LOGF(info, "Mixed event globalIndex: (%d, %d) , ngpcm: (%d, %d), ngphos: (%d, %d), ngemc: (%d, %d)", collision1.globalIndex(), collision2.globalIndex(), collision1.ngpcm(), collision2.ngpcm(), collision1.ngphos(), collision2.ngphos(), collision1.ngemc(), collision2.ngemc());

      const float centralities1[3] = {collision1.centFT0M(), collision1.centFT0A(), collision1.centFT0C()};
      const float centralities2[3] = {collision2.centFT0M(), collision2.centFT0A(), collision2.centFT0C()};

      if (centralities1[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities1[cfgCentEstimator]) {
        continue;
      }
      if (centralities2[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities2[cfgCentEstimator]) {
        continue;
      }
      if (!fEMEventCut.IsSelected(collision1) || !fEMEventCut.IsSelected(collision2)) {
        continue;
      }

      auto photons_coll1 = photons1.sliceBy(perCollision1, collision1.globalIndex());
      auto photons_coll2 = photons2.sliceBy(perCollision2, collision2.globalIndex());
      // LOGF(info, "collision1: posZ = %f, numContrib = %d , sel8 = %d | collision2: posZ = %f, numContrib = %d , sel8 = %d", collision1.posZ(), collision1.numContrib(), collision1.sel8(), collision2.posZ(), collision2.numContrib(), collision2.sel8());

      for (auto& cut1 : cuts1) {
        for (auto& cut2 : cuts2) {
          for (auto& paircut : paircuts) {
            for (auto& [g1, g2] : combinations(soa::CombinationsFullIndexPolicy(photons_coll1, photons_coll2))) {
              // LOGF(info, "Mixed event photon pair: (%d, %d) from events (%d, %d), photon event: (%d, %d)", g1.index(), g2.index(), collision1.index(), collision2.index(), g1.globalIndex(), g2.globalIndex());

              if ((pairtype == PairType::kPCMPCM || pairtype == PairType::kPHOSPHOS || pairtype == PairType::kEMCEMC) && (TString(cut1.GetName()) != TString(cut2.GetName()))) {
                continue;
              }

              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) {
                auto pos_pv = g2.template posTrack_as<MyPrimaryElectrons>();
                auto ele_pv = g2.template negTrack_as<MyPrimaryElectrons>();
                std::tuple<MyPrimaryElectron, MyPrimaryElectron, float> pair2 = std::make_tuple(pos_pv, ele_pv, collision1.bz());
                if (!IsSelectedPair<pairtype>(g1, pair2, cut1, cut2)) {
                  continue;
                }
              } else {
                if (!IsSelectedPair<pairtype>(g1, g2, cut1, cut2)) {
                  continue;
                }
              }

              if (!paircut.IsSelected(g1, g2)) {
                continue;
              }

              ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
              ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
              if constexpr (pairtype == PairType::kPCMDalitzEE || pairtype == PairType::kPCMDalitzMuMu) {
                v2.SetM(g2.mass());
              }
              ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
              if (abs(v12.Rapidity()) > maxY) {
                continue;
              }
              reinterpret_cast<TH2F*>(list_pair_ss->FindObject(Form("%s_%s", cut1.GetName(), cut2.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Mixed"))->Fill(v12.M(), v12.Pt());

            } // end of different photon combinations
          }   // end of pair cut loop
        }     // end of cut2 loop
      }       // end of cut1 loop
    }         // end of different collision combinations
  }

  /// \brief Calculate background (using rotation background method only for EMCal!)
  template <typename TPhotons>
  void RotationBackground(const ROOT::Math::PtEtaPhiMVector& meson, ROOT::Math::PtEtaPhiMVector photon1, ROOT::Math::PtEtaPhiMVector photon2, TPhotons const& photons_coll, unsigned int ig1, unsigned int ig2, EMCPhotonCut const& cut, PairCut const& paircut, SkimEMCMTs const& /*emcmatchedtracks*/)
  {
    // if less than 3 clusters are present skip event since we need at least 3 clusters
    if (photons_coll.size() < 3) {
      return;
    }
    const double rotationAngle = M_PI / 2.0; // rotaion angle 90°
    // ROOT::Math::XYZVector meson3D(meson.Px(), meson.Py(), meson.Pz());
    ROOT::Math::AxisAngle rotationAxis(meson.Vect(), rotationAngle);
    // LOG(info) << "rotationAxis.Angle() = " << rotationAxis.Angle();
    ROOT::Math::Rotation3D rotationMatrix(rotationAxis);
    // ROOT::Math::XYZVector test(photon1.Px(), photon1.Py(), photon1.Pz());
    // ROOT::Math::XYZVector photon1_3D(photon1.Px(), photon1.Py(), photon1.Pz());
    // ROOT::Math::XYZVector photon2_3D(photon2.Px(), photon2.Py(), photon2.Pz());

    // float openingAngleTest = std::acos(photon1_3D.Dot(test) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(test.Mag2())));
    // LOG(info) << "openingAngleTest before rotation = " << openingAngleTest;
    photon1 = rotationMatrix * photon1;
    photon2 = rotationMatrix * photon2;

    // openingAngleTest = std::acos(photon1_3D.Dot(test) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(test.Mag2())));
    // LOG(info) << "openingAngleTest = " << openingAngleTest;

    for (auto& photon : photons_coll) {
      if (photon.globalIndex() == ig1 || photon.globalIndex() == ig2) {
        // only combine rotated photons with other photons
        continue;
      }
      if (!cut.template IsSelected<aod::SkimEMCMTs>(photon)) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector photon3;
      photon3.SetPt(photon.pt());
      photon3.SetEta(photon.eta());
      photon3.SetPhi(photon.phi());
      photon3.SetM(0.);
      // ROOT::Math::XYZVector photon3_3D(photon3.Px(), photon3.Py(), photon3.Pz());
      ROOT::Math::PtEtaPhiMVector mother1 = photon1 + photon3;
      ROOT::Math::PtEtaPhiMVector mother2 = photon2 + photon3;

      float openingAngle1 = std::acos(photon1.Vect().Dot(photon3.Vect()) / (photon1.P() * photon3.P()));
      float openingAngle2 = std::acos(photon2.Vect().Dot(photon3.Vect()) / (photon2.P() * photon3.P()));
      // float openingAngle1_1 = std::acos(photon1_3D.Dot(photon3_3D) / (std::sqrt(photon1_3D.Mag2()) * std::sqrt(photon3_3D.Mag2())));
      // float openingAngle2_2 = std::acos(photon2_3D.Dot(photon3_3D) / (std::sqrt(photon2_3D.Mag2()) * std::sqrt(photon3_3D.Mag2())));
      // LOG(info) << "openingAngle1 = " << openingAngle1;
      // LOG(info) << "openingAngle2 = " << openingAngle2;
      // LOG(info) << "openingAngle1_1 = " << openingAngle1_1;
      // LOG(info) << "openingAngle2_2 = " << openingAngle2_2;

      // Fill histograms
      if (openingAngle1 > minOpenAngle) {
        reinterpret_cast<TH2F*>(fMainList->FindObject("Pair")->FindObject("EMCEMC")->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same_RotatedBkg"))->Fill(mother1.M(), mother1.Pt());
      }
      if (openingAngle2 > minOpenAngle) {
        reinterpret_cast<TH2F*>(fMainList->FindObject("Pair")->FindObject("EMCEMC")->FindObject(Form("%s_%s", cut.GetName(), cut.GetName()))->FindObject(paircut.GetName())->FindObject("hMggPt_Same_RotatedBkg"))->Fill(mother2.M(), mother2.Pt());
      }
    }
  }

  Filter DalitzEEFilter = o2::aod::dalitzee::sign == 0; // analyze only uls
  using MyFilteredDalitzEEs = soa::Filtered<MyDalitzEEs>;

  Filter DalitzMuMuFilter = o2::aod::dalitzmumu::sign == 0; // analyze only uls
  using MyFilteredDalitzMuMus = soa::Filtered<MyDalitzMuMus>;

  Partition<MyCollisions> grouped_collisions = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax); // this goes to same event.
  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  Filter collisionFilter_subsys = (o2::aod::emevent::ngpcm >= 1) || (o2::aod::emevent::ngphos >= 1) || (o2::aod::emevent::ngemc >= 1) || (aod::emevent::neeuls >= 1);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to mixed event.

  void processPCMPCM(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMPCM>(grouped_collisions, v0photons, v0photons, perCollision, perCollision, fPCMCuts, fPCMCuts, fPairCuts, legs, nullptr, nullptr, nullptr);
    if (cfgCentEstimator == 0) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fPCMCuts, fPCMCuts, fPairCuts, legs, nullptr, nullptr, nullptr, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fPCMCuts, fPCMCuts, fPairCuts, legs, nullptr, nullptr, nullptr, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing<PairType::kPCMPCM>(filtered_collisions, v0photons, v0photons, perCollision, perCollision, fPCMCuts, fPCMCuts, fPairCuts, legs, nullptr, nullptr, nullptr, colBinning_C);
    }
  }

  void processPHOSPHOS(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters)
  {
    SameEventPairing<PairType::kPHOSPHOS>(grouped_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fPHOSCuts, fPHOSCuts, fPairCuts, nullptr, nullptr, nullptr, nullptr);
    MixedEventPairing<PairType::kPHOSPHOS>(filtered_collisions, phosclusters, phosclusters, perCollision_phos, perCollision_phos, fPHOSCuts, fPHOSCuts, fPairCuts, nullptr, nullptr, nullptr, nullptr, colBinning_C);
  }

  void processEMCEMC(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, aod::SkimEMCClusters const& emcclusters, aod::SkimEMCMTs const& emcmatchedtracks)
  {
    SameEventPairing<PairType::kEMCEMC>(grouped_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fEMCCuts, fEMCCuts, fPairCuts, nullptr, nullptr, nullptr, emcmatchedtracks);
    MixedEventPairing<PairType::kEMCEMC>(filtered_collisions, emcclusters, emcclusters, perCollision_emc, perCollision_emc, fEMCCuts, fEMCCuts, fPairCuts, nullptr, nullptr, nullptr, emcmatchedtracks, colBinning_C);
  }

  void processPCMDalitzEE(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs, MyFilteredDalitzEEs const& dileptons, MyPrimaryElectrons const& emprimaryelectrons)
  {
    SameEventPairing<PairType::kPCMDalitzEE>(grouped_collisions, v0photons, dileptons, perCollision, perCollision_dalitzee, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, nullptr, nullptr);
    if (cfgCentEstimator == 0) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzee, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, nullptr, nullptr, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzee, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, nullptr, nullptr, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing<PairType::kPCMDalitzEE>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzee, fPCMCuts, fDalitzEECuts, fPairCuts, legs, emprimaryelectrons, nullptr, nullptr, colBinning_C);
    }
  }

  void processPCMDalitzMuMu(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::V0Legs const& legs, MyFilteredDalitzMuMus const& dileptons, MyPrimaryMuons const& emprimarymuons)
  {
    SameEventPairing<PairType::kPCMDalitzMuMu>(grouped_collisions, v0photons, dileptons, perCollision, perCollision_dalitzmumu, fPCMCuts, fDalitzMuMuCuts, fPairCuts, legs, nullptr, emprimarymuons, nullptr);
    if (cfgCentEstimator == 0) {
      MixedEventPairing<PairType::kPCMDalitzMuMu>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzmumu, fPCMCuts, fDalitzMuMuCuts, fPairCuts, legs, nullptr, emprimarymuons, nullptr, colBinning_M);
    } else if (cfgCentEstimator == 1) {
      MixedEventPairing<PairType::kPCMDalitzMuMu>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzmumu, fPCMCuts, fDalitzMuMuCuts, fPairCuts, legs, nullptr, emprimarymuons, nullptr, colBinning_A);
    } else if (cfgCentEstimator == 2) {
      MixedEventPairing<PairType::kPCMDalitzMuMu>(filtered_collisions, v0photons, dileptons, perCollision, perCollision_dalitzmumu, fPCMCuts, fDalitzMuMuCuts, fPairCuts, legs, nullptr, emprimarymuons, nullptr, colBinning_C);
    }
  }

  void processPCMPHOS(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::PHOSClusters const& phosclusters, aod::V0Legs const& legs)
  {
    SameEventPairing<PairType::kPCMPHOS>(grouped_collisions, v0photons, phosclusters, perCollision, perCollision_phos, fPCMCuts, fPHOSCuts, fPairCuts, legs, nullptr, nullptr, nullptr);
    MixedEventPairing<PairType::kPCMPHOS>(filtered_collisions, v0photons, phosclusters, perCollision, perCollision_phos, fPCMCuts, fPHOSCuts, fPairCuts, legs, nullptr, nullptr, nullptr, colBinning_C);
  }

  void processPCMEMC(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, MyV0Photons const& v0photons, aod::SkimEMCClusters const& emcclusters, aod::V0Legs const& legs, aod::SkimEMCMTs const& emcmatchedtracks)
  {
    SameEventPairing<PairType::kPCMEMC>(grouped_collisions, v0photons, emcclusters, perCollision, perCollision_emc, fPCMCuts, fEMCCuts, fPairCuts, legs, nullptr, nullptr, emcmatchedtracks);
    MixedEventPairing<PairType::kPCMEMC>(filtered_collisions, v0photons, emcclusters, perCollision, perCollision_emc, fPCMCuts, fEMCCuts, fPairCuts, legs, nullptr, nullptr, emcmatchedtracks, colBinning_C);
  }

  void processPHOSEMC(MyCollisions const&, MyFilteredCollisions const& filtered_collisions, aod::PHOSClusters const& phosclusters, aod::SkimEMCClusters const& emcclusters, aod::SkimEMCMTs const& emcmatchedtracks)
  {
    SameEventPairing<PairType::kPHOSEMC>(grouped_collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc, fPHOSCuts, fEMCCuts, fPairCuts, nullptr, nullptr, nullptr, emcmatchedtracks);
    MixedEventPairing<PairType::kPHOSEMC>(filtered_collisions, phosclusters, emcclusters, perCollision_phos, perCollision_emc, fPHOSCuts, fEMCCuts, fPairCuts, nullptr, nullptr, nullptr, emcmatchedtracks, colBinning_C);
  }

  void processDummy(MyCollisions const&) {}

  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMPCM, "pairing PCM-PCM", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPHOSPHOS, "pairing PHOS-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processEMCEMC, "pairing EMCal-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMPHOS, "pairing PCM-PHOS", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMEMC, "pairing PCM-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMDalitzEE, "pairing PCM-DalitzEE", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPCMDalitzMuMu, "pairing PCM-DalitzMuMu", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processPHOSEMC, "pairing PHOS-EMCal", false);
  PROCESS_SWITCH(Pi0EtaToGammaGamma, processDummy, "Dummy function", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaToGammaGamma>(cfgc, TaskName{"pi0eta-to-gammagamma"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPCMPCM, MyV0Photons, aod::V0Legs>>(cfgc, TaskName{"pi0eta-to-gammagamma-pcmpcm"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPCMDalitzEE, MyV0Photons, aod::V0Legs, MyPrimaryElectrons>>(cfgc, TaskName{"pi0eta-to-gammagamma-pcmdalitzee"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPCMDalitzMuMu, MyV0Photons, aod::V0Legs, MyPrimaryMuons>>(cfgc, TaskName{"pi0eta-to-gammagamma-pcmdalitzmumu"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kEMCEMC, MyEMCClusters, aod::SkimEMCMTs>>(cfgc, TaskName{"pi0eta-to-gammagamma-emcemc"}),
    adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPHOSPHOS, MyPHOSClusters>>(cfgc, TaskName{"pi0eta-to-gammagamma-phosphos"}),
    // adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPCMEMC, MyV0Photons, aod::V0Legs, MyEMCClusters, aod::SkimEMCMTs>>(cfgc, TaskName{"pi0eta-to-gammagamma-pcmemc"}),
    // adaptAnalysisTask<Pi0EtaToGammaGamma_v2<PairType::kPCMPHOS, MyV0Photons, aod::V0Legs, MyPHOSClusters>>(cfgc, TaskName{"pi0eta-to-gammagamma-pcmphos"}),
  };
}
