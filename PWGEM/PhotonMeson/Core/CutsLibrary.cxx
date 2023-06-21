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
// Contact: daiki.sekihata@cern.ch
//
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"

V0PhotonCut* o2::aod::pcmcuts::GetCut(const char* cutName)
{
  V0PhotonCut* cut = new V0PhotonCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("analysis")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.998);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    cut->SetMaxMeePsiPairDep([](float psipair) { return psipair < 0.4 ? 0.06 : 0.015; });
    return cut;
  }
  if (!nameStr.compare("analysis_eta12")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-1.2, +1.2);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-1.2, +1.2);
    cut->SetMinCosPA(0.998);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    cut->SetMaxMeePsiPairDep([](float psipair) { return psipair < 0.4 ? 0.06 : 0.015; });
    return cut;
  }
  if (!nameStr.compare("analysis_wo_mee")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.998);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }
  if (!nameStr.compare("analysis_wo_mee_eta12")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-1.2, +1.2);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-1.2, +1.2);
    cut->SetMinCosPA(0.998);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }
  if (!nameStr.compare("qc_w_mee")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    cut->SetMaxMeePsiPairDep([](float psipair) { return psipair < 0.4 ? 0.06 : 0.015; });
    return cut;
  }
  if (!nameStr.compare("qc")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_eta12")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-1.2, +1.2);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-1.2, +1.2);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_ITS")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITS(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_TPConly")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPConly(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_AntiTPConly")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireAntiTPConly(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("wwire")) { // conversion only on tungstate wire
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetOnWwireIB(true);
    return cut;
  }
  if (!nameStr.compare("nopid")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("nocut")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(false);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("tag_track")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.6);
    cut->SetMaxChi2PerClusterTPC(4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITS(true);
    // for v0
    cut->SetV0PtRange(0.01f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.998);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

PHOSPhotonCut* o2::aod::phoscuts::GetCut(const char* cutName)
{
  PHOSPhotonCut* cut = new PHOSPhotonCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("test01")) {
    cut->SetEnergyRange(0.1f, 1e10f);
    return cut;
  }
  if (!nameStr.compare("test02")) {
    cut->SetEnergyRange(0.2f, 1e10f);
    return cut;
  }
  if (!nameStr.compare("test03")) {
    cut->SetEnergyRange(0.3f, 1e10f);
    return cut;
  }
  if (!nameStr.compare("test05")) {
    cut->SetEnergyRange(0.5f, 1e10f);
    return cut;
  }
  if (!nameStr.compare("test10")) {
    cut->SetEnergyRange(1.0f, 1e10f);
    return cut;
  }
  if (!nameStr.compare("tag")) {
    cut->SetEnergyRange(1.0f, 1e10f);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

EMCPhotonCut* o2::aod::emccuts::GetCut(const char* cutName)
{
  EMCPhotonCut* cut = new EMCPhotonCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("standard")) {
    cut->SetMinE(0.7f);
    cut->SetMinNCell(1);
    cut->SetM02Range(0.1f, 0.7f);
    cut->SetTimeRange(-20.f, 25.f);

    cut->SetTrackMatchingEta([](float pT) {
      return 0.01f + pow(pT + 4.07f, -2.5f);
    });
    cut->SetTrackMatchingPhi([](float pT) {
      return 0.015f + pow(pT + 3.65f, -2.f);
    });
    cut->SetMinEoverP(1.75f);
    cut->SetUseExoticCut(true);
    return cut;
  }
  if (!nameStr.compare("nocut")) {
    cut->SetMinE(0.f);
    cut->SetMinNCell(1);
    cut->SetM02Range(0.0f, 1000.f);
    cut->SetTimeRange(-500.f, 500.f);

    cut->SetTrackMatchingEta([](float pT) {
      return -1.f;
    });
    cut->SetTrackMatchingPhi([](float pT) {
      return -1.f;
    });
    cut->SetMinEoverP(0.f);
    cut->SetUseExoticCut(false);
    return cut;
  }
  if (!nameStr.compare("tag")) {
    cut->SetMinE(0.7f);
    cut->SetMinNCell(1);
    cut->SetM02Range(0.1f, 0.7f);
    cut->SetTimeRange(-20.f, 25.f);

    cut->SetTrackMatchingEta([](float pT) {
      return 0.01f + pow(pT + 4.07f, -2.5f);
    });
    cut->SetTrackMatchingPhi([](float pT) {
      return 0.015f + pow(pT + 3.65f, -2.f);
    });
    cut->SetMinEoverP(1.75f);
    cut->SetUseExoticCut(true);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

PairCut* o2::aod::paircuts::GetCut(const char* cutName)
{
  PairCut* cut = new PairCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("nocut")) {
    cut->SetAsymRange(-1e+10f, +1e+10f);
    return cut;
  }
  if (!nameStr.compare("asym08")) {
    cut->SetAsymRange(0, 0.8);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}
