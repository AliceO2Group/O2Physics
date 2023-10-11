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
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }
  if (!nameStr.compare("qc")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_ITSTPC")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSTPC(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 180);
    cut->SetAPRange(0.95, 0.02);
    return cut;
  }
  if (!nameStr.compare("qc_ITSonly")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSonly(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    cut->SetAPRange(0.95, 0.02);
    return cut;
  }
  if (!nameStr.compare("qc_TPConly")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPConly(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_TPCTRD")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPCTRD(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_TPCTOF")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPCTOF(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.5);
    cut->SetRxyRange(1, 180);
    return cut;
  }
  if (!nameStr.compare("qc_TPCTRDTOF")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPCTRDTOF(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
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
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-4, +4);
    cut->SetNClustersITS(2, 4);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetOnWwireIB(true);
    cut->SetOnWwireOB(true);
    cut->SetAPRange(0.95, 0.02);
    return cut;
  }
  if (!nameStr.compare("wwire_ib")) { // conversion only on tungstate wire outside of ITSib
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetNClustersITS(2, 4);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetOnWwireIB(true);
    cut->SetOnWwireOB(false);
    cut->SetAPRange(0.95, 0.02);
    return cut;
  }
  if (!nameStr.compare("wwire_ob")) { // conversion only on tungstate wire outside of ITSob (middle layer)
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetNClustersITS(0, 2);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetOnWwireIB(false);
    cut->SetOnWwireOB(true);
    cut->SetAPRange(0.95, 0.02);
    return cut;
  }
  if (!nameStr.compare("nopid")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
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
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(false);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
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
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSTPC(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

DalitzEECut* o2::aod::dalitzeecuts::GetCut(const char* cutName)
{
  DalitzEECut* cut = new DalitzEECut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("mee_all_tpchadrejortofreq_lowB")) {
    // for pair
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaMuRange(-2, +2);
    cut->SetTPCNsigmaPiRange(-3, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetMuonExclusionTPC(true);
    return cut;
  }

  if (!nameStr.compare("mee_all_tpchadrej_lowB")) {
    // for pair
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrej);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaMuRange(-2, +2);
    cut->SetTPCNsigmaPiRange(-3, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetMuonExclusionTPC(true);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    return cut;
  }

  if (!nameStr.compare("mee_all_tofreq_lowB")) {
    // for pair
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaPiRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_all_nopid_lowB")) {
    // for pair
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kUnDef);
    return cut;
  }
  if (!nameStr.compare("mee_0_120_tpconly")) {
    // for pair
    cut->SetMeeRange(0., 0.12);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    return cut;
  }
  if (!nameStr.compare("mee_120_500_tpconly")) {
    // for pair
    cut->SetMeeRange(0.12, 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_500_tpconly")) {
    // for pair
    cut->SetMeeRange(0., 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_120_tpconly_lowB")) {
    // for pair
    cut->SetMeeRange(0., 0.12);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_120_500_tpconly_lowB")) {
    // for pair
    cut->SetMeeRange(0.12, 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_500_tpconly_lowB")) {
    // for pair
    cut->SetMeeRange(0., 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_120_tpchadrejortofreq_lowB")) {
    // for pair
    cut->SetMeeRange(0., 0.12);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaMuRange(-2, +2);
    cut->SetTPCNsigmaPiRange(-3, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetMuonExclusionTPC(true);
    return cut;
  }
  if (!nameStr.compare("mee_120_500_tpchadrejortofreq_lowB")) {
    // for pair
    cut->SetMeeRange(0.12, 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaMuRange(-2, +2);
    cut->SetTPCNsigmaPiRange(-3, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetMuonExclusionTPC(true);
    return cut;
  }
  if (!nameStr.compare("mee_0_500_tpchadrejortofreq_lowB")) {
    // for pair
    cut->SetMeeRange(0., 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaMuRange(-2, +2);
    cut->SetTPCNsigmaPiRange(-3, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetMuonExclusionTPC(true);
    return cut;
  }

  if (!nameStr.compare("mee_all_tpchadrejortofreq")) {
    // for pair
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_120_tpchadrejortofreq")) {
    // for pair
    cut->SetMeeRange(0, 0.12);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_120_tpchadrejortofreq_wo_phiv")) {
    // for pair
    cut->SetMeeRange(0, 0.12);

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_120_500_tpchadrejortofreq")) {
    // for pair
    cut->SetMeeRange(0.12, 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    return cut;
  }
  if (!nameStr.compare("mee_0_500_tpchadrejortofreq")) {
    // for pair
    cut->SetMeeRange(0., 0.5);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return mee < 0.01 ? 1.5 : (mee < 0.02 ? 2.0 : (mee < 0.03 ? 2.5 : 3.2));
    });

    // for track
    cut->SetTrackPtRange(0.2f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-2, +3);
    cut->SetTPCNsigmaPiRange(-1e+10, +3);
    cut->SetTPCNsigmaKaRange(-3, +3);
    cut->SetTPCNsigmaPrRange(-3, +3);
    cut->SetTOFNsigmaElRange(-3, +3);
    return cut;
  }

  if (!nameStr.compare("nocut")) {
    // for track
    cut->SetTrackPtRange(0.01f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(80);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);
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
