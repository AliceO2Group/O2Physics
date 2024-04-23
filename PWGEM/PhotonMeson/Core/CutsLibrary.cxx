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
#include <string>
#include <vector>
#include <regex>
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"

//_______________________________________________
int customAtoi(const std::string& str)
{
  std::regex pattern(R"(\d+)"); // extract only numbers
  std::smatch match;
  if (std::regex_search(str, match, pattern)) {
    return std::stoi(match.str());
  }
  return -1;
}
//_______________________________________________
std::vector<std::string> splitString(const std::string& str, char delimiter)
{
  std::istringstream iss(str);
  std::vector<std::string> ret;
  for (std::string temp; std::getline(iss, temp, delimiter); ret.push_back(temp))
    ;
  return ret;
}
//_______________________________________________
EMEventCut* o2::aod::pwgem::photon::eventcuts::GetCut(const char* cutName)
{
  EMEventCut* cut = new EMEventCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("nocut")) {
    cut->SetRequireSel8(false);
    cut->SetRequireFT0AND(false);
    cut->SetZvtxRange(-1e+10, +1e+10);
    cut->SetRequireNoTFB(false);
    cut->SetRequireNoITSROFB(false);
    cut->SetRequireNoSameBunchPileup(false);
    cut->SetRequireVertexITSTPC(false);
    cut->SetRequireIsGoodZvtxFT0vsPV(false);
    return cut;
  }

  if (!nameStr.compare("ft0and")) {
    cut->SetRequireSel8(false);
    cut->SetRequireFT0AND(true);
    cut->SetZvtxRange(-10.f, +10.f);
    cut->SetRequireNoTFB(false);
    cut->SetRequireNoITSROFB(false);
    cut->SetRequireNoSameBunchPileup(false);
    cut->SetRequireVertexITSTPC(false);
    cut->SetRequireIsGoodZvtxFT0vsPV(false);
    return cut;
  }

  if (!nameStr.compare("minbias")) {
    cut->SetRequireSel8(true);
    cut->SetRequireFT0AND(true);
    cut->SetZvtxRange(-10.f, +10.f);
    cut->SetRequireNoTFB(false);     // included in sel8
    cut->SetRequireNoITSROFB(false); // included in sel8
    cut->SetRequireNoSameBunchPileup(false);
    cut->SetRequireVertexITSTPC(false);
    cut->SetRequireIsGoodZvtxFT0vsPV(false);

    if (nameStr.find("notfb") != std::string::npos) {
      cut->SetRequireNoTFB(true);
    }
    if (nameStr.find("noitsrofb") != std::string::npos) {
      cut->SetRequireNoITSROFB(true);
    }
    if (nameStr.find("nosbp") != std::string::npos) {
      cut->SetRequireNoSameBunchPileup(true);
    }
    if (nameStr.find("vtxitstpc") != std::string::npos) {
      cut->SetRequireVertexITSTPC(true);
    }
    if (nameStr.find("goodvtx") != std::string::npos) {
      cut->SetRequireIsGoodZvtxFT0vsPV(true);
    }

    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

V0PhotonCut* o2::aod::pwgem::photon::pcmcuts::GetCut(const char* cutName)
{
  V0PhotonCut* cut = new V0PhotonCut(cutName, cutName);
  std::string nameStr = cutName;

  if (!nameStr.compare("analysis")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.995);
    cut->SetMaxPCA(0.5);
    cut->SetRxyRange(1, 90);
    return cut;
  }
  if (!nameStr.compare("qc")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_ITSTPC")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetRequireITSTPC(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_ITSonly")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSonly(true);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(4, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_TPConly")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireTPConly(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(36, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("wwire")) { // conversion only on tungstate wire
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetNClustersITS(2, 4);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(0.3);
    cut->SetOnWwireIB(true);
    cut->SetOnWwireOB(true);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("wwire_ib")) { // conversion only on tungstate wire outside of ITSib
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetNClustersITS(2, 4);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(0.3);
    cut->SetOnWwireIB(true);
    cut->SetOnWwireOB(false);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("wwire_ob")) { // conversion only on tungstate wire outside of ITSob (middle layer)
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetNClustersITS(0, 2);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 2);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    // cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(0.3);
    cut->SetOnWwireIB(false);
    cut->SetOnWwireOB(true);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }

  // for low B=0.2T
  if (!nameStr.compare("qc_lowB")) {
    // for track
    cut->SetTrackPtRange(0.02f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_ITSTPC_lowB")) {
    // for track
    cut->SetTrackPtRange(0.02f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSTPC(true);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_ITSonly_lowB")) {
    // for track
    cut->SetTrackPtRange(0.02f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSonly(true);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(4, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    // for v0
    cut->SetV0PtRange(0.05f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("qc_TPConly_lowB")) {
    // for track
    cut->SetTrackPtRange(0.02f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
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
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(36, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }

  if (!nameStr.compare("nopid")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    cut->RejectITSib(true);
    return cut;
  }
  if (!nameStr.compare("nocut")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(false);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.95);
    cut->SetMaxPCA(3.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    return cut;
  }
  if (!nameStr.compare("tag_track")) {
    // for track
    cut->SetTrackPtRange(0.04f, 1e10f);
    // cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(20);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetChi2PerClusterITS(-1e+10, 5.0);
    cut->SetNClustersITS(2, 4);
    cut->SetMeanClusterSizeITSob(0.0, 16.0);
    cut->SetIsWithinBeamPipe(true);
    cut->SetRequireITSTPC(true);
    // for v0
    cut->SetV0PtRange(0.1f, 1e10f);
    cut->SetV0EtaRange(-0.9, +0.9);
    cut->SetMinCosPA(0.99);
    cut->SetMaxPCA(1.0);
    cut->SetRxyRange(1, 90);
    cut->SetAPRange(0.95, 0.01);
    return cut;
  }

  delete cut;
  LOGF(info, Form("Did not find cut %s", cutName));
  return nullptr;
}

DalitzEECut* o2::aod::pwgem::photon::dalitzeecuts::GetCut(const char* cutName)
{
  DalitzEECut* cut = new DalitzEECut(cutName, cutName);
  std::string nameStr = cutName;
  // cut name should be like this mee0_120_minpt200_maxeta09_dca20_100_tpchadronbandrej_lowB in unit of MeV.

  if (!nameStr.compare("nocut")) {
    // apply kinetic cuts
    cut->SetTrackPtRange(0.1, 1e+10f);
    cut->SetTrackEtaRange(-0.9, +0.9);

    // for pair
    cut->SetMeeRange(0, 1e+10);
    cut->SetMaxPhivPairMeeDep([](float mee) { return (mee - -0.028) / 0.0185; });
    cut->ApplyPhiV(false);
    cut->ApplyPrefilter(false);

    // for track cuts
    cut->SetMinNCrossedRowsTPC(40);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);
    return cut;
  }

  if (!nameStr.compare("pc_itsib")) {
    // for pair
    cut->ApplyPhiV(true);
    cut->SelectPhotonConversion(true);
    cut->SetMeeRange(0.f, 0.03f);
    cut->SetMaxPhivPairMeeDep([](float mee) {
      return (mee - -0.028) / 0.0185;
    });

    // for track
    cut->SetTrackPtRange(0.05f, 1e10f);
    cut->SetTrackEtaRange(-0.9, +0.9);
    cut->SetMinNCrossedRowsTPC(40);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(2, 7);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);

    // for PID
    cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
    cut->SetTOFbetaRange(true, 0.0, 0.95);
    cut->SetTPCNsigmaElRange(-3, +3);
    cut->SetTPCNsigmaPiRange(0, 0);
    return cut;
  }

  float min_mass = 0.0;
  float max_mass = 1e+10;
  float min_pt = 0.05;
  float max_eta = 0.9;
  float min_dca3d_pair = 0.0;
  float max_dca3d_pair = 1e+10;
  std::vector<std::string> tmp = splitString(nameStr, '_');
  for (size_t i = 0; i < tmp.size(); i++) {
    // printf("string = %s , num = %d\n", tmp[i].data(), customAtoi(tmp[i]));
    if (tmp[i].find("mee") != std::string::npos || tmp[i].find("mmumu") != std::string::npos) {
      min_mass = static_cast<float>(customAtoi(tmp[i])) * 1e-3;     // convert MeV to GeV
      max_mass = static_cast<float>(customAtoi(tmp[i + 1])) * 1e-3; // convert MeV to GeV
    }
    if (tmp[i].find("minpt") != std::string::npos) {
      min_pt = static_cast<float>(customAtoi(tmp[i])) * 1e-3; // convert MeV to GeV
    }
    if (tmp[i].find("maxeta") != std::string::npos) {
      max_eta = static_cast<float>(customAtoi(tmp[i])) * 0.1;
    }
    if (tmp[i].find("dca") != std::string::npos) {
      min_dca3d_pair = static_cast<float>(customAtoi(tmp[i])) * 0.1;     // 3d dca in sigma
      max_dca3d_pair = static_cast<float>(customAtoi(tmp[i + 1])) * 0.1; // 3d dca in sigma
    }
  } // end of split string loop

  // apply kinetic cuts
  cut->SetTrackPtRange(min_pt, 1e+10f);
  cut->SetTrackEtaRange(-max_eta, +max_eta);

  // for pair
  cut->SetMeeRange(min_mass, max_mass);
  cut->SetMaxPhivPairMeeDep([](float mee) { return (mee - -0.028) / 0.0185; });
  cut->SetPairDCARange(min_dca3d_pair, max_dca3d_pair); // in sigma

  cut->ApplyPhiV(true);
  if (nameStr.find("wophiv") != std::string::npos) {
    cut->ApplyPhiV(false);
  } else if (nameStr.find("wphiv") != std::string::npos) {
    cut->ApplyPhiV(true);
  }

  cut->ApplyPrefilter(false);
  if (nameStr.find("wpf") != std::string::npos) {
    cut->ApplyPrefilter(true);
  } else if (nameStr.find("wopf") != std::string::npos) {
    cut->ApplyPrefilter(false);
  }

  // for track cuts
  cut->SetMinNCrossedRowsTPC(100);
  cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
  cut->SetChi2PerClusterTPC(0.0, 4.0);
  cut->SetChi2PerClusterITS(0.0, 5.0);
  cut->SetNClustersITS(5, 7);
  cut->SetMeanClusterSizeITSob(0, 16);
  cut->SetMaxDcaXY(1.0);
  cut->SetMaxDcaZ(1.0);

  if (nameStr.find("lowB") != std::string::npos) {
    cut->SetMinNCrossedRowsTPC(40);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(4, 7);
    cut->SetMeanClusterSizeITSob(0, 16);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);
  } else if (nameStr.find("nominalB") != std::string::npos) {
    cut->SetMinNCrossedRowsTPC(100);
    cut->SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    cut->SetChi2PerClusterTPC(0.0, 4.0);
    cut->SetChi2PerClusterITS(0.0, 5.0);
    cut->SetNClustersITS(5, 7);
    cut->SetMeanClusterSizeITSob(0, 16);
    cut->SetMaxDcaXY(1.0);
    cut->SetMaxDcaZ(1.0);
  }

  // long pid name should be top, because of if--else if conditions.

  if (nameStr.find("mee") != std::string::npos) { // for electron
    if (nameStr.find("tpcmupikaprrejortofreq") != std::string::npos) {
      // for PID
      cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
      cut->SetTOFbetaRange(true, 0.0, 0.95);
      cut->SetTPCNsigmaElRange(-2, +3);
      cut->SetTPCNsigmaPiRange(-3, +3);
      cut->SetTPCNsigmaKaRange(-3, +3);
      cut->SetTPCNsigmaPrRange(-3, +3);
      cut->SetTOFNsigmaElRange(-3, +3);
      cut->SetTPCNsigmaMuRange(-2, +2);
      cut->SetMuonExclusionTPC(true);
      return cut;
    } else if (nameStr.find("tpcpikaprrejortofreq") != std::string::npos) {
      // for PID
      cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPChadrejORTOFreq);
      cut->SetTOFbetaRange(true, 0.0, 0.95);
      cut->SetTPCNsigmaElRange(-2, +3);
      cut->SetTPCNsigmaPiRange(-3, +3);
      cut->SetTPCNsigmaKaRange(-3, +3);
      cut->SetTPCNsigmaPrRange(-3, +3);
      cut->SetTOFNsigmaElRange(-3, +3);
      return cut;
    } else if (nameStr.find("tpconly") != std::string::npos) {
      // for PID
      cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
      cut->SetTOFbetaRange(true, 0.0, 0.95);
      cut->SetTPCNsigmaElRange(-2, +3);
      cut->SetTPCNsigmaPiRange(-3, +3);
      return cut;
    } else if (nameStr.find("tpcelonly") != std::string::npos) {
      // for PID
      cut->SetPIDScheme(DalitzEECut::PIDSchemes::kTPConly);
      cut->SetTOFbetaRange(true, 0.0, 0.95);
      cut->SetTPCNsigmaElRange(-2, +3);
      cut->SetTPCNsigmaPiRange(0, 0);
      return cut;
    } else { // not match electron cut
      LOGF(info, Form("Did not find electron ID cut %s", cutName));
      return cut;
    }
  } else if (nameStr.find("mmumu") != std::string::npos) { // for muon
    if (nameStr.find("tpctof") != std::string::npos) {
      // for PID
      cut->SetPIDScheme(DalitzEECut::PIDSchemes::kMuon_lowB);
      cut->SetTPCNsigmaElRange(-2, +2); // exclusion
      cut->SetTPCNsigmaMuRange(-3, +3);
      cut->SetTPCNsigmaPiRange(-3, +1e+10);
      cut->SetTOFNsigmaMuRange(-3, +3);
      cut->SetTOFNsigmaPiRange(-3, +1e+10);
      return cut;
    } else { // not match muon cut
      LOGF(info, Form("Did not find muon ID cut %s", cutName));
      return cut;
    }
  } else { // match neither electron nor electron
    LOGF(info, Form("Did not find any pid cut %s", cutName));
    return cut;
  }
}

PHOSPhotonCut* o2::aod::pwgem::photon::phoscuts::GetCut(const char* cutName)
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

EMCPhotonCut* o2::aod::pwgem::photon::emccuts::GetCut(const char* cutName)
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

    cut->SetTrackMatchingEta([](float /*pT*/) {
      return -1.f;
    });
    cut->SetTrackMatchingPhi([](float /*pT*/) {
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

PairCut* o2::aod::pwgem::photon::paircuts::GetCut(const char* cutName)
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
