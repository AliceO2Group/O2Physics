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
// Class for track selection
//

#include "Framework/Logger.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"

ClassImp(V0PhotonCut);

const char* V0PhotonCut::mCutNames[static_cast<int>(V0PhotonCut::V0PhotonCuts::kNCuts)] = {"Mee", "V0PtRange", "V0EtaRange", "AP", " PsiPair", "PhivPair", "Rxy", "CosPA", "PCA", "RZLine", "OnWwireIB", "OnWwireOB", "TrackPtRange", "TrackEtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCNsigmaEl", "TPCNsigmaPi", "DCAxy", "DCAz", "ITSNCls", "ITSChi2NDF", "IsWithinBeamPipe", "RequireITSTPC", "RequireITSonly", "RequireTPConly", "RequireTPCTRD", "RequireTPCTOF", "RequireTPCTRDTOF"};

const std::pair<int8_t, std::set<uint8_t>> V0PhotonCut::its_ib_Requirement = {0, {0, 1, 2}};           // no hit on 3 ITS ib layers.
const std::pair<int8_t, std::set<uint8_t>> V0PhotonCut::its_ob_Requirement = {4, {3, 4, 5, 6}};        // all hits on 4 ITS ob layers.
const std::pair<int8_t, std::set<uint8_t>> V0PhotonCut::its_ob_Requirement_ITSTPC = {2, {3, 4, 5, 6}}; // at least 2 hits on 4 ITS ob layers.

void V0PhotonCut::SetV0PtRange(float minPt, float maxPt)
{
  mMinV0Pt = minPt;
  mMaxV0Pt = maxPt;
  LOG(info) << "V0 Photon Cut, set v0 pt range: " << mMinV0Pt << " - " << mMaxV0Pt;
}
void V0PhotonCut::SetV0EtaRange(float minEta, float maxEta)
{
  mMinV0Eta = minEta;
  mMaxV0Eta = maxEta;
  LOG(info) << "V0 Photon Cut, set v0 eta range: " << mMinV0Eta << " - " << mMaxV0Eta;
}
void V0PhotonCut::SetMeeRange(float min, float max)
{
  mMinMee = min;
  mMaxMee = max;
  LOG(info) << "V0 Photon selection, set mee range: " << mMinMee << " - " << mMaxMee;
}
void V0PhotonCut::SetPsiPairRange(float min, float max)
{
  mMinPsiPair = min;
  mMaxPsiPair = max;
  LOG(info) << "V0 Photon selection, set psi pair range: " << mMinPsiPair << " - " << mMaxPsiPair;
}
void V0PhotonCut::SetPhivPairRange(float min, float max)
{
  mMinPhivPair = min;
  mMaxPhivPair = max;
  LOG(info) << "V0 Photon selection, set phiv pair range: " << mMinPhivPair << " - " << mMaxPhivPair;
}
void V0PhotonCut::SetAPRange(float max_alpha, float max_qt)
{
  mMaxAlpha = max_alpha;
  mMaxQt = max_qt;
  LOG(info) << "V0 Photon selection, set Armenteroz-Podolanski range: " << mMaxAlpha << " - " << mMaxQt;
}
void V0PhotonCut::SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut)
{
  mMaxMeePsiPairDep = psiDepCut;
  LOG(info) << "V0 Photon Cut, set max mee psi pair dep: " << mMaxMeePsiPairDep(0.1);
}
void V0PhotonCut::SetRxyRange(float min, float max)
{
  mMinRxy = min;
  mMaxRxy = max;
  LOG(info) << "V0 Photon selection, set Rxy range: " << mMinRxy << " - " << mMaxRxy;
}
void V0PhotonCut::SetMinCosPA(float min)
{
  mMinCosPA = min;
  LOG(info) << "V0 Photon Cut, set min cosine pointing angle: " << mMinCosPA;
}
void V0PhotonCut::SetMaxPCA(float max)
{
  mMaxPCA = max;
  LOG(info) << "V0 Photon Cut, set max distance between 2 legs: " << mMaxPCA;
}
void V0PhotonCut::SetMaxMarginZ(float max)
{
  mMaxMarginZ = max;
  LOG(info) << "V0 Photon Cut, set max margin z: " << mMaxMarginZ;
}
void V0PhotonCut::SetOnWwireIB(bool flag)
{
  mIsOnWwireIB = flag;
  LOG(info) << "V0 Photon Cut, select photon on Tungsten wire IB: " << mIsOnWwireIB;
}
void V0PhotonCut::SetOnWwireOB(bool flag)
{
  mIsOnWwireOB = flag;
  LOG(info) << "V0 Photon Cut, select photon on Tungsten wire OB: " << mIsOnWwireOB;
}
void V0PhotonCut::SetTPCNsigmaElRange(float min, float max)
{
  mMinTPCNsigmaEl = min;
  mMaxTPCNsigmaEl = max;
  LOG(info) << "V0 Photon selection, set TPC n sigma El range: " << mMinTPCNsigmaEl << " - " << mMaxTPCNsigmaEl;
}
void V0PhotonCut::SetTPCNsigmaPiRange(float min, float max)
{
  mMinTPCNsigmaPi = min;
  mMaxTPCNsigmaPi = max;
  LOG(info) << "V0 Photon selection, set TPC n sigma Pi range: " << mMinTPCNsigmaPi << " - " << mMaxTPCNsigmaPi;
}
void V0PhotonCut::SetTrackPtRange(float minPt, float maxPt)
{
  mMinTrackPt = minPt;
  mMaxTrackPt = maxPt;
  LOG(info) << "V0 Photon Cut, set track pt range: " << mMinTrackPt << " - " << mMaxTrackPt;
}
void V0PhotonCut::SetTrackEtaRange(float minEta, float maxEta)
{
  mMinTrackEta = minEta;
  mMaxTrackEta = maxEta;
  LOG(info) << "V0 Photon Cut, set track eta range: " << mMinTrackEta << " - " << mMaxTrackEta;
}
void V0PhotonCut::SetMinNClustersTPC(int minNClustersTPC)
{
  mMinNClustersTPC = minNClustersTPC;
  LOG(info) << "V0 Photon Cut, set min N clusters TPC: " << mMinNClustersTPC;
}
void V0PhotonCut::SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
{
  mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  LOG(info) << "V0 Photon Cut, set min N crossed rows TPC: " << mMinNCrossedRowsTPC;
}
void V0PhotonCut::SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
{
  mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  LOG(info) << "V0 Photon Cut, set min N crossed rows over findable clusters TPC: " << mMinNCrossedRowsOverFindableClustersTPC;
}
void V0PhotonCut::SetChi2PerClusterTPC(float min, float max)
{
  mMinChi2PerClusterTPC = min;
  mMaxChi2PerClusterTPC = max;
  LOG(info) << "V0 Photon Cut, set chi2 per cluster TPC range: " << mMinChi2PerClusterTPC << " - " << mMaxChi2PerClusterTPC;
}
void V0PhotonCut::SetMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "V0 Photon Cut, set max DCA xy: " << mMaxDcaXY;
}
void V0PhotonCut::SetMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "V0 Photon Cut, set max DCA z: " << mMaxDcaZ;
}

void V0PhotonCut::SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
{
  mMaxDcaXYPtDep = ptDepCut;
  LOG(info) << "V0 Photon Cut, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
}

void V0PhotonCut::SetNClustersITS(int min, int max)
{
  mMinNClustersITS = min;
  mMaxNClustersITS = max;
  LOG(info) << "V0 Photon Cut, set N clusters ITS range: " << mMinNClustersITS << " - " << mMaxNClustersITS;
}
void V0PhotonCut::SetChi2PerClusterITS(float min, float max)
{
  mMinChi2PerClusterITS = min;
  mMaxChi2PerClusterITS = max;
  LOG(info) << "V0 Photon Cut, set chi2 per cluster ITS range: " << mMinChi2PerClusterITS << " - " << mMaxChi2PerClusterITS;
}

void V0PhotonCut::SetIsWithinBeamPipe(bool flag)
{
  mIsWithinBP = flag;
  LOG(info) << "V0 Photon Cut, propagated to within beam pipe: " << mIsWithinBP;
}

void V0PhotonCut::SetRequireITSTPC(bool flag)
{
  mRequireITSTPC = flag;
  LOG(info) << "V0 Photon Cut, require ITS-TPC matched track: " << mRequireITSTPC;
}

void V0PhotonCut::SetRequireITSonly(bool flag)
{
  mRequireITSonly = flag;
  LOG(info) << "V0 Photon Cut, require ITS only track: " << mRequireITSonly;
}

void V0PhotonCut::SetRequireTPConly(bool flag)
{
  mRequireTPConly = flag;
  LOG(info) << "V0 Photon Cut, require TPConly track: " << mRequireTPConly;
}

void V0PhotonCut::SetRequireTPCTRD(bool flag)
{
  mRequireTPCTRD = flag;
  LOG(info) << "V0 Photon Cut, require TPC-TRD track: " << mRequireTPCTRD;
}

void V0PhotonCut::SetRequireTPCTOF(bool flag)
{
  mRequireTPCTOF = flag;
  LOG(info) << "V0 Photon Cut, require TPC-TOF track: " << mRequireTPCTOF;
}

void V0PhotonCut::SetRequireTPCTRDTOF(bool flag)
{
  mRequireTPCTRDTOF = flag;
  LOG(info) << "V0 Photon Cut, require TPC-TOF track: " << mRequireTPCTRDTOF;
}

void V0PhotonCut::print() const
{
  LOG(info) << "V0 Photon Cut:";
  for (int i = 0; i < static_cast<int>(V0PhotonCuts::kNCuts); i++) {
    switch (static_cast<V0PhotonCuts>(i)) {
      case V0PhotonCuts::kTrackPtRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackPt << ", " << mMaxTrackPt << "]";
        break;
      case V0PhotonCuts::kTrackEtaRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackEta << ", " << mMaxTrackEta << "]";
        break;
      case V0PhotonCuts::kTPCNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNClustersTPC;
        break;
      case V0PhotonCuts::kTPCCrossedRows:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsTPC;
        break;
      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsOverFindableClustersTPC;
        break;
      case V0PhotonCuts::kTPCChi2NDF:
        LOG(info) << mCutNames[i] << " < " << mMaxChi2PerClusterTPC;
        break;
      case V0PhotonCuts::kDCAxy:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaXY;
        break;
      case V0PhotonCuts::kDCAz:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaZ;
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
