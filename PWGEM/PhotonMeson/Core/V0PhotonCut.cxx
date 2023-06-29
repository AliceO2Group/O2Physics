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

const char* V0PhotonCut::mCutNames[static_cast<int>(V0PhotonCut::V0PhotonCuts::kNCuts)] = {"Mee", "V0PtRange", "V0EtaRange", "PsiPair", "Rxy", "CosPA", "PCA", "RZLine", "OnWwireIB", "OnWwireOB", "TrackPtRange", "TrackEtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCNsigmaEl", "TPCNsigmaPi", "DCAxy", "DCAz", "IsWithinBeamPipe", "RequireITS", "TPConly", "AntiTPConly"};

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
  LOG(info) << "V0 Photon Cut, select photon on Tungstate wire IB: " << mIsOnWwireIB;
}
void V0PhotonCut::SetOnWwireOB(bool flag)
{
  mIsOnWwireOB = flag;
  LOG(info) << "V0 Photon Cut, select photon on Tungstate wire OB: " << mIsOnWwireOB;
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
void V0PhotonCut::SetMaxChi2PerClusterTPC(float maxChi2PerClusterTPC)
{
  mMaxChi2PerClusterTPC = maxChi2PerClusterTPC;
  LOG(info) << "V0 Photon Cut, set max chi2 per cluster TPC: " << mMaxChi2PerClusterTPC;
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

void V0PhotonCut::SetIsWithinBeamPipe(bool flag)
{
  mIsWithinBP = flag;
  LOG(info) << "V0 Photon Cut, propagated to within beam pipe: " << mIsWithinBP;
}

void V0PhotonCut::SetRequireITS(bool flag)
{
  mRequireITS = flag;
  LOG(info) << "V0 Photon Cut, require ITS hit: " << mRequireITS;
}

void V0PhotonCut::SetRequireTPConly(bool flag)
{
  mRequireTPConly = flag;
  LOG(info) << "V0 Photon Cut, require TPConly: " << mRequireTPConly;
}

void V0PhotonCut::SetRequireAntiTPConly(bool flag)
{
  mRequireAntiTPConly = flag;
  LOG(info) << "V0 Photon Cut, require AntiTPConly: " << mRequireAntiTPConly;
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
