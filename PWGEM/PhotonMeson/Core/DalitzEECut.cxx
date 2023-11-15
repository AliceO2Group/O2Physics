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
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"

ClassImp(DalitzEECut);

const char* DalitzEECut::mCutNames[static_cast<int>(DalitzEECut::DalitzEECuts::kNCuts)] = {"Mee", "PairPtRange", "PairEtaRange", "PhivPair", "TrackPtRange", "TrackEtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCNsigmaEl", "TPCNsigmaMu", "TPCNsigmaPi", "TPCNsigmaKa", "TPCNsigmaPr", "TOFNsigmaEl", "TOFNsigmaMu", "TOFNsigmaPi", "TOFNsigmaKa", "TOFNsigmaPr", "DCAxy", "DCAz", "ITSNCls", "ITSChi2NDF"};

void DalitzEECut::SetPairPtRange(float minPt, float maxPt)
{
  mMinPairPt = minPt;
  mMaxPairPt = maxPt;
  LOG(info) << "DalitzEE Cut, set pair pt range: " << mMinPairPt << " - " << mMaxPairPt;
}
void DalitzEECut::SetPairEtaRange(float minEta, float maxEta)
{
  mMinPairEta = minEta;
  mMaxPairEta = maxEta;
  LOG(info) << "DalitzEE Cut, set pair eta range: " << mMinPairEta << " - " << mMaxPairEta;
}
void DalitzEECut::SetMeeRange(float min, float max)
{
  mMinMee = min;
  mMaxMee = max;
  LOG(info) << "DalitzEE selection, set mee range: " << mMinMee << " - " << mMaxMee;
}
void DalitzEECut::SetMaxPhivPairMeeDep(std::function<float(float)> meeDepCut)
{
  mMaxPhivPairMeeDep = meeDepCut;
  LOG(info) << "DalitzEE Cut, set max phiv pair mee dep: " << mMaxPhivPairMeeDep(0.02);
}
void DalitzEECut::SetTrackPtRange(float minPt, float maxPt)
{
  mMinTrackPt = minPt;
  mMaxTrackPt = maxPt;
  LOG(info) << "DalitzEE Cut, set track pt range: " << mMinTrackPt << " - " << mMaxTrackPt;
}
void DalitzEECut::SetTrackEtaRange(float minEta, float maxEta)
{
  mMinTrackEta = minEta;
  mMaxTrackEta = maxEta;
  LOG(info) << "DalitzEE Cut, set track eta range: " << mMinTrackEta << " - " << mMaxTrackEta;
}
void DalitzEECut::SetMinNClustersTPC(int minNClustersTPC)
{
  mMinNClustersTPC = minNClustersTPC;
  LOG(info) << "DalitzEE Cut, set min N clusters TPC: " << mMinNClustersTPC;
}
void DalitzEECut::SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
{
  mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  LOG(info) << "DalitzEE Cut, set min N crossed rows TPC: " << mMinNCrossedRowsTPC;
}
void DalitzEECut::SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
{
  mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  LOG(info) << "DalitzEE Cut, set min N crossed rows over findable clusters TPC: " << mMinNCrossedRowsOverFindableClustersTPC;
}
void DalitzEECut::SetChi2PerClusterTPC(float min, float max)
{
  mMinChi2PerClusterTPC = min;
  mMaxChi2PerClusterTPC = max;
  LOG(info) << "DalitzEE Cut, set chi2 per cluster TPC range: " << mMinChi2PerClusterTPC << " - " << mMaxChi2PerClusterTPC;
}

void DalitzEECut::SetNClustersITS(int min, int max)
{
  mMinNClustersITS = min;
  mMaxNClustersITS = max;
  LOG(info) << "DalitzEE Cut, set N clusters ITS range: " << mMinNClustersITS << " - " << mMaxNClustersITS;
}
void DalitzEECut::SetChi2PerClusterITS(float min, float max)
{
  mMinChi2PerClusterITS = min;
  mMaxChi2PerClusterITS = max;
  LOG(info) << "DalitzEE Cut, set chi2 per cluster ITS range: " << mMinChi2PerClusterITS << " - " << mMaxChi2PerClusterITS;
}

void DalitzEECut::SetMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "DalitzEE Cut, set max DCA xy: " << mMaxDcaXY;
}
void DalitzEECut::SetMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "DalitzEE Cut, set max DCA z: " << mMaxDcaZ;
}

void DalitzEECut::SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
{
  mMaxDcaXYPtDep = ptDepCut;
  LOG(info) << "DalitzEE Cut, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
}

void DalitzEECut::SetPIDScheme(PIDSchemes scheme)
{
  mPIDScheme = scheme;
  LOG(info) << "DalitzEE Cut, PID scheme: " << static_cast<int>(mPIDScheme);
}
void DalitzEECut::SetMinPinTOF(float min)
{
  mMinPinTOF = min;
  LOG(info) << "DalitzEE Cut, set min pin for TOF: " << mMinPinTOF;
}

void DalitzEECut::SetMuonExclusionTPC(bool flag)
{
  mMuonExclusionTPC = flag;
  LOG(info) << "DalitzEE Cut, set flag for muon exclusion in TPC: " << mMuonExclusionTPC;
}

void DalitzEECut::SetTOFbetaRange(bool flag, float min, float max)
{
  mApplyTOFbeta = flag;
  mMinTOFbeta = min;
  mMaxTOFbeta = max;
  LOG(info) << "DalitzEE selection, set TOF beta rejection range: " << mMinTOFbeta << " - " << mMaxTOFbeta;
}

void DalitzEECut::SetTPCNsigmaElRange(float min, float max)
{
  mMinTPCNsigmaEl = min;
  mMaxTPCNsigmaEl = max;
  LOG(info) << "DalitzEE selection, set TPC n sigma El range: " << mMinTPCNsigmaEl << " - " << mMaxTPCNsigmaEl;
}
void DalitzEECut::SetTPCNsigmaMuRange(float min, float max)
{
  mMinTPCNsigmaMu = min;
  mMaxTPCNsigmaMu = max;
  LOG(info) << "DalitzEE selection, set TPC n sigma Mu range: " << mMinTPCNsigmaMu << " - " << mMaxTPCNsigmaMu;
}
void DalitzEECut::SetTPCNsigmaPiRange(float min, float max)
{
  mMinTPCNsigmaPi = min;
  mMaxTPCNsigmaPi = max;
  LOG(info) << "DalitzEE selection, set TPC n sigma Pi range: " << mMinTPCNsigmaPi << " - " << mMaxTPCNsigmaPi;
}
void DalitzEECut::SetTPCNsigmaKaRange(float min, float max)
{
  mMinTPCNsigmaKa = min;
  mMaxTPCNsigmaKa = max;
  LOG(info) << "DalitzEE selection, set TPC n sigma Ka range: " << mMinTPCNsigmaKa << " - " << mMaxTPCNsigmaKa;
}
void DalitzEECut::SetTPCNsigmaPrRange(float min, float max)
{
  mMinTPCNsigmaPr = min;
  mMaxTPCNsigmaPr = max;
  LOG(info) << "DalitzEE selection, set TPC n sigma Pr range: " << mMinTPCNsigmaPr << " - " << mMaxTPCNsigmaPr;
}

void DalitzEECut::SetTOFNsigmaElRange(float min, float max)
{
  mMinTOFNsigmaEl = min;
  mMaxTOFNsigmaEl = max;
  LOG(info) << "DalitzEE selection, set TOF n sigma El range: " << mMinTOFNsigmaEl << " - " << mMaxTOFNsigmaEl;
}
void DalitzEECut::SetTOFNsigmaMuRange(float min, float max)
{
  mMinTOFNsigmaMu = min;
  mMaxTOFNsigmaMu = max;
  LOG(info) << "DalitzEE selection, set TOF n sigma Mu range: " << mMinTOFNsigmaMu << " - " << mMaxTOFNsigmaMu;
}
void DalitzEECut::SetTOFNsigmaPiRange(float min, float max)
{
  mMinTOFNsigmaPi = min;
  mMaxTOFNsigmaPi = max;
  LOG(info) << "DalitzEE selection, set TOF n sigma Pi range: " << mMinTOFNsigmaPi << " - " << mMaxTOFNsigmaPi;
}
void DalitzEECut::SetTOFNsigmaKaRange(float min, float max)
{
  mMinTOFNsigmaKa = min;
  mMaxTOFNsigmaKa = max;
  LOG(info) << "DalitzEE selection, set TOF n sigma Ka range: " << mMinTOFNsigmaKa << " - " << mMaxTOFNsigmaKa;
}
void DalitzEECut::SetTOFNsigmaPrRange(float min, float max)
{
  mMinTOFNsigmaPr = min;
  mMaxTOFNsigmaPr = max;
  LOG(info) << "DalitzEE selection, set TOF n sigma Pr range: " << mMinTOFNsigmaPr << " - " << mMaxTOFNsigmaPr;
}

void DalitzEECut::print() const
{
  LOG(info) << "Dalitz EE Cut:";
  for (int i = 0; i < static_cast<int>(DalitzEECuts::kNCuts); i++) {
    switch (static_cast<DalitzEECuts>(i)) {
      case DalitzEECuts::kTrackPtRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackPt << ", " << mMaxTrackPt << "]";
        break;
      case DalitzEECuts::kTrackEtaRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackEta << ", " << mMaxTrackEta << "]";
        break;
      case DalitzEECuts::kTPCNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNClustersTPC;
        break;
      case DalitzEECuts::kTPCCrossedRows:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsTPC;
        break;
      case DalitzEECuts::kTPCCrossedRowsOverNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsOverFindableClustersTPC;
        break;
      case DalitzEECuts::kTPCChi2NDF:
        LOG(info) << mCutNames[i] << " < " << mMaxChi2PerClusterTPC;
        break;
      case DalitzEECuts::kDCAxy:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaXY;
        break;
      case DalitzEECuts::kDCAz:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaZ;
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
