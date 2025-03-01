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
// Class for dilepton Cut
//

#include <utility>
#include <set>

#include "Framework/Logger.h"
#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"

ClassImp(DalitzEECut);

// const char* DalitzEECut::mCutNames[static_cast<int>(DalitzEECut::DalitzEECuts::kNCuts)] = {"Mee", "PairPtRange", "PairRapidityRange", "PairDCARange", "PhivPair", "TrackPtRange", "TrackEtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCNsigmaEl", "TPCNsigmaMu", "TPCNsigmaPi", "TPCNsigmaKa", "TPCNsigmaPr", "TOFNsigmaEl", "TOFNsigmaMu", "TOFNsigmaPi", "TOFNsigmaKa", "TOFNsigmaPr", "DCA3Dsigma", "DCAxy", "DCAz", "ITSNCls", "ITSChi2NDF", "ITSClusterSize", "Prefilter"};

const std::pair<int8_t, std::set<uint8_t>> DalitzEECut::its_ib_any_Requirement = {1, {0, 1, 2}}; // hits on any ITS ib layers.
const std::pair<int8_t, std::set<uint8_t>> DalitzEECut::its_ib_1st_Requirement = {1, {0}};       // hit on 1st ITS ib layers.

void DalitzEECut::SetPairPtRange(float minPt, float maxPt)
{
  mMinPairPt = minPt;
  mMaxPairPt = maxPt;
  LOG(info) << "DalitzEE Cut, set pair pt range: " << mMinPairPt << " - " << mMaxPairPt;
}
void DalitzEECut::SetPairYRange(float minY, float maxY)
{
  mMinPairY = minY;
  mMaxPairY = maxY;
  LOG(info) << "DalitzEE Cut, set pair eta range: " << mMinPairY << " - " << mMaxPairY;
}
void DalitzEECut::SetMeeRange(float min, float max)
{
  mMinMee = min;
  mMaxMee = max;
  LOG(info) << "DalitzEE Cut, set mee range: " << mMinMee << " - " << mMaxMee;
}
void DalitzEECut::SetMaxPhivPairMeeDep(std::function<float(float)> meeDepCut)
{
  mMaxPhivPairMeeDep = meeDepCut;
  LOG(info) << "DalitzEE Cut, set max phiv pair mee dep: " << mMaxPhivPairMeeDep(0.02);
}
void DalitzEECut::SelectPhotonConversion(bool flag)
{
  mSelectPC = flag;
  LOG(info) << "DalitzEE Cut, select photon conversion: " << mSelectPC;
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
void DalitzEECut::SetMeanClusterSizeITS(float min, float max)
{
  mMinMeanClusterSizeITS = min;
  mMaxMeanClusterSizeITS = max;
  LOG(info) << "DalitzEE Cut, set mean cluster size ITS range: " << mMinMeanClusterSizeITS << " - " << mMaxMeanClusterSizeITS;
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
void DalitzEECut::ApplyPhiV(bool flag)
{
  mApplyPhiV = flag;
  LOG(info) << "DalitzEE Cut, apply phiv cut: " << mApplyPhiV;
}

void DalitzEECut::SetTPCNsigmaElRange(float min, float max)
{
  mMinTPCNsigmaEl = min;
  mMaxTPCNsigmaEl = max;
  LOG(info) << "DalitzEE Cut, set TPC n sigma El range: " << mMinTPCNsigmaEl << " - " << mMaxTPCNsigmaEl;
}
void DalitzEECut::SetTPCNsigmaPiRange(float min, float max)
{
  mMinTPCNsigmaPi = min;
  mMaxTPCNsigmaPi = max;
  LOG(info) << "DalitzEE Cut, set TPC n sigma Pi range: " << mMinTPCNsigmaPi << " - " << mMaxTPCNsigmaPi;
}

void DalitzEECut::SetTOFNsigmaElRange(float min, float max)
{
  mMinTOFNsigmaEl = min;
  mMaxTOFNsigmaEl = max;
  LOG(info) << "DalitzEE Cut, set TOF n sigma El range: " << mMinTOFNsigmaEl << " - " << mMaxTOFNsigmaEl;
}
void DalitzEECut::RequireITSibAny(bool flag)
{
  mRequireITSibAny = flag;
  LOG(info) << "DalitzEE Cut, require ITS ib any: " << mRequireITSibAny;
}
void DalitzEECut::RequireITSib1st(bool flag)
{
  mRequireITSib1st = flag;
  LOG(info) << "DalitzEE Cut, require ITS ib 1st: " << mRequireITSib1st;
}
void DalitzEECut::SetChi2TOF(float min, float max)
{
  mMinChi2TOF = min;
  mMaxChi2TOF = max;
  LOG(info) << "Dielectron Cut, set chi2 TOF range: " << mMinChi2TOF << " - " << mMaxChi2TOF;
}
void DalitzEECut::SetPIDScheme(int scheme)
{
  mPIDScheme = scheme;
  LOG(info) << "DalitzEE Cut, PID scheme: " << static_cast<int>(mPIDScheme);
}
