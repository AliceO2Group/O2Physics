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
// Class for dimuon Cut
//

#include "PWGEM/Dilepton/Core/DimuonCut.h"

#include "Framework/Logger.h"

#include <vector>

ClassImp(DimuonCut);

void DimuonCut::SetMassRange(float min, float max)
{
  mMinMass = min;
  mMaxMass = max;
  LOG(info) << "Dimuon Cut, set mass range: " << mMinMass << " - " << mMaxMass;
}
void DimuonCut::SetPairPtRange(float minPt, float maxPt)
{
  mMinPairPt = minPt;
  mMaxPairPt = maxPt;
  LOG(info) << "Dimuon Cut, set pair pt range: " << mMinPairPt << " - " << mMaxPairPt;
}
void DimuonCut::SetPairYRange(float minY, float maxY)
{
  mMinPairY = minY;
  mMaxPairY = maxY;
  LOG(info) << "Dimuon Cut, set pair rapidity range: " << mMinPairY << " - " << mMaxPairY;
}
void DimuonCut::SetPairDCAxyRange(float min, float max)
{
  mMinPairDCAxy = min;
  mMaxPairDCAxy = max;
  LOG(info) << "Dimuon Cut, set pair 3d dca range: " << mMinPairDCAxy << " - " << mMaxPairDCAxy;
}
void DimuonCut::SetMindEtadPhi(bool flag, float min_deta, float min_dphi)
{
  mApplydEtadPhi = flag;
  mMinDeltaEta = min_deta;
  mMinDeltaPhi = min_dphi;
  LOG(info) << "Dimuon Cut, set apply deta-dphi cut: " << mApplydEtadPhi << " min_deta: " << mMinDeltaEta << " min_dphi: " << mMinDeltaPhi;
}
void DimuonCut::SetTrackType(int track_type)
{
  mTrackType = track_type;
  LOG(info) << "Dimuon Cut, set track type: " << mTrackType;
}
void DimuonCut::SetTrackPtRange(float minPt, float maxPt)
{
  mMinTrackPt = minPt;
  mMaxTrackPt = maxPt;
  LOG(info) << "Dimuon Cut, set track pt range: " << mMinTrackPt << " - " << mMaxTrackPt;
}
void DimuonCut::SetTrackEtaRange(float minEta, float maxEta)
{
  mMinTrackEta = minEta;
  mMaxTrackEta = maxEta;
  LOG(info) << "Dimuon Cut, set track eta range: " << mMinTrackEta << " - " << mMaxTrackEta;
}
void DimuonCut::SetTrackPhiRange(float minPhi, float maxPhi)
{
  mMinTrackPhi = minPhi;
  mMaxTrackPhi = maxPhi;
  LOG(info) << "Dimuon Cut, set track phi range (rad.): " << mMinTrackPhi << " - " << mMaxTrackPhi;
}
void DimuonCut::SetChi2(float min, float max)
{
  mMinChi2 = min;
  mMaxChi2 = max;
  LOG(info) << "Dimuon Cut, set chi2 range: " << mMinChi2 << " - " << mMaxChi2;
}
void DimuonCut::SetChi2MFT(float min, float max)
{
  mMinChi2MFT = min;
  mMaxChi2MFT = max;
  LOG(info) << "Dimuon Cut, set chi2mft range: " << mMinChi2MFT << " - " << mMaxChi2MFT;
}
void DimuonCut::SetMatchingChi2MCHMFT(float min, float max)
{
  mMinMatchingChi2MCHMFT = min;
  mMaxMatchingChi2MCHMFT = max;
  LOG(info) << "Dimuon Cut, set matching chi2 MFT-MCH range: " << mMinMatchingChi2MCHMFT << " - " << mMaxMatchingChi2MCHMFT;
}
void DimuonCut::SetMatchingChi2MCHMID(float min, float max)
{
  mMinMatchingChi2MCHMID = min;
  mMaxMatchingChi2MCHMID = max;
  LOG(info) << "Dimuon Cut, set matching chi2 MCH-MID range: " << mMinMatchingChi2MCHMID << " - " << mMaxMatchingChi2MCHMID;
}
void DimuonCut::SetNClustersMFT(int min, int max)
{
  mMinNClustersMFT = min;
  mMaxNClustersMFT = max;
  LOG(info) << "Dimuon Cut, set N clusters MFT range: " << mMinNClustersMFT << " - " << mMaxNClustersMFT;
}
void DimuonCut::SetNClustersMCHMID(int min, int max)
{
  mMinNClustersMCHMID = min;
  mMaxNClustersMCHMID = max;
  LOG(info) << "Dimuon Cut, set N clusters MCHMID range: " << mMinNClustersMCHMID << " - " << mMaxNClustersMCHMID;
}
void DimuonCut::SetRabs(float min, float max)
{
  mMinRabs = min;
  mMaxRabs = max;
  LOG(info) << "Dimuon Cut, set Rabs range: " << mMinRabs << " - " << mMaxRabs;
}
void DimuonCut::SetDCAxy(float min, float max)
{
  mMinDcaXY = min;
  mMaxDcaXY = max;
  LOG(info) << "Dimuon Cut, set DCAxy range: " << mMinDcaXY << " - " << mMaxDcaXY;
}
void DimuonCut::SetMaxPDCARabsDep(std::function<float(float)> RabsDepCut)
{
  mMaxPDCARabsDep = RabsDepCut;
  LOG(info) << "Dimuon Cut, set max pDCA as a function of Rabs: " << mMaxPDCARabsDep(10.0);
}
void DimuonCut::SetMFTHitMap(bool flag, std::vector<int> hitMap)
{
  mApplyMFTHitMap = flag;
  mRequiredMFTDisks = hitMap;
  if (mApplyMFTHitMap) {
    for (const auto& iDisk : mRequiredMFTDisks) {
      LOG(info) << "Dimuon Cut, require MFT hit on Disk: " << iDisk;
    }
  }
}
void DimuonCut::SetMaxdPtdEtadPhiwrtMCHMID(float reldPtMax, float dEtaMax, float dPhiMax)
{
  mMaxReldPtwrtMCHMID = reldPtMax;
  mMaxdEtawrtMCHMID = dEtaMax;
  mMaxdPhiwrtMCHMID = dPhiMax;
  LOG(info) << "Dimuon Cut, set max rel. dpt between MFT-MCH-MID and associated MCH-MID: " << mMaxReldPtwrtMCHMID;
  LOG(info) << "Dimuon Cut, set max deta between MFT-MCH-MID and associated MCH-MID: " << mMaxdEtawrtMCHMID;
  LOG(info) << "Dimuon Cut, set max dphi between MFT-MCH-MID and associated MCH-MID: " << mMaxdPhiwrtMCHMID;
}
