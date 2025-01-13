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
// Class for dielectron Cut
//

#include <utility>
#include <set>

#include "Framework/Logger.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"

ClassImp(DielectronCut);

const std::pair<int8_t, std::set<uint8_t>> DielectronCut::its_ib_any_Requirement = {1, {0, 1, 2}}; // hits on any ITS ib layers.
const std::pair<int8_t, std::set<uint8_t>> DielectronCut::its_ib_1st_Requirement = {1, {0}};       // hit on 1st ITS ib layers.

void DielectronCut::SetPairPtRange(float minPt, float maxPt)
{
  mMinPairPt = minPt;
  mMaxPairPt = maxPt;
  LOG(info) << "Dielectron Cut, set pair pt range: " << mMinPairPt << " - " << mMaxPairPt;
}
void DielectronCut::SetPairYRange(float minY, float maxY)
{
  mMinPairY = minY;
  mMaxPairY = maxY;
  LOG(info) << "Dielectron Cut, set pair y range: " << mMinPairY << " - " << mMaxPairY;
}
void DielectronCut::SetPairDCARange(float min, float max)
{
  mMinPairDCA3D = min;
  mMaxPairDCA3D = max;
  LOG(info) << "Dielectron Cut, set pair 3d dca range: " << mMinPairDCA3D << " - " << mMaxPairDCA3D;
}
void DielectronCut::SetMeeRange(float min, float max)
{
  mMinMee = min;
  mMaxMee = max;
  LOG(info) << "Dielectron Cut, set mee range: " << mMinMee << " - " << mMaxMee;
}
void DielectronCut::SetPairOpAng(float minOpAng, float maxOpAng)
{
  mMinOpAng = minOpAng;
  mMaxOpAng = maxOpAng;
  LOG(info) << "Dielectron Cut, set pair opening angle range: " << mMinOpAng << " - " << mMaxOpAng;
}
void DielectronCut::SetMaxMeePhiVDep(std::function<float(float)> phivDepCut, float min_phiv, float max_phiv)
{
  mMaxMeePhiVDep = phivDepCut;
  mMinPhivPair = min_phiv;
  mMaxPhivPair = max_phiv;
  LOG(info) << "Dielectron Cut, set max mee phiv dep: " << mMaxMeePhiVDep(2.5);
}
void DielectronCut::SelectPhotonConversion(bool flag)
{
  mSelectPC = flag;
  LOG(info) << "Dielectron Cut, select photon conversion: " << mSelectPC;
}
void DielectronCut::SetMindEtadPhi(bool flag, float min_deta, float min_dphi)
{
  mApplydEtadPhi = flag;
  mMinDeltaEta = min_deta;
  mMinDeltaPhi = min_dphi;
  LOG(info) << "Dielectron Cut, set apply deta-dphi cut: " << mApplydEtadPhi << " min_deta: " << mMinDeltaEta << " min_dphi: " << mMinDeltaPhi;
}
void DielectronCut::SetRequireDifferentSides(bool flag)
{
  mRequireDiffSides = flag;
  LOG(info) << "Dielectron Cut, require 2 tracks to be from different sides: " << mRequireDiffSides;
}
void DielectronCut::SetTrackPtRange(float minPt, float maxPt)
{
  mMinTrackPt = minPt;
  mMaxTrackPt = maxPt;
  LOG(info) << "Dielectron Cut, set track pt range: " << mMinTrackPt << " - " << mMaxTrackPt;
}
void DielectronCut::SetTrackEtaRange(float minEta, float maxEta)
{
  mMinTrackEta = minEta;
  mMaxTrackEta = maxEta;
  LOG(info) << "Dielectron Cut, set track eta range: " << mMinTrackEta << " - " << mMaxTrackEta;
}
void DielectronCut::SetTrackPhiRange(float minPhi, float maxPhi)
{
  mMinTrackPhi = minPhi;
  mMaxTrackPhi = maxPhi;
  LOG(info) << "Dielectron Cut, set track phi range (rad.): " << mMinTrackPhi << " - " << mMaxTrackPhi;
}
void DielectronCut::SetMinNClustersTPC(int minNClustersTPC)
{
  mMinNClustersTPC = minNClustersTPC;
  LOG(info) << "Dielectron Cut, set min N clusters TPC: " << mMinNClustersTPC;
}
void DielectronCut::SetMinNCrossedRowsTPC(int minNCrossedRowsTPC)
{
  mMinNCrossedRowsTPC = minNCrossedRowsTPC;
  LOG(info) << "Dielectron Cut, set min N crossed rows TPC: " << mMinNCrossedRowsTPC;
}
void DielectronCut::SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC)
{
  mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC;
  LOG(info) << "Dielectron Cut, set min N crossed rows over findable clusters TPC: " << mMinNCrossedRowsOverFindableClustersTPC;
}
void DielectronCut::SetMaxFracSharedClustersTPC(float max)
{
  mMaxFracSharedClustersTPC = max;
  LOG(info) << "Dielectron Cut, set max fraction of shared clusters in  TPC: " << mMaxFracSharedClustersTPC;
}
void DielectronCut::SetRelDiffPin(float min, float max)
{
  mMinRelDiffPin = min;
  mMaxRelDiffPin = max;
  LOG(info) << "Dielectron Cut, set rel. diff. between Pin and Ppv range: " << mMinRelDiffPin << " - " << mMaxRelDiffPin;
}
void DielectronCut::SetChi2PerClusterTPC(float min, float max)
{
  mMinChi2PerClusterTPC = min;
  mMaxChi2PerClusterTPC = max;
  LOG(info) << "Dielectron Cut, set chi2 per cluster TPC range: " << mMinChi2PerClusterTPC << " - " << mMaxChi2PerClusterTPC;
}

void DielectronCut::SetNClustersITS(int min, int max)
{
  mMinNClustersITS = min;
  mMaxNClustersITS = max;
  LOG(info) << "Dielectron Cut, set N clusters ITS range: " << mMinNClustersITS << " - " << mMaxNClustersITS;
}
void DielectronCut::SetChi2PerClusterITS(float min, float max)
{
  mMinChi2PerClusterITS = min;
  mMaxChi2PerClusterITS = max;
  LOG(info) << "Dielectron Cut, set chi2 per cluster ITS range: " << mMinChi2PerClusterITS << " - " << mMaxChi2PerClusterITS;
}
void DielectronCut::SetMeanClusterSizeITS(float min, float max, float minP, float maxP)
{
  mMinMeanClusterSizeITS = min;
  mMaxMeanClusterSizeITS = max;
  mMinP_ITSClusterSize = minP;
  mMaxP_ITSClusterSize = maxP;
  LOG(info) << "Dielectron Cut, set mean cluster size ITS range: " << mMinMeanClusterSizeITS << " - " << mMaxMeanClusterSizeITS;
}
void DielectronCut::SetChi2TOF(float min, float max)
{
  mMinChi2TOF = min;
  mMaxChi2TOF = max;
  LOG(info) << "Dielectron Cut, set chi2 TOF range: " << mMinChi2TOF << " - " << mMaxChi2TOF;
}

void DielectronCut::SetTrackDca3DRange(float min, float max)
{
  mMinDca3D = min;
  mMaxDca3D = max;
  LOG(info) << "Dielectron Cut, set DCA 3D range in sigma: " << mMinDca3D << " - " << mMaxDca3D;
}
void DielectronCut::SetTrackMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "Dielectron Cut, set max DCA xy: " << mMaxDcaXY;
}
void DielectronCut::SetTrackMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "Dielectron Cut, set max DCA z: " << mMaxDcaZ;
}

void DielectronCut::SetTrackMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
{
  mMaxDcaXYPtDep = ptDepCut;
  LOG(info) << "Dielectron Cut, set max DCA xy pt dep: " << mMaxDcaXYPtDep(1.0);
}
void DielectronCut::ApplyPhiV(bool flag)
{
  mApplyPhiV = flag;
  LOG(info) << "Dielectron Cut, apply phiv cut: " << mApplyPhiV;
}
void DielectronCut::ApplyPrefilter(bool flag)
{
  mApplyPF = flag;
  LOG(info) << "Dielectron Cut, apply prefilter: " << mApplyPF;
}

void DielectronCut::SetPIDScheme(int scheme)
{
  mPIDScheme = scheme;
  LOG(info) << "Dielectron Cut, PID scheme: " << static_cast<int>(mPIDScheme);
}
void DielectronCut::SetMinPinTOF(float min)
{
  mMinPinTOF = min;
  LOG(info) << "Dielectron Cut, set min pin for TOF: " << mMinPinTOF;
}

void DielectronCut::SetMuonExclusionTPC(bool flag)
{
  mMuonExclusionTPC = flag;
  LOG(info) << "Dielectron Cut, set flag for muon exclusion in TPC: " << mMuonExclusionTPC;
}

void DielectronCut::SetTOFbetaRange(float min, float max)
{
  mMinTOFbeta = min;
  mMaxTOFbeta = max;
  LOG(info) << "Dielectron Cut, set TOF beta range (TOFif): " << mMinTOFbeta << " - " << mMaxTOFbeta;
}

void DielectronCut::SetTPCNsigmaElRange(float min, float max)
{
  mMinTPCNsigmaEl = min;
  mMaxTPCNsigmaEl = max;
  LOG(info) << "Dielectron Cut, set TPC n sigma El range: " << mMinTPCNsigmaEl << " - " << mMaxTPCNsigmaEl;
}
void DielectronCut::SetTPCNsigmaMuRange(float min, float max)
{
  mMinTPCNsigmaMu = min;
  mMaxTPCNsigmaMu = max;
  LOG(info) << "Dielectron Cut, set TPC n sigma Mu range: " << mMinTPCNsigmaMu << " - " << mMaxTPCNsigmaMu;
}
void DielectronCut::SetTPCNsigmaPiRange(float min, float max)
{
  mMinTPCNsigmaPi = min;
  mMaxTPCNsigmaPi = max;
  LOG(info) << "Dielectron Cut, set TPC n sigma Pi range: " << mMinTPCNsigmaPi << " - " << mMaxTPCNsigmaPi;
}
void DielectronCut::SetTPCNsigmaKaRange(float min, float max)
{
  mMinTPCNsigmaKa = min;
  mMaxTPCNsigmaKa = max;
  LOG(info) << "Dielectron Cut, set TPC n sigma Ka range: " << mMinTPCNsigmaKa << " - " << mMaxTPCNsigmaKa;
}
void DielectronCut::SetTPCNsigmaPrRange(float min, float max)
{
  mMinTPCNsigmaPr = min;
  mMaxTPCNsigmaPr = max;
  LOG(info) << "Dielectron Cut, set TPC n sigma Pr range: " << mMinTPCNsigmaPr << " - " << mMaxTPCNsigmaPr;
}

void DielectronCut::SetTOFNsigmaElRange(float min, float max)
{
  mMinTOFNsigmaEl = min;
  mMaxTOFNsigmaEl = max;
  LOG(info) << "Dielectron Cut, set TOF n sigma El range: " << mMinTOFNsigmaEl << " - " << mMaxTOFNsigmaEl;
}
void DielectronCut::SetTOFNsigmaMuRange(float min, float max)
{
  mMinTOFNsigmaMu = min;
  mMaxTOFNsigmaMu = max;
  LOG(info) << "Dielectron Cut, set TOF n sigma Mu range: " << mMinTOFNsigmaMu << " - " << mMaxTOFNsigmaMu;
}
void DielectronCut::SetTOFNsigmaPiRange(float min, float max)
{
  mMinTOFNsigmaPi = min;
  mMaxTOFNsigmaPi = max;
  LOG(info) << "Dielectron Cut, set TOF n sigma Pi range: " << mMinTOFNsigmaPi << " - " << mMaxTOFNsigmaPi;
}
void DielectronCut::SetTOFNsigmaKaRange(float min, float max)
{
  mMinTOFNsigmaKa = min;
  mMaxTOFNsigmaKa = max;
  LOG(info) << "Dielectron Cut, set TOF n sigma Ka range: " << mMinTOFNsigmaKa << " - " << mMaxTOFNsigmaKa;
}
void DielectronCut::SetTOFNsigmaPrRange(float min, float max)
{
  mMinTOFNsigmaPr = min;
  mMaxTOFNsigmaPr = max;
  LOG(info) << "Dielectron Cut, set TOF n sigma Pr range: " << mMinTOFNsigmaPr << " - " << mMaxTOFNsigmaPr;
}

void DielectronCut::SetITSNsigmaElRange(float min, float max)
{
  mMinITSNsigmaEl = min;
  mMaxITSNsigmaEl = max;
  LOG(info) << "Dielectron Cut, set ITS n sigma El range: " << mMinITSNsigmaEl << " - " << mMaxITSNsigmaEl;
}
void DielectronCut::SetITSNsigmaMuRange(float min, float max)
{
  mMinITSNsigmaMu = min;
  mMaxITSNsigmaMu = max;
  LOG(info) << "Dielectron Cut, set ITS n sigma Mu range: " << mMinITSNsigmaMu << " - " << mMaxITSNsigmaMu;
}
void DielectronCut::SetITSNsigmaPiRange(float min, float max)
{
  mMinITSNsigmaPi = min;
  mMaxITSNsigmaPi = max;
  LOG(info) << "Dielectron Cut, set ITS n sigma Pi range: " << mMinITSNsigmaPi << " - " << mMaxITSNsigmaPi;
}
void DielectronCut::SetITSNsigmaKaRange(float min, float max)
{
  mMinITSNsigmaKa = min;
  mMaxITSNsigmaKa = max;
  LOG(info) << "Dielectron Cut, set ITS n sigma Ka range: " << mMinITSNsigmaKa << " - " << mMaxITSNsigmaKa;
}
void DielectronCut::SetITSNsigmaPrRange(float min, float max)
{
  mMinITSNsigmaPr = min;
  mMaxITSNsigmaPr = max;
  LOG(info) << "Dielectron Cut, set ITS n sigma Pr range: " << mMinITSNsigmaPr << " - " << mMaxITSNsigmaPr;
}

void DielectronCut::SetPRangeForITSNsigmaKa(float min, float max)
{
  mMinP_ITSNsigmaKa = min;
  mMaxP_ITSNsigmaKa = max;
  LOG(info) << "Dielectron Cut, set p range for ITS n sigma Ka: " << mMinP_ITSNsigmaKa << " - " << mMaxP_ITSNsigmaKa;
}

void DielectronCut::SetPRangeForITSNsigmaPr(float min, float max)
{
  mMinP_ITSNsigmaPr = min;
  mMaxP_ITSNsigmaPr = max;
  LOG(info) << "Dielectron Cut, set p range for ITS n sigma Pr: " << mMinP_ITSNsigmaPr << " - " << mMaxP_ITSNsigmaPr;
}

void DielectronCut::SetMaxPinMuonTPConly(float max)
{
  mMaxPinMuonTPConly = max;
  LOG(info) << "Dielectron Cut, set max pin for Muon ID with TPC only: " << mMaxPinMuonTPConly;
}
void DielectronCut::SetMaxPinForPionRejectionTPC(float max)
{
  mMaxPinForPionRejectionTPC = max;
  LOG(info) << "Dielectron Cut, set max pin for pion rejection in TPC: " << mMaxPinForPionRejectionTPC;
}
void DielectronCut::RequireITSibAny(bool flag)
{
  mRequireITSibAny = flag;
  LOG(info) << "Dielectron Cut, require ITS ib any: " << mRequireITSibAny;
}
void DielectronCut::RequireITSib1st(bool flag)
{
  mRequireITSib1st = flag;
  LOG(info) << "Dielectron Cut, require ITS ib 1st: " << mRequireITSib1st;
}
