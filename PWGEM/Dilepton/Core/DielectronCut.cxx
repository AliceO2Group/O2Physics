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

#include "Framework/Logger.h"
#include "PWGEM/Dilepton/Core/DielectronCut.h"

ClassImp(DielectronCut);

const char* DielectronCut::mCutNames[static_cast<int>(DielectronCut::DielectronCuts::kNCuts)] = {"Mee", "PairPtRange", "PairRapidityRange", "PairDCARange", "PhivPair", "TrackPtRange", "TrackEtaRange", "TPCNCls", "TPCCrossedRows", "TPCCrossedRowsOverNCls", "TPCChi2NDF", "TPCNsigmaEl", "TPCNsigmaMu", "TPCNsigmaPi", "TPCNsigmaKa", "TPCNsigmaPr", "TOFNsigmaEl", "TOFNsigmaMu", "TOFNsigmaPi", "TOFNsigmaKa", "TOFNsigmaPr", "DCA3Dsigma", "DCAxy", "DCAz", "ITSNCls", "ITSChi2NDF", "ITSClusterSize", "Prefilter"};

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
  LOG(info) << "Dielectron Cut, set pair eta range: " << mMinPairY << " - " << mMaxPairY;
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
void DielectronCut::SetMaxPhivPairMeeDep(std::function<float(float)> meeDepCut)
{
  mMaxPhivPairMeeDep = meeDepCut;
  LOG(info) << "Dielectron Cut, set max phiv pair mee dep: " << mMaxPhivPairMeeDep(0.02);
}
void DielectronCut::SelectPhotonConversion(bool flag)
{
  mSelectPC = flag;
  LOG(info) << "Dielectron Cut, select photon conversion: " << mSelectPC;
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
void DielectronCut::SetMeanClusterSizeITS(float min, float max, float maxP)
{
  mMinMeanClusterSizeITS = min;
  mMaxMeanClusterSizeITS = max;
  mMaxP_ITSClusterSize = maxP;
  LOG(info) << "Dielectron Cut, set mean cluster size ITS range: " << mMinMeanClusterSizeITS << " - " << mMaxMeanClusterSizeITS;
}
void DielectronCut::SetDca3DRange(float min, float max)
{
  mMinDca3D = min;
  mMaxDca3D = max;
  LOG(info) << "Dielectron Cut, set DCA 3D range in sigma: " << mMinDca3D << " - " << mMaxDca3D;
}
void DielectronCut::SetMaxDcaXY(float maxDcaXY)
{
  mMaxDcaXY = maxDcaXY;
  LOG(info) << "Dielectron Cut, set max DCA xy: " << mMaxDcaXY;
}
void DielectronCut::SetMaxDcaZ(float maxDcaZ)
{
  mMaxDcaZ = maxDcaZ;
  LOG(info) << "Dielectron Cut, set max DCA z: " << mMaxDcaZ;
}

void DielectronCut::SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut)
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
void DielectronCut::SetMaxPinMuonTPConly(float max)
{
  mMaxPinMuonTPConly = max;
  LOG(info) << "Dielectron Cut, set max pin for Muon ID with TPC only: " << mMaxPinMuonTPConly;
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

void DielectronCut::print() const
{
  LOG(info) << "Dalitz EE Cut:";
  for (int i = 0; i < static_cast<int>(DielectronCuts::kNCuts); i++) {
    switch (static_cast<DielectronCuts>(i)) {
      case DielectronCuts::kTrackPtRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackPt << ", " << mMaxTrackPt << "]";
        break;
      case DielectronCuts::kTrackEtaRange:
        LOG(info) << mCutNames[i] << " in [" << mMinTrackEta << ", " << mMaxTrackEta << "]";
        break;
      case DielectronCuts::kTPCNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNClustersTPC;
        break;
      case DielectronCuts::kTPCCrossedRows:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsTPC;
        break;
      case DielectronCuts::kTPCCrossedRowsOverNCls:
        LOG(info) << mCutNames[i] << " > " << mMinNCrossedRowsOverFindableClustersTPC;
        break;
      case DielectronCuts::kTPCChi2NDF:
        LOG(info) << mCutNames[i] << " < " << mMaxChi2PerClusterTPC;
        break;
      case DielectronCuts::kDCAxy:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaXY;
        break;
      case DielectronCuts::kDCAz:
        LOG(info) << mCutNames[i] << " < " << mMaxDcaZ;
        break;
      default:
        LOG(fatal) << "Cut unknown!";
    }
  }
}
