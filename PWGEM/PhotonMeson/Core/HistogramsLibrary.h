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

/// \file HistogramsLibrary.h
/// \brief Small histogram library for photon and meson analysis.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
#define PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"

#include <TH2.h>
#include <THashList.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstddef>

enum EMHistType {
  kEvent = 0,
  kEvent_Cent = 1,
  kEvent_Cent_Qvec = 2,
  kV0 = 3,
  kV0Leg = 4,
  kDalitzEE = 5,
  kTrack = 6,
  kPHOSCluster = 7,
  kEMCCluster = 8,
};

const float maxZ = 10.f;

namespace o2::aod
{
namespace pwgem::photon::histogram
{
void DefineHistograms(THashList* list, const char* histClass, const char* subGroup = "");
THashList* AddHistClass(THashList* list, const char* histClass);

template <EMHistType htype1, typename T1>
void FillHistClass(THashList* list, const char* subGroup, T1 const& obj1 /*, const float weight = 1.f*/)
{
  if constexpr (htype1 == EMHistType::kEvent || htype1 == EMHistType::kEvent_Cent || htype1 == EMHistType::kEvent_Cent_Qvec) {
    dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("all", 1.f);
    if (obj1.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("No TF border", 1.f);
    }
    if (obj1.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("No ITS ROF border", 1.f);
    }
    if (obj1.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("No Same Bunch Pileup", 1.f);
    }
    if (obj1.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("Is Vertex ITSTPC", 1.f);
    }
    if (obj1.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("Is Good Zvtx FT0vsPV", 1.f);
    }
    if (obj1.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("FT0AND", 1.f);
    }
    if (obj1.sel8()) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("sel8", 1.f);
    }
    if (std::abs(obj1.posZ()) < maxZ) {
      dynamic_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("|Z_{vtx}| < 10 cm", 1.f);
    }

    dynamic_cast<TH1F*>(list->FindObject("hZvtx"))->Fill(obj1.posZ());
    dynamic_cast<TH1F*>(list->FindObject("hMultNTracksPV"))->Fill(obj1.multNTracksPV());
    dynamic_cast<TH1F*>(list->FindObject("hMultNTracksPVeta1"))->Fill(obj1.multNTracksPVeta1());
    dynamic_cast<TH2F*>(list->FindObject("hMultFT0"))->Fill(obj1.multFT0A(), obj1.multFT0C());
    dynamic_cast<TH1F*>(list->FindObject("hCentFT0A"))->Fill(obj1.centFT0A());
    dynamic_cast<TH1F*>(list->FindObject("hCentFT0C"))->Fill(obj1.centFT0C());
    dynamic_cast<TH1F*>(list->FindObject("hCentFT0M"))->Fill(obj1.centFT0M());
    dynamic_cast<TH2F*>(list->FindObject("hCentFT0MvsMultNTracksPV"))->Fill(obj1.centFT0M(), obj1.multNTracksPV());
    dynamic_cast<TH2F*>(list->FindObject("hMultFT0MvsMultNTracksPV"))->Fill(obj1.multFT0A() + obj1.multFT0C(), obj1.multNTracksPV());

    if constexpr (htype1 == EMHistType::kEvent_Cent_Qvec) {
      dynamic_cast<TH2F*>(list->FindObject("hQ2xFT0M_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xft0m());
      dynamic_cast<TH2F*>(list->FindObject("hQ2yFT0M_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2yft0m());
      dynamic_cast<TH2F*>(list->FindObject("hQ2xFT0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xft0a());
      dynamic_cast<TH2F*>(list->FindObject("hQ2yFT0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2yft0a());
      dynamic_cast<TH2F*>(list->FindObject("hQ2xFT0C_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xft0c());
      dynamic_cast<TH2F*>(list->FindObject("hQ2yFT0C_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2yft0c());
      // static_cast<TH2F*>(list->FindObject("hQ2xFV0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xfv0a());
      // static_cast<TH2F*>(list->FindObject("hQ2yFV0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2yfv0a());
      dynamic_cast<TH2F*>(list->FindObject("hQ2xBPos_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xbpos());
      dynamic_cast<TH2F*>(list->FindObject("hQ2yBPos_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2ybpos());
      dynamic_cast<TH2F*>(list->FindObject("hQ2xBNeg_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2xbneg());
      dynamic_cast<TH2F*>(list->FindObject("hQ2yBNeg_CentFT0C"))->Fill(obj1.centFT0C(), obj1.q2ybneg());

      dynamic_cast<TH2F*>(list->FindObject("hEP2FT0M_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2ft0m());
      dynamic_cast<TH2F*>(list->FindObject("hEP2FT0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2ft0a());
      dynamic_cast<TH2F*>(list->FindObject("hEP2FT0C_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2ft0c());
      // static_cast<TH2F*>(list->FindObject("hEP2FV0A_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2fv0a());
      dynamic_cast<TH2F*>(list->FindObject("hEP2BPos_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2bpos());
      dynamic_cast<TH2F*>(list->FindObject("hEP2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), obj1.ep2bneg());

      std::array<float, 2> q2ft0m = {obj1.q2xft0m(), obj1.q2yft0m()};
      std::array<float, 2> q2ft0a = {obj1.q2xft0a(), obj1.q2yft0a()};
      std::array<float, 2> q2ft0c = {obj1.q2xft0c(), obj1.q2yft0c()};
      // std::array<float, 2> q2fv0a = {obj1.q2xfv0a(), obj1.q2yfv0a()};
      std::array<float, 2> q2bpos = {obj1.q2xbpos(), obj1.q2ybpos()};
      std::array<float, 2> q2bneg = {obj1.q2xbneg(), obj1.q2ybneg()};

      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0MQ2BPos_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bpos));
      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0MQ2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0m, q2bneg));
      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0AQ2BPos_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bpos));
      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0AQ2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0a, q2bneg));
      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0CQ2BPos_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bpos));
      dynamic_cast<TH2F*>(list->FindObject("hQ2FT0CQ2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2ft0c, q2bneg));
      // static_cast<TH2F*>(list->FindObject("hQ2FV0AQ2BPos_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bpos));
      // static_cast<TH2F*>(list->FindObject("hQ2FV0AQ2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2fv0a, q2bneg));
      dynamic_cast<TH2F*>(list->FindObject("hQ2BPosQ2BNeg_CentFT0C"))->Fill(obj1.centFT0C(), RecoDecay::dotProd(q2bpos, q2bneg));
    }
  } else if constexpr (htype1 == EMHistType::kV0) {
    dynamic_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj1.pt());
    dynamic_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj1.phi(), obj1.eta());
    dynamic_cast<TH2F*>(list->FindObject("hRadius"))->Fill(obj1.vz(), obj1.v0radius());
    dynamic_cast<TH1F*>(list->FindObject("hCosPA"))->Fill(obj1.cospa());
    dynamic_cast<TH2F*>(list->FindObject("hCosPA_Rxy"))->Fill(obj1.v0radius(), obj1.cospa());
    dynamic_cast<TH1F*>(list->FindObject("hPCA"))->Fill(obj1.pca());
    dynamic_cast<TH1F*>(list->FindObject("hPCA_CosPA"))->Fill(obj1.cospa(), obj1.pca());
    dynamic_cast<TH2F*>(list->FindObject("hPCA_Rxy"))->Fill(obj1.v0radius(), obj1.pca());
    dynamic_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj1.dcaXYtopv(), obj1.dcaZtopv());
    dynamic_cast<TH2F*>(list->FindObject("hAPplot"))->Fill(obj1.alpha(), obj1.qtarm());
    dynamic_cast<TH2F*>(list->FindObject("hMassGamma"))->Fill(obj1.v0radius(), obj1.mGamma());
    dynamic_cast<TH2F*>(list->FindObject("hGammaRxy"))->Fill(obj1.vx(), obj1.vy());
    dynamic_cast<TH2F*>(list->FindObject("hKFChi2vsM"))->Fill(obj1.mGamma(), obj1.chiSquareNDF());
    dynamic_cast<TH2F*>(list->FindObject("hKFChi2vsR"))->Fill(obj1.v0radius(), obj1.chiSquareNDF());
    dynamic_cast<TH2F*>(list->FindObject("hKFChi2vsX"))->Fill(obj1.vx(), obj1.chiSquareNDF());
    dynamic_cast<TH2F*>(list->FindObject("hKFChi2vsY"))->Fill(obj1.vy(), obj1.chiSquareNDF());
    dynamic_cast<TH2F*>(list->FindObject("hKFChi2vsZ"))->Fill(obj1.vz(), obj1.chiSquareNDF());
  } else if constexpr (htype1 == EMHistType::kV0Leg) {
    dynamic_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj1.pt());
    dynamic_cast<TH1F*>(list->FindObject("hQoverPt"))->Fill(obj1.sign() / obj1.pt());
    dynamic_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj1.phi(), obj1.eta());
    dynamic_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj1.dcaXY(), obj1.dcaZ());
    dynamic_cast<TH1F*>(list->FindObject("hNclsITS"))->Fill(obj1.itsNCls());
    dynamic_cast<TH1F*>(list->FindObject("hNclsTPC"))->Fill(obj1.tpcNClsFound());
    dynamic_cast<TH1F*>(list->FindObject("hNcrTPC"))->Fill(obj1.tpcNClsCrossedRows());
    dynamic_cast<TH1F*>(list->FindObject("hTPCNcr2Nf"))->Fill(obj1.tpcCrossedRowsOverFindableCls());
    dynamic_cast<TH1F*>(list->FindObject("hTPCNcls2Nf"))->Fill(obj1.tpcFoundOverFindableCls());
    dynamic_cast<TH1F*>(list->FindObject("hChi2TPC"))->Fill(obj1.tpcChi2NCl());
    dynamic_cast<TH1F*>(list->FindObject("hChi2ITS"))->Fill(obj1.itsChi2NCl());
    dynamic_cast<TH1F*>(list->FindObject("hITSClusterMap"))->Fill(obj1.itsClusterMap());
    dynamic_cast<TH1F*>(list->FindObject("hMeanClusterSizeITS"))->Fill(obj1.meanClusterSizeITS() * std::cos(std::atan(obj1.tgl())));
    dynamic_cast<TH2F*>(list->FindObject("hTPCdEdx"))->Fill(obj1.tpcInnerParam(), obj1.tpcSignal());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaEl"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaEl());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaPi"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaPi());
    dynamic_cast<TH2F*>(list->FindObject("hXY"))->Fill(obj1.x(), obj1.y());
    dynamic_cast<TH2F*>(list->FindObject("hZX"))->Fill(obj1.z(), obj1.x());
    dynamic_cast<TH2F*>(list->FindObject("hZY"))->Fill(obj1.z(), obj1.y());
  } else if constexpr (htype1 == EMHistType::kTrack) {
    dynamic_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj1.pt());
    dynamic_cast<TH1F*>(list->FindObject("hQoverPt"))->Fill(obj1.sign() / obj1.pt());
    dynamic_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj1.phi(), obj1.eta());
    dynamic_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj1.dcaXY(), obj1.dcaZ());
    dynamic_cast<TH2F*>(list->FindObject("hDCAxyzSigma"))->Fill(obj1.dcaXY() / std::sqrt(obj1.cYY()), obj1.dcaZ() / std::sqrt(obj1.cZZ()));
    dynamic_cast<TH2F*>(list->FindObject("hDCAxyRes_Pt"))->Fill(obj1.pt(), std::sqrt(obj1.cYY()) * 1e+4); // convert cm to um
    dynamic_cast<TH2F*>(list->FindObject("hDCAzRes_Pt"))->Fill(obj1.pt(), std::sqrt(obj1.cZZ()) * 1e+4);  // convert cm to um
    dynamic_cast<TH1F*>(list->FindObject("hNclsITS"))->Fill(obj1.itsNCls());
    dynamic_cast<TH1F*>(list->FindObject("hNclsTPC"))->Fill(obj1.tpcNClsFound());
    dynamic_cast<TH1F*>(list->FindObject("hNcrTPC"))->Fill(obj1.tpcNClsCrossedRows());
    dynamic_cast<TH1F*>(list->FindObject("hTPCNcr2Nf"))->Fill(obj1.tpcCrossedRowsOverFindableCls());
    dynamic_cast<TH1F*>(list->FindObject("hTPCNcls2Nf"))->Fill(obj1.tpcFoundOverFindableCls());
    dynamic_cast<TH1F*>(list->FindObject("hChi2TPC"))->Fill(obj1.tpcChi2NCl());
    dynamic_cast<TH1F*>(list->FindObject("hChi2ITS"))->Fill(obj1.itsChi2NCl());
    dynamic_cast<TH1F*>(list->FindObject("hITSClusterMap"))->Fill(obj1.itsClusterMap());
    dynamic_cast<TH1F*>(list->FindObject("hMeanClusterSizeITS"))->Fill(obj1.meanClusterSizeITS() * std::cos(std::atan(obj1.tgl())));
    dynamic_cast<TH2F*>(list->FindObject("hTPCdEdx"))->Fill(obj1.tpcInnerParam(), obj1.tpcSignal());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaEl"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaEl());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaMu"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaMu());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaPi"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaPi());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaKa"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaKa());
    dynamic_cast<TH2F*>(list->FindObject("hTPCNsigmaPr"))->Fill(obj1.tpcInnerParam(), obj1.tpcNSigmaPr());
    dynamic_cast<TH2F*>(list->FindObject("hTOFbeta"))->Fill(obj1.tpcInnerParam(), obj1.beta());
    dynamic_cast<TH2F*>(list->FindObject("hTOFNsigmaEl"))->Fill(obj1.tpcInnerParam(), obj1.tofNSigmaEl());
    dynamic_cast<TH2F*>(list->FindObject("hTOFNsigmaMu"))->Fill(obj1.tpcInnerParam(), obj1.tofNSigmaMu());
    dynamic_cast<TH2F*>(list->FindObject("hTOFNsigmaPi"))->Fill(obj1.tpcInnerParam(), obj1.tofNSigmaPi());
    dynamic_cast<TH2F*>(list->FindObject("hTOFNsigmaKa"))->Fill(obj1.tpcInnerParam(), obj1.tofNSigmaKa());
    dynamic_cast<TH2F*>(list->FindObject("hTOFNsigmaPr"))->Fill(obj1.tpcInnerParam(), obj1.tofNSigmaPr());
  } else if constexpr (htype1 == EMHistType::kPHOSCluster) {
    dynamic_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj1.pt());
    dynamic_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj1.phi(), obj1.eta());
    dynamic_cast<TH2F*>(list->FindObject("hEvsNcell"))->Fill(obj1.e(), obj1.nCells());
    dynamic_cast<TH2F*>(list->FindObject("hEvsM02"))->Fill(obj1.e(), obj1.m02());
    dynamic_cast<TH2F*>(list->FindObject("hEvsM20"))->Fill(obj1.e(), obj1.m20());
    dynamic_cast<TH1F*>(list->FindObject("hDistToBC"))->Fill(obj1.distanceToBadChannel());
    static_cast<TH2F*>(list->FindObject(Form("hClusterXZM%d", obj1.mod())))->Fill(obj1.cellx(), obj1.cellz());
  } else if constexpr (htype1 == EMHistType::kEMCCluster) {
    if (TString(subGroup) == "2D") {
      dynamic_cast<TH2F*>(list->FindObject("hNCell"))->Fill(obj1.nCells(), obj1.e());
      dynamic_cast<TH2F*>(list->FindObject("hM02"))->Fill(obj1.m02(), obj1.e());
      dynamic_cast<TH2F*>(list->FindObject("hM20"))->Fill(obj1.m20(), obj1.e());
      dynamic_cast<TH2F*>(list->FindObject("hTime"))->Fill(obj1.time(), obj1.e());
      dynamic_cast<TH2F*>(list->FindObject("hDistToBC"))->Fill(obj1.distanceToBadChannel(), obj1.e());
    } else {
      dynamic_cast<TH1F*>(list->FindObject("hNCell"))->Fill(obj1.nCells());
      dynamic_cast<TH1F*>(list->FindObject("hM02"))->Fill(obj1.m02());
      dynamic_cast<TH1F*>(list->FindObject("hM20"))->Fill(obj1.m20());
      dynamic_cast<TH1F*>(list->FindObject("hTime"))->Fill(obj1.time());
      dynamic_cast<TH1F*>(list->FindObject("hDistToBC"))->Fill(obj1.distanceToBadChannel());
    }
    dynamic_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj1.pt());
    dynamic_cast<TH1F*>(list->FindObject("hE"))->Fill(obj1.e());
    dynamic_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj1.phi(), obj1.eta());
    for (size_t itrack = 0; itrack < obj1.deltaEta().size(); itrack++) { // Fill TrackEtaPhi histogram with delta phi and delta eta of all tracks saved in the vectors in skimmerGammaCalo.cxx
      dynamic_cast<TH2F*>(list->FindObject("hTrackEtaPhi"))->Fill(obj1.deltaPhi()[itrack], obj1.deltaEta()[itrack]);
    }
  }
}
} // namespace pwgem::photon::histogram
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
