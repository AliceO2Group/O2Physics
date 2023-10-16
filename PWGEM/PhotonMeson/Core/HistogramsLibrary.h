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

#ifndef PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
#define PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_

#include <iostream>
#include <array>
#include <TString.h>
#include <THashList.h>
#include <TObject.h>
#include <TObjArray.h>
#include <THashList.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THn.h>
#include <THnSparse.h>
#include <TIterator.h>
#include <TClass.h>

enum class EMHistType : int {
  kEvent = 0,
  kV0 = 1,
  kV0Leg = 2,
  kDalitzEE = 3,
  kTrack = 4,
  kPHOSCluster = 5,
  kEMCCluster = 6,
  kPhoton = 7, // photon candidates
};

namespace o2::aod
{
namespace emphotonhistograms
{
void DefineHistograms(THashList* list, const char* histClass, const char* subGroup = "");
THashList* AddHistClass(THashList* list, const char* histClass);

template <EMHistType htype, typename T>
void FillHistClass(THashList* list, const char* subGroup, T const& obj)
{
  if constexpr (htype == EMHistType::kEvent) {
    reinterpret_cast<TH1F*>(list->FindObject("hMultNTracksPV"))->Fill(obj.multNTracksPV());
    reinterpret_cast<TH1F*>(list->FindObject("hMultNTracksPVeta1"))->Fill(obj.multNTracksPVeta1());
    reinterpret_cast<TH2F*>(list->FindObject("hMultFT0"))->Fill(obj.multFT0A(), obj.multFT0C());
    reinterpret_cast<TH1F*>(list->FindObject("hCentFT0M"))->Fill(obj.centFT0M());
    reinterpret_cast<TH2F*>(list->FindObject("hCentFT0MvsMultNTracksPV"))->Fill(obj.centFT0M(), obj.multNTracksPV());
  } else if constexpr (htype == EMHistType::kPhoton) { // ROOT::Math::PtEtaPhiMVector
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.Pt());
    reinterpret_cast<TH1F*>(list->FindObject("hY"))->Fill(obj.Rapidity());
    reinterpret_cast<TH1F*>(list->FindObject("hPhi"))->Fill(obj.Phi() < 0.f ? obj.Phi() + TMath::TwoPi() : obj.Phi());
  } else if constexpr (htype == EMHistType::kV0) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hRadius"))->Fill(obj.vz(), obj.v0radius());
    reinterpret_cast<TH2F*>(list->FindObject("hRadius_recalc"))->Fill(obj.recalculatedVtxZ(), obj.recalculatedVtxR());
    reinterpret_cast<TH1F*>(list->FindObject("hCosPA"))->Fill(obj.cospa());
    reinterpret_cast<TH1F*>(list->FindObject("hPCA"))->Fill(obj.pca());
    reinterpret_cast<TH2F*>(list->FindObject("hAPplot"))->Fill(obj.alpha(), obj.qtarm());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGamma"))->Fill(obj.v0radius(), obj.mGamma());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGamma_recalc"))->Fill(obj.recalculatedVtxR(), obj.mGamma());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGammaKF_SV_Rxy"))->Fill(obj.recalculatedVtxR(), obj.mGammaKFSV());
    reinterpret_cast<TH2F*>(list->FindObject("hGammaRxy"))->Fill(obj.vx(), obj.vy());
    reinterpret_cast<TH2F*>(list->FindObject("hGammaRxy_recalc"))->Fill(obj.recalculatedVtxX(), obj.recalculatedVtxY());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsR"))->Fill(obj.recalculatedVtxR(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsX"))->Fill(obj.recalculatedVtxX(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsY"))->Fill(obj.recalculatedVtxY(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsZ"))->Fill(obj.recalculatedVtxZ(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGammaKF_PV_SV"))->Fill(obj.mGammaKFPV(), obj.mGammaKFSV());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGammaKF_SV_PsiPair"))->Fill(abs(obj.psipair()), obj.mGammaKFSV());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGammaKF_SV_PhiV"))->Fill(obj.phiv(), obj.mGammaKFSV());
  } else if constexpr (htype == EMHistType::kV0Leg) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH1F*>(list->FindObject("hQoverPt"))->Fill(obj.sign() / obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj.dcaXY(), obj.dcaZ());
    reinterpret_cast<TH1F*>(list->FindObject("hNclsITS"))->Fill(obj.itsNCls());
    reinterpret_cast<TH1F*>(list->FindObject("hNclsTPC"))->Fill(obj.tpcNClsFound());
    reinterpret_cast<TH1F*>(list->FindObject("hNcrTPC"))->Fill(obj.tpcNClsCrossedRows());
    reinterpret_cast<TH1F*>(list->FindObject("hTPCNcr2Nf"))->Fill(obj.tpcCrossedRowsOverFindableCls());
    reinterpret_cast<TH1F*>(list->FindObject("hTPCNcls2Nf"))->Fill(obj.tpcFoundOverFindableCls());
    reinterpret_cast<TH1F*>(list->FindObject("hChi2TPC"))->Fill(obj.tpcChi2NCl());
    reinterpret_cast<TH1F*>(list->FindObject("hChi2ITS"))->Fill(obj.itsChi2NCl());
    reinterpret_cast<TH1F*>(list->FindObject("hITSClusterMap"))->Fill(obj.itsClusterMap());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCdEdx"))->Fill(obj.tpcInnerParam(), obj.tpcSignal());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaEl"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaEl());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaPi"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaPi());
    reinterpret_cast<TH2F*>(list->FindObject("hXY"))->Fill(obj.x(), obj.y());
    reinterpret_cast<TH2F*>(list->FindObject("hZX"))->Fill(obj.z(), obj.x());
    reinterpret_cast<TH2F*>(list->FindObject("hZY"))->Fill(obj.z(), obj.y());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyZ"))->Fill(obj.z(), obj.dcaXY());
  } else if constexpr (htype == EMHistType::kTrack) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH1F*>(list->FindObject("hQoverPt"))->Fill(obj.sign() / obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj.dcaXY(), obj.dcaZ());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyRes_Pt"))->Fill(obj.pt(), sqrt(obj.cYY()) * 1e+4); // convert cm to um
    reinterpret_cast<TH2F*>(list->FindObject("hDCAzRes_Pt"))->Fill(obj.pt(), sqrt(obj.cZZ()) * 1e+4);  // convert cm to um
    reinterpret_cast<TH1F*>(list->FindObject("hNclsITS"))->Fill(obj.itsNCls());
    reinterpret_cast<TH1F*>(list->FindObject("hNclsTPC"))->Fill(obj.tpcNClsFound());
    reinterpret_cast<TH1F*>(list->FindObject("hNcrTPC"))->Fill(obj.tpcNClsCrossedRows());
    reinterpret_cast<TH1F*>(list->FindObject("hTPCNcr2Nf"))->Fill(obj.tpcCrossedRowsOverFindableCls());
    reinterpret_cast<TH1F*>(list->FindObject("hTPCNcls2Nf"))->Fill(obj.tpcFoundOverFindableCls());
    reinterpret_cast<TH1F*>(list->FindObject("hChi2TPC"))->Fill(obj.tpcChi2NCl());
    reinterpret_cast<TH1F*>(list->FindObject("hChi2ITS"))->Fill(obj.itsChi2NCl());
    reinterpret_cast<TH1F*>(list->FindObject("hITSClusterMap"))->Fill(obj.itsClusterMap());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCdEdx"))->Fill(obj.tpcInnerParam(), obj.tpcSignal());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaEl"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaEl());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaMu"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaMu());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaPi"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaPi());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaKa"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaKa());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaPr"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaPr());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFbeta"))->Fill(obj.tpcInnerParam(), obj.beta());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFNsigmaEl"))->Fill(obj.tpcInnerParam(), obj.tofNSigmaEl());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFNsigmaMu"))->Fill(obj.tpcInnerParam(), obj.tofNSigmaMu());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFNsigmaPi"))->Fill(obj.tpcInnerParam(), obj.tofNSigmaPi());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFNsigmaKa"))->Fill(obj.tpcInnerParam(), obj.tofNSigmaKa());
    reinterpret_cast<TH2F*>(list->FindObject("hTOFNsigmaPr"))->Fill(obj.tpcInnerParam(), obj.tofNSigmaPr());
  } else if constexpr (htype == EMHistType::kPHOSCluster) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hEvsNcell"))->Fill(obj.e(), obj.nCells());
    reinterpret_cast<TH2F*>(list->FindObject("hEvsM02"))->Fill(obj.e(), obj.m02());
    reinterpret_cast<TH2F*>(list->FindObject("hEvsM20"))->Fill(obj.e(), obj.m20());
    reinterpret_cast<TH1F*>(list->FindObject("hDistToBC"))->Fill(obj.distanceToBadChannel());
    reinterpret_cast<TH2F*>(list->FindObject(Form("hClusterXZM%d", obj.mod())))->Fill(obj.cellx(), obj.cellz());
  }
}
} // namespace emphotonhistograms
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
