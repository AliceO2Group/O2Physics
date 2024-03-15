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
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"

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
namespace pwgem::photon::histogram
{
void DefineHistograms(THashList* list, const char* histClass, const char* subGroup = "");
THashList* AddHistClass(THashList* list, const char* histClass);

template <EMHistType htype, typename T>
void FillHistClass(THashList* list, const char* subGroup, T const& obj, const float weight = 1.f)
{
  if constexpr (htype == EMHistType::kEvent) {
    reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("all", 1.f);
    if (obj.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("No TF border", 1.f);
    }
    if (obj.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("No ITS ROF border", 1.f);
    }
    if (obj.sel8()) {
      reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("FT0AND", 1.f);
    }
    if (obj.numContrib() > 0.5) {
      reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("N_{contrib}^{PV} > 0", 1.f);
    }
    if (abs(obj.posZ()) < 10.0) {
      reinterpret_cast<TH1F*>(list->FindObject("hCollisionCounter"))->Fill("|Z_{vtx}| < 10 cm", 1.f);
    }

    reinterpret_cast<TH1F*>(list->FindObject("hZvtx"))->Fill(obj.posZ());
    reinterpret_cast<TH1F*>(list->FindObject("hMultNTracksPV"))->Fill(obj.multNTracksPV());
    reinterpret_cast<TH1F*>(list->FindObject("hMultNTracksPVeta1"))->Fill(obj.multNTracksPVeta1());
    reinterpret_cast<TH2F*>(list->FindObject("hMultFT0"))->Fill(obj.multFT0A(), obj.multFT0C());
    reinterpret_cast<TH1F*>(list->FindObject("hCentFT0A"))->Fill(obj.centFT0A());
    reinterpret_cast<TH1F*>(list->FindObject("hCentFT0C"))->Fill(obj.centFT0C());
    reinterpret_cast<TH1F*>(list->FindObject("hCentFT0M"))->Fill(obj.centFT0M());
    reinterpret_cast<TH2F*>(list->FindObject("hCentFT0MvsMultNTracksPV"))->Fill(obj.centFT0M(), obj.multNTracksPV());
    reinterpret_cast<TH2F*>(list->FindObject("hMultFT0MvsMultNTracksPV"))->Fill(obj.multFT0A() + obj.multFT0C(), obj.multNTracksPV());
  } else if constexpr (htype == EMHistType::kPhoton) { // ROOT::Math::PtEtaPhiMVector
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.Pt());
    reinterpret_cast<TH1F*>(list->FindObject("hY"))->Fill(obj.Rapidity());
    reinterpret_cast<TH1F*>(list->FindObject("hPhi"))->Fill(obj.Phi() < 0.f ? obj.Phi() + TMath::TwoPi() : obj.Phi());
  } else if constexpr (htype == EMHistType::kV0) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hRadius"))->Fill(obj.vz(), obj.v0radius());
    reinterpret_cast<TH1F*>(list->FindObject("hCosPA"))->Fill(obj.cospa());
    reinterpret_cast<TH2F*>(list->FindObject("hCosPA_Rxy"))->Fill(obj.v0radius(), obj.cospa());
    reinterpret_cast<TH1F*>(list->FindObject("hPCA"))->Fill(obj.pca());
    reinterpret_cast<TH1F*>(list->FindObject("hPCA_CosPA"))->Fill(obj.cospa(), obj.pca());
    reinterpret_cast<TH2F*>(list->FindObject("hPCA_Rxy"))->Fill(obj.v0radius(), obj.pca());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj.dcaXYtopv(), obj.dcaZtopv());
    reinterpret_cast<TH2F*>(list->FindObject("hAPplot"))->Fill(obj.alpha(), obj.qtarm());
    reinterpret_cast<TH2F*>(list->FindObject("hMassGamma"))->Fill(obj.v0radius(), obj.mGamma());
    reinterpret_cast<TH2F*>(list->FindObject("hGammaRxy"))->Fill(obj.vx(), obj.vy());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsM"))->Fill(obj.mGamma(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsR"))->Fill(obj.v0radius(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsX"))->Fill(obj.vx(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsY"))->Fill(obj.vy(), obj.chiSquareNDF());
    reinterpret_cast<TH2F*>(list->FindObject("hKFChi2vsZ"))->Fill(obj.vz(), obj.chiSquareNDF());
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
    reinterpret_cast<TH1F*>(list->FindObject("hMeanClusterSizeITS"))->Fill(obj.meanClusterSizeITS() * std::cos(std::atan(obj.tgl())));
    reinterpret_cast<TH2F*>(list->FindObject("hTPCdEdx"))->Fill(obj.tpcInnerParam(), obj.tpcSignal());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaEl"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaEl());
    reinterpret_cast<TH2F*>(list->FindObject("hTPCNsigmaPi"))->Fill(obj.tpcInnerParam(), obj.tpcNSigmaPi());
    reinterpret_cast<TH2F*>(list->FindObject("hXY"))->Fill(obj.x(), obj.y());
    reinterpret_cast<TH2F*>(list->FindObject("hZX"))->Fill(obj.z(), obj.x());
    reinterpret_cast<TH2F*>(list->FindObject("hZY"))->Fill(obj.z(), obj.y());
  } else if constexpr (htype == EMHistType::kTrack) {
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH1F*>(list->FindObject("hQoverPt"))->Fill(obj.sign() / obj.pt());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyz"))->Fill(obj.dcaXY(), obj.dcaZ());
    reinterpret_cast<TH2F*>(list->FindObject("hDCAxyzSigma"))->Fill(obj.dcaXY() / sqrt(obj.cYY()), obj.dcaZ() / sqrt(obj.cZZ()));
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
    reinterpret_cast<TH1F*>(list->FindObject("hMeanClusterSizeITS"))->Fill(obj.meanClusterSizeITS() * std::cos(std::atan(obj.tgl())));
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
  } else if constexpr (htype == EMHistType::kEMCCluster) {
    if (TString(subGroup) == "2D") {
      reinterpret_cast<TH2F*>(list->FindObject("hNCell"))->Fill(obj.nCells(), obj.e());
      reinterpret_cast<TH2F*>(list->FindObject("hM02"))->Fill(obj.m02(), obj.e());
      reinterpret_cast<TH2F*>(list->FindObject("hM20"))->Fill(obj.m20(), obj.e());
      reinterpret_cast<TH2F*>(list->FindObject("hTime"))->Fill(obj.time(), obj.e());
      reinterpret_cast<TH2F*>(list->FindObject("hDistToBC"))->Fill(obj.distanceToBadChannel(), obj.e());
    } else {
      reinterpret_cast<TH1F*>(list->FindObject("hNCell"))->Fill(obj.nCells());
      reinterpret_cast<TH1F*>(list->FindObject("hM02"))->Fill(obj.m02());
      reinterpret_cast<TH1F*>(list->FindObject("hM20"))->Fill(obj.m20());
      reinterpret_cast<TH1F*>(list->FindObject("hTime"))->Fill(obj.time());
      reinterpret_cast<TH1F*>(list->FindObject("hDistToBC"))->Fill(obj.distanceToBadChannel());
    }
    reinterpret_cast<TH1F*>(list->FindObject("hPt"))->Fill(obj.pt());
    reinterpret_cast<TH1F*>(list->FindObject("hE"))->Fill(obj.e());
    reinterpret_cast<TH2F*>(list->FindObject("hEtaPhi"))->Fill(obj.phi(), obj.eta());
    for (size_t itrack = 0; itrack < obj.tracketa().size(); itrack++) { // Fill TrackEtaPhi histogram with delta phi and delta eta of all tracks saved in the vectors in skimmerGammaCalo.cxx
      reinterpret_cast<TH2F*>(list->FindObject("hTrackEtaPhi"))->Fill(obj.trackphi()[itrack] - obj.phi(), obj.tracketa()[itrack] - obj.eta());
    }
  }
}
} // namespace pwgem::photon::histogram
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
