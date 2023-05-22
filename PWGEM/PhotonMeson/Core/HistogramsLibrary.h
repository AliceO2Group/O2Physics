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
using namespace std;
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

enum EMHistType {
  kEvent = 0,
  kV0 = 1,
  kTrack = 2,
  kCluster = 3,
  kPhoton = 4, // photon candidates
};

namespace o2::aod
{
namespace emphotonhistograms
{
void DefineHistograms(THashList* list, const char* histClass, const char* subGroup = "");
void AddHistClass(THashList* list, const char* histClass);

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
  }
}
} // namespace emphotonhistograms
} // namespace o2::aod

#endif // PWGEM_PHOTONMESON_CORE_HISTOGRAMSLIBRARY_H_
