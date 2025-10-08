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

/// \file gammaConvDefinitions.h
/// \brief commonly used definitions for gammaConv tasks
/// \author stephan.friedrich.stiefelmaier@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_GAMMACONVDEFINITIONS_H_
#define PWGEM_PHOTONMESON_UTILS_GAMMACONVDEFINITIONS_H_

#include <CommonConstants/MathConstants.h>
#include <Framework/HistogramSpec.h>

using namespace o2::framework;

AxisSpec const gAxis_zColl{800, -50.f, 50.f};
AxisSpec const gAxis_pT{800, 0.f, 25.f};
AxisSpec const gAxis_pT_armenteros{400, 0.f, 1.f};
AxisSpec const gAxis_pT2d{400, 0.f, 25.f};
AxisSpec const gAxis_eta{800, -o2::constants::math::PIHalf, o2::constants::math::PIHalf};
AxisSpec const gAxis_eta2d{400, -o2::constants::math::PIHalf, o2::constants::math::PIHalf};
AxisSpec const gAxis_phi{800, 0.f, o2::constants::math::TwoPI};
AxisSpec const gAxis_r{800, 0.f, 200.f};
AxisSpec const gAxis_r_extended{800, 0.f, 500.f};
AxisSpec const gAxis_dr{200, -100.f, 100.f};
AxisSpec const gAxis_r2d{400, 0.f, 250.f};
AxisSpec const gAxis_z2d{1000, -250.f, 250.f};
AxisSpec const gAxis_TPCdEdxSig{401, -10.025f, 10.025f};
AxisSpec const gAxis_radRes{800, -o2::constants::math::PI, o2::constants::math::PI};
AxisSpec const gAxis_xyz{2400, -300.f, 300.f};
AxisSpec const gAxis_chi2{501, -1.f, 500.f};
AxisSpec gAxis_pT_log{800, 0.01f, 25.f};

HistogramSpec const gHistoSpec_hCollisionZ_all_MCTrue{"hCollisionZ_all_MCTrue", "hCollisionZ_all_MCTrue;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};
HistogramSpec const gHistoSpec_hCollisionZ_MCTrue{"hCollisionZ_MCTrue", "hCollisionZ_MCTrue;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};
HistogramSpec const gHistoSpec_hCollisionZ_MCRec{"hCollisionZ_MCRec", "hCollisionZ_MCRec;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};

#endif // PWGEM_PHOTONMESON_UTILS_GAMMACONVDEFINITIONS_H_
