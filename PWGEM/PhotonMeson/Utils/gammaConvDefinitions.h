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

/// \brief commonly used definitions for gammaConv tasks
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"

using namespace o2::framework;

AxisSpec const gAxis_zColl{800, -50.f, 50.f};
AxisSpec const gAxis_pT{800, 0.f, 25.f};
AxisSpec const gAxis_pT2d{400, 0.f, 25.f};
AxisSpec const gAxis_eta{800, -PIHalf, PIHalf};
AxisSpec const gAxis_eta2d{400, -PIHalf, PIHalf};
AxisSpec const gAxis_phi{800, 0.f, TwoPI};
AxisSpec const gAxis_r{800, 0.f, 250.f};
AxisSpec const gAxis_r2d{400, 0.f, 250.f};
AxisSpec const gAxis_z2d{400, -250.f, 250.f};
AxisSpec const gAxis_TPCdEdxSig{800, -10.f, 10.f};
AxisSpec const gAxis_radRes{800, -PI, PI};

HistogramSpec const gHistoSpec_hCollisionZ_all_MCTrue{"hCollisionZ_all_MCTrue", "hCollisionZ_all_MCTrue;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};
HistogramSpec const gHistoSpec_hCollisionZ_MCTrue{"hCollisionZ_MCTrue", "hCollisionZ_MCTrue;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};
HistogramSpec const gHistoSpec_hCollisionZ_MCRec{"hCollisionZ_MCRec", "hCollisionZ_MCRec;z (cm);counts", {HistType::kTH1F, {gAxis_zColl}}};
