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

/// \file qaKFParticle.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>

#ifndef TOOLS_KFPARTICLE_QAKFPARTICLE_H_
#define TOOLS_KFPARTICLE_QAKFPARTICLE_H_

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

enum Source {
  kPrompt = BIT(0),
  kNonPrompt = BIT(1),
  kReflection = BIT(2)
};

namespace o2::aod
{
namespace kfparticle
{
DECLARE_SOA_COLUMN(RUNNUMBER, Runnumber, int);
DECLARE_SOA_COLUMN(PVCONTRIB, PVContrib, int);
DECLARE_SOA_COLUMN(PtPi, ptpi, float);
DECLARE_SOA_COLUMN(PtKa, ptka, float);
DECLARE_SOA_COLUMN(EtaPi, etapi, float);
DECLARE_SOA_COLUMN(EtaKa, etaka, float);
DECLARE_SOA_COLUMN(RapPi, rappi, float);
DECLARE_SOA_COLUMN(RapKa, rapka, float);
DECLARE_SOA_COLUMN(TPCNclsPi, tpcnclspi, float);
DECLARE_SOA_COLUMN(TPCNclsKa, tpcnclska, float);
DECLARE_SOA_COLUMN(PhiPi, phipi, float);
DECLARE_SOA_COLUMN(PhiKa, phika, float);
DECLARE_SOA_COLUMN(DCAXYPiToPV, dcaxypvpi, float);
DECLARE_SOA_COLUMN(DCAXYKaToPV, dcaxypvka, float);
DECLARE_SOA_COLUMN(DCAPiToPV, dcapvpi, float);
DECLARE_SOA_COLUMN(DCAKaToPV, dcapvka, float);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcnsigmapi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKA, tpcnsigmaka, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofnsigmapi, float);
DECLARE_SOA_COLUMN(TOFNSigmaKA, tofnsigmaka, float);
DECLARE_SOA_COLUMN(DCAXYPiToSV, distpisv, float);
DECLARE_SOA_COLUMN(DCAXYKaToSV, distkasv, float);
DECLARE_SOA_COLUMN(CosThetaStar, costhetastar, float);
DECLARE_SOA_COLUMN(DevDau, deviationdau, float);
DECLARE_SOA_COLUMN(DCADau, distdau, float);
DECLARE_SOA_COLUMN(ImpParPiKa, d0pid0ka, float);

DECLARE_SOA_COLUMN(PtDGeo, ptdgeo, float);
DECLARE_SOA_COLUMN(EtaDGeo, etadgeo, float);
DECLARE_SOA_COLUMN(PhiDGeo, phidgeo, float);
DECLARE_SOA_COLUMN(MassGeo, massgeo, float);
DECLARE_SOA_COLUMN(CosPaGeo, cospageo, float);
DECLARE_SOA_COLUMN(CosPaXYGeo, cospaxygeo, float);
DECLARE_SOA_COLUMN(DCADPVGeo, distdpvgeo, float);
DECLARE_SOA_COLUMN(DCADPVXYGeo, distdpvxygeo, float);
DECLARE_SOA_COLUMN(Chi2Geo, chi2geo, float);

DECLARE_SOA_COLUMN(PtD, ptd, float);
DECLARE_SOA_COLUMN(EtaD, etad, float);
DECLARE_SOA_COLUMN(PhiD, phid, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(DecayLengthD, decaylengthd, float);
DECLARE_SOA_COLUMN(DecayLengthDXY, decaylengthdxy, float);
DECLARE_SOA_COLUMN(CosPa, cospa, float);
DECLARE_SOA_COLUMN(Lifetime, lifetime, float);
DECLARE_SOA_COLUMN(NormDecayLength, normdecaylength, float);
DECLARE_SOA_COLUMN(DCADPV, distdpv, float);
DECLARE_SOA_COLUMN(DCADPVXY, distdpvxy, float);
DECLARE_SOA_COLUMN(Chi2Topo, chi2topo, float);
DECLARE_SOA_COLUMN(PartSource, source, float);

} // namespace kfparticle

DECLARE_SOA_TABLE(TreeKF, "AOD", "TREEKF",
                  kfparticle::RUNNUMBER,
                  kfparticle::PVCONTRIB,
                  kfparticle::PtPi,
                  kfparticle::PtKa,
                  kfparticle::EtaPi,
                  kfparticle::EtaKa,
                  kfparticle::RapPi,
                  kfparticle::RapKa,
                  kfparticle::PhiPi,
                  kfparticle::PhiKa,
                  kfparticle::DCAXYPiToPV,
                  kfparticle::DCAXYKaToPV,
                  kfparticle::DCAPiToPV,
                  kfparticle::DCAKaToPV,
                  kfparticle::TPCNclsPi,
                  kfparticle::TPCNclsKa,
                  kfparticle::TPCNSigmaPi,
                  kfparticle::TPCNSigmaKA,
                  kfparticle::TOFNSigmaPi,
                  kfparticle::TOFNSigmaKA,
                  kfparticle::DCAXYPiToSV,
                  kfparticle::DCAXYKaToSV,
                  kfparticle::CosThetaStar,
                  kfparticle::DCADau,
                  kfparticle::ImpParPiKa,
                  kfparticle::PtDGeo,
                  kfparticle::EtaDGeo,
                  kfparticle::PhiDGeo,
                  kfparticle::MassGeo,
                  kfparticle::CosPaGeo,
                  kfparticle::CosPaXYGeo,
                  kfparticle::DCADPVGeo,
                  kfparticle::DCADPVXYGeo,
                  kfparticle::Chi2Geo,
                  kfparticle::PtD,
                  kfparticle::EtaD,
                  kfparticle::PhiD,
                  kfparticle::Mass,
                  kfparticle::DecayLengthD,
                  kfparticle::DecayLengthDXY,
                  kfparticle::CosPa,
                  kfparticle::Lifetime,
                  kfparticle::NormDecayLength,
                  kfparticle::DCADPV,
                  kfparticle::DCADPVXY,
                  kfparticle::Chi2Topo,
                  kfparticle::PartSource);

} // namespace o2::aod

#endif // TOOLS_KFPARTICLE_QAKFPARTICLE_H_
