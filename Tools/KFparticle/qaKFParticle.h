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
using namespace o2;
using namespace o2::framework;
using namespace o2::track;

enum Source {
  kPrompt = BIT(0),
  kNonPrompt = BIT(1),
  kReflection = BIT(2)
};

namespace o2::aod
{
namespace kfparticle
{
DECLARE_SOA_COLUMN(PtPi, ptpi, float);
DECLARE_SOA_COLUMN(PtKa, ptka, float);
DECLARE_SOA_COLUMN(DCAXYPiToPV, dcaxypvpi, float);
DECLARE_SOA_COLUMN(DCAXYKaToPV, dcaxypvka, float);
DECLARE_SOA_COLUMN(NSigmaPosPi, tpcnsigmapospi, float);
DECLARE_SOA_COLUMN(NSigmaNegPi, tpcnsigmanegpi, float);
DECLARE_SOA_COLUMN(NSigmaPosKA, tpcnsigmaposka, float);
DECLARE_SOA_COLUMN(NSigmaNegKA, tpcnsigmanegka, float);
DECLARE_SOA_COLUMN(DCAXYPiToSV, distpisv, float);
DECLARE_SOA_COLUMN(DCAXYKaToSV, distkasv, float);
DECLARE_SOA_COLUMN(CosThetaStarPi, costhetastarpi, float);
DECLARE_SOA_COLUMN(CosThetaStarKa, costhetastarka, float);
DECLARE_SOA_COLUMN(DevDau, deviationdau, float);
DECLARE_SOA_COLUMN(DCADau, distdau, float);
DECLARE_SOA_COLUMN(ImpParPiKa, d0pid0ka, float);

DECLARE_SOA_COLUMN(PtDGeo, ptdgeo, float);
DECLARE_SOA_COLUMN(MassGeo, massgeo, float);
DECLARE_SOA_COLUMN(CosPaGeo, cospageo, float);
DECLARE_SOA_COLUMN(DCADPVGeo, distdpvgeo, float);
DECLARE_SOA_COLUMN(DCADPVXYGeo, distdpvxygeo, float);
DECLARE_SOA_COLUMN(Chi2Geo, chi2geo, float);

DECLARE_SOA_COLUMN(PtD, ptd, float);
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

DECLARE_SOA_COLUMN(PtDTopo, ptdtopo, float);
DECLARE_SOA_COLUMN(MassTopo, masstopo, float);
DECLARE_SOA_COLUMN(CosPaTopo, cospatopo, float);
DECLARE_SOA_COLUMN(DCADPVTopo, distdpvtopo, float);
DECLARE_SOA_COLUMN(DCADPVXYTopo, distdpvxytopo, float);
DECLARE_SOA_COLUMN(Chi2GeoTopo, chi2geotopo, float);

} // namespace kfparticle

DECLARE_SOA_TABLE(TreeKF, "AOD", "TREEKF",
                  kfparticle::PtPi,
                  kfparticle::PtKa,
                  kfparticle::DCAXYPiToPV,
                  kfparticle::DCAXYKaToPV,
                  kfparticle::NSigmaPosPi,
                  kfparticle::NSigmaNegPi,
                  kfparticle::NSigmaPosKA,
                  kfparticle::NSigmaNegKA,
                  kfparticle::DCAXYPiToSV,
                  kfparticle::DCAXYKaToSV,
                  kfparticle::CosThetaStarPi,
                  kfparticle::CosThetaStarKa,
                  kfparticle::DCADau,
                  kfparticle::ImpParPiKa,
                  kfparticle::PtDGeo,
                  kfparticle::MassGeo,
                  kfparticle::CosPaGeo,
                  kfparticle::DCADPVGeo,
                  kfparticle::DCADPVXYGeo,
                  kfparticle::Chi2Geo,
                  kfparticle::PtD,
                  kfparticle::Mass,
                  kfparticle::DecayLengthD,
                  kfparticle::DecayLengthDXY,
                  kfparticle::CosPa,
                  kfparticle::Lifetime,
                  kfparticle::NormDecayLength,
                  kfparticle::DCADPV,
                  kfparticle::DCADPVXY,
                  kfparticle::Chi2Topo,
                  kfparticle::PartSource,
                  kfparticle::PtDTopo,
                  kfparticle::MassTopo,
                  kfparticle::DCADPVTopo,
                  kfparticle::DCADPVXYTopo,
                  kfparticle::Chi2GeoTopo);

} // namespace o2::aod

#endif // TOOLS_KFPARTICLE_QAKFPARTICLE_H_
