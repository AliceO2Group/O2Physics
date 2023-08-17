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

/// \file qaKFParticleLc.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>

#ifndef TOOLS_KFPARTICLE_QAKFPARTICLELC_H_
#define TOOLS_KFPARTICLE_QAKFPARTICLELC_H_

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
namespace kfparticleLc
{
DECLARE_SOA_COLUMN(RUNNUMBER, Runnumber, int);
DECLARE_SOA_COLUMN(PtPi, ptpi, float);
DECLARE_SOA_COLUMN(PtKa, ptka, float);
DECLARE_SOA_COLUMN(PtPr, ptpr, float);
DECLARE_SOA_COLUMN(DCAXYPiToPV, dcaxypvpi, float);
DECLARE_SOA_COLUMN(DCAXYKaToPV, dcaxypvka, float);
DECLARE_SOA_COLUMN(DCAXYPrToPV, dcaxypvpr, float);
DECLARE_SOA_COLUMN(DCAPiToPV, dcapvpi, float);
DECLARE_SOA_COLUMN(DCAKaToPV, dcapvka, float);
DECLARE_SOA_COLUMN(DCAPrToPV, dcapvpr, float);
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcnsigmapi, float);
DECLARE_SOA_COLUMN(TPCNSigmaKA, tpcnsigmaka, float);
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcnsigmapr, float);
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofnsigmapi, float);
DECLARE_SOA_COLUMN(TOFNSigmaKA, tofnsigmaka, float);
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofnsigmapr, float);
DECLARE_SOA_COLUMN(DCAXYPiToSV, distpisv, float);
DECLARE_SOA_COLUMN(DCAXYKaToSV, distkasv, float);
DECLARE_SOA_COLUMN(DCAXYPrToSV, distprsv, float);
DECLARE_SOA_COLUMN(PtGeo, ptgeo, float);
DECLARE_SOA_COLUMN(MassGeo, massgeo, float);
DECLARE_SOA_COLUMN(CosPaGeo, cospageo, float);
DECLARE_SOA_COLUMN(CosPaXYGeo, cospaxygeo, float);
DECLARE_SOA_COLUMN(DCAPVGeo, distpvgeo, float);
DECLARE_SOA_COLUMN(DCAPVXYGeo, distpvxygeo, float);
DECLARE_SOA_COLUMN(Chi2Geo, chi2geo, float);
DECLARE_SOA_COLUMN(PtTopo, pttopo, float);
DECLARE_SOA_COLUMN(MassTopo, masstopo, float);
DECLARE_SOA_COLUMN(DecayLength, decaylength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decaylengthxy, float);
DECLARE_SOA_COLUMN(CosPaTopo, cospatopo, float);
DECLARE_SOA_COLUMN(Lifetime, lifetime, float);
DECLARE_SOA_COLUMN(NormDecayLength, normdecaylength, float);
DECLARE_SOA_COLUMN(Chi2Topo, chi2topo, float);
DECLARE_SOA_COLUMN(PartSource, source, float);

} // namespace kfparticleLc

DECLARE_SOA_TABLE(TreeKFLc, "AOD", "TREEKFLc",
                  kfparticleLc::RUNNUMBER,
                  kfparticleLc::PtPi,
                  kfparticleLc::PtKa,
                  kfparticleLc::PtPr,
                  kfparticleLc::DCAXYPiToPV,
                  kfparticleLc::DCAXYKaToPV,
                  kfparticleLc::DCAXYPrToPV,
                  kfparticleLc::DCAPiToPV,
                  kfparticleLc::DCAKaToPV,
                  kfparticleLc::DCAPrToPV,
                  kfparticleLc::TPCNSigmaPi,
                  kfparticleLc::TPCNSigmaKA,
                  kfparticleLc::TPCNSigmaPr,
                  kfparticleLc::TOFNSigmaPi,
                  kfparticleLc::TOFNSigmaKA,
                  kfparticleLc::TOFNSigmaPr,
                  kfparticleLc::DCAXYPiToSV,
                  kfparticleLc::DCAXYKaToSV,
                  kfparticleLc::DCAXYPrToSV,
                  kfparticleLc::PtGeo,
                  kfparticleLc::MassGeo,
                  kfparticleLc::CosPaGeo,
                  kfparticleLc::CosPaXYGeo,
                  kfparticleLc::DCAPVGeo,
                  kfparticleLc::DCAPVXYGeo,
                  kfparticleLc::Chi2Geo,
                  kfparticleLc::PtTopo,
                  kfparticleLc::MassTopo,
                  kfparticleLc::DecayLength,
                  kfparticleLc::DecayLengthXY,
                  kfparticleLc::CosPaTopo,
                  kfparticleLc::Lifetime,
                  kfparticleLc::NormDecayLength,
                  kfparticleLc::Chi2Topo,
                  kfparticleLc::PartSource);

} // namespace o2::aod

#endif // TOOLS_KFPARTICLE_QAKFPARTICLELC_H_
