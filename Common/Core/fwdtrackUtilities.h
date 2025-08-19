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

/// \file fwdtrackUtilities.h
/// \brief Utilities for manipulating parameters of fwdtracks
/// \author Maurice Coquet <maurice.louis.coquet@cern.ch>
/// \author Luca Micheletti <luca.micheletti@cern.ch>
/// \author Daiki Sekihata <daiki.sekihata@cern.ch>

#ifndef COMMON_CORE_FWDTRACKUTILITIES_H_
#define COMMON_CORE_FWDTRACKUTILITIES_H_

#include <DetectorsBase/GeometryManager.h>
#include <Field/MagneticField.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <TGeoGlobalMagField.h>

#include <vector>

namespace o2::aod
{
namespace fwdtrackutils
{
// Index used to set different options for muon propagation
enum class propagationPoint : int {
  kToVertex = 0,
  kToDCA = 1,
  kToRabs = 2,
};
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

/// propagate fwdtrack to a certain point.
template <typename TFwdTrack, typename TCollision>
o2::dataformats::GlobalFwdTrack propagateMuon(TFwdTrack const& muon, TCollision const& collision, const propagationPoint endPoint)
{
  double chi2 = muon.chi2();
  SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
  std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                         muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                         muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
  SMatrix55 tcovs(v1.begin(), v1.end());
  o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
  o2::dataformats::GlobalFwdTrack propmuon;
  o2::globaltracking::MatchGlobalFwd mMatching;

  if (static_cast<int>(muon.trackType()) > 2) { // MCH-MID or MCH standalone
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = mMatching.FwdtoMCH(track);

    if (endPoint == propagationPoint::kToVertex) {
      o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());
    } else if (endPoint == propagationPoint::kToDCA) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.posZ());
    } else if (endPoint == propagationPoint::kToRabs) {
      o2::mch::TrackExtrap::extrapToZ(mchTrack, -505.);
    }

    auto proptrack = mMatching.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
  } else if (static_cast<int>(muon.trackType()) < 2) { // MFT-MCH-MID
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
    auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.posX(), collision.posY(), collision.posZ());
    auto x2x0 = static_cast<float>(geoMan.meanX2X0);
    if (endPoint == propagationPoint::kToVertex) {
      fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
    } else if (endPoint == propagationPoint::kToDCA) {
      fwdtrack.propagateToZhelix(collision.posZ(), Bz);
    }
    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());
  }

  v1.clear();
  v1.shrink_to_fit();

  return propmuon;
}
} // namespace fwdtrackutils
} // namespace o2::aod

#endif // COMMON_CORE_FWDTRACKUTILITIES_H_
