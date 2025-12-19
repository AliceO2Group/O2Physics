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
/// \author Nazar Burmasov, nazar.burmasov@cern.ch
/// \author Diana Krupova, diana.krupova@cern.ch
/// \since 04.06.2024

#include "PWGUD/DataModel/UDTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Math/SMatrix.h"
#include "TGeoGlobalMagField.h"

using namespace o2::framework;
using namespace o2::framework::expressions;
using o2::aod::fwdtrack::ForwardTrackTypeEnum;

struct FwdTrackPropagation {
  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
  using SMatrix5 = ROOT::Math::SVector<double, 5>;

  Produces<o2::aod::UDFwdTracksProp> propFwdTracks;
  Produces<o2::aod::UDFwdTracksCovProp> propFwdTracksCov;

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::globaltracking::MatchGlobalFwd fMatching;

  int fRun = 0;

  void init(InitContext&)
  {
    fCCDB->setURL("http://alice-ccdb.cern.ch");
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();
    fCCDBApi.init("http://alice-ccdb.cern.ch");
  }

  template <typename T>
  auto propagateFwdToVtx(const T& muon, const std::array<float, 3>& vtx, const std::array<float, 2>& vtxCov)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                           muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                           muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::dataformats::GlobalFwdTrack propmuon;

    if (muon.trackType() > ForwardTrackTypeEnum::GlobalForwardTrack) { // tracks without MFT
      o2::dataformats::GlobalFwdTrack track;
      track.setParameters(tpars);
      track.setZ(muon.z());
      track.setCovariances(tcovs);
      auto mchTrack = fMatching.FwdtoMCH(track);
      o2::mch::TrackExtrap::extrapToVertex(mchTrack, vtx[0], vtx[1], vtx[2], vtxCov[0], vtxCov[1]);
      auto proptrack = fMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());
    } else if (muon.trackType() <= ForwardTrackTypeEnum::GlobalForwardTrack) { // tracks with MFT
      double centerMFT[3] = {0, 0, -61.4};
      o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
      auto* field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
      auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), vtx[0], vtx[1], vtx[2]);
      auto x2x0 = static_cast<float>(geoMan.meanX2X0);
      fwdtrack.propagateToVtxhelixWithMCS(vtx[2], {vtx[0], vtx[1]}, {vtxCov[0], vtxCov[1]}, Bz, x2x0);
      propmuon.setParameters(fwdtrack.getParameters());
      propmuon.setZ(fwdtrack.getZ());
      propmuon.setCovariances(fwdtrack.getCovariances());
    }
    return propmuon;
  }

  void process(o2::aod::BCs const& bcs, ForwardTracks const& fwdTracks, o2::aod::Collisions const& cols)
  {
    int run = bcs.begin().runNumber();

    if (run != fRun) {
      fRun = run;
      std::map<std::string, std::string> metadata;
      auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, run);
      auto ts = soreor.first;
      auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (!o2::base::GeometryManager::isGeometryLoaded())
        fCCDB->get<TGeoManager>("GLO/Config/GeometryAligned");
      o2::mch::TrackExtrap::setField();
    }

    propFwdTracks.reserve(fwdTracks.size());
    propFwdTracksCov.reserve(fwdTracks.size());

    constexpr float ZAbsorberFront = -90.f;
    for (const auto& t : fwdTracks) {
      if (t.z() < ZAbsorberFront) {
        std::array<float, 3> vtx = {0.f, 0.f, 0.f};
        std::array<float, 2> vtxCov = {0.f, 0.f};
        if (t.has_collision() && t.trackType() < ForwardTrackTypeEnum::GlobalForwardTrack) { // propagate only global muon tracks to collision vtx
          auto col = cols.iteratorAt(t.collisionId());
          vtx[0] = col.posX();
          vtx[1] = col.posY();
          vtx[2] = col.posZ();
          vtxCov[0] = col.covXX();
          vtxCov[1] = col.covYY();
        }
        auto pft = propagateFwdToVtx(t, vtx, vtxCov);
        propFwdTracks(t.collisionId(), t.trackType(),
                      pft.getX(), pft.getY(), pft.getZ(), pft.getPhi(), pft.getTgl(), pft.getInvQPt(),
                      pft.getEta(), pft.getPt(), pft.getP(),
                      t.nClusters(), t.pDca(), t.rAtAbsorberEnd(),
                      t.chi2(), t.chi2MatchMCHMID(), t.chi2MatchMCHMFT(),
                      t.matchScoreMCHMFT(), t.matchMFTTrackId(), t.matchMCHTrackId(),
                      t.mchBitMap(), t.midBoards(), t.midBitMap(),
                      t.trackTime(), t.trackTimeRes());
        // debug
        // LOGP(info, "track {}, before: {} {} {} {} {} {}", t.globalIndex(), t.x(), t.y(), t.z(), t.phi(), t.tgl(), t.signed1Pt());
        // LOGP(info, "track {}, after: {} {} {} {} {} {}", t.globalIndex(), pft.getX(), pft.getY(), pft.getZ(), pft.getPhi(), pft.getTgl(), pft.getInvQPt());
        SMatrix55 cov = pft.getCovariances();
        float sigX = TMath::Sqrt(cov(0, 0));
        float sigY = TMath::Sqrt(cov(1, 1));
        float sigPhi = TMath::Sqrt(cov(2, 2));
        float sigTgl = TMath::Sqrt(cov(3, 3));
        float sig1Pt = TMath::Sqrt(cov(4, 4));
        auto rhoXY = static_cast<int8_t>(128. * cov(0, 1) / (sigX * sigY));
        auto rhoPhiX = static_cast<int8_t>(128. * cov(0, 2) / (sigPhi * sigX));
        auto rhoPhiY = static_cast<int8_t>(128. * cov(1, 2) / (sigPhi * sigY));
        auto rhoTglX = static_cast<int8_t>(128. * cov(0, 3) / (sigTgl * sigX));
        auto rhoTglY = static_cast<int8_t>(128. * cov(1, 3) / (sigTgl * sigY));
        auto rhoTglPhi = static_cast<int8_t>(128. * cov(2, 3) / (sigTgl * sigPhi));
        auto rho1PtX = static_cast<int8_t>(128. * cov(0, 4) / (sig1Pt * sigX));
        auto rho1PtY = static_cast<int8_t>(128. * cov(1, 4) / (sig1Pt * sigY));
        auto rho1PtPhi = static_cast<int8_t>(128. * cov(2, 4) / (sig1Pt * sigPhi));
        auto rho1PtTgl = static_cast<int8_t>(128. * cov(3, 4) / (sig1Pt * sigTgl));
        propFwdTracksCov(sigX, sigY, sigTgl, sigPhi, sig1Pt,
                         rhoXY, rhoPhiX, rhoPhiY, rhoTglX,
                         rhoTglY, rhoTglPhi, rho1PtX, rho1PtY,
                         rho1PtPhi, rho1PtTgl);
      } else {
        propFwdTracks(t.collisionId(), t.trackType(),
                      t.x(), t.y(), t.z(), t.phi(), t.tgl(), t.signed1Pt(),
                      t.eta(), t.pt(), t.p(),
                      t.nClusters(), t.pDca(), t.rAtAbsorberEnd(),
                      t.chi2(), t.chi2MatchMCHMID(), t.chi2MatchMCHMFT(),
                      t.matchScoreMCHMFT(), t.matchMFTTrackId(), t.matchMCHTrackId(),
                      t.mchBitMap(), t.midBoards(), t.midBitMap(),
                      t.trackTime(), t.trackTimeRes());
        propFwdTracksCov(t.sigmaX(), t.sigmaY(), t.sigmaTgl(), t.sigmaPhi(), t.sigma1Pt(),
                         t.rhoXY(), t.rhoPhiX(), t.rhoPhiY(), t.rhoTglX(),
                         t.rhoTglY(), t.rhoTglPhi(), t.rho1PtX(), t.rho1PtY(),
                         t.rho1PtPhi(), t.rho1PtTgl());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FwdTrackPropagation>(cfgc)};
}
