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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

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

    o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = fMatching.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertex(mchTrack, vtx[0], vtx[1], vtx[2], vtxCov[0], vtxCov[1]);
    auto proptrack = fMatching.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
    return propmuon;
  }

  void process(o2::aod::BCs const& bcs, ForwardTracks const& fwdTracks, o2::aod::Collisions const& cols)
  {
    int run = bcs.begin().runNumber();

    if (run != fRun) {
      fRun = run;
      std::map<string, string> metadata, headers;
      headers = fCCDBApi.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadata, -1);
      int64_t ts = std::atol(headers["SOR"].c_str());
      auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (!o2::base::GeometryManager::isGeometryLoaded())
        fCCDB->get<TGeoManager>("GLO/Config/GeometryAligned");
      o2::mch::TrackExtrap::setField();
    }

    propFwdTracks.reserve(fwdTracks.size());
    propFwdTracksCov.reserve(fwdTracks.size());

    for (const auto& t : fwdTracks) {
      std::array<float, 3> vtx = {0.f, 0.f, 0.f};
      std::array<float, 2> vtxCov = {0.f, 0.f};
      if (t.has_collision()) {
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
                    pft.getTrackChi2(), t.chi2MatchMCHMID(), t.chi2MatchMCHMFT(),
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
      auto rhoXY = static_cast<Char_t>(128. * cov(0, 1) / (sigX * sigY));
      auto rhoPhiX = static_cast<Char_t>(128. * cov(0, 2) / (sigPhi * sigX));
      auto rhoPhiY = static_cast<Char_t>(128. * cov(1, 2) / (sigPhi * sigY));
      auto rhoTglX = static_cast<Char_t>(128. * cov(0, 3) / (sigTgl * sigX));
      auto rhoTglY = static_cast<Char_t>(128. * cov(1, 3) / (sigTgl * sigY));
      auto rhoTglPhi = static_cast<Char_t>(128. * cov(2, 3) / (sigTgl * sigPhi));
      auto rho1PtX = static_cast<Char_t>(128. * cov(0, 4) / (sig1Pt * sigX));
      auto rho1PtY = static_cast<Char_t>(128. * cov(1, 4) / (sig1Pt * sigY));
      auto rho1PtPhi = static_cast<Char_t>(128. * cov(2, 4) / (sig1Pt * sigPhi));
      auto rho1PtTgl = static_cast<Char_t>(128. * cov(3, 4) / (sig1Pt * sigTgl));
      propFwdTracksCov(sigX, sigY, sigTgl, sigPhi, sig1Pt,
                       rhoXY, rhoPhiX, rhoPhiY, rhoTglX,
                       rhoTglY, rhoTglPhi, rho1PtX, rho1PtY,
                       rho1PtPhi, rho1PtTgl);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FwdTrackPropagation>(cfgc)};
}
