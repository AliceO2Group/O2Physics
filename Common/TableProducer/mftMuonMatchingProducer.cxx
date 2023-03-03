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
// Made by mooya

#include <iostream>
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CommonUtils/NameConf.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/EventSelectionParams.h"

#include "DetectorsBase/Propagator.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/SMatrix.h"
#include "MFTTracking/Tracker.h"

using namespace std;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using o2::track::TrackParCovFwd;
using o2::track::TrackParFwd;
using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

namespace o2::aod
{

namespace variables
{

DECLARE_SOA_COLUMN(FilteringFlags, filteringFlags, uint64_t);

DECLARE_SOA_COLUMN(PRIMZ, prim_z, float);
DECLARE_SOA_COLUMN(BCID, bcid, int);
DECLARE_SOA_COLUMN(COLL, coll, int);
DECLARE_SOA_COLUMN(NMUONS, nmuons, float);
DECLARE_SOA_COLUMN(NMFTS, nmfts, float);

DECLARE_SOA_COLUMN(MFT_X, mft_x, float);
DECLARE_SOA_COLUMN(MFT_Y, mft_y, float);
DECLARE_SOA_COLUMN(MFT_Z, mft_z, float);
DECLARE_SOA_COLUMN(MFT_PHI, mft_phi, float);
DECLARE_SOA_COLUMN(MFT_PT, mft_pt, float);
DECLARE_SOA_COLUMN(MFT_ETA, mft_eta, float);
DECLARE_SOA_COLUMN(MFT_Z_AT_PRIM, mft_z_at_prim, float);
DECLARE_SOA_COLUMN(MFT_PHI_AT_PRIM, mft_phi_at_prim, float);
DECLARE_SOA_COLUMN(MFT_ETA_AT_PRIM, mft_eta_at_prim, float);
DECLARE_SOA_COLUMN(MFT_Q, mft_q, int);

DECLARE_SOA_COLUMN(MCH_X, mch_x, float);
DECLARE_SOA_COLUMN(MCH_Y, mch_y, float);
DECLARE_SOA_COLUMN(MCH_PHI, mch_phi, float);
DECLARE_SOA_COLUMN(MCH_PT, mch_pt, float);
DECLARE_SOA_COLUMN(MCH_ETA, mch_eta, float);
DECLARE_SOA_COLUMN(MCH_Q, mch_q, int);
DECLARE_SOA_COLUMN(MCH_TYPE, mch_type, int);
} // namespace variables

DECLARE_SOA_TABLE(Events, "AOD", "Events",
                  variables::PRIMZ,
                  variables::BCID,
                  variables::COLL,
                  variables::NMUONS,
                  variables::NMFTS);

DECLARE_SOA_TABLE(Muons, "AOD", "Muons",
                  variables::FilteringFlags,
                  variables::MCH_X,
                  variables::MCH_Y,
                  variables::MCH_PHI,
                  variables::MCH_PT,
                  variables::MCH_ETA,
                  variables::MCH_Q,
                  variables::MCH_TYPE,
                  variables::BCID,
                  variables::COLL);

DECLARE_SOA_TABLE(MFTs, "AOD", "MFTs",
                  variables::MFT_X,
                  variables::MFT_Y,
                  variables::MFT_Z,
                  variables::MFT_PHI,
                  variables::MFT_PT,
                  variables::MFT_ETA,
                  variables::MFT_Z_AT_PRIM,
                  variables::MFT_PHI_AT_PRIM,
                  variables::MFT_ETA_AT_PRIM,
                  variables::MFT_Q,
                  variables::BCID,
                  variables::COLL);

} // namespace o2::aod

struct mftMuonMatchingProducer {

  float mAbsEnd = -505.;
  float mBz = 0.5;
  int iColl = 0;

  Produces<o2::aod::Muons> mchs;
  Produces<o2::aod::MFTs> mfts;
  Produces<o2::aod::Events> events;

  void init(o2::framework::InitContext&)
  {
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::FwdTracks const& muons,
               aod::MFTTracks const& mfttracks)
  {

    if (muons.size() < 1)
      return;
    if (!collision.selection()[kNoPileupFromSPD])
      return;

    int nmchtrack = 0;
    int nmfttrack = 0;

    vector<float> mft_x;
    vector<float> mft_y;
    vector<float> mft_z;
    vector<float> mft_phi;
    vector<float> mft_pt;
    vector<float> mft_eta;
    vector<float> mft_z_prim;
    vector<float> mft_phi_prim;
    vector<float> mft_eta_prim;
    vector<int> mft_q;

    vector<uint8_t> mch_cut;
    vector<float> mch_x;
    vector<float> mch_y;
    vector<float> mch_phi;
    vector<float> mch_pt;
    vector<float> mch_eta;
    vector<int> mch_q;
    vector<int> mch_type;

    float primZ = collision.posZ();

    for (auto const& muon : muons) {
      if (muon.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
        // if (muon.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack)
        continue;
      if (!muon.has_collision())
        continue;
      uint8_t filler = 0;
      float rAbs = muon.rAtAbsorberEnd();
      float eta = muon.eta();
      float phi = muon.phi();
      float x = rAbs * TMath::Cos(phi);
      float y = rAbs * TMath::Sin(phi);
      float pt = muon.pt();
      int q = muon.sign();

      if (eta < -4.0 || eta > -2.5)
        continue;
      if (17.6 < rAbs && rAbs < 26.5) {
        if (muon.pDca() > 594.0)
          continue;
      } else if (26.5 < rAbs && rAbs < 89.5) {
        if (muon.pDca() > 324.0)
          continue;
      }
      if (muon.chi2() > 1e6)
        continue;
      if (muon.chi2MatchMCHMID() > 1e6)
        continue;

      mch_cut.push_back(filler);
      mch_x.push_back(x);
      mch_y.push_back(y);
      mch_phi.push_back(phi);
      mch_pt.push_back(pt);
      mch_eta.push_back(eta);
      mch_q.push_back(q);
      mch_type.push_back(muon.trackType());
      ++nmchtrack;
    } // end of loop tracks

    for (auto const& mfttrack : mfttracks) {
      if (!mfttrack.has_collision())
        continue;

      int q = mfttrack.sign();
      float mftchi2 = mfttrack.chi2();
      SMatrix5 mftpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
      std::vector<double> mftv1;
      SMatrix55 mftcovs(mftv1.begin(), mftv1.end());
      o2::track::TrackParCovFwd mftpars1{mfttrack.z(), mftpars, mftcovs, mftchi2};
      mftpars1.propagateToZ(primZ, mBz);

      if (mftpars1.getEta() < -3.6 || mftpars1.getEta() > -2.5)
        continue;
      mft_z_prim.push_back(mftpars1.getZ());
      mft_phi_prim.push_back(mftpars1.getPhi());
      mft_eta_prim.push_back(mftpars1.getEta());

      mftpars1.propagateToZ(mAbsEnd, mBz);
      float x = mftpars1.getX();
      float y = mftpars1.getY();
      float z = mftpars1.getZ();
      mft_x.push_back(x);
      mft_y.push_back(y);
      mft_z.push_back(z);
      mft_phi.push_back(mftpars1.getPhi());
      mft_pt.push_back(mftpars1.getPt());
      mft_eta.push_back(mftpars1.getEta());
      mft_q.push_back(q);
      nmfttrack++;
    } // end of loop mfttracks

    events(collision.posZ(), collision.bcId(), iColl, nmchtrack, nmfttrack);

    if (nmchtrack > 0) {
      for (int itrack = 0; itrack < nmchtrack; ++itrack) {
        mchs(mch_cut[itrack], mch_x[itrack], mch_y[itrack], mch_phi[itrack], mch_pt[itrack], mch_eta[itrack],
             mch_q[itrack], mch_type[itrack], collision.bcId(), iColl);
      }
      for (int itrack = 0; itrack < nmfttrack; ++itrack) {
        mfts(mft_x[itrack], mft_y[itrack], mft_z[itrack], mft_phi[itrack], mft_pt[itrack], mft_eta[itrack],
             mft_z_prim[itrack], mft_phi_prim[itrack], mft_eta_prim[itrack], mft_q[itrack], collision.bcId(), iColl);
      }
    }
    ++iColl;
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftMuonMatchingProducer>(cfgc),
  };
}
