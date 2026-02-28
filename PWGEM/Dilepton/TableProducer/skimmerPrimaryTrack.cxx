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

/// \brief write relevant information about primary tracks.
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/PIDResponseITS.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

struct skimmerPrimaryTrack {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::EMEoIs>;
  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
  using MyTrack = MyTracks::iterator;
  using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;
  using MyTrackMC = MyTracksMC::iterator;

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryTracks> emprimarytracks;
  Produces<aod::EMPrimaryTrackEMEventIdsTMP> prmtrackeventidtmp;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};

  // Operation and minimisation criteria
  Configurable<bool> fillQAHistogram{"fillQAHistogram", false, "flag to fill QA histograms"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};

  Configurable<float> minpt{"minpt", 0.2, "min pt for ITS-TPC track"};
  Configurable<float> maxpt{"maxpt", 5.0, "max pt for ITS-TPC track"};
  Configurable<float> maxeta{"maxeta", 0.8, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0, "max DCAz in cm"};
  Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};

  // Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  // Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  // Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  // Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  // Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  // Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  // Configurable<float> maxchi2its{"maxchi2its", 36.0, "max. chi2/NclsITS"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (fillQAHistogram) {
      fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
      fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{4000, -20, 20}}, false);
      fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{90, 0, 2 * M_PI}, {80, -2.0f, 2.0f}}, false);
      fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{800, 0, 40}}, false);
      fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    }
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());

      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsChi2NCl() > 36.f) {
      return false;
    }
    if (track.itsNCls() < 4) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < 1) {
      return false;
    }

    if (track.tpcChi2NCl() > 5.f) {
      return false;
    }

    if (track.tpcNClsFound() < 0) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < 50) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }

    if (track.tpcFractionSharedCls() > max_frac_shared_clusters_tpc) {
      return false;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > dca_xy_max || std::fabs(dcaZ) > dca_z_max) {
      return false;
    }
    if (std::fabs(dcaZ) > 3.f) {
      return false;
    }

    if (std::fabs(trackParCov.getEta()) > maxeta || trackParCov.getPt() < minpt || maxpt < trackParCov.getPt()) {
      return false;
    }
    if (trackParCov.getPt() > 5.f) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    float pt = trackParCov.getPt();
    float eta = trackParCov.getEta();
    float phi = trackParCov.getPhi();
    o2::math_utils::bringTo02Pi(phi);
    uint16_t trackBit = 0;

    // As minimal cuts, following cuts are applied. The cut values are hardcoded on the purpose for consistent bit operation.
    // has info on ITS and TPC
    // a hit on ITSib any
    // Ncls ITS >= 4
    // chi2/Ncls ITS < 36
    // Ncr TPC >= 50
    // chi2/Ncls TPC < 5
    // Ncr/Nf ratio in TPC > 0.8

    if (track.itsNCls() >= 5) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNclsITS5);
    }
    if (track.itsNCls() >= 6) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNclsITS6);
    }

    if (track.tpcNClsCrossedRows() >= 70) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNcrTPC70);
    }
    if (track.tpcNClsCrossedRows() >= 90) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNcrTPC90);
    }
    if (track.tpcNClsFound() >= 50) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNclsTPC50);
    }
    if (track.tpcNClsFound() >= 70) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNclsTPC70);
    }
    if (track.tpcNClsFound() >= 90) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kNclsTPC90);
    }
    if (track.tpcChi2NCl() < 4.f) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kChi2TPC4);
    }
    if (track.tpcChi2NCl() < 3.f) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kChi2TPC3);
    }
    if (track.tpcFractionSharedCls() < 0.7) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kFracSharedTPC07);
    }

    if (std::fabs(dcaZ) < 0.5) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kDCAz05cm);
    }
    if (std::fabs(dcaZ) < 0.3) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kDCAz03cm);
    }

    if (std::fabs(dcaXY) < 0.5) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kDCAxy05cm);
    }
    if (std::fabs(dcaXY) < 0.3) {
      trackBit |= static_cast<uint16_t>(RefTrackBit::kDCAxy03cm);
    }

    emprimarytracks(/*collision.globalIndex(),*/ /*track.globalIndex(),*/ track.sign() / pt, eta, phi, trackBit);
    prmtrackeventidtmp(collision.globalIndex());

    if (fillQAHistogram) {
      fRegistry.fill(HIST("Track/hPt"), pt);
      fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / pt);
      fRegistry.fill(HIST("Track/hEtaPhi"), phi, eta);
      fRegistry.fill(HIST("Track/hDCAxyz"), dcaXY, dcaZ);
      fRegistry.fill(HIST("Track/hDCAxyzSigma"), dcaXY / std::sqrt(trackParCov.getSigmaY2()), dcaZ / std::sqrt(trackParCov.getSigmaZ2()));
      fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), pt, std::sqrt(trackParCov.getSigmaY2()) * 1e+4); // convert cm to um
      fRegistry.fill(HIST("Track/hDCAzRes_Pt"), pt, std::sqrt(trackParCov.getSigmaZ2()) * 1e+4);  // convert cm to um
      fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
      fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
      fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
      fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
      fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
      fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
      fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
      fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
      fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    }
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Filter trackFilter = o2::aod::track::itsChi2NCl < 36.f && o2::aod::track::tpcChi2NCl < 5.f && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  // ---------- for data ----------

  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.isEoI()) { // events with at least 1 lepton for data reduction.
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!checkTrack<false>(collision, track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }

    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryTrack, processRec, "process reconstructed info only", true); // standalone

  void processRec_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.isEoI()) { // events with at least 1 lepton for data reduction.
        continue;
      }
      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!checkTrack<false>(collision, track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }

    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryTrack, processRec_SWT, "process reconstructed info only", false); // standalone with swt

  // ---------- for MC ----------

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  void processMC(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.isEoI()) { // events with at least 1 lepton for data reduction.
        continue;
      }

      auto tracks_per_coll = tracks.sliceBy(perCol, collision.globalIndex());
      for (const auto& track : tracks_per_coll) {
        if (!checkTrack<true>(collision, track)) {
          continue;
        }
        fillTrackTable(collision, track);
      }
    } // end of collision loop
  }
  PROCESS_SWITCH(skimmerPrimaryTrack, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryTrack>(cfgc, TaskName{"skimmer-primary-track"})};
}
