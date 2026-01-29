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

/// \brief write relevant information about MFTsa tracks for reference flow.
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/fwdtrackUtilities.h"
// #include "Common/DataModel/CollisionAssociationTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::fwdtrackutils;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;

struct skimmerPrimaryMFTTrack {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels, aod::EMEoIs>;
  using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;
  using MyCollisionsMC = soa::Join<MyCollisions, aod::McCollisionLabels>;

  using MyTracks = aod::MFTTracks;
  using MyTrack = MyTracks::iterator;
  using MyTracksMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;
  using MyTrackMC = MyTracksMC::iterator;

  SliceCache cache;
  Preslice<aod::MFTTracks> perCol = o2::aod::fwdtrack::collisionId;
  Produces<aod::EMPrimaryTracks> emprimarytracks;

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<bool> fillQAHistogram{"fillQAHistogram", true, "flag to fill QA histograms"};

  Configurable<float> cfgPtMin{"cfgPtMin", 0.2, "min pt for MFTsa track"};
  Configurable<float> cfgPtMax{"cfgPtMax", 1e+10, "max pt for MFTsa track"};
  Configurable<float> cfgEtaMin{"cfgEtaMin", -4, "min eta acceptance"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", -2, "max eta acceptance"};

  // for z shift for propagation
  Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift from CCDB"};
  Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  Configurable<float> cfgManualZShift{"cfgManualZShift", 0, "manual z-shift for propagation of global muon to PV"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int mRunNumber = 0;
  float mBz = 0;
  float mZShift = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  void init(InitContext&)
  {
    mRunNumber = 0;
    mBz = 0;
    mZShift = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    if (fillQAHistogram) {
      fRegistry.add("MFT/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{100, 0.0f, 10}}, false);
      fRegistry.add("MFT/hPtEta", "pT vs. #eta;#eta;p_{T} (GeV/c)", kTH2F, {{200, -4.f, -2.f}, {100, 0.0f, 10}}, false);
      fRegistry.add("MFT/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{2000, -10, 10}}, false);
      fRegistry.add("MFT/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {200, -4.0f, -2.0f}}, false);
      fRegistry.add("MFT/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", kTH2F, {{200, -0.1f, 0.1f}, {200, -0.1f, 0.1f}}, false);
      fRegistry.add("MFT/hDCAxy", "DCA xy;DCA_{xy} (cm)", kTH1F, {{100, 0, 0.1f}}, false);
      fRegistry.add("MFT/hDCAz", "DCA z;DCA_{z} (cm)", kTH1F, {{100, -0.05f, 0.05f}}, false);
      fRegistry.add("MFT/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{100, 0, 0.1f}, {100, -0.05f, 0.05f}}, false);
      fRegistry.add("MFT/hNclsMFT", "number of MFT clusters", kTH1F, {{11, -0.5, 10.5}}, false);
      fRegistry.add("MFT/hChi2MFT", "chi2/ndf;chi2/ndf", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("MFT/hDCAx_PosZ", "DCAx vs. posZ;Z_{vtx} (cm);DCA_{x} (cm)", kTH2F, {{200, -10, +10}, {200, -0.1, +0.1}}, false);
      fRegistry.add("MFT/hDCAy_PosZ", "DCAy vs. posZ;Z_{vtx} (cm);DCA_{y} (cm)", kTH2F, {{200, -10, +10}, {200, -0.1, +0.1}}, false);
      fRegistry.add("MFT/hDCAx_Phi", "DCAx vs. #varphi;#varphi (rad.);DCA_{x} (cm)", kTH2F, {{180, 0, 2 * M_PI}, {200, -0.1, +0.1}}, false);
      fRegistry.add("MFT/hDCAy_Phi", "DCAy vs. #varphi;#varphi (rad.);DCA_{y} (cm)", kTH2F, {{180, 0, 2 * M_PI}, {200, -0.1, +0.1}}, false);
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOGF(info, "mRunNumber = %d", mRunNumber);
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    mBz = field->getBz(centerMFT); // Get field at centre of MFT
    LOGF(info, "Bz at center of MFT = %f kZG", mBz);

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShift = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShift = 0;
      }
    } else {
      LOGF(info, "z shift is manually set to %f cm", cfgManualZShift.value);
      mZShift = cfgManualZShift;
    }
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const&, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    return true;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& mfttrack)
  {
    // propagate MFTsa track to PV
    std::array<double, 3> dcaInfOrig{999.f, 999.f, 999.f};
    std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
    SMatrix55 tcovs(v1.begin(), v1.end());
    SMatrix5 tpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
    o2::track::TrackParCovFwd trackPar{mfttrack.z() + mZShift, tpars, tcovs, mfttrack.chi2()};
    trackPar.propagateToDCAhelix(mBz, {collision.posX(), collision.posY(), collision.posZ()}, dcaInfOrig);
    v1.clear();
    v1.shrink_to_fit();
    float dcaX = dcaInfOrig[0];
    float dcaY = dcaInfOrig[1];
    float dcaZ = dcaInfOrig[2];
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

    float pt = trackPar.getPt();
    float eta = trackPar.getEta();
    float phi = trackPar.getPhi();
    o2::math_utils::bringTo02Pi(phi);

    if (pt < cfgPtMin || cfgPtMax < pt) {
      return;
    }
    if (eta < cfgEtaMin || cfgEtaMax < eta) {
      return;
    }

    uint16_t trackBit = 0;
    int ndf = 2 * mfttrack.nClusters() - 5;

    // As minimal cuts, following cuts are applied. The cut values are hardcoded on the purpose for consistent bit operation.
    // Ncls MFT >= 5
    // chi2/ndf MFT < 5
    // |dcaXY| < 0.06 cm

    if (mfttrack.nClusters() < 6 || mfttrack.chi2() / ndf > 5.f || std::fabs(dcaXY) > 0.05) {
      return;
    }

    if (mfttrack.nClusters() >= 7) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kNclsMFT7);
    }
    if (mfttrack.nClusters() >= 8) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kNclsMFT8);
    }

    if (mfttrack.chi2() / ndf < 4.f) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kChi2MFT4);
    }
    if (mfttrack.chi2() / ndf < 3.f) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kChi2MFT3);
    }

    if (std::fabs(dcaXY) < 0.04) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kDCAxy004cm);
    }
    if (std::fabs(dcaXY) < 0.03) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kDCAxy003cm);
    }
    if (std::fabs(dcaXY) < 0.02) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kDCAxy002cm);
    }
    if (std::fabs(dcaXY) < 0.01) {
      trackBit |= static_cast<uint16_t>(RefMFTTrackBit::kDCAxy001cm);
    }

    emprimarytracks(collision.globalIndex(), mfttrack.globalIndex(), mfttrack.sign() / pt, eta, phi, trackBit);

    if (fillQAHistogram) {
      fRegistry.fill(HIST("MFT/hPt"), pt);
      fRegistry.fill(HIST("MFT/hPtEta"), eta, pt);
      fRegistry.fill(HIST("MFT/hQoverPt"), mfttrack.sign() / pt);
      fRegistry.fill(HIST("MFT/hEtaPhi"), phi, eta);
      fRegistry.fill(HIST("MFT/hDCAxy2D"), dcaX, dcaY);
      fRegistry.fill(HIST("MFT/hDCAxy"), dcaXY);
      fRegistry.fill(HIST("MFT/hDCAz"), dcaZ);
      fRegistry.fill(HIST("MFT/hDCAxyz"), dcaXY, dcaZ);
      fRegistry.fill(HIST("MFT/hNclsMFT"), mfttrack.nClusters());
      fRegistry.fill(HIST("MFT/hChi2MFT"), mfttrack.chi2() / ndf);
      fRegistry.fill(HIST("MFT/hDCAx_PosZ"), collision.posZ(), dcaX);
      fRegistry.fill(HIST("MFT/hDCAy_PosZ"), collision.posZ(), dcaY);
      fRegistry.fill(HIST("MFT/hDCAx_Phi"), phi, dcaX);
      fRegistry.fill(HIST("MFT/hDCAy_Phi"), phi, dcaY);
    }
  }

  // ---------- for data ----------

  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks)
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
  PROCESS_SWITCH(skimmerPrimaryMFTTrack, processRec, "process reconstructed info only", true); // standalone

  void processRec_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks)
  {

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }
      if (!collision.isEoI()) {
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
  PROCESS_SWITCH(skimmerPrimaryMFTTrack, processRec_SWT, "process reconstructed info only", false); // standalone with swt

  // ---------- for MC ----------

  void processMC(MyCollisionsMC const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyTracksMC const& tracks)
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
  PROCESS_SWITCH(skimmerPrimaryMFTTrack, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryMFTTrack>(cfgc, TaskName{"skimmer-primary-mfttrack"})};
}
