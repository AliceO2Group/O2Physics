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

/// \brief write relevant information about primary electrons.
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/ElectronModule.h"
#include "PWGEM/Dilepton/Utils/MlResponseO2Track.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Utils.h>
#include <PID/PIDTOFParamService.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <math.h>

// using namespace o2;
// using namespace o2::soa;
// using namespace o2::framework;
// using namespace o2::framework::expressions;
// using namespace o2::constants::physics;
// using namespace o2::common::core;

struct skimmerPrimaryElectronSCT {

  using MyBCs = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::BcSels>;
  using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::EMEvSels>;
  using MyCollisionsWithSWT = o2::soa::Join<MyCollisions, o2::aod::EMSWTriggerBitsTMP>;

  using MyTracks = o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksExtra, o2::aod::TracksCovIU,
                                 o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                 o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr, o2::aod::pidTOFbeta, o2::aod::TOFSignal, o2::aod::TOFEvTime>;
  using MyTrack = MyTracks::iterator;
  using MyTracksMC = o2::soa::Join<MyTracks, o2::aod::McTrackLabels, o2::aod::mcTPCTuneOnData>;
  using MyTrackMC = MyTracksMC::iterator;

  using MyV0s = o2::soa::Join<o2::aod::V0Datas, o2::aod::V0Covs>;
  using MyCascades = o2::soa::Join<o2::aod::CascDatas, o2::aod::CascCovs>;

  struct : o2::framework::ConfigurableGroup {
    o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    o2::framework::Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    o2::framework::Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    o2::framework::Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfig;

  o2::framework::Configurable<bool> doSCTwithTracks{"doSCTwithTracks", true, "flag to tag electrons with tracks"};
  o2::framework::Configurable<bool> doSCTwithV0s{"doSCTwithV0s", false, "flag to tag electrons with v0s"};
  o2::framework::Configurable<bool> doSCTwithCascades{"doSCTwithCascades", false, "flag to tag electrons with cascades"};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::pid::tof::TOFResponse> mTOFResponse;

  o2::framework::SliceCache cache;
  o2::framework::Preslice<o2::aod::TracksIU> perCol_track = o2::aod::track::collisionId;
  o2::framework::Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;
  o2::framework::Preslice<o2::aod::V0Datas> perCol_v0 = o2::aod::v0data::collisionId;
  o2::framework::Preslice<o2::aod::CascDatas> perCol_casc = o2::aod::cascdata::collisionId;

  void init(o2::framework::InitContext& initContext)
  {
    mRunNumber = 0;
    d_bz = 0;

    LOGF(info, "initializing CCDB");
    ccdb->setURL(ccdbConfig.ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    // ccdbApi.init(ccdbConfig.ccdburl);
    // mTOFResponse->initSetup(ccdb, initContext);

    LOGF(info, "initializing electronModule");
    electronModule.init(cfgEelectronCut, cfgEelectronPFCut, cfgHadronCut, cfgV0Cut, cfgCascadeCut, cfgDFeT, cfgDFeV0, cfgDFeC, initContext, ccdb, mTOFResponse, ccdbConfig.ccdburl.value);
    // electronModule.setTOFResponse(mTOFResponse);
    electronModule.addHistograms(mRegistry);
    // electronModule.doPFB(doPF.value);

    electronModule.doSCTwithTracks(doSCTwithTracks.value);
    electronModule.doSCTwithV0s(doSCTwithV0s.value);
    electronModule.doSCTwithCascades(doSCTwithCascades.value);
  }

  o2::pwgem::dilepton::utils::ElectronModule electronModule;
  o2::pwgem::dilepton::utils::ElectronProducts products;
  o2::pwgem::dilepton::utils::electronCut cfgEelectronCut;
  o2::pwgem::dilepton::utils::electronPFCut cfgEelectronPFCut;
  o2::pwgem::dilepton::utils::hadronCut cfgHadronCut;
  o2::pwgem::dilepton::utils::v0Cut cfgV0Cut;
  o2::pwgem::dilepton::utils::cascadeCut cfgCascadeCut;
  o2::pwgem::dilepton::utils::cfgDFeT cfgDFeT;
  o2::pwgem::dilepton::utils::cfgDFeV0 cfgDFeV0;
  o2::pwgem::dilepton::utils::cfgDFeC cfgDFeC;

  o2::framework::HistogramRegistry mRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  // ---------- for data ----------

  int mRunNumber{0};
  float d_bz{0};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  // const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // load matLUT for this timestamp
    if (!lut) {
      LOG(info) << "Loading material look-up table for timestamp: " << bc.timestamp();
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(ccdbConfig.lutPath, bc.timestamp()));
    } else {
      LOG(info) << "Material look-up table already in place. Not reloading.";
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbConfig.grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      // mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(ccdbConfig.mVtxPath, bc.timestamp());
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfig.grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfig.grpmagPath << " of object GRPMagField and " << ccdbConfig.grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      o2::base::Propagator::Instance()->setMatLUT(lut);
      // mMeanVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(ccdbConfig.mVtxPath, bc.timestamp());

      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }

    mRunNumber = bc.runNumber();
  }

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cut), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  o2::framework::expressions::Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > cfgV0Cut.cfg_min_cospa&& o2::aod::v0data::dcaV0daughters<cfgV0Cut.cfg_max_dca2legs && nabs(o2::aod::v0data::dcanegtopv)> cfgV0Cut.cfg_min_dcaxy&& nabs(o2::aod::v0data::dcanegtopv) > cfgV0Cut.cfg_min_dcaxy;
  using filteredMyV0s = o2::soa::Filtered<MyV0s>;

  o2::framework::expressions::Filter cascadeFilter = nabs(o2::aod::cascdata::dcanegtopv) > cfgCascadeCut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcanegtopv) > cfgCascadeCut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcabachtopv) > cfgCascadeCut.cfg_min_dcaxy_bachelor&& o2::aod::cascdata::dcacascdaughters < cfgCascadeCut.cfg_max_dcadau&& o2::aod::cascdata::dcaV0daughters < cfgCascadeCut.cfg_max_dcadau_v0;
  using filteredMyCascades = o2::soa::Filtered<MyCascades>;

  void processRec_SA(MyCollisions const& collisions, MyBCs const& bcs, MyTracks const& tracks, filteredMyV0s const& v0s, filteredMyCascades const& cascades)
  {
    initCCDB(bcs.begin());
    electronModule.processWithoutTTCA<false, false>(bcs, collisions, tracks, v0s, cascades, nullptr, products, mRegistry);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processRec_SA, "process reconstructed info only", true); // standalone

  void processRec_TTCA(MyCollisions const& collisions, MyBCs const& bcs, MyTracks const& tracks, o2::aod::TrackAssoc const& trackIndices, filteredMyV0s const& v0s, filteredMyCascades const& cascades)
  {
    initCCDB(bcs.begin());
    electronModule.processWithTTCA<false, false>(bcs, collisions, tracks, v0s, cascades, trackIndices, nullptr, products, mRegistry, cache, perCol_track, trackIndicesPerCollision, perCol_v0, perCol_casc);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processRec_TTCA, "process reconstructed info only", false); // with TTCA

  void processRec_SA_SWT(MyCollisionsWithSWT const& collisions, MyBCs const& bcs, MyTracks const& tracks, filteredMyV0s const& v0s, filteredMyCascades const& cascades)
  {
    initCCDB(bcs.begin());
    electronModule.processWithoutTTCA<false, true>(bcs, collisions, tracks, v0s, cascades, nullptr, products, mRegistry);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processRec_SA_SWT, "process reconstructed info only", false); // standalone with swt

  void processRec_TTCA_SWT(MyCollisionsWithSWT const& collisions, MyBCs const& bcs, MyTracks const& tracks, o2::aod::TrackAssoc const& trackIndices, filteredMyV0s const& v0s, filteredMyCascades const& cascades)
  {
    initCCDB(bcs.begin());
    electronModule.processWithTTCA<false, true>(bcs, collisions, tracks, v0s, cascades, trackIndices, nullptr, products, mRegistry, cache, perCol_track, trackIndicesPerCollision, perCol_v0, perCol_casc);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processRec_TTCA_SWT, "process reconstructed info only", false); // with TTCA with swt

  // ---------- for MC ----------

  void processMC_SA(o2::soa::Join<MyCollisions, o2::aod::McCollisionLabels> const& collisions, MyBCs const& bcs, MyTracksMC const& tracks, filteredMyV0s const& v0s, filteredMyCascades const& cascades, o2::aod::McParticles const& mcParticles)
  {
    initCCDB(bcs.begin());
    electronModule.processWithoutTTCA<true, false>(bcs, collisions, tracks, v0s, cascades, mcParticles, products, mRegistry);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processMC_SA, "process reconstructed and MC info ", false); // without TTCA in MC

  void processMC_TTCA(o2::soa::Join<MyCollisions, o2::aod::McCollisionLabels> const& collisions, MyBCs const& bcs, MyTracksMC const& tracks, o2::aod::TrackAssoc const& trackIndices, filteredMyV0s const& v0s, filteredMyCascades const& cascades, o2::aod::McParticles const& mcParticles)
  {
    initCCDB(bcs.begin());
    electronModule.processWithTTCA<true, false>(bcs, collisions, tracks, v0s, cascades, trackIndices, mcParticles, products, mRegistry, cache, perCol_track, trackIndicesPerCollision, perCol_v0, perCol_casc);
  }
  PROCESS_SWITCH(skimmerPrimaryElectronSCT, processMC_TTCA, "process reconstructed info only", false); // with TTCA in MC
};

struct associateAmbiguousElectron {
  o2::framework::Produces<o2::aod::EMAmbiguousElectronSelfIds> em_amb_ele_ids;

  o2::framework::SliceCache cache;
  o2::framework::PresliceUnsorted<o2::aod::EMPrimaryElectrons> perTrack = o2::aod::emprimaryelectron::trackId;
  std::vector<int> ambele_self_Ids;

  void process(o2::aod::EMPrimaryElectrons const& electrons)
  {
    for (const auto& electron : electrons) {
      auto electrons_with_same_trackId = electrons.sliceBy(perTrack, electron.trackId());
      ambele_self_Ids.reserve(electrons_with_same_trackId.size());
      for (const auto& amb_ele : electrons_with_same_trackId) {
        if (amb_ele.globalIndex() == electron.globalIndex()) { // don't store myself.
          continue;
        }
        ambele_self_Ids.emplace_back(amb_ele.globalIndex());
      }
      em_amb_ele_ids(ambele_self_Ids);
      ambele_self_Ids.clear();
      ambele_self_Ids.shrink_to_fit();
    }
  }
};
o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  return o2::framework::WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectronSCT>(cfgc, o2::framework::TaskName{"skimmer-primary-electron-sct"}),
    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, o2::framework::TaskName{"associate-ambiguous-electron"})};
}
