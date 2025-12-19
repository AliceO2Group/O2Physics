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
#include "PWGEM/Dilepton/Utils/MlResponseO2Track.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
// #include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/ML/MlResponse.h"

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

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EMEvSels>;
using MyCollisionsWithSWT = soa::Join<MyCollisions, aod::EMSWTriggerBitsTMP>;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                           aod::pidTPCFullEl, /*aod::pidTPCFullMu,*/ aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, /*aod::pidTOFFullMu,*/ aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels, aod::mcTPCTuneOnData>;
using MyTrackMC = MyTracksMC::iterator;

struct skimmerPrimaryElectronQC {
  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryElectrons> emprimaryelectrons;
  Produces<aod::EMPrimaryElectronsCov> emprimaryelectronscov;

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

  struct : ConfigurableGroup {
    std::string prefix = "trackcut";
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> mincrossedrows{"mincrossedrows", 40, "min. crossed rows"};
    Configurable<int> min_ncluster_its{"min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
    Configurable<float> maxchi2its{"maxchi2its", 36.0, "max. chi2/NclsITS"};
    Configurable<float> minchi2its{"minchi2its", -1e+10, "min. chi2/NclsITS"};
    Configurable<float> minpt{"minpt", 0.05, "min pt for ITS-TPC track"};
    Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
    Configurable<float> dca_xy_max{"dca_xy_max", 1.0, "max DCAxy in cm"};
    Configurable<float> dca_z_max{"dca_z_max", 1.0, "max DCAz in cm"};
    Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to include ITSsa tracks"};
    Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -3.5, "min. TPC n sigma for electron inclusion"};
    Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", +3.5, "max. TPC n sigma for electron inclusion"};
  } trackcut;

  struct : ConfigurableGroup {
    std::string prefix = "tighttrackcut";
    Configurable<int> min_ncluster_tpc_pid{"min_ncluster_tpc_pid", 60, "min ncluster tpc used for PID"};
    Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> mincrossedrows{"mincrossedrows", 100, "min. crossed rows"};
    Configurable<int> min_ncluster_its{"min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 3, "min ncluster itsib"};
    Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max. chi2/NclsTPC"};
    Configurable<float> maxchi2its{"maxchi2its", 5.0, "max. chi2/NclsITS"};
    Configurable<float> dca_xy_max{"dca_xy_max", 0.2, "max DCAxy in cm"};
    Configurable<float> dca_z_max{"dca_z_max", 0.2, "max DCAz in cm"};
    Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"}; // Don't change.
    Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", +2.0, "max. TPC n sigma for electron inclusion"}; // Don't change.
    Configurable<float> minTOFNsigmaEl{"minTOFNsigmaEl", -2.0, "min. TOF n sigma for electron inclusion"}; // Don't change.
    Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", +2.0, "max. TOF n sigma for electron inclusion"}; // Don't change.
  } tighttrackcut;

  Configurable<bool> storeOnlyTrueElectronMC{"storeOnlyTrueElectronMC", false, "Flag to store only true electron in MC"};
  Configurable<float> maxMee{"maxMee", 0.005, "max mee for pi0 -> ee"};
  Configurable<float> maxPhiV{"maxPhiV", M_PI / 2, "max phiv for pi0 -> ee"};

  // configuration for PID ML
  Configurable<bool> usePIDML{"usePIDML", false, "Flag to use PID ML"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
  Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
  Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{-999999., 999999.}, "Bin limits for ML application"};
  Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.95}, "ML cuts per bin"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature"}, "Names of ML model input features"};
  Configurable<std::string> nameBinningFeature{"nameBinningFeature", "pt", "Names of ML model binning feature"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB.  Exceptions: > 0 for the specific timestamp, 0 gets the run dependent timestamp"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};
  o2::analysis::MlResponseO2Track<float> mlResponseSingleTrack;

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  o2::dataformats::VertexBase mVtx;
  o2::dataformats::DCA mDcaInfoCov;

  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::ccdb::CcdbApi ccdbApi;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);

    if (fillQAHistogram) {
      fRegistry.add("Track/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{1000, 0.0f, 10}}, false);
      fRegistry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", kTH1F, {{4000, -20, 20}}, false);
      fRegistry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", kTH2F, {{180, 0, 2 * M_PI}, {20, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
      fRegistry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
      fRegistry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
      fRegistry.add("Track/hNclsTPC", "number of TPC clusters", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hNcrTPC", "number of TPC crossed rows", kTH1F, {{161, -0.5, 160.5}}, false);
      fRegistry.add("Track/hChi2TPC", "chi2/number of TPC clusters", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hChi2TOF", "chi2 of TOF", kTH1F, {{100, 0, 10}}, false);
      fRegistry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", kTH1F, {{200, 0, 2}}, false);
      fRegistry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Track/hNclsITS", "number of ITS clusters", kTH1F, {{8, -0.5, 7.5}}, false);
      fRegistry.add("Track/hChi2ITS", "chi2/number of ITS clusters", kTH1F, {{400, 0, 40}}, false);
      fRegistry.add("Track/hITSClusterMap", "ITS cluster map", kTH1F, {{128, -0.5, 127.5}}, false);
      fRegistry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTPCdEdxMC", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
      fRegistry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
      fRegistry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
      fRegistry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<ITS cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
      fRegistry.add("Track/hMeanClusterSizeITSib", "mean cluster size ITSib;p_{pv} (GeV/c);<ITSib cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
      fRegistry.add("Track/hMeanClusterSizeITSob", "mean cluster size ITSob;p_{pv} (GeV/c);<ITSob cluster size> #times cos(#lambda)", kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
      fRegistry.add("Track/hProbElBDT", "probability to be e from BDT;p_{in} (GeV/c);BDT score;", kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
      fRegistry.add("Pair/hMvsPhiV", "m_{ee} vs. #varphi_{V} ULS;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", kTH2F, {{90, 0.f, M_PI}, {100, 0, 0.1}});
    }

    if (usePIDML) {
      static constexpr int nClassesMl = 2;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      const std::vector<std::string> labelsClasses = {"Background", "Signal"};
      const uint32_t nBinsMl = binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = 0.0;
        cutsMlArr[i][1] = cutsMl.value[i];
      }
      o2::framework::LabeledArray<double> cutsMl = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponseSingleTrack.configure(binsMl.value, cutsMl, cutDirMl, nClassesMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdburl);
        mlResponseSingleTrack.setModelPathsCCDB(onnxFileNames.value, ccdbApi, onnxPathsCCDB.value, timestampCCDB.value);
      } else {
        mlResponseSingleTrack.setModelPathsLocal(onnxFileNames.value);
      }
      mlResponseSingleTrack.cacheInputFeaturesIndices(namesInputFeatures);
      mlResponseSingleTrack.cacheBinningIndex(nameBinningFeature);
      mlResponseSingleTrack.init(enableOptimizations.value);
    } // end of PID ML
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
      if (storeOnlyTrueElectronMC) {
        const auto& mcParticle = track.template mcParticle_as<aod::McParticles>();
        if (std::abs(mcParticle.pdgCode()) != 11) {
          return false;
        }
      }
    }

    if (!track.hasITS()) {
      return false;
    }

    if (track.itsChi2NCl() < trackcut.minchi2its || trackcut.maxchi2its < track.itsChi2NCl()) { // accept ITS afterburner (itsChi2NCl = -999)
      return false;
    }
    if (track.itsNCls() < trackcut.min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < trackcut.min_ncluster_itsib) {
      return false;
    }

    if (!trackcut.includeITSsa && (!track.hasITS() || !track.hasTPC())) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcChi2NCl() < 0.f || trackcut.maxchi2tpc < track.tpcChi2NCl()) {
        return false;
      }

      if (track.tpcNClsFound() < trackcut.min_ncluster_tpc) {
        return false;
      }

      if (track.tpcNClsCrossedRows() < trackcut.mincrossedrows) {
        return false;
      }

      if (track.tpcCrossedRowsOverFindableCls() < trackcut.min_tpc_cr_findable_ratio) {
        return false;
      }

      if (track.tpcFractionSharedCls() > trackcut.max_frac_shared_clusters_tpc) {
        return false;
      }
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > trackcut.dca_xy_max || std::fabs(dcaZ) > trackcut.dca_z_max) {
      return false;
    }

    if (trackParCov.getPt() < trackcut.minpt || std::fabs(trackParCov.getEta()) > trackcut.maxeta) {
      return false;
    }

    return true;
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrackTight(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
      if (storeOnlyTrueElectronMC) {
        const auto& mcParticle = track.template mcParticle_as<aod::McParticles>();
        if (std::abs(mcParticle.pdgCode()) != 11) {
          return false;
        }
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || tighttrackcut.maxchi2its < track.itsChi2NCl()) {
      return false;
    }
    if (track.itsNCls() < tighttrackcut.min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < tighttrackcut.min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || tighttrackcut.maxchi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < tighttrackcut.min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < tighttrackcut.mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < tighttrackcut.min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > tighttrackcut.max_frac_shared_clusters_tpc) {
      return false;
    }

    if (track.tpcNClsPID() < tighttrackcut.min_ncluster_tpc_pid) {
      return false;
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > tighttrackcut.dca_xy_max || std::fabs(dcaZ) > tighttrackcut.dca_z_max) {
      return false;
    }

    if (trackParCov.getPt() < trackcut.minpt || std::fabs(trackParCov.getEta()) > trackcut.maxeta) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    if (track.hasTPC()) {
      if (trackcut.minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < trackcut.maxTPCNsigmaEl) {
        return true;
      } else {
        return false;
      }
    } else { // accept ITSsa too
      return true;
    }
  }

  template <typename TTrack>
  bool isElectronTight(TTrack const& track)
  {
    bool is_El_TPC = tighttrackcut.minTPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < tighttrackcut.maxTPCNsigmaEl;
    bool is_El_TOF = tighttrackcut.minTOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < tighttrackcut.maxTOFNsigmaEl;
    return is_El_TPC && is_El_TOF;
  }

  template <bool isMC, typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), std::pair<int, int>{collision.globalIndex(), track.globalIndex()}) == stored_trackIds.end()) {
      mDcaInfoCov.set(999, 999, 999, 999, 999);
      auto trackParCov = getTrackParCov(track);
      trackParCov.setPID(o2::track::PID::Electron);
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
      if (!isPropOK) {
        return;
      }
      float dcaXY = mDcaInfoCov.getY();
      float dcaZ = mDcaInfoCov.getZ();

      float pt_recalc = trackParCov.getPt();
      float eta_recalc = trackParCov.getEta();
      float phi_recalc = trackParCov.getPhi();
      o2::math_utils::bringTo02Pi(phi_recalc);

      bool isAssociatedToMPC = collision.globalIndex() == track.collisionId();
      float mcTunedTPCSignal = 0.f;
      if constexpr (isMC) {
        mcTunedTPCSignal = track.mcTunedTPCSignal();
      }

      float probaEl = 1.0;
      if (usePIDML) {
        std::vector<float> inputFeatures = mlResponseSingleTrack.getInputFeatures(track, trackParCov, collision);
        float binningFeature = mlResponseSingleTrack.getBinningFeature(track, trackParCov, collision);

        int pbin = lower_bound(binsMl.value.begin(), binsMl.value.end(), binningFeature) - binsMl.value.begin() - 1;
        if (pbin < 0) {
          pbin = 0;
        } else if (static_cast<int>(binsMl.value.size()) - 2 < pbin) {
          pbin = static_cast<int>(binsMl.value.size()) - 2;
        }
        probaEl = mlResponseSingleTrack.getModelOutput(inputFeatures, pbin)[1]; // 0: hadron, 1:electron
      }

      emprimaryelectrons(collision.globalIndex(), track.globalIndex(), track.sign(),
                         pt_recalc, eta_recalc, phi_recalc,
                         dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(),
                         track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusPID(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                         track.tpcChi2NCl(), track.tpcInnerParam(),
                         track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                         track.beta(), track.tofNSigmaEl(), /*track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),*/
                         track.itsClusterSizes(),
                         track.itsChi2NCl(), track.tofChi2(), track.detectorMap(),
                         // trackParCov.getTgl(),
                         isAssociatedToMPC, false, probaEl, mcTunedTPCSignal);

      emprimaryelectronscov(
        trackParCov.getX(),
        trackParCov.getAlpha(),
        trackParCov.getY(),
        trackParCov.getZ(),
        trackParCov.getSnp(),
        // trackParCov.getTgl(),
        // trackParCov.getSigmaY2(),
        // trackParCov.getSigmaZY(),
        // trackParCov.getSigmaZ2(),
        trackParCov.getSigmaSnpY(),
        trackParCov.getSigmaSnpZ(),
        trackParCov.getSigmaSnp2(),
        trackParCov.getSigmaTglY(),
        trackParCov.getSigmaTglZ(),
        trackParCov.getSigmaTglSnp(),
        trackParCov.getSigmaTgl2(),
        trackParCov.getSigma1PtY(),
        trackParCov.getSigma1PtZ(),
        trackParCov.getSigma1PtSnp(),
        trackParCov.getSigma1PtTgl(),
        trackParCov.getSigma1Pt2());

      stored_trackIds.emplace_back(std::pair<int, int>{collision.globalIndex(), track.globalIndex()});

      if (fillQAHistogram) {
        int total_cluster_size = 0, nl = 0;
        for (unsigned int layer = 0; layer < 7; layer++) {
          int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
          if (cluster_size_per_layer > 0) {
            nl++;
          }
          total_cluster_size += cluster_size_per_layer;
        }

        int total_cluster_size_ib = 0, nl_ib = 0;
        for (unsigned int layer = 0; layer < 3; layer++) {
          int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
          if (cluster_size_per_layer > 0) {
            nl_ib++;
          }
          total_cluster_size_ib += cluster_size_per_layer;
        }

        int total_cluster_size_ob = 0, nl_ob = 0;
        for (unsigned int layer = 3; layer < 7; layer++) {
          int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
          if (cluster_size_per_layer > 0) {
            nl_ob++;
          }
          total_cluster_size_ob += cluster_size_per_layer;
        }

        fRegistry.fill(HIST("Track/hPt"), pt_recalc);
        fRegistry.fill(HIST("Track/hQoverPt"), track.sign() / pt_recalc);
        fRegistry.fill(HIST("Track/hEtaPhi"), phi_recalc, eta_recalc);
        fRegistry.fill(HIST("Track/hDCAxyz"), dcaXY, dcaZ);
        fRegistry.fill(HIST("Track/hDCAxyzSigma"), dcaXY / std::sqrt(trackParCov.getSigmaY2()), dcaZ / std::sqrt(trackParCov.getSigmaZ2()));
        fRegistry.fill(HIST("Track/hDCAxyRes_Pt"), pt_recalc, std::sqrt(trackParCov.getSigmaY2()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Track/hDCAzRes_Pt"), pt_recalc, std::sqrt(trackParCov.getSigmaZ2()) * 1e+4);  // convert cm to um
        fRegistry.fill(HIST("Track/hNclsITS"), track.itsNCls());
        fRegistry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
        fRegistry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
        fRegistry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
        fRegistry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
        fRegistry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
        fRegistry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
        fRegistry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
        fRegistry.fill(HIST("Track/hChi2TOF"), track.tofChi2());
        fRegistry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
        fRegistry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
        fRegistry.fill(HIST("Track/hTPCdEdxMC"), track.tpcInnerParam(), mcTunedTPCSignal);
        fRegistry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
        fRegistry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        fRegistry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
        fRegistry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
        fRegistry.fill(HIST("Track/hTOFbeta"), trackParCov.getP(), track.beta());
        fRegistry.fill(HIST("Track/hTOFNsigmaEl"), track.tpcInnerParam(), track.tofNSigmaEl());
        fRegistry.fill(HIST("Track/hTOFNsigmaPi"), track.tpcInnerParam(), track.tofNSigmaPi());
        fRegistry.fill(HIST("Track/hTOFNsigmaKa"), track.tpcInnerParam(), track.tofNSigmaKa());
        fRegistry.fill(HIST("Track/hTOFNsigmaPr"), track.tpcInnerParam(), track.tofNSigmaPr());
        fRegistry.fill(HIST("Track/hMeanClusterSizeITS"), trackParCov.getP(), static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(trackParCov.getTgl())));
        fRegistry.fill(HIST("Track/hMeanClusterSizeITSib"), trackParCov.getP(), static_cast<float>(total_cluster_size_ib) / static_cast<float>(nl_ib) * std::cos(std::atan(trackParCov.getTgl())));
        fRegistry.fill(HIST("Track/hMeanClusterSizeITSob"), trackParCov.getP(), static_cast<float>(total_cluster_size_ob) / static_cast<float>(nl_ob) * std::cos(std::atan(trackParCov.getTgl())));
        fRegistry.fill(HIST("Track/hProbElBDT"), track.tpcInnerParam(), probaEl);
      }
    }
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool isDielectronFromPi0(TCollision const& collision, TTrack const& t1, TTrack const& t2)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    std::array<float, 3> pVec1 = {0, 0, 0}; // px, py, pz
    std::array<float, 3> pVec2 = {0, 0, 0}; // px, py, pz
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

    auto t1ParCov = getTrackParCov(t1);
    t1ParCov.setPID(o2::track::PID::Electron);
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, t1ParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    getPxPyPz(t1ParCov, pVec1);

    auto t2ParCov = getTrackParCov(t2);
    t2ParCov.setPID(o2::track::PID::Electron);
    isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, t2ParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    getPxPyPz(t2ParCov, pVec2);

    ROOT::Math::PtEtaPhiMVector v1(t1ParCov.getPt(), t1ParCov.getEta(), t1ParCov.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(t2ParCov.getPt(), t2ParCov.getEta(), t2ParCov.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float mee = v12.M();
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pVec1[0], pVec1[1], pVec1[2], pVec2[0], pVec2[1], pVec2[2], t1.sign(), t2.sign(), d_bz);

    if (fillQAHistogram) {
      fRegistry.fill(HIST("Pair/hMvsPhiV"), phiv, mee);
    }
    if (mee < maxMee && phiv < maxPhiV) {
      return true;
    } else {
      return false;
    }
  }

  std::vector<std::pair<int, int>> stored_trackIds;
  Filter trackFilter = trackcut.minpt < o2::aod::track::pt && nabs(o2::aod::track::eta) < trackcut.maxeta && o2::aod::track::itsChi2NCl < trackcut.maxchi2its && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  // ---------- for data ----------

  void processRec(MyCollisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      const auto& posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (const auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if ((checkTrackTight<false>(collision, pos) && isElectronTight(pos)) && (checkTrack<false>(collision, ele) && isElectron(ele))) {
          if (isDielectronFromPi0<false>(collision, pos, ele)) {
            fillTrackTable<false>(collision, ele);
          }
        }
        if ((checkTrackTight<false>(collision, ele) && isElectronTight(ele)) && (checkTrack<false>(collision, pos) && isElectron(pos))) {
          if (isDielectronFromPi0<false>(collision, pos, ele)) {
            fillTrackTable<false>(collision, pos);
          }
        }

      } // end of ULS pairing
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronQC, processRec, "process reconstructed info only", true); // standalone

  void processRec_SWT(MyCollisionsWithSWT const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }

      if (collision.swtaliastmp_raw() == 0) {
        continue;
      }

      const auto& posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (const auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if ((checkTrackTight<false>(collision, pos) && isElectronTight(pos)) && (checkTrack<false>(collision, ele) && isElectron(ele))) {
          if (isDielectronFromPi0<false>(collision, pos, ele)) {
            fillTrackTable<false>(collision, ele);
          }
        }
        if ((checkTrackTight<false>(collision, ele) && isElectronTight(ele)) && (checkTrack<false>(collision, pos) && isElectron(pos))) {
          if (isDielectronFromPi0<false>(collision, pos, ele)) {
            fillTrackTable<false>(collision, pos);
          }
        }

      } // end of ULS pairing

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronQC, processRec_SWT, "process reconstructed info only", false); // standalone with swt

  // ---------- for MC ----------
  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC(soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const&)
  {
    stored_trackIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collision.isSelected()) {
        continue;
      }
      const auto& posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      const auto& negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);

      for (const auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if ((checkTrackTight<true>(collision, pos) && isElectronTight(pos)) && (checkTrack<true>(collision, ele) && isElectron(ele))) {
          if (isDielectronFromPi0<true>(collision, pos, ele)) {
            fillTrackTable<true>(collision, ele);
          }
        }
        if ((checkTrackTight<true>(collision, ele) && isElectronTight(ele)) && (checkTrack<true>(collision, pos) && isElectron(pos))) {
          if (isDielectronFromPi0<true>(collision, pos, ele)) {
            fillTrackTable<true>(collision, pos);
          }
        }

      } // end of ULS pairing

    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryElectronQC, processMC, "process reconstructed and MC info ", false);
};

struct prefilterPrimaryElectron {
  Produces<aod::EMPrimaryElectronsPrefilterBit> ele_pfb;

  void init(InitContext&) {}

  void process(aod::EMPrimaryElectrons const& primaryelectrons)
  {
    for (int i = 0; i < primaryelectrons.size(); i++) {
      ele_pfb(0);
    }
  }
};

struct associateAmbiguousElectron {
  Produces<aod::EMAmbiguousElectronSelfIds> em_amb_ele_ids;

  SliceCache cache;
  PresliceUnsorted<aod::EMPrimaryElectrons> perTrack = o2::aod::emprimaryelectron::trackId;
  std::vector<int> ambele_self_Ids;

  void process(aod::EMPrimaryElectrons const& electrons)
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
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<skimmerPrimaryElectronQC>(cfgc, TaskName{"skimmer-primary-electron-qc"}),
    adaptAnalysisTask<prefilterPrimaryElectron>(cfgc, TaskName{"prefilter-primary-electron"}),
    adaptAnalysisTask<associateAmbiguousElectron>(cfgc, TaskName{"associate-ambiguous-electron"})};
}
