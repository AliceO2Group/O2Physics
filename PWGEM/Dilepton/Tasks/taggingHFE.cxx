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

/// \file taggingHFE.cxx
/// \brief a task to study tagging e from charm hadron decays in MC
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

struct taggingHFE {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using MyCollisionsWithMCLabel = soa::Join<MyCollisions, aod::McCollisionLabels>;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
  using MyTracksWithMCLabel = soa::Join<MyTracks, aod::McTrackLabels>;

  using MyV0s = soa::Join<aod::V0Datas, aod::V0Covs>;
  using MyCascades = soa::Join<aod::CascDatas, aod::CascCovs>;

  struct EBPair { // electron-baryon pair
    float mass{-999.f};
    float dca2legs{-999.f};
    float cospa{-999.f};
    float lxy{-999.f};
    float lz{-999.f};
    float ptepv{-999.f};
    float dca3dinsigma{-999.f};
  };

  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};

  struct : ConfigurableGroup {
    std::string prefix = "electroncut";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 3, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2, "min n sigma el in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3, "max n sigma el in TPC"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3, "min n sigma el in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3, "max n sigma el in TOF"};
  } electroncut;

  struct : ConfigurableGroup {
    std::string prefix = "loose_electroncut";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -1.2, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +1.2, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 40, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
  } loose_electroncut;

  struct : ConfigurableGroup {
    std::string prefix = "kaoncut";
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -3, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +3, "max n sigma ka in TOF"};
  } kaoncut;

  struct : ConfigurableGroup {
    std::string prefix = "v0cut";
    Configurable<float> cfg_min_mass_k0s_veto{"cfg_min_mass_k0s_veto", 0.47, "min mass for K0S veto"};
    Configurable<float> cfg_max_mass_k0s_veto{"cfg_max_mass_k0s_veto", 0.52, "max mass for K0S veto"};
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.113, "min mass for Lambda"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.118, "max mass for Lambda"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min cospa for v0hadron"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 0.1, "max distance between 2 legs for v0hadron"};
    // Configurable<float> cfg_min_radius{"cfg_min_radius", 0.1, "min rxy for v0hadron"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_min_dcaxy{"cfg_min_dcaxy", 0.1, "min dca XY for v0 legs in cm"};

    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -2, "min n sigma pi in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +2, "max n sigma pi in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -2, "min n sigma pr in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +2, "max n sigma pr in TPC"};
    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -2, "min n sigma pi in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +2, "max n sigma pi in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -2, "min n sigma pr in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +2, "max n sigma pr in TOF"};
  } v0cut;

  struct : ConfigurableGroup {
    std::string prefix = "cascadecut";
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for lambda in cascade"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for lambda in cascade"};
    Configurable<float> cfg_min_mass_Xi{"cfg_min_mass_Xi", 1.316, "min mass for Xi"};
    Configurable<float> cfg_max_mass_Xi{"cfg_max_mass_Xi", 1.326, "max mass for Xi"};
    Configurable<float> cfg_min_mass_Xi_veto{"cfg_min_mass_Xi_veto", 1.31, "min mass for Xi veto"};
    Configurable<float> cfg_max_mass_Xi_veto{"cfg_max_mass_Xi_veto", 1.33, "max mass for Xi veto"};
    Configurable<float> cfg_min_mass_Omega{"cfg_min_mass_Omega", 1.669, "min mass for Omega"};
    Configurable<float> cfg_max_mass_Omega{"cfg_max_mass_Omega", 1.675, "max mass for Omega"};
    Configurable<float> cfg_min_mass_Omega_veto{"cfg_min_mass_Omega_veto", 1.66, "min mass for Omega veto"};
    Configurable<float> cfg_max_mass_Omega_veto{"cfg_max_mass_Omega_veto", 1.68, "max mass for Omega veto"};
    Configurable<float> cfg_min_cospa_v0{"cfg_min_cospa_v0", 0.995, "minimum V0 CosPA in cascade"};
    Configurable<float> cfg_max_dcadau_v0{"cfg_max_dcadau_v0", 0.1, "max distance between V0 Daughters in cascade"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.9998, "minimum cascade CosPA"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.1, "max distance between bachelor and V0"};
    Configurable<float> cfg_min_rxy_v0{"cfg_min_rxy_v0", 1.2, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_rxy{"cfg_min_rxy", 0.5, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};
    Configurable<float> cfg_min_dcaxy_bachelor{"cfg_min_dcaxy_bachelor", 0.05, "min dca XY for bachelor in cm"};
    Configurable<float> cfg_min_dcaxy_v0{"cfg_min_dcaxy_v0", 0.05, "min dca XY for V0 in cm"};
  } cascadecut;

  struct : ConfigurableGroup {
    std::string prefix = "eventcut";
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  } eventcut;

  Configurable<float> cfgMeeMaxPF{"cfgMeeMaxPF", 0.04, "max mee for prefilter to reject pi0->ee and gamma->ee in LMR"};

  HistogramRegistry fRegistry{"fRegistry"};
  static constexpr std::string_view hadron_names[6] = {"LF/", "Jpsi/", "D0/", "Dpm/", "Ds/", "Lc/"};
  static constexpr std::string_view pair_names[3] = {"e_Kpm/", "e_K0S/", "e_Lambda/"};
  static constexpr std::string_view hTypes[4] = {"findable/", "correct/", "fake/", "miss/"};
  static constexpr std::string_view promptTypes[2] = {"prompt/", "nonprompt/"};

  void init(o2::framework::InitContext&)
  {
    // if (doprocessSA && doprocessTTCA) {
    //   LOGF(fatal, "Cannot enable doprocessWithoutFTTCA and doprocessWithFTTCA at the same time. Please choose one.");
    // }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(10.f);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);
    fitter.setMatCorrType(matCorr);

    addHistograms();
  }

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;
  o2::dataformats::DCA mDcaInfoCov;
  o2::dataformats::VertexBase mVtx;

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
    fitter.setBz(d_bz);
  }

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("Event/hCollisionCounter", "collision counter", kTH1D, {{5, -0.5f, 4.5f}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    fRegistry.add("Event/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Event/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    fRegistry.add("Event/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{200, 0, 200000}, {60, 0, 60000}}, false);
    fRegistry.add("Event/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    fRegistry.add("Event/hCentFT0CvsMultNTracksPV", "hCentFT0CvsMultNTracksPV;centrality FT0C (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    fRegistry.add("Event/hMultFT0CvsMultNTracksPV", "hMultFT0CvsMultNTracksPV;mult. FT0C;N_{track} to PV", kTH2F, {{60, 0, 60000}, {600, 0, 6000}}, false);

    fRegistry.add("Prefilter/before/uls/hMee", "hMee;m_{ee} (GeV/c^{2});", kTH1F, {{500, 0, 5}}, false);
    fRegistry.addClone("Prefilter/before/uls/", "Prefilter/before/lspp/");
    fRegistry.addClone("Prefilter/before/uls/", "Prefilter/before/lsmm/");
    fRegistry.addClone("Prefilter/before/", "Prefilter/after/");

    // electron-related histograms
    fRegistry.add("Data/electron/hs", "hs;p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);DCA_{e}^{3D} (#sigma);", kTHnSparseF, {{100, 0, 10}, {20, -1, +1}, {90, 0, 2 * M_PI}, {100, 0, 10}}, false);
    fRegistry.addClone("Data/electron/", "MC/eFromPromptLF/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptLF/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptJpsi/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptJpsi/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptD0/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptDpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptDs/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptLcpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptXic0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptXicpm/"); // cannot be detected
    fRegistry.addClone("Data/electron/", "MC/eFromPromptOmegac0/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptD0/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptDpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptDs/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptLcpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptXic0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptXicpm/"); // cannot be detected
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptOmegac0/");
    fRegistry.addClone("Data/electron/", "MC/eFromB0/");
    fRegistry.addClone("Data/electron/", "MC/eFromBpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromBs/");
    fRegistry.addClone("Data/electron/", "MC/eFromBc/");
    fRegistry.addClone("Data/electron/", "MC/eFromLb0/");

    // for V0 (Lambda)
    fRegistry.add("Data/V0/hPt", "pT of V0;p_{T} (GeV/c)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/V0/hYPhi", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{#Lambda}", kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    fRegistry.add("Data/V0/hAP", "Ap plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1, 1}, {250, 0, 0.25}}, false);
    fRegistry.add("Data/V0/hLxy", "decay length from PV;L_{xy} (cm)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/V0/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Data/V0/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("Data/V0/hMassK0S", "K0S mass;m_{#pi#pi} (GeV/c^{2})", kTH1F, {{100, 0.45, 0.55}}, false);
    fRegistry.add("Data/V0/hMassLambda", "Lambda mass;m_{p#pi^{-}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);
    fRegistry.add("Data/V0/hMassAntiLambda", "Anti-Lambda mass;m_{#bar{p}#pi^{+}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);

    // for cascade
    fRegistry.add("Data/Cascade/hPt", "pT of V0;p_{T} (GeV/c)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/Cascade/hYPhi", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{#Lambda}", kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    fRegistry.add("Data/Cascade/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Data/Cascade/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("Data/Cascade/hV0CosPA", "cosPA of V0 in cascade;cosine of pointing angle", kTH1F, {{100, 0.99, 1}}, false);
    fRegistry.add("Data/Cascade/hV0DCA2Legs", "distance between 2 legs at PCA of V0 in cascade;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);

    fRegistry.add("Data/Cascade/hMassLambda", "Lambda mass;m_{p#pi^{-}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);
    fRegistry.add("Data/Cascade/hMassXi", "#Xi mass;m_{#Lambda#pi} (GeV/c^{2})", kTH1F, {{100, 1.27, 1.37}}, false);
    fRegistry.add("Data/Cascade/hMassOmega", "#Omega mass;m_{#LambdaK} (GeV/c^{2})", kTH1F, {{100, 1.62, 1.72}}, false);

    // for e-L pair
    fRegistry.add("Data/eL/RS/hs", "hs;m_{e#Lambda} (GeV/c^{2});p_{T,e} (GeV/c);DCA_{e}^{3D} (#sigma);L_{xy} (cm);", kTHnSparseF, {{20, 1.1, 3.1}, {100, 0, 10}, {100, 0, 10}, {100, 0, 1.0}}, false);
    fRegistry.add("Data/eL/RS/hCosPA", "cos PA;cosPA", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Data/eL/RS/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs at PCA (cm)", kTH1F, {{500, 0.0, 0.5}}, false);
    fRegistry.add("Data/eL/RS/hLxy", "distance between PV and SV in XY;L_{xy} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.add("Data/eL/RS/hLz", "distance between PV and SV in Z;L_{z} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.addClone("Data/eL/RS/", "Data/eL/WS/"); // right and wrong sign
    fRegistry.addClone("Data/eL/RS/", "MC/eLfromPromptLcpm/");
    fRegistry.addClone("Data/eL/RS/", "MC/eLfromNonPromptLcpm/");

    // for e-Xi pair
    fRegistry.add("Data/eXi/RS/hs", "hs;m_{e#Xi} (GeV/c^{2});p_{T,e} (GeV/c);DCA_{e}^{3D} (#sigma);L_{xy} (cm);", kTHnSparseF, {{20, 1.3, 3.3}, {100, 0, 10}, {100, 0, 10}, {100, 0, 1.0}}, false);
    fRegistry.add("Data/eXi/RS/hCosPA", "cos PA;cosPA", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Data/eXi/RS/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs at PCA (cm)", kTH1F, {{500, 0.0, 0.5}}, false);
    fRegistry.add("Data/eXi/RS/hLxy", "distance between PV and SV in XY;L_{xy} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.add("Data/eXi/RS/hLz", "distance between PV and SV in Z;L_{z} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.addClone("Data/eXi/RS/", "Data/eXi/WS/"); // right and wrong sign
    fRegistry.addClone("Data/eXi/RS/", "MC/eXifromPromptXic0/");
    fRegistry.addClone("Data/eXi/RS/", "MC/eXifromNonPromptXic0/");

    // for e-Omega pair
    fRegistry.add("Data/eOmega/RS/hs", "hs;m_{e#Omega} (GeV/c^{2});p_{T,e} (GeV/c);DCA_{e}^{3D} (#sigma);L_{xy} (cm);", kTHnSparseF, {{20, 1.6, 3.6}, {100, 0, 10}, {100, 0, 10}, {100, 0, 1.0}}, false);
    fRegistry.add("Data/eOmega/RS/hCosPA", "cos PA;cosPA", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Data/eOmega/RS/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs at PCA (cm)", kTH1F, {{500, 0.0, 0.5}}, false);
    fRegistry.add("Data/eOmega/RS/hLxy", "distance between PV and SV in XY;L_{xy} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.add("Data/eOmega/RS/hLz", "distance between PV and SV in Z;L_{z} (cm)", kTH1F, {{100, 0.0, 1}}, false);
    fRegistry.addClone("Data/eOmega/RS/", "Data/eOmega/WS/"); // right and wrong sign
    fRegistry.addClone("Data/eOmega/RS/", "MC/eOmegafromPromptOmegac0/");
    fRegistry.addClone("Data/eOmega/RS/", "MC/eOmegafromNonPromptOmegac0/");
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    // TOFif
    bool is_ka_included_TPC = kaoncut.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < kaoncut.cfg_max_TPCNsigmaKa;
    bool is_ka_included_TOF = track.hasTOF() ? (kaoncut.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < kaoncut.cfg_max_TOFNsigmaKa) : true;
    return is_ka_included_TPC && is_ka_included_TOF;
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    // TOFif
    bool is_pi_included_TPC = v0cut.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0cut.cfg_max_TPCNsigmaPi;
    bool is_pi_included_TOF = track.hasTOF() ? (v0cut.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < v0cut.cfg_max_TOFNsigmaPi) : true;
    return is_pi_included_TPC && is_pi_included_TOF;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    // TOFif
    bool is_pr_included_TPC = v0cut.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0cut.cfg_max_TPCNsigmaPr;
    bool is_pr_included_TOF = track.hasTOF() ? (v0cut.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < v0cut.cfg_max_TOFNsigmaPr) : true;
    return is_pr_included_TPC && is_pr_included_TOF;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    // TOFif
    bool is_el_included_TPC = electroncut.cfg_min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < electroncut.cfg_max_TPCNsigmaEl;
    bool is_el_included_TOF = track.hasTOF() ? (electroncut.cfg_min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < electroncut.cfg_max_TOFNsigmaEl) : true;
    return is_el_included_TPC && is_el_included_TOF;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedElectron(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < electroncut.cfg_min_pt_track || electroncut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < electroncut.cfg_min_eta_track || electroncut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > electroncut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > electroncut.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || electroncut.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < electroncut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < electroncut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || electroncut.cfg_max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < electroncut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electroncut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electroncut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electroncut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    if (!isElectron(track)) {
      return false;
    }

    return true;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedElectronLoose(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS()) {
      return false;
    }

    if (trackParCov.getPt() < loose_electroncut.cfg_min_pt_track || loose_electroncut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < loose_electroncut.cfg_min_eta_track || loose_electroncut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > loose_electroncut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > loose_electroncut.cfg_max_dcaz) {
      return false;
    }

    if (loose_electroncut.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < loose_electroncut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < loose_electroncut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.hasTPC()) {
      if (loose_electroncut.cfg_max_chi2tpc < track.tpcChi2NCl()) {
        return false;
      }

      if (track.tpcNClsFound() < loose_electroncut.cfg_min_ncluster_tpc) {
        return false;
      }

      if (track.tpcNClsCrossedRows() < loose_electroncut.cfg_min_ncrossedrows_tpc) {
        return false;
      }

      if (track.tpcCrossedRowsOverFindableCls() < loose_electroncut.cfg_min_cr2findable_ratio_tpc) {
        return false;
      }

      if (track.tpcFractionSharedCls() > loose_electroncut.cfg_max_frac_shared_clusters_tpc) {
        return false;
      }

      if (!isElectron(track)) {
        return false;
      }
    }

    return true;
  }

  template <typename TV0>
  bool isLambda(TV0 const& v0)
  {
    return (v0cut.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < v0cut.cfg_max_mass_lambda) && (v0.mK0Short() < v0cut.cfg_min_mass_k0s_veto || v0cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TV0>
  bool isAntiLambda(TV0 const& v0)
  {
    return (v0cut.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < v0cut.cfg_max_mass_lambda) && (v0.mK0Short() < v0cut.cfg_min_mass_k0s_veto || v0cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TCascade>
  bool isXi(TCascade const& cascade)
  {
    return (cascadecut.cfg_min_mass_Xi < cascade.mXi() && cascade.mXi() < cascadecut.cfg_max_mass_Xi) && (cascade.mOmega() < cascadecut.cfg_min_mass_Omega_veto || cascadecut.cfg_max_mass_Omega_veto < cascade.mOmega());
  }

  template <typename TCascade>
  bool isOmega(TCascade const& cascade)
  {
    return (cascadecut.cfg_min_mass_Omega < cascade.mOmega() && cascade.mOmega() < cascadecut.cfg_max_mass_Omega) && (cascade.mXi() < cascadecut.cfg_min_mass_Xi || cascadecut.cfg_max_mass_Xi < cascade.mXi());
  }

  template <bool isMC, typename TTrack>
  bool isSelectedV0Leg(TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsChi2NCl() > v0cut.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < v0cut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < v0cut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > v0cut.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < v0cut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < v0cut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < v0cut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > v0cut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  void fillEventHistograms(TCollision const& collision)
  {
    fRegistry.fill(HIST("Event/hZvtx"), collision.posZ());
    fRegistry.fill(HIST("Event/hMultNTracksPV"), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    fRegistry.fill(HIST("Event/hMultFT0"), collision.multFT0A(), collision.multFT0C());
    fRegistry.fill(HIST("Event/hCentFT0A"), collision.centFT0A());
    fRegistry.fill(HIST("Event/hCentFT0C"), collision.centFT0C());
    fRegistry.fill(HIST("Event/hCentFT0M"), collision.centFT0M());
    fRegistry.fill(HIST("Event/hCentFT0CvsMultNTracksPV"), collision.centFT0C(), collision.multNTracksPV());
    fRegistry.fill(HIST("Event/hMultFT0CvsMultNTracksPV"), collision.multFT0C(), collision.multNTracksPV());
  }

  template <typename TV0>
  void fillV0Histograms(TV0 const& v0)
  {
    fRegistry.fill(HIST("Data/V0/hPt"), v0.pt());
    fRegistry.fill(HIST("Data/V0/hYPhi"), v0.phi(), v0.yLambda());
    fRegistry.fill(HIST("Data/V0/hAP"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("Data/V0/hCosPA"), v0.v0cosPA());
    fRegistry.fill(HIST("Data/V0/hLxy"), v0.v0radius());
    fRegistry.fill(HIST("Data/V0/hDCA2Legs"), v0.dcaV0daughters());
    fRegistry.fill(HIST("Data/V0/hMassK0S"), v0.mK0Short());
    fRegistry.fill(HIST("Data/V0/hMassLambda"), v0.mLambda());
    fRegistry.fill(HIST("Data/V0/hMassAntiLambda"), v0.mAntiLambda());
  }

  template <typename TCollision, typename TCascade>
  void fillCascadeHistograms(TCollision const& collision, TCascade const& cascade)
  {
    fRegistry.fill(HIST("Data/Cascade/hPt"), cascade.pt());
    fRegistry.fill(HIST("Data/Cascade/hMassLambda"), cascade.mLambda());
    fRegistry.fill(HIST("Data/Cascade/hCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
    fRegistry.fill(HIST("Data/Cascade/hDCA2Legs"), cascade.dcacascdaughters());
    fRegistry.fill(HIST("Data/Cascade/hV0CosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
    fRegistry.fill(HIST("Data/Cascade/hV0DCA2Legs"), cascade.dcaV0daughters());
    fRegistry.fill(HIST("Data/Cascade/hMassXi"), cascade.mXi());
    fRegistry.fill(HIST("Data/Cascade/hMassOmega"), cascade.mOmega());
  }

  template <bool isMC, typename TTrack, typename TMCParticles>
  void fillElectronHistograms(TTrack const& track, TMCParticles const& mcParticles)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();
    float dca3DinSigma = dca3DinSigmaOTF(dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());

    if (!isSelectedElectron(track, trackParCov, dcaXY, dcaZ)) {
      return;
    }
    fRegistry.fill(HIST("Data/electron/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);

    if constexpr (isMC) {
      const auto& mctrack = track.template mcParticle_as<aod::McParticles>();
      if (std::abs(mctrack.pdgCode()) != 11) {
        return;
      }
      if (!(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
        return;
      }
      const auto& mcmother = mctrack.template mothers_first_as<aod::McParticles>(); // mother particle of electron
      int pdg_mother = std::abs(mcmother.pdgCode());

      if (pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 331 || pdg_mother == 113 || pdg_mother == 223 || pdg_mother == 333) { // LF
        if (IsFromCharm(mcmother, mcParticles) < 0 && IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptLF/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptLF/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 443) { // Jpsi
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptJpsi/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptJpsi/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 411) { // D+/-
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptDpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptDpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 421) { // D0
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptD0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptD0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 431) { // Ds+/-
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptDs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptDs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 4122) { // Lc+/-
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptLcpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptLcpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 4132) { // Xic0
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptXic0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptXic0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 4332) { // Omegac0
        if (IsFromBeauty(mcmother, mcParticles) < 0) {
          fRegistry.fill(HIST("MC/eFromPromptOmegac0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        } else {
          fRegistry.fill(HIST("MC/eFromNonPromptOmegac0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        }
      } else if (pdg_mother == 511) { // B0
        fRegistry.fill(HIST("MC/eFromB0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      } else if (pdg_mother == 521) { // B+/-
        fRegistry.fill(HIST("MC/eFromBpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      } else if (pdg_mother == 531) { // Bs0
        fRegistry.fill(HIST("MC/eFromBs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      } else if (pdg_mother == 541) { // Bc+/-
        fRegistry.fill(HIST("MC/eFromBc/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      } else if (pdg_mother == 5122) { // Lb0
        fRegistry.fill(HIST("MC/eFromLb0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      }
    }
  }

  float dca3DinSigmaOTF(const float dcaXY, const float dcaZ, const float cYY, const float cZZ, const float cZY)
  {
    float det = cYY * cZZ - cZY * cZY; // determinant
    if (det < 0) {
      return 999.f;
    } else {
      return std::sqrt(std::fabs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
    }
  }

  template <int pdgLepton, int pdgNeutrino, typename TMCParticle, typename TMCParticles>
  bool isSemiLeptonic(TMCParticle const& mcParticle, TMCParticles const& mcParticles)
  {
    if (!mcParticle.has_daughters()) {
      return false;
    }
    bool is_lepton_involved = false;
    bool is_neutrino_involved = false;
    for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
      if (d < mcParticles.size()) { // protect against bad daughter indices
        const auto& daughter = mcParticles.rawIteratorAt(d);
        if (daughter.pdgCode() == pdgLepton) {
          is_lepton_involved = true;
        } else if (daughter.pdgCode() == pdgNeutrino) {
          is_neutrino_involved = true;
        }
      } else {
        std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")" << std::endl;
        std::cout << " Check the MC generator" << std::endl;
        return false;
      }
    }

    if (is_lepton_involved && is_neutrino_involved) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TTrack, typename TV0>
  EBPair makeELPair(TCollision const& collision, TTrack const& track, TV0 const& v0)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();
    float dca3DinSigma = dca3DinSigmaOTF(dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());

    EBPair eLpair;
    eLpair.ptepv = trackParCov.getPt();
    eLpair.dca3dinsigma = dca3DinSigma;

    const std::array<float, 3> vertex = {collision.posX(), collision.posY(), collision.posZ()};
    const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
    const std::array<float, 3> momV0 = {v0.px(), v0.py(), v0.pz()};
    std::array<float, 21> covV0 = {0.f};

    constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      covV0[MomInd[i]] = v0.momentumCovMat()[i];
      covV0[i] = v0.positionCovMat()[i];
    }

    auto v0ParCov = o2::track::TrackParCov(vertexV0, momV0, covV0, 0, true);
    v0ParCov.setAbsCharge(0);
    v0ParCov.setPID(o2::track::PID::Lambda);

    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    int nCand = 0;
    try {
      nCand = fitter.process(trackParCov, v0ParCov);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return eLpair;
    }
    if (nCand == 0) {
      return eLpair;
    }

    fitter.propagateTracksToVertex(); // propagate e and K to D vertex
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
    fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
    std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

    float cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
    float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
    float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
    float lz = std::fabs(svpos[2] - collision.posZ());
    ROOT::Math::PxPyPzMVector v1(pvec0[0], pvec0[1], pvec0[2], o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector v2(pvec1[0], pvec1[1], pvec1[2], o2::constants::physics::MassLambda);
    ROOT::Math::PxPyPzMVector v12 = v1 + v2;

    eLpair.mass = v12.M();
    eLpair.dca2legs = dca2legs;
    eLpair.cospa = cospa;
    eLpair.lxy = lxy;
    eLpair.lz = lz;
    return eLpair;
  }

  template <uint8_t cascType, typename TCollision, typename TTrack, typename TCascade>
  EBPair makeECascadePair(TCollision const& collision, TTrack const& track, TCascade const& cascade)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();
    float dca3DinSigma = dca3DinSigmaOTF(dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());

    EBPair eCascPair;
    eCascPair.ptepv = trackParCov.getPt();
    eCascPair.dca3dinsigma = dca3DinSigma;

    const std::array<float, 3> vertex = {collision.posX(), collision.posY(), collision.posZ()};
    const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
    const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};

    std::array<float, 21> covCasc = {0.};
    constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
      covCasc[i] = cascade.positionCovMat()[i];
    }

    auto cascParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
    cascParCov.setAbsCharge(1);
    if constexpr (cascType == 0) {
      cascParCov.setPID(o2::track::PID::XiMinus);
    } else if constexpr (cascType == 1) {
      cascParCov.setPID(o2::track::PID::OmegaMinus);
    }
    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    int nCand = 0;
    try {
      nCand = fitter.process(trackParCov, cascParCov);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return eCascPair;
    }
    if (nCand == 0) {
      return eCascPair;
    }

    fitter.propagateTracksToVertex(); // propagate e and Xi/Omega to decay vertex of charm baryon
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
    fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
    std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

    float cospa = RecoDecay::cpa(vertex, svpos, pvecSum);
    float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
    float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
    float lz = std::fabs(svpos[2] - collision.posZ());
    ROOT::Math::PxPyPzMVector v1(pvec0[0], pvec0[1], pvec0[2], o2::constants::physics::MassElectron);
    ROOT::Math::PxPyPzMVector v2(pvec1[0], pvec1[1], pvec1[2], o2::constants::physics::MassXiMinus);
    if constexpr (cascType == 0) {
      v2.SetM(o2::constants::physics::MassXiMinus);
    } else if constexpr (cascType == 1) {
      v2.SetM(o2::constants::physics::MassOmegaMinus);
    }
    ROOT::Math::PxPyPzMVector v12 = v1 + v2;

    eCascPair.mass = v12.M();
    eCascPair.dca2legs = dca2legs;
    eCascPair.cospa = cospa;
    eCascPair.lxy = lxy;
    eCascPair.lz = lz;
    return eCascPair;
  }

  template <bool isMC, typename TBCs, typename TCollisions, typename TTracks, typename TTrackAssoc, typename TV0s, typename TCascades, typename TMCCollisions, typename TMCParticles>
  void runPairing(TBCs const&, TCollisions const& collisions, TTracks const& tracks, TTrackAssoc const& trackIndices, TV0s const& v0s, TCascades const& cascades, TMCCollisions const&, TMCParticles const& mcParticles)
  {
    used_electronIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      const auto& bc = collision.template foundBC_as<TBCs>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[eventcut.cfgCentEstimator] < eventcut.cfgCentMin || eventcut.cfgCentMax < centralities[eventcut.cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      const auto& mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (eventcut.cfgEventGeneratorType < 0 || mcCollision.getSubGeneratorId() == eventcut.cfgEventGeneratorType) {
        fillEventHistograms(collision);
      }
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      electronIds.reserve(trackIdsThisCollision.size());
      positronIds.reserve(trackIdsThisCollision.size());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = trackId.template track_as<TTracks>();
        if (!track.hasITS() || !track.hasTPC()) {
          continue;
        }

        if constexpr (isMC) {
          if (!track.has_mcParticle()) {
            continue;
          }
          const auto& mctrack = track.template mcParticle_as<aod::McParticles>();
          const auto& mcCollision = mctrack.template mcCollision_as<aod::McCollisions>();
          if (eventcut.cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventcut.cfgEventGeneratorType) {
            continue;
          }
          if (!mctrack.has_mothers() || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
            continue;
          }

          fillElectronHistograms<isMC>(track, mcParticles);
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto trackParCov = getTrackParCov(track);
        trackParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();

        if (isSelectedElectron(track, trackParCov, dcaXY, dcaZ)) {
          if (track.sign() > 0) { // positron
            positronIds.emplace_back(trackId.trackId());
          } else { // electron
            electronIds.emplace_back(trackId.trackId());
          }
        }

        if (isSelectedElectronLoose(track, trackParCov, dcaXY, dcaZ)) {
          if (track.sign() > 0) { // positron
            positronIdsLoose.emplace_back(trackId.trackId());
          } else { // electron
            electronIdsLoose.emplace_back(trackId.trackId());
          }
        }
      } // end of track loop for electron selection

      // First, apply pi0 prefilter to imporove S/B
      std::vector<int> vec_eFromPi0;
      vec_eFromPi0.reserve(electronIds.size() + positronIds.size());

      for (const auto& positronId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(positronId);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto posParCov = getTrackParCov(pos);
        posParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, posParCov, 2.f, matCorr, &mDcaInfoCov);
        ROOT::Math::PtEtaPhiMVector v1(posParCov.getPt(), posParCov.getEta(), posParCov.getPhi(), o2::constants::physics::MassElectron);

        for (const auto& electronId : electronIdsLoose) {
          const auto& ele = tracks.rawIteratorAt(electronId);
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto eleParCov = getTrackParCov(ele);
          eleParCov.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, eleParCov, 2.f, matCorr, &mDcaInfoCov);
          ROOT::Math::PtEtaPhiMVector v2(eleParCov.getPt(), eleParCov.getEta(), eleParCov.getPhi(), o2::constants::physics::MassElectron);

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          fRegistry.fill(HIST("Prefilter/before/uls/hMee"), mee);
          if (mee < cfgMeeMaxPF) {
            vec_eFromPi0.emplace_back(positronId);
          }
        } // end of loose electron sample
      } // end of main positron sample

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto eleParCov = getTrackParCov(ele);
        eleParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, eleParCov, 2.f, matCorr, &mDcaInfoCov);
        ROOT::Math::PtEtaPhiMVector v1(eleParCov.getPt(), eleParCov.getEta(), eleParCov.getPhi(), o2::constants::physics::MassElectron);

        for (const auto& positronId : positronIdsLoose) {
          const auto& pos = tracks.rawIteratorAt(positronId);
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto posParCov = getTrackParCov(pos);
          posParCov.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, posParCov, 2.f, matCorr, &mDcaInfoCov);
          ROOT::Math::PtEtaPhiMVector v2(posParCov.getPt(), posParCov.getEta(), posParCov.getPhi(), o2::constants::physics::MassElectron);

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          fRegistry.fill(HIST("Prefilter/before/uls/hMee"), mee);
          if (mee < cfgMeeMaxPF) {
            vec_eFromPi0.emplace_back(electronId);
          }
        } // end of loose positron sample
      } // end of main electron sample

      for (const auto& positronId1 : positronIds) {
        const auto& pos1 = tracks.rawIteratorAt(positronId1);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto pos1ParCov = getTrackParCov(pos1);
        pos1ParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, pos1ParCov, 2.f, matCorr, &mDcaInfoCov);
        ROOT::Math::PtEtaPhiMVector v1(pos1ParCov.getPt(), pos1ParCov.getEta(), pos1ParCov.getPhi(), o2::constants::physics::MassElectron);

        for (const auto& positronId2 : positronIdsLoose) {
          const auto& pos2 = tracks.rawIteratorAt(positronId2);
          if (positronId1 == positronId2) {
            continue;
          }
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto pos2ParCov = getTrackParCov(pos2);
          pos2ParCov.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, pos2ParCov, 2.f, matCorr, &mDcaInfoCov);
          ROOT::Math::PtEtaPhiMVector v2(pos2ParCov.getPt(), pos2ParCov.getEta(), pos2ParCov.getPhi(), o2::constants::physics::MassElectron);

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          fRegistry.fill(HIST("Prefilter/before/lspp/hMee"), mee);
        } // end of loose positron sample
      } // end of main positron sample

      for (const auto& electronId1 : electronIds) {
        const auto& ele1 = tracks.rawIteratorAt(electronId1);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto ele1ParCov = getTrackParCov(ele1);
        ele1ParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, ele1ParCov, 2.f, matCorr, &mDcaInfoCov);
        ROOT::Math::PtEtaPhiMVector v1(ele1ParCov.getPt(), ele1ParCov.getEta(), ele1ParCov.getPhi(), o2::constants::physics::MassElectron);

        for (const auto& electronId2 : electronIdsLoose) {
          const auto& ele2 = tracks.rawIteratorAt(electronId2);
          if (electronId1 == electronId2) {
            continue;
          }
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto ele2ParCov = getTrackParCov(ele2);
          ele2ParCov.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, ele2ParCov, 2.f, matCorr, &mDcaInfoCov);
          ROOT::Math::PtEtaPhiMVector v2(ele2ParCov.getPt(), ele2ParCov.getEta(), ele2ParCov.getPhi(), o2::constants::physics::MassElectron);

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          fRegistry.fill(HIST("Prefilter/before/lsmm/hMee"), mee);
        } // end of loose electron sample
      } // end of main electron sample

      std::vector<int> vec_diff_pos;
      std::set_difference(positronIds.begin(), positronIds.end(), vec_eFromPi0.begin(), vec_eFromPi0.end(), std::back_inserter(vec_diff_pos));
      positronIds = vec_diff_pos;

      std::vector<int> vec_diff_ele;
      std::set_difference(electronIds.begin(), electronIds.end(), vec_eFromPi0.begin(), vec_eFromPi0.end(), std::back_inserter(vec_diff_ele));
      electronIds = vec_diff_ele;

      vec_eFromPi0.clear();
      vec_eFromPi0.shrink_to_fit();

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto eleParCov = getTrackParCov(ele);
        eleParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, eleParCov, 2.f, matCorr, &mDcaInfoCov);
        ROOT::Math::PtEtaPhiMVector v1(eleParCov.getPt(), eleParCov.getEta(), eleParCov.getPhi(), o2::constants::physics::MassElectron);

        for (const auto& positronId : positronIds) {
          const auto& pos = tracks.rawIteratorAt(positronId);
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto posParCov = getTrackParCov(pos);
          posParCov.setPID(o2::track::PID::Electron);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, posParCov, 2.f, matCorr, &mDcaInfoCov);
          ROOT::Math::PtEtaPhiMVector v2(posParCov.getPt(), posParCov.getEta(), posParCov.getPhi(), o2::constants::physics::MassElectron);

          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float mee = v12.M();
          fRegistry.fill(HIST("Prefilter/after/uls/hMee"), mee);
        } // end of main positron sample
      } // end of main electron sample

      const auto& v0s_per_coll = v0s.sliceBy(perCol_v0, collision.globalIndex());
      lambdaIds.reserve(v0s_per_coll.size());
      lambdaIds.reserve(v0s_per_coll.size());
      antilambdaIds.reserve(v0s_per_coll.size());
      antilambdaIds.reserve(v0s_per_coll.size());
      for (const auto& v0 : v0s_per_coll) {
        auto pos = v0.template posTrack_as<TTracks>();
        auto neg = v0.template negTrack_as<TTracks>();
        if (!isSelectedV0Leg<isMC>(pos) || !isSelectedV0Leg<isMC>(neg)) {
          continue;
        }
        fillV0Histograms(v0);
        if (isLambda(v0) && isProton(pos) && isPion(neg)) {
          lambdaIds.emplace_back(v0.globalIndex());
        } else if (isAntiLambda(v0) && isProton(neg) && isPion(pos)) {
          antilambdaIds.emplace_back(v0.globalIndex());
        }
      } // end of V0 loop

      const auto& cascades_per_coll = cascades.sliceBy(perCol_casc, collision.globalIndex());
      xiPlusIds.reserve(cascades_per_coll.size());
      xiPlusIds.reserve(cascades_per_coll.size());
      xiMinusIds.reserve(cascades_per_coll.size());
      xiMinusIds.reserve(cascades_per_coll.size());
      omegaPlusIds.reserve(cascades_per_coll.size());
      omegaPlusIds.reserve(cascades_per_coll.size());
      omegaMinusIds.reserve(cascades_per_coll.size());
      omegaMinusIds.reserve(cascades_per_coll.size());
      for (const auto& cascade : cascades_per_coll) {
        auto pos = cascade.template posTrack_as<TTracks>();
        auto neg = cascade.template negTrack_as<TTracks>();
        auto bachelor = cascade.template bachelor_as<TTracks>();
        if (pos.sign() * neg.sign() > 0) {
          continue;
        }
        if (cascade.mLambda() < cascadecut.cfg_min_mass_lambda || cascadecut.cfg_max_mass_lambda < cascade.mLambda()) {
          continue;
        }

        if (!isSelectedV0Leg<isMC>(pos) || !isSelectedV0Leg<isMC>(neg) || !isSelectedV0Leg<isMC>(bachelor)) {
          continue;
        }

        if (cascade.sign() < 0) { // L-> p pi-
          if (!isProton(pos) || !isPion(neg)) {
            continue;
          }
        } else { // Lbar-> pbar pi+
          if (!isProton(neg) || !isPion(pos)) {
            continue;
          }
        }

        if (cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecut.cfg_min_cospa) {
          continue;
        }
        if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecut.cfg_min_cospa_v0) {
          continue;
        }

        fillCascadeHistograms(collision, cascade);
        if (cascade.sign() < 0) { // Xi- or Omega-
          if (isXi(cascade) && isPion(bachelor)) {
            xiMinusIds.emplace_back(cascade.globalIndex());
          } else if (isOmega(cascade) && isKaon(bachelor)) {
            omegaMinusIds.emplace_back(cascade.globalIndex());
          }
        } else { // Xi+ or Omega+
          if (isXi(cascade) && isPion(bachelor)) {
            xiPlusIds.emplace_back(cascade.globalIndex());
          } else if (isOmega(cascade) && isKaon(bachelor)) {
            omegaPlusIds.emplace_back(cascade.globalIndex());
          }
        }
      } // end of cascade loop

      // Lc+ -> e+ Lambda nu_e, br = 0.0356, ctau = 60.75 um, m = 2286 MeV/c2
      for (const auto& positronId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(positronId);

        for (const auto& lambdaId : lambdaIds) {
          const auto& lambda = v0s.rawIteratorAt(lambdaId);
          const auto& eLpair = makeELPair(collision, pos, lambda); // RS
          fRegistry.fill(HIST("Data/eL/RS/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/RS/hCosPA"), eLpair.cospa);
          fRegistry.fill(HIST("Data/eL/RS/hDCA2Legs"), eLpair.dca2legs);
          fRegistry.fill(HIST("Data/eL/RS/hLxy"), eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/RS/hLz"), eLpair.lz);

          if constexpr (isMC) {
            const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
            auto posLeg = lambda.template posTrack_as<TTracks>();
            auto negLeg = lambda.template negTrack_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 2212, -211, 3122, mcParticles);
            if (mcLambdaId > 0) { // true lambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcLambdacId = FindCommonMotherFrom2Prongs(mcpos, mcLambda, -11, 3122, 4122, mcParticles);
              if (mcLambdacId > 0) { // true Lc0
                const auto& mcLambdac0 = mcParticles.rawIteratorAt(mcLambdacId);
                if (IsFromBeauty(mcLambdac0, mcParticles) < 0) {
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hCosPA"), eLpair.cospa);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hDCA2Legs"), eLpair.dca2legs);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hLxy"), eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hLz"), eLpair.lz);
                } else {
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hCosPA"), eLpair.cospa);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hDCA2Legs"), eLpair.dca2legs);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hLxy"), eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hLz"), eLpair.lz);
                }
              }
            }
          } // end of MC truth
        } // end of Lambda loop

        for (const auto& antilambdaId : antilambdaIds) {
          const auto& antilambda = v0s.rawIteratorAt(antilambdaId);
          const auto& eLpair = makeELPair(collision, pos, antilambda); // WS
          fRegistry.fill(HIST("Data/eL/WS/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/WS/hCosPA"), eLpair.cospa);
          fRegistry.fill(HIST("Data/eL/WS/hDCA2Legs"), eLpair.dca2legs);
          fRegistry.fill(HIST("Data/eL/WS/hLxy"), eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/WS/hLz"), eLpair.lz);
        } // end of AntiLambda loop

      } // end of main positron sample

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);

        for (const auto& lambdaId : lambdaIds) {
          const auto& lambda = v0s.rawIteratorAt(lambdaId);
          const auto& eLpair = makeELPair(collision, ele, lambda); // WS
          fRegistry.fill(HIST("Data/eL/WS/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/WS/hCosPA"), eLpair.cospa);
          fRegistry.fill(HIST("Data/eL/WS/hDCA2Legs"), eLpair.dca2legs);
          fRegistry.fill(HIST("Data/eL/WS/hLxy"), eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/WS/hLz"), eLpair.lz);
        } // end of Lambda loop

        for (const auto& antilambdaId : antilambdaIds) {
          const auto& antilambda = v0s.rawIteratorAt(antilambdaId);
          const auto& eLpair = makeELPair(collision, ele, antilambda); // RS
          fRegistry.fill(HIST("Data/eL/RS/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/RS/hCosPA"), eLpair.cospa);
          fRegistry.fill(HIST("Data/eL/RS/hDCA2Legs"), eLpair.dca2legs);
          fRegistry.fill(HIST("Data/eL/RS/hLxy"), eLpair.lxy);
          fRegistry.fill(HIST("Data/eL/RS/hLz"), eLpair.lz);

          if constexpr (isMC) {
            const auto& mcele = ele.template mcParticle_as<aod::McParticles>();
            auto posLeg = antilambda.template posTrack_as<TTracks>();
            auto negLeg = antilambda.template negTrack_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 221, -2212, -3122, mcParticles);
            if (mcLambdaId > 0) { // true lambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcLambdacId = FindCommonMotherFrom2Prongs(mcele, mcLambda, 11, -3122, -4122, mcParticles);
              if (mcLambdacId > 0) { // true Lc0
                const auto& mcLambdac0 = mcParticles.rawIteratorAt(mcLambdacId);
                if (IsFromBeauty(mcLambdac0, mcParticles) < 0) {
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hCosPA"), eLpair.cospa);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hDCA2Legs"), eLpair.dca2legs);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hLxy"), eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromPromptLcpm/hLz"), eLpair.lz);
                } else {
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hs"), eLpair.mass, eLpair.ptepv, eLpair.dca3dinsigma, eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hCosPA"), eLpair.cospa);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hDCA2Legs"), eLpair.dca2legs);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hLxy"), eLpair.lxy);
                  fRegistry.fill(HIST("MC/eLfromNonPromptLcpm/hLz"), eLpair.lz);
                }
              }
            }
          } // end of MC truth

        } // end of AntiLambda loop

      } // end of main electron sample

      // Xic0 -> e+ Xi- nu_e, br = 0.0105, ctau = 45.1 um, m = 2470 MeV/c2
      for (const auto& positronId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(positronId);

        for (const auto& xiMinusId : xiMinusIds) {
          const auto& xiMinus = cascades.rawIteratorAt(xiMinusId);
          const auto& eXipair = makeECascadePair<0>(collision, pos, xiMinus); // RS
          fRegistry.fill(HIST("Data/eXi/RS/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/RS/hCosPA"), eXipair.cospa);
          fRegistry.fill(HIST("Data/eXi/RS/hDCA2Legs"), eXipair.dca2legs);
          fRegistry.fill(HIST("Data/eXi/RS/hLxy"), eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/RS/hLz"), eXipair.lz);

          if constexpr (isMC) {
            const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
            auto posLeg = xiMinus.template posTrack_as<TTracks>();
            auto negLeg = xiMinus.template negTrack_as<TTracks>();
            auto bachelor = xiMinus.template bachelor_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 2212, -211, 3122, mcParticles);
            if (mcLambdaId > 0) { // true Lambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcXiId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, 3122, -211, 3312, mcParticles);
              if (mcXiId > 0) { // true xiMinus
                const auto& mcXi = mcParticles.rawIteratorAt(mcXiId);
                int mcXic0Id = FindCommonMotherFrom2Prongs(mcpos, mcXi, -11, 3312, 4132, mcParticles);
                if (mcXic0Id > 0) { // true Xic0
                  const auto& mcXic0 = mcParticles.rawIteratorAt(mcXic0Id);
                  if (IsFromBeauty(mcXic0, mcParticles) < 0) {
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hCosPA"), eXipair.cospa);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hDCA2Legs"), eXipair.dca2legs);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hLxy"), eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hLz"), eXipair.lz);
                  } else {
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hCosPA"), eXipair.cospa);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hDCA2Legs"), eXipair.dca2legs);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hLxy"), eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hLz"), eXipair.lz);
                  }
                }
              }
            }
          } // end of MC truth
        } // end of Lambda loop

        for (const auto& xiPlusId : xiPlusIds) {
          const auto& xiPlus = cascades.rawIteratorAt(xiPlusId);
          const auto& eXipair = makeECascadePair<0>(collision, pos, xiPlus); // WS
          fRegistry.fill(HIST("Data/eXi/WS/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/WS/hCosPA"), eXipair.cospa);
          fRegistry.fill(HIST("Data/eXi/WS/hDCA2Legs"), eXipair.dca2legs);
          fRegistry.fill(HIST("Data/eXi/WS/hLxy"), eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/WS/hLz"), eXipair.lz);
        } // end of AntiLambda loop
      } // end of main positron sample

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);

        for (const auto& xiMinusId : xiMinusIds) {
          const auto& xiMinus = cascades.rawIteratorAt(xiMinusId);
          const auto& eXipair = makeECascadePair<0>(collision, ele, xiMinus); // WS
          fRegistry.fill(HIST("Data/eXi/WS/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/WS/hCosPA"), eXipair.cospa);
          fRegistry.fill(HIST("Data/eXi/WS/hDCA2Legs"), eXipair.dca2legs);
          fRegistry.fill(HIST("Data/eXi/WS/hLxy"), eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/WS/hLz"), eXipair.lz);

        } // end of Xi- loop

        for (const auto& xiPlusId : xiPlusIds) {
          const auto& xiPlus = cascades.rawIteratorAt(xiPlusId);
          const auto& eXipair = makeECascadePair<0>(collision, ele, xiPlus); // RS
          fRegistry.fill(HIST("Data/eXi/RS/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/RS/hCosPA"), eXipair.cospa);
          fRegistry.fill(HIST("Data/eXi/RS/hDCA2Legs"), eXipair.dca2legs);
          fRegistry.fill(HIST("Data/eXi/RS/hLxy"), eXipair.lxy);
          fRegistry.fill(HIST("Data/eXi/RS/hLz"), eXipair.lz);

          if constexpr (isMC) {
            const auto& mcele = ele.template mcParticle_as<aod::McParticles>();
            auto posLeg = xiPlus.template posTrack_as<TTracks>();
            auto negLeg = xiPlus.template negTrack_as<TTracks>();
            auto bachelor = xiPlus.template bachelor_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -2212, -3122, mcParticles);
            if (mcLambdaId > 0) { // true AntiLambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcXiId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, -3122, 211, -3312, mcParticles);
              if (mcXiId > 0) { // true xiPlus
                const auto& mcXi = mcParticles.rawIteratorAt(mcXiId);
                int mcXic0Id = FindCommonMotherFrom2Prongs(mcele, mcXi, 11, -3312, 4132, mcParticles);
                if (mcXic0Id > 0) { // true Xic0
                  const auto& mcXic0 = mcParticles.rawIteratorAt(mcXic0Id);
                  if (IsFromBeauty(mcXic0, mcParticles) < 0) {
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hCosPA"), eXipair.cospa);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hDCA2Legs"), eXipair.dca2legs);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hLxy"), eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromPromptXic0/hLz"), eXipair.lz);
                  } else {
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hs"), eXipair.mass, eXipair.ptepv, eXipair.dca3dinsigma, eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hCosPA"), eXipair.cospa);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hDCA2Legs"), eXipair.dca2legs);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hLxy"), eXipair.lxy);
                    fRegistry.fill(HIST("MC/eXifromNonPromptXic0/hLz"), eXipair.lz);
                  }
                }
              }
            }
          } // end of MC truth

        } // end of Xi+ loop
      } // end of main electron sample

      // Omegac0 -> e+ Omega- nu_e, br(Omegac0 -> e+ Omega- nu_e) / br(Omegac0 -> Omega- pi+) = 1.98, ctau = 82 um, m = 2695 MeV/c2
      for (const auto& positronId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(positronId);

        for (const auto& omegaMinusId : omegaMinusIds) {
          const auto& omegaMinus = cascades.rawIteratorAt(omegaMinusId);
          const auto& eOmegapair = makeECascadePair<1>(collision, pos, omegaMinus); // RS
          fRegistry.fill(HIST("Data/eOmega/RS/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/RS/hCosPA"), eOmegapair.cospa);
          fRegistry.fill(HIST("Data/eOmega/RS/hDCA2Legs"), eOmegapair.dca2legs);
          fRegistry.fill(HIST("Data/eOmega/RS/hLxy"), eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/RS/hLz"), eOmegapair.lz);

          if constexpr (isMC) {
            const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
            auto posLeg = omegaMinus.template posTrack_as<TTracks>();
            auto negLeg = omegaMinus.template negTrack_as<TTracks>();
            auto bachelor = omegaMinus.template bachelor_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 2212, -211, 3122, mcParticles);
            if (mcLambdaId > 0) { // true Lambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcOmegaId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, 3122, -321, 3334, mcParticles);
              if (mcOmegaId > 0) { // true omegaMinus
                const auto& mcOmega = mcParticles.rawIteratorAt(mcOmegaId);
                int mcOmegac0Id = FindCommonMotherFrom2Prongs(mcpos, mcOmega, -11, 3334, 4332, mcParticles);
                if (mcOmegac0Id > 0) { // true Omegac0
                  const auto& mcOmegac0 = mcParticles.rawIteratorAt(mcOmegac0Id);
                  if (IsFromBeauty(mcOmegac0, mcParticles) < 0) {
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hCosPA"), eOmegapair.cospa);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hDCA2Legs"), eOmegapair.dca2legs);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hLxy"), eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hLz"), eOmegapair.lz);
                  } else {
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hCosPA"), eOmegapair.cospa);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hDCA2Legs"), eOmegapair.dca2legs);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hLxy"), eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hLz"), eOmegapair.lz);
                  }
                }
              }
            }
          } // end of MC truth
        } // end of Lambda loop

        for (const auto& omegaPlusId : omegaPlusIds) {
          const auto& omegaPlus = cascades.rawIteratorAt(omegaPlusId);
          const auto& eOmegapair = makeECascadePair<1>(collision, pos, omegaPlus); // WS
          fRegistry.fill(HIST("Data/eOmega/WS/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/WS/hCosPA"), eOmegapair.cospa);
          fRegistry.fill(HIST("Data/eOmega/WS/hDCA2Legs"), eOmegapair.dca2legs);
          fRegistry.fill(HIST("Data/eOmega/WS/hLxy"), eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/WS/hLz"), eOmegapair.lz);
        } // end of AntiLambda loop
      } // end of main positron sample

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);

        for (const auto& omegaMinusId : omegaMinusIds) {
          const auto& omegaMinus = cascades.rawIteratorAt(omegaMinusId);
          const auto& eOmegapair = makeECascadePair<1>(collision, ele, omegaMinus); // WS
          fRegistry.fill(HIST("Data/eOmega/WS/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/WS/hCosPA"), eOmegapair.cospa);
          fRegistry.fill(HIST("Data/eOmega/WS/hDCA2Legs"), eOmegapair.dca2legs);
          fRegistry.fill(HIST("Data/eOmega/WS/hLxy"), eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/WS/hLz"), eOmegapair.lz);

        } // end of Omega- loop

        for (const auto& omegaPlusId : omegaPlusIds) {
          const auto& omegaPlus = cascades.rawIteratorAt(omegaPlusId);
          const auto& eOmegapair = makeECascadePair<1>(collision, ele, omegaPlus); // RS
          fRegistry.fill(HIST("Data/eOmega/RS/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/RS/hCosPA"), eOmegapair.cospa);
          fRegistry.fill(HIST("Data/eOmega/RS/hDCA2Legs"), eOmegapair.dca2legs);
          fRegistry.fill(HIST("Data/eOmega/RS/hLxy"), eOmegapair.lxy);
          fRegistry.fill(HIST("Data/eOmega/RS/hLz"), eOmegapair.lz);

          if constexpr (isMC) {
            const auto& mcele = ele.template mcParticle_as<aod::McParticles>();
            auto posLeg = omegaPlus.template posTrack_as<TTracks>();
            auto negLeg = omegaPlus.template negTrack_as<TTracks>();
            auto bachelor = omegaPlus.template bachelor_as<TTracks>();
            const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
            const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
            int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -2212, -3122, mcParticles);
            if (mcLambdaId > 0) { // true AntiLambda
              const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
              int mcOmegaId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, -3122, 321, -3334, mcParticles);
              if (mcOmegaId > 0) { // true omegaPlus
                const auto& mcOmega = mcParticles.rawIteratorAt(mcOmegaId);
                int mcOmegac0Id = FindCommonMotherFrom2Prongs(mcele, mcOmega, 11, -3334, 4332, mcParticles);
                if (mcOmegac0Id > 0) { // true Omegac0
                  const auto& mcOmegac0 = mcParticles.rawIteratorAt(mcOmegac0Id);
                  if (IsFromBeauty(mcOmegac0, mcParticles) < 0) {
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hCosPA"), eOmegapair.cospa);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hDCA2Legs"), eOmegapair.dca2legs);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hLxy"), eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromPromptOmegac0/hLz"), eOmegapair.lz);
                  } else {
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hs"), eOmegapair.mass, eOmegapair.ptepv, eOmegapair.dca3dinsigma, eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hCosPA"), eOmegapair.cospa);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hDCA2Legs"), eOmegapair.dca2legs);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hLxy"), eOmegapair.lxy);
                    fRegistry.fill(HIST("MC/eOmegafromNonPromptOmegac0/hLz"), eOmegapair.lz);
                  }
                }
              }
            }
          } // end of MC truth

        } // end of Omega+ loop
      } // end of main electron sample

      electronIdsLoose.clear();
      electronIdsLoose.shrink_to_fit();
      positronIdsLoose.clear();
      positronIdsLoose.shrink_to_fit();
      electronIds.clear();
      electronIds.shrink_to_fit();
      positronIds.clear();
      positronIds.shrink_to_fit();

      lambdaIds.clear();
      lambdaIds.shrink_to_fit();
      antilambdaIds.clear();
      antilambdaIds.shrink_to_fit();

      xiPlusIds.clear();
      xiPlusIds.shrink_to_fit();
      xiMinusIds.clear();
      xiMinusIds.shrink_to_fit();

      omegaPlusIds.clear();
      omegaPlusIds.shrink_to_fit();
      omegaMinusIds.clear();
      omegaMinusIds.shrink_to_fit();
    } // end of collision loop

    used_electronIds.clear();
    used_electronIds.shrink_to_fit();
  }

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Preslice<aod::V0Datas> perCol_v0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> perCol_casc = o2::aod::cascdata::collisionId;

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && (eventcut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventcut.cfgZvtxMax);
  Filter collisionFilter_centrality = (eventcut.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventcut.cfgCentMax) || (eventcut.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventcut.cfgCentMax) || (eventcut.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventcut.cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollisionsWithMCLabel = soa::Filtered<MyCollisionsWithMCLabel>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  std::vector<std::pair<int, int>> stored_trackIds;

  // Filter trackFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  // using MyFilteredTracks = soa::Filtered<MyTracks>;
  // using MyFilteredTracksWithMCLabel = soa::Filtered<MyTracksWithMCLabel>;
  // Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  // Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cut), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > v0cut.cfg_min_cospa&& o2::aod::v0data::dcaV0daughters<v0cut.cfg_max_dca2legs && nabs(o2::aod::v0data::dcanegtopv)> v0cut.cfg_min_dcaxy&& nabs(o2::aod::v0data::dcanegtopv) > v0cut.cfg_min_dcaxy;
  using filteredV0s = soa::Filtered<MyV0s>;

  Filter cascadeFilter = nabs(o2::aod::cascdata::dcanegtopv) > cascadecut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcanegtopv) > cascadecut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcabachtopv) > cascadecut.cfg_min_dcaxy_bachelor;
  Filter cascadeFilter_dca = o2::aod::cascdata::dcacascdaughters < cascadecut.cfg_max_dcadau && o2::aod::cascdata::dcaV0daughters < cascadecut.cfg_max_dcadau_v0;
  using filteredMyCascades = soa::Filtered<MyCascades>;

  std::vector<int> electronIdsLoose;
  std::vector<int> positronIdsLoose;
  std::vector<int> electronIds;
  std::vector<int> positronIds;

  std::vector<int> lambdaIds;
  std::vector<int> antilambdaIds;

  std::vector<int> xiPlusIds;
  std::vector<int> xiMinusIds;
  std::vector<int> omegaPlusIds;
  std::vector<int> omegaMinusIds;

  std::vector<std::pair<int, int>> used_electronIds; // pair of hTypeId and electronId

  void processMC(FilteredMyCollisionsWithMCLabel const& collisions, aod::BCsWithTimestamps const& bcs, MyTracksWithMCLabel const& tracks, aod::TrackAssoc const& trackIndices, filteredV0s const& v0s, filteredMyCascades const& cascades, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    runPairing<true>(bcs, collisions, tracks, trackIndices, v0s, cascades, mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(taggingHFE, processMC, "process with TTCA", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taggingHFE>(cfgc, TaskName{"tagging-hfe"})};
}
