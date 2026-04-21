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
#include "PWGEM/Dilepton/Utils/SemiCharmTag.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <math.h>

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
  using MyTracksWithMCLabel = soa::Join<MyTracks, aod::McTrackLabels, aod::mcTPCTuneOnData>;

  // using MyV0s = soa::Join<aod::V0Datas, aod::V0Covs, aod::V0TOFNSigmas, aod::V0CoreMCLabels>;
  // using MyCascades = soa::Join<aod::CascDatas, aod::CascCovs, aod::CascTOFNSigmas, aod::CascCoreMCLabels>;
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
  Configurable<int> cfgPdgLepton{"cfgPdgLepton", 11, "pdg code of desired lepton: 11 or 13"};

  struct : ConfigurableGroup {
    std::string prefix = "dcaFitterGroup_eK";
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  } dcaFitterGroup_eK;

  struct : ConfigurableGroup {
    std::string prefix = "dcaFitterGroup_eV0";
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  } dcaFitterGroup_eV0;

  struct : ConfigurableGroup {
    std::string prefix = "dcaFitterGroup_eCascade";
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  } dcaFitterGroup_eCascade;

  struct : ConfigurableGroup {
    std::string prefix = "electronCut";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.4, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.5, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.5, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 3, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    // Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2, "min n sigma el in TPC"};
    // Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3, "max n sigma el in TPC"};
    // Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3, "min n sigma el in TOF"};
    // Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3, "max n sigma el in TOF"};
  } electronCut;

  struct : ConfigurableGroup {
    std::string prefix = "kaonCut";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 4, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -2, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +2, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -2, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +2, "max n sigma ka in TOF"};
    Configurable<bool> requireTOF{"requireTOF", true, "require TOF hit"};
    Configurable<float> cfg_min_pin_TOFreq{"cfg_min_pin_TOFreq", 0.4, "min pin for TOFreq"};
  } kaonCut;

  struct : ConfigurableGroup {
    std::string prefix = "v0Cut";
    Configurable<float> cfg_min_mass_k0s{"cfg_min_mass_k0s", 0.49, "min mass for K0S"};
    Configurable<float> cfg_max_mass_k0s{"cfg_max_mass_k0s", 0.51, "max mass for K0S"};
    Configurable<float> cfg_min_mass_k0s_veto{"cfg_min_mass_k0s_veto", 0.47, "min mass for K0S veto for Lambda"};
    Configurable<float> cfg_max_mass_k0s_veto{"cfg_max_mass_k0s_veto", 0.52, "max mass for K0S veto for Lambda"};
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.113, "min mass for Lambda"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.118, "max mass for Lambda"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.9, "min cospa for v0"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 0.1, "max distance between 2 legs for v0"};
    // Configurable<float> cfg_min_radius{"cfg_min_radius", 0.1, "min rxy for v"};
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
  } v0Cut;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeCut";
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
    Configurable<float> cfg_min_cospa_v0{"cfg_min_cospa_v0", 0.9, "minimum V0 CosPA in cascade"};
    Configurable<float> cfg_max_dcadau_v0{"cfg_max_dcadau_v0", 0.1, "max distance between V0 Daughters in cascade"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.9, "minimum cascade CosPA"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.1, "max distance between bachelor and V0"};
    Configurable<float> cfg_min_rxy_v0{"cfg_min_rxy_v0", 1.2, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_rxy{"cfg_min_rxy", 0.5, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};
    Configurable<float> cfg_min_dcaxy_bachelor{"cfg_min_dcaxy_bachelor", 0.05, "min dca XY for bachelor in cm"};
    Configurable<float> cfg_min_dcaxy_v0{"cfg_min_dcaxy_v0", 0.05, "min dca XY for V0 in cm"};
  } cascadeCut;

  struct : ConfigurableGroup {
    std::string prefix = "eventCut";
    Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
    Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
    Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  } eventCut;

  struct : ConfigurableGroup {
    std::string prefix = "lKPairCut";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0, "min mass at SV"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass at SV"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", -1e+10, "min cospa"};
    Configurable<float> cfg_max_lxyz{"cfg_max_lxyz", 1e+10, "min rxy for v0hadron"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 1e+10, "max distance between 2 legs"};
  } lKPairCut;

  struct : ConfigurableGroup {
    std::string prefix = "lV0PairCut";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0, "min mass at SV"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass at SV"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", -1e+10, "min cospa"};
    Configurable<float> cfg_max_lxyz{"cfg_max_lxyz", 1e+10, "min rxy for v0hadron"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 1e+10, "max distance between 2 legs"};
  } lV0PairCut;

  struct : ConfigurableGroup {
    std::string prefix = "lCPairCut";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0, "min mass at SV"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 1e+10, "max mass at SV"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", -1e+10, "min cospa"};
    Configurable<float> cfg_max_lxyz{"cfg_max_lxyz", 1e+10, "min rxy for v0hadron"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 1e+10, "max distance between 2 legs"};
  } lCPairCut;

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

    fitter_eK.setPropagateToPCA(true);
    fitter_eK.setMaxR(20.f);
    fitter_eK.setMinParamChange(1e-3);
    fitter_eK.setMinRelChi2Change(0.9);
    fitter_eK.setMaxDZIni(1e9);
    fitter_eK.setMaxChi2(1e9);
    fitter_eK.setUseAbsDCA(dcaFitterGroup_eK.useAbsDCA);
    fitter_eK.setWeightedFinalPCA(dcaFitterGroup_eK.useWeightedFinalPCA);
    fitter_eK.setMatCorrType(matCorr);

    fitter_eV0.setPropagateToPCA(true);
    fitter_eV0.setMaxR(20.f);
    fitter_eV0.setMinParamChange(1e-3);
    fitter_eV0.setMinRelChi2Change(0.9);
    fitter_eV0.setMaxDZIni(1e9);
    fitter_eV0.setMaxChi2(1e9);
    fitter_eV0.setUseAbsDCA(dcaFitterGroup_eV0.useAbsDCA);
    fitter_eV0.setWeightedFinalPCA(dcaFitterGroup_eV0.useWeightedFinalPCA);
    fitter_eV0.setMatCorrType(matCorr);

    fitter_eCascade.setPropagateToPCA(true);
    fitter_eCascade.setMaxR(20.f);
    fitter_eCascade.setMinParamChange(1e-3);
    fitter_eCascade.setMinRelChi2Change(0.9);
    fitter_eCascade.setMaxDZIni(1e9);
    fitter_eCascade.setMaxChi2(1e9);
    fitter_eCascade.setUseAbsDCA(dcaFitterGroup_eCascade.useAbsDCA);
    fitter_eCascade.setWeightedFinalPCA(dcaFitterGroup_eCascade.useWeightedFinalPCA);
    fitter_eCascade.setMatCorrType(matCorr);

    addHistograms();
  }

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter_eK;
  o2::vertexing::DCAFitterN<2> fitter_eV0;
  o2::vertexing::DCAFitterN<2> fitter_eCascade;
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
    fitter_eK.setBz(d_bz);
    fitter_eV0.setBz(d_bz);
    fitter_eCascade.setBz(d_bz);
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

    fRegistry.add("Generated/Dpm/hsAcc", "pT-#eta acc.;p_{T,l} (GeV/c);p_{T,K} (GeV/c);#eta_{l};#eta_{K};", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, -5, +5}, {100, -5, +5}}, false);
    fRegistry.add("Generated/D0/hsAcc", "pT-#eta acc.;p_{T,l} (GeV/c);p_{T,K} (GeV/c);#eta_{l};#eta_{K};", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, -5, +5}, {100, -5, +5}}, false);
    fRegistry.add("Generated/Lc/hsAcc", "pT-#eta acc.;p_{T,l} (GeV/c);p_{T,#Lambda} (GeV/c);#eta_{l};#eta_{#Lambda};", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, -5, +5}, {100, -5, +5}}, false);

    fRegistry.add("Track/Electron/hTPCdEdx", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/Electron/hTOFbeta", "TOF #beta vs. p;p_{pv} (GeV/c);TOF #beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    fRegistry.add("Track/Kaon/hTPCdEdx", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("Track/Kaon/hTOFbeta", "TOF #beta vs. p;p_{pv} (GeV/c);TOF #beta", kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);

    // electron-related histograms
    fRegistry.add("Data/electron/hs", "hs;p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);DCA_{e}^{3D} (#sigma);", kTHnSparseF, {{100, 0, 10}, {20, -1, +1}, {90, 0, 2 * M_PI}, {100, 0, 10}}, false);
    fRegistry.addClone("Data/electron/", "MC/eFromPromptLF/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptLF/");
    fRegistry.addClone("Data/electron/", "MC/eFromPromptJpsi/");
    fRegistry.addClone("Data/electron/", "MC/eFromNonPromptJpsi/");

    fRegistry.addClone("Data/electron/", "MC/eFromD0/");
    fRegistry.addClone("Data/electron/", "MC/eFromDpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromDs/");
    fRegistry.addClone("Data/electron/", "MC/eFromLcpm/");
    fRegistry.addClone("Data/electron/", "MC/eFromXic0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromXicpm/"); // cannot be detected
    fRegistry.addClone("Data/electron/", "MC/eFromOmegac0/");

    // fRegistry.addClone("Data/electron/", "MC/eFromPromptD0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptDpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptDs/");
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptLcpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptXic0/");
    // // fRegistry.addClone("Data/electron/", "MC/eFromPromptXicpm/"); // cannot be detected
    // fRegistry.addClone("Data/electron/", "MC/eFromPromptOmegac0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptD0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptDpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptDs/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptLcpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptXic0/");
    // // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptXicpm/"); // cannot be detected
    // fRegistry.addClone("Data/electron/", "MC/eFromNonPromptOmegac0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromB0/");
    // fRegistry.addClone("Data/electron/", "MC/eFromBpm/");
    // fRegistry.addClone("Data/electron/", "MC/eFromBs/");
    // fRegistry.addClone("Data/electron/", "MC/eFromBc/");
    // fRegistry.addClone("Data/electron/", "MC/eFromLb0/");

    // for V0 (Lambda)
    fRegistry.add("Data/V0/hPt", "pT of V0;p_{T} (GeV/c)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/V0/hYPhi", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{#Lambda}", kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    fRegistry.add("Data/V0/hAP", "Ap plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1, 1}, {250, 0, 0.25}}, false);
    fRegistry.add("Data/V0/hLxy", "decay length from PV;L_{xy} (cm)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/V0/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{100, 0.9, 1}}, false);
    fRegistry.add("Data/V0/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("Data/V0/hMassK0S", "K0S mass;m_{#pi#pi} (GeV/c^{2})", kTH1F, {{100, 0.45, 0.55}}, false);
    fRegistry.add("Data/V0/hMassLambda", "Lambda mass;m_{p#pi^{#minus}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);
    fRegistry.add("Data/V0/hMassAntiLambda", "Anti-Lambda mass;m_{#bar{p}#pi^{+}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);

    // for cascade
    fRegistry.add("Data/Cascade/hPt", "pT of V0;p_{T} (GeV/c)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Data/Cascade/hYPhi", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{#Lambda}", kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    fRegistry.add("Data/Cascade/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{100, 0.9, 1}}, false);
    fRegistry.add("Data/Cascade/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("Data/Cascade/hV0CosPA", "cosPA of V0 in cascade;cosine of pointing angle", kTH1F, {{100, 0.9, 1}}, false);
    fRegistry.add("Data/Cascade/hV0DCA2Legs", "distance between 2 legs at PCA of V0 in cascade;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);

    fRegistry.add("Data/Cascade/hMassLambda", "Lambda mass;m_{p#pi^{-}} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);
    fRegistry.add("Data/Cascade/hMassXi", "#Xi mass;m_{#Lambda#pi} (GeV/c^{2})", kTH1F, {{100, 1.27, 1.37}}, false);
    fRegistry.add("Data/Cascade/hMassOmega", "#Omega mass;m_{#LambdaK} (GeV/c^{2})", kTH1F, {{100, 1.62, 1.72}}, false);

    // for e-K pair
    fRegistry.add("Data/eK/hs", "hs;p_{T,l} (GeV/c);DCA_{l}^{3D} (#sigma);p_{T,K} (GeV/c);DCA_{K}^{3D} (#sigma);m_{eK} (GeV/c^{2});L_{xyz} (#sigma);cosPA;DCA 2 legs (cm);", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {90, 0.5, 5.0}, {100, 0, 10}, {200, -1, 1}, {500, 0.0, 0.5}}, false);
    fRegistry.addClone("Data/eK/", "MC/eKfromD0/");
    fRegistry.addClone("Data/eK/", "MC/eKfromDpm/");

    // for e-K0 pair
    fRegistry.add("Data/eK0/hs", "hs;p_{T,l} (GeV/c);DCA_{l}^{3D} (#sigma);p_{T,K0} (GeV/c);L_{xyz}^{K0} (cm);m_{eK0} (GeV/c^{2});L_{xyz} (#sigma);cosPA;DCA 2 legs (cm);", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {90, 0.5, 5.0}, {100, 0, 10}, {200, -1, 1}, {500, 0.0, 0.5}}, false);
    fRegistry.addClone("Data/eK0/", "MC/eK0fromD0/");
    fRegistry.addClone("Data/eK0/", "MC/eK0fromDpm/");
    // fRegistry.addClone("Data/eK0/", "MC/eK0fromDspm/");

    // for e-L pair
    fRegistry.add("Data/eL/hs", "hs;p_{T,l} (GeV/c);DCA_{l}^{3D} (#sigma);p_{T,#Lambda} (GeV/c);L_{xyz}^{#Lambda} (cm);m_{e#Lambda} (GeV/c^{2});L_{xyz} (#sigma);cosPA;DCA 2 legs (cm);", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {90, 0.5, 5.0}, {100, 0, 10}, {200, -1, 1}, {500, 0.0, 0.5}}, false);
    fRegistry.addClone("Data/eL/", "MC/eLfromLcpm/");

    // for e-Xi pair
    fRegistry.add("Data/eXi/hs", "hs;p_{T,l} (GeV/c);DCA_{l}^{3D} (#sigma);p_{T,#Xi} (GeV/c);L_{xyz}^{#Xi} (cm);m_{e#Xi} (GeV/c^{2});L_{xyz} (#sigma);cosPA;DCA 2 legs (cm);", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {90, 0.5, 5.0}, {100, 0, 10}, {200, -1, 1}, {500, 0.0, 0.5}}, false);
    fRegistry.addClone("Data/eXi/", "MC/eXifromXic0/");

    // for e-Omega pair
    fRegistry.add("Data/eOmega/hs", "hs;p_{T,l} (GeV/c);DCA_{l}^{3D} (#sigma);p_{T,#Omega} (GeV/c);L_{xyz}^{#Omega} (cm);m_{e#Omega} (GeV/c^{2});L_{#sigma} (cm);cosPA;DCA 2 legs (cm);", kTHnSparseF, {{100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {100, 0, 10}, {90, 0.5, 5.0}, {100, 0, 10}, {200, -1, 1}, {500, 0.0, 0.5}}, false);
    fRegistry.addClone("Data/eOmega/", "MC/eOmegafromOmegac0/");
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    // TOFif
    bool is_ka_included_TPC = kaonCut.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < kaonCut.cfg_max_TPCNsigmaKa;
    bool is_ka_included_TOF = track.hasTOF() ? (kaonCut.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < kaonCut.cfg_max_TOFNsigmaKa) : true;
    if (kaonCut.requireTOF) {
      if (track.tpcInnerParam() < kaonCut.cfg_min_pin_TOFreq) {
        return is_ka_included_TPC;
      } else {
        is_ka_included_TOF = kaonCut.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < kaonCut.cfg_max_TOFNsigmaKa;
        return is_ka_included_TPC && is_ka_included_TOF;
      }
    } else {
      return is_ka_included_TPC && is_ka_included_TOF; // TOFif
    }
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    // TOFif
    bool is_pi_included_TPC = v0Cut.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0Cut.cfg_max_TPCNsigmaPi;
    bool is_pi_included_TOF = track.hasTOF() ? (v0Cut.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < v0Cut.cfg_max_TOFNsigmaPi) : true;
    return is_pi_included_TPC && is_pi_included_TOF;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    // TOFif
    bool is_pr_included_TPC = v0Cut.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0Cut.cfg_max_TPCNsigmaPr;
    bool is_pr_included_TOF = track.hasTOF() ? (v0Cut.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < v0Cut.cfg_max_TOFNsigmaPr) : true;
    return is_pr_included_TPC && is_pr_included_TOF;
  }

  // template <typename TTrack>
  // bool isElectron(TTrack const& track)
  // {
  //   // TOFif
  //   bool is_el_included_TPC = electronCut.cfg_min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < electronCut.cfg_max_TPCNsigmaEl;
  //   bool is_el_included_TOF = track.hasTOF() ? (electronCut.cfg_min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < electronCut.cfg_max_TOFNsigmaEl) : true;
  //   return is_el_included_TPC && is_el_included_TOF;
  // }

  template <typename TTrack>
  bool isKaonBachelor(TTrack const& track)
  {
    // TOFif
    bool is_ka_included_TPC = kaonCut.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < kaonCut.cfg_max_TPCNsigmaKa;
    bool is_ka_included_TOF = track.hasTOF() ? (kaonCut.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < kaonCut.cfg_max_TOFNsigmaKa) : true;
    return is_ka_included_TPC && is_ka_included_TOF;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedTrack(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < electronCut.cfg_min_pt_track || electronCut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < electronCut.cfg_min_eta_track || electronCut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > electronCut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > electronCut.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || electronCut.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < electronCut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < electronCut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || electronCut.cfg_max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < electronCut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < electronCut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < electronCut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > electronCut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    // if (!isElectron(track)) {
    //   return false;
    // }

    return true;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedKaon(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < kaonCut.cfg_min_pt_track || kaonCut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < kaonCut.cfg_min_eta_track || kaonCut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > kaonCut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > kaonCut.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || kaonCut.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < kaonCut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < kaonCut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || kaonCut.cfg_max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < kaonCut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < kaonCut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < kaonCut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > kaonCut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    if (!isKaon(track)) {
      return false;
    }

    return true;
  }

  template <typename TV0>
  bool isK0S(TV0 const& v0)
  {
    return (v0Cut.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < v0Cut.cfg_max_mass_k0s);
  }

  template <typename TV0>
  bool isLambda(TV0 const& v0)
  {
    return (v0Cut.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < v0Cut.cfg_max_mass_lambda) && (v0.mK0Short() < v0Cut.cfg_min_mass_k0s_veto || v0Cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TV0>
  bool isAntiLambda(TV0 const& v0)
  {
    return (v0Cut.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < v0Cut.cfg_max_mass_lambda) && (v0.mK0Short() < v0Cut.cfg_min_mass_k0s_veto || v0Cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TCascade>
  bool isXi(TCascade const& cascade)
  {
    return (cascadeCut.cfg_min_mass_Xi < cascade.mXi() && cascade.mXi() < cascadeCut.cfg_max_mass_Xi) && (cascade.mOmega() < cascadeCut.cfg_min_mass_Omega_veto || cascadeCut.cfg_max_mass_Omega_veto < cascade.mOmega());
  }

  template <typename TCascade>
  bool isOmega(TCascade const& cascade)
  {
    return (cascadeCut.cfg_min_mass_Omega < cascade.mOmega() && cascade.mOmega() < cascadeCut.cfg_max_mass_Omega) && (cascade.mXi() < cascadeCut.cfg_min_mass_Xi_veto || cascadeCut.cfg_max_mass_Xi_veto < cascade.mXi());
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

    if (track.itsChi2NCl() > v0Cut.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < v0Cut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < v0Cut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > v0Cut.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < v0Cut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < v0Cut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < v0Cut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > v0Cut.cfg_max_frac_shared_clusters_tpc) {
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
    auto pos = v0.template posTrack_as<MyTracksWithMCLabel>();
    auto neg = v0.template negTrack_as<MyTracksWithMCLabel>();
    fRegistry.fill(HIST("Data/V0/hPt"), v0.pt());
    fRegistry.fill(HIST("Data/V0/hYPhi"), v0.phi(), v0.yLambda());
    fRegistry.fill(HIST("Data/V0/hAP"), v0.alpha(), v0.qtarm());
    fRegistry.fill(HIST("Data/V0/hCosPA"), v0.v0cosPA());
    fRegistry.fill(HIST("Data/V0/hLxy"), v0.v0radius());
    fRegistry.fill(HIST("Data/V0/hDCA2Legs"), v0.dcaV0daughters());

    if (isPion(pos) && isPion(neg)) {
      fRegistry.fill(HIST("Data/V0/hMassK0S"), v0.mK0Short());
    }
    if (isProton(pos) && isPion(neg)) {
      fRegistry.fill(HIST("Data/V0/hMassLambda"), v0.mLambda());
    }
    if (isProton(neg) && isPion(pos)) {
      fRegistry.fill(HIST("Data/V0/hMassAntiLambda"), v0.mAntiLambda());
    }
  }

  template <typename TCollision, typename TCascade>
  void fillCascadeHistograms(TCollision const& collision, TCascade const& cascade)
  {
    auto bachelor = cascade.template bachelor_as<MyTracksWithMCLabel>();
    fRegistry.fill(HIST("Data/Cascade/hPt"), cascade.pt());
    fRegistry.fill(HIST("Data/Cascade/hMassLambda"), cascade.mLambda());
    fRegistry.fill(HIST("Data/Cascade/hCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
    fRegistry.fill(HIST("Data/Cascade/hDCA2Legs"), cascade.dcacascdaughters());
    fRegistry.fill(HIST("Data/Cascade/hV0CosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
    fRegistry.fill(HIST("Data/Cascade/hV0DCA2Legs"), cascade.dcaV0daughters());

    if (isPion(bachelor)) {
      fRegistry.fill(HIST("Data/Cascade/hMassXi"), cascade.mXi());
    }
    if (isKaonBachelor(bachelor)) {
      fRegistry.fill(HIST("Data/Cascade/hMassOmega"), cascade.mOmega());
    }
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

    if (!isSelectedTrack(track, trackParCov, dcaXY, dcaZ)) {
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
        fRegistry.fill(HIST("MC/eFromDpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptDpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptDpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      } else if (pdg_mother == 421) { // D0
        fRegistry.fill(HIST("MC/eFromD0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptD0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptD0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      } else if (pdg_mother == 431) { // Ds+/-
        // fRegistry.fill(HIST("MC/eFromDs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptDs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptDs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      } else if (pdg_mother == 4122) { // Lc+/-
        fRegistry.fill(HIST("MC/eFromLcpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptLcpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptLcpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      } else if (pdg_mother == 4132) { // Xic0
        fRegistry.fill(HIST("MC/eFromXic0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptXic0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptXic0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      } else if (pdg_mother == 4332) { // Omegac0
        fRegistry.fill(HIST("MC/eFromOmegac0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // if (IsFromBeauty(mcmother, mcParticles) < 0) {
        //   fRegistry.fill(HIST("MC/eFromPromptOmegac0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // } else {
        //   fRegistry.fill(HIST("MC/eFromNonPromptOmegac0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
        // }
      }
      // else if (pdg_mother == 511) { // B0
      //   fRegistry.fill(HIST("MC/eFromB0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      // } else if (pdg_mother == 521) { // B+/-
      //   fRegistry.fill(HIST("MC/eFromBpm/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      // } else if (pdg_mother == 531) { // Bs0
      //   fRegistry.fill(HIST("MC/eFromBs/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      // } else if (pdg_mother == 541) { // Bc+/-
      //   fRegistry.fill(HIST("MC/eFromBc/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      // } else if (pdg_mother == 5122) { // Lb0
      //   fRegistry.fill(HIST("MC/eFromLb0/hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      // }
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

  template <typename TMCParticle, typename TMCParticles>
  bool isSemiLeptonic(TMCParticle const& mcParticle, TMCParticles const& mcParticles, const int pdgLepton, const int pdgNeutrino, const int pdgStrHad)
  {
    if (!mcParticle.has_daughters()) {
      return false;
    }
    bool is_lepton_involved = false;
    bool is_neutrino_involved = false;
    bool is_strhad_involved = false;
    for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
      if (d < mcParticles.size()) { // protect against bad daughter indices
        const auto& daughter = mcParticles.rawIteratorAt(d);
        if (daughter.pdgCode() == pdgLepton) {
          is_lepton_involved = true;
        } else if (daughter.pdgCode() == pdgNeutrino) {
          is_neutrino_involved = true;
        } else if (daughter.pdgCode() == pdgStrHad) {
          is_strhad_involved = true;
        }

      } else {
        std::cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")" << std::endl;
        std::cout << " Check the MC generator" << std::endl;
        return false;
      }
    }

    if (is_lepton_involved && is_neutrino_involved && is_strhad_involved) {
      return true;
    } else {
      return false;
    }
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
      if (centralities[eventCut.cfgCentEstimator] < eventCut.cfgCentMin || eventCut.cfgCentMax < centralities[eventCut.cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      const auto& mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (eventCut.cfgEventGeneratorType < 0 || mcCollision.getSubGeneratorId() == eventCut.cfgEventGeneratorType) {
        fillEventHistograms(collision);
      }
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      electronIds.reserve(trackIdsThisCollision.size());
      positronIds.reserve(trackIdsThisCollision.size());
      kaonPlusIds.reserve(trackIdsThisCollision.size());
      kaonMinusIds.reserve(trackIdsThisCollision.size());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = trackId.template track_as<TTracks>();
        if (!track.hasITS() || !track.hasTPC()) {
          continue;
        }

        if (!track.has_mcParticle()) {
          continue;
        }
        const auto& mcParticle = track.template mcParticle_as<aod::McParticles>();
        const auto& mcCollision = mcParticle.template mcCollision_as<aod::McCollisions>();
        if (eventCut.cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventCut.cfgEventGeneratorType) {
          continue;
        }
        if (!mcParticle.has_mothers() || !(mcParticle.isPhysicalPrimary() || mcParticle.producedByGenerator())) {
          continue;
        }

        fillElectronHistograms<isMC>(track, mcParticles);

        auto trackParCov = getTrackParCov(track);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        trackParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();

        if (isSelectedTrack(track, trackParCov, dcaXY, dcaZ) && std::abs(mcParticle.pdgCode()) == cfgPdgLepton) {
          fRegistry.fill(HIST("Track/Electron/hTPCdEdx"), track.tpcInnerParam(), track.mcTunedTPCSignal());
          fRegistry.fill(HIST("Track/Electron/hTOFbeta"), track.p(), track.beta());
          if (track.sign() > 0) { // positron
            positronIds.emplace_back(trackId.trackId());
          } else { // electron
            electronIds.emplace_back(trackId.trackId());
          }
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        trackParCov.setPID(o2::track::PID::Kaon);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        dcaXY = mDcaInfoCov.getY();
        dcaZ = mDcaInfoCov.getZ();

        if (isSelectedKaon(track, trackParCov, dcaXY, dcaZ)) {
          fRegistry.fill(HIST("Track/Kaon/hTPCdEdx"), track.tpcInnerParam(), track.mcTunedTPCSignal());
          fRegistry.fill(HIST("Track/Kaon/hTOFbeta"), track.p(), track.beta());
          if (track.sign() > 0) { // K+
            kaonPlusIds.emplace_back(trackId.trackId());
          } else { // K-
            kaonMinusIds.emplace_back(trackId.trackId());
          }
        }

      } // end of track loop for electron selection

      const auto& v0s_per_coll = v0s.sliceBy(perCol_v0, collision.globalIndex());
      k0Ids.reserve(v0s_per_coll.size());
      k0Ids.reserve(v0s_per_coll.size());
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
        fillV0Histograms<>(v0);

        if (isK0S(v0) && isPion(pos) && isPion(neg)) {
          k0Ids.emplace_back(v0.globalIndex());
        }

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
        if (cascade.mLambda() < cascadeCut.cfg_min_mass_lambda || cascadeCut.cfg_max_mass_lambda < cascade.mLambda()) {
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

        if (cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeCut.cfg_min_cospa) {
          continue;
        }
        if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadeCut.cfg_min_cospa_v0) {
          continue;
        }

        fillCascadeHistograms(collision, cascade);
        if (cascade.sign() < 0) { // Xi- or Omega-
          if (isXi(cascade) && isPion(bachelor)) {
            xiMinusIds.emplace_back(cascade.globalIndex());
          }
          if (isOmega(cascade) && isKaonBachelor(bachelor)) {
            omegaMinusIds.emplace_back(cascade.globalIndex());
          }
        } else { // Xi+ or Omega+
          if (isXi(cascade) && isPion(bachelor)) {
            xiPlusIds.emplace_back(cascade.globalIndex());
          }
          if (isOmega(cascade) && isKaonBachelor(bachelor)) {
            omegaPlusIds.emplace_back(cascade.globalIndex());
          }
        }
      } // end of cascade loop

      // // if (electronIds.size() > 0 || positronIds.size() > 0) {
      // if ((electronIds.size() > 0 || positronIds.size() > 0) && (xiMinusIds.size() > 0 || xiPlusIds.size() > 0)) {
      // LOGF(info, "collision.globalIndex() = %d, electronIds.size() = %d, positronIds.size() = %d, kaonMinusIds.size() = %d, kaonPlusIds.size() = %d, k0Ids.size() = %d, lambdaIds.size() = %d, antilambdaIds.size() = %d, xiMinusIds.size() = %d, xiPlusIds.size() = %d, omegaMinusIds.size() = %d, omegaPlusIds.size() = %d",
      //     collision.globalIndex(), electronIds.size(), positronIds.size(), kaonMinusIds.size(), kaonPlusIds.size(), k0Ids.size(), lambdaIds.size(), antilambdaIds.size(), xiMinusIds.size(), xiPlusIds.size(), omegaMinusIds.size(), omegaPlusIds.size());
      // }

      for (const auto& positronId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(positronId);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto leptonParCov = getTrackParCov(pos);
        leptonParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, leptonParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY_lepton = mDcaInfoCov.getY();
        float dcaZ_lepton = mDcaInfoCov.getZ();
        float dca3DinSigma_lepton = dca3DinSigmaOTF(dcaXY_lepton, dcaZ_lepton, leptonParCov.getSigmaY2(), leptonParCov.getSigmaZ2(), leptonParCov.getSigmaZY());
        const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();

        // D0 -> e+ nu_e K-, br = 0.03538, ctau = 123.01 um, m = 1864 MeV/c2
        for (const auto& kaonId : kaonMinusIds) {
          const auto& kaon = tracks.rawIteratorAt(kaonId);
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto trackParCov = getTrackParCov(kaon);
          trackParCov.setPID(o2::track::PID::Kaon);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
          float dcaXY_kaon = mDcaInfoCov.getY();
          float dcaZ_kaon = mDcaInfoCov.getZ();
          float dca3DinSigma_kaon = dca3DinSigmaOTF(dcaXY_kaon, dcaZ_kaon, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());

          if (positronId == kaonId) {
            continue;
          }

          const auto& eKpair = o2::aod::pwgem::dilepton::utils::makePairLeptonTrack(fitter_eK, collision, pos, kaon, o2::track::PID::Electron, o2::track::PID::Kaon);
          if (!eKpair.isOK) {
            continue;
          }

          if (!(lKPairCut.cfg_min_mass < eKpair.mass && eKpair.mass < lKPairCut.cfg_max_mass) || eKpair.cospa < lKPairCut.cfg_min_cospa || lKPairCut.cfg_max_lxyz < eKpair.lxyz / eKpair.lxyzErr || lKPairCut.cfg_max_dca2legs < eKpair.dca2legs) {
            continue;
          }

          fRegistry.fill(HIST("Data/eK/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);

          const auto& mckaon = kaon.template mcParticle_as<aod::McParticles>();
          int mcD0Id = FindCommonMotherFrom2Prongs(mcpos, mckaon, -11, -321, 421, mcParticles);
          int mcDpmId = FindCommonMotherFrom2Prongs(mcpos, mckaon, -11, -321, 411, mcParticles);
          if (mcD0Id > 0) { // true D0
            // const auto& mcD0 = mcParticles.rawIteratorAt(mcD0Id);
            fRegistry.fill(HIST("MC/eKfromD0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);
          } else if (mcDpmId > 0) { // true D+
            // const auto& mcD0 = mcParticles.rawIteratorAt(mcD0Id);
            fRegistry.fill(HIST("MC/eKfromDpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);
          }
        } // end of kaon loop

        // D+ -> e+ K0S nu_e
        for (const auto& k0Id : k0Ids) {
          const auto& v0 = v0s.rawIteratorAt(k0Id);
          float lxyz_v0 = std::sqrt(std::pow(v0.x() - collision.posX(), 2) + std::pow(v0.y() - collision.posY(), 2) + std::pow(v0.z() - collision.posZ(), 2));
          const auto& eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(fitter_eV0, collision, pos, v0, o2::track::PID::Electron, o2::track::PID::K0);

          if (!eV0pair.isOK) {
            continue;
          }
          if (!(lV0PairCut.cfg_min_mass < eV0pair.mass && eV0pair.mass < lV0PairCut.cfg_max_mass) || eV0pair.cospa < lV0PairCut.cfg_min_cospa || lV0PairCut.cfg_max_lxyz < eV0pair.lxyz / eV0pair.lxyzErr || lV0PairCut.cfg_max_dca2legs < eV0pair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eK0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);

          auto posLeg = v0.template posTrack_as<TTracks>();
          auto negLeg = v0.template negTrack_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          int mcK0Id = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -211, 310, mcParticles);
          if (mcK0Id > 0) { // true K0S
            const auto& mcK0 = mcParticles.rawIteratorAt(mcK0Id);
            int mcDpmId = FindCommonMotherFrom2Prongs(mcpos, mcK0, -11, 310, 411, mcParticles);
            if (mcDpmId > 0) { // true D+
              // const auto& mcDpm = mcParticles.rawIteratorAt(mcDpmId);
              fRegistry.fill(HIST("MC/eK0fromDpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);
            }
          }
        } // end of K0S loop

        // Lc+ -> e+ Lambda nu_e, br = 0.0356, ctau = 60.75 um, m = 2286 MeV/c2
        for (const auto& lambdaId : lambdaIds) {
          const auto& v0 = v0s.rawIteratorAt(lambdaId);
          float lxyz_v0 = std::sqrt(std::pow(v0.x() - collision.posX(), 2) + std::pow(v0.y() - collision.posY(), 2) + std::pow(v0.z() - collision.posZ(), 2));
          const auto& eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(fitter_eV0, collision, pos, v0, o2::track::PID::Electron, o2::track::PID::Lambda);

          if (!eV0pair.isOK) {
            continue;
          }
          if (!(lV0PairCut.cfg_min_mass < eV0pair.mass && eV0pair.mass < lV0PairCut.cfg_max_mass) || eV0pair.cospa < lV0PairCut.cfg_min_cospa || lV0PairCut.cfg_max_lxyz < eV0pair.lxyz / eV0pair.lxyzErr || lV0PairCut.cfg_max_dca2legs < eV0pair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eL/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);

          auto posLeg = v0.template posTrack_as<TTracks>();
          auto negLeg = v0.template negTrack_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 2212, -211, 3122, mcParticles);
          if (mcLambdaId > 0) { // true v0
            const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
            int mcLambdacId = FindCommonMotherFrom2Prongs(mcpos, mcLambda, -11, 3122, 4122, mcParticles);
            if (mcLambdacId > 0) { // true Lc0
              // const auto& mcLambdac0 = mcParticles.rawIteratorAt(mcLambdacId);
              fRegistry.fill(HIST("MC/eLfromLcpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);
            }
          }
        } // end of Lambda loop

        for (const auto& cascadeId : xiMinusIds) {
          const auto& cascade = cascades.rawIteratorAt(cascadeId);
          float lxyz_cascade = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
          const auto& eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(fitter_eCascade, collision, pos, cascade, o2::track::PID::Electron, o2::track::PID::XiMinus);

          if (!eCpair.isOK) {
            continue;
          }
          if (!(lCPairCut.cfg_min_mass < eCpair.mass && eCpair.mass < lCPairCut.cfg_max_mass) || eCpair.cospa < lCPairCut.cfg_min_cospa || lCPairCut.cfg_max_lxyz < eCpair.lxyz / eCpair.lxyzErr || lCPairCut.cfg_max_dca2legs < eCpair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eXi/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);

          auto posLeg = cascade.template posTrack_as<TTracks>();
          auto negLeg = cascade.template negTrack_as<TTracks>();
          auto bachelor = cascade.template bachelor_as<TTracks>();
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
                // const auto& mcXic0 = mcParticles.rawIteratorAt(mcXic0Id);
                fRegistry.fill(HIST("MC/eXifromXic0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);
              }
            }
          }
        } // end of Xi- loop

        for (const auto& cascadeId : omegaMinusIds) {
          const auto& cascade = cascades.rawIteratorAt(cascadeId);
          float lxyz_cascade = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
          const auto& eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(fitter_eCascade, collision, pos, cascade, o2::track::PID::Electron, o2::track::PID::OmegaMinus);

          if (!eCpair.isOK) {
            continue;
          }
          if (!(lCPairCut.cfg_min_mass < eCpair.mass && eCpair.mass < lCPairCut.cfg_max_mass) || eCpair.cospa < lCPairCut.cfg_min_cospa || lCPairCut.cfg_max_lxyz < eCpair.lxyz / eCpair.lxyzErr || lCPairCut.cfg_max_dca2legs < eCpair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eOmega/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);

          auto posLeg = cascade.template posTrack_as<TTracks>();
          auto negLeg = cascade.template negTrack_as<TTracks>();
          auto bachelor = cascade.template bachelor_as<TTracks>();
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
                // const auto& mcOmegac0 = mcParticles.rawIteratorAt(mcOmegac0Id);
                fRegistry.fill(HIST("MC/eOmegafromOmegac0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);
              }
            }
          }
        } // end of Omega- loop

      } // end of main positron sample

      for (const auto& electronId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(electronId);
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto leptonParCov = getTrackParCov(ele);
        leptonParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, leptonParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY_lepton = mDcaInfoCov.getY();
        float dcaZ_lepton = mDcaInfoCov.getZ();
        float dca3DinSigma_lepton = dca3DinSigmaOTF(dcaXY_lepton, dcaZ_lepton, leptonParCov.getSigmaY2(), leptonParCov.getSigmaZ2(), leptonParCov.getSigmaZY());
        const auto& mcele = ele.template mcParticle_as<aod::McParticles>();

        // D0bar -> e- anti-nu_e K+, br = 0.03538, ctau = 123.01 um, m = 1864 MeV/c2
        for (const auto& kaonId : kaonPlusIds) {
          const auto& kaon = tracks.rawIteratorAt(kaonId);
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto trackParCov = getTrackParCov(kaon);
          trackParCov.setPID(o2::track::PID::Kaon);
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
          float dcaXY_kaon = mDcaInfoCov.getY();
          float dcaZ_kaon = mDcaInfoCov.getZ();
          float dca3DinSigma_kaon = dca3DinSigmaOTF(dcaXY_kaon, dcaZ_kaon, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());

          if (electronId == kaonId) {
            continue;
          }

          const auto& eKpair = o2::aod::pwgem::dilepton::utils::makePairLeptonTrack(fitter_eK, collision, ele, kaon, o2::track::PID::Electron, o2::track::PID::Kaon);
          if (!eKpair.isOK) {
            continue;
          }
          if (!(lKPairCut.cfg_min_mass < eKpair.mass && eKpair.mass < lKPairCut.cfg_max_mass) || eKpair.cospa < lKPairCut.cfg_min_cospa || lKPairCut.cfg_max_lxyz < eKpair.lxyz / eKpair.lxyzErr || lKPairCut.cfg_max_dca2legs < eKpair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eK/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);

          const auto& mckaon = kaon.template mcParticle_as<aod::McParticles>();
          int mcD0Id = FindCommonMotherFrom2Prongs(mcele, mckaon, 11, 321, -421, mcParticles);
          int mcDpmId = FindCommonMotherFrom2Prongs(mcele, mckaon, 11, 321, -411, mcParticles);
          if (mcD0Id > 0) { // true D0
            // const auto& mcD0 = mcParticles.rawIteratorAt(mcD0Id);
            fRegistry.fill(HIST("MC/eKfromD0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);
          } else if (mcDpmId > 0) { // true D-
            // const auto& mcD0 = mcParticles.rawIteratorAt(mcD0Id);
            fRegistry.fill(HIST("MC/eKfromDpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, trackParCov.getPt(), dca3DinSigma_kaon, eKpair.mass, eKpair.lxyz / eKpair.lxyzErr, eKpair.cospa, eKpair.dca2legs);
          }
        } // end of kaon loop

        // D- -> e0 anti-K0S anti-nu_e
        for (const auto& k0Id : k0Ids) {
          const auto& v0 = v0s.rawIteratorAt(k0Id);
          float lxyz_v0 = std::sqrt(std::pow(v0.x() - collision.posX(), 2) + std::pow(v0.y() - collision.posY(), 2) + std::pow(v0.z() - collision.posZ(), 2));
          const auto& eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(fitter_eV0, collision, ele, v0, o2::track::PID::Electron, o2::track::PID::K0);

          if (!eV0pair.isOK) {
            continue;
          }
          if (!(lV0PairCut.cfg_min_mass < eV0pair.mass && eV0pair.mass < lV0PairCut.cfg_max_mass) || eV0pair.cospa < lV0PairCut.cfg_min_cospa || lV0PairCut.cfg_max_lxyz < eV0pair.lxyz / eV0pair.lxyzErr || lV0PairCut.cfg_max_dca2legs < eV0pair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eK0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);

          auto posLeg = v0.template posTrack_as<TTracks>();
          auto negLeg = v0.template negTrack_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          int mcK0Id = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -211, 310, mcParticles);
          if (mcK0Id > 0) { // true K0S
            const auto& mcK0 = mcParticles.rawIteratorAt(mcK0Id);
            int mcDpmId = FindCommonMotherFrom2Prongs(mcele, mcK0, 11, 310, -411, mcParticles);
            if (mcDpmId > 0) { // true D+
              // const auto& mcDpm = mcParticles.rawIteratorAt(mcDpmId);
              fRegistry.fill(HIST("MC/eK0fromDpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);
            }
          }
        } // end of K0S loop

        // Lc- -> e- anti-Lambda anti-nu_e, br = 0.0356, ctau = 60.75 um, m = 2286 MeV/c2
        for (const auto& lambdaId : antilambdaIds) {
          const auto& v0 = v0s.rawIteratorAt(lambdaId);
          float lxyz_v0 = std::sqrt(std::pow(v0.x() - collision.posX(), 2) + std::pow(v0.y() - collision.posY(), 2) + std::pow(v0.z() - collision.posZ(), 2));
          const auto& eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(fitter_eV0, collision, ele, v0, o2::track::PID::Electron, o2::track::PID::Lambda);

          if (!eV0pair.isOK) {
            continue;
          }
          if (!(lV0PairCut.cfg_min_mass < eV0pair.mass && eV0pair.mass < lV0PairCut.cfg_max_mass) || eV0pair.cospa < lV0PairCut.cfg_min_cospa || lV0PairCut.cfg_max_lxyz < eV0pair.lxyz / eV0pair.lxyzErr || lV0PairCut.cfg_max_dca2legs < eV0pair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eL/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);

          auto posLeg = v0.template posTrack_as<TTracks>();
          auto negLeg = v0.template negTrack_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -2212, -3122, mcParticles);
          if (mcLambdaId > 0) { // true v0
            const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
            int mcLambdacId = FindCommonMotherFrom2Prongs(mcele, mcLambda, 11, -3122, -4122, mcParticles);
            if (mcLambdacId > 0) { // true Lc0
              // const auto& mcLambdac0 = mcParticles.rawIteratorAt(mcLambdacId);
              fRegistry.fill(HIST("MC/eLfromLcpm/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, v0.pt(), lxyz_v0, eV0pair.mass, eV0pair.lxyz / eV0pair.lxyzErr, eV0pair.cospa, eV0pair.dca2legs);
            }
          }
        } // end of Anti-Lambda loop

        for (const auto& cascadeId : xiPlusIds) {
          const auto& cascade = cascades.rawIteratorAt(cascadeId);
          float lxyz_cascade = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
          const auto& eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(fitter_eCascade, collision, ele, cascade, o2::track::PID::Electron, o2::track::PID::XiMinus);

          if (!eCpair.isOK) {
            continue;
          }
          if (!(lCPairCut.cfg_min_mass < eCpair.mass && eCpair.mass < lCPairCut.cfg_max_mass) || eCpair.cospa < lCPairCut.cfg_min_cospa || lCPairCut.cfg_max_lxyz < eCpair.lxyz / eCpair.lxyzErr || lCPairCut.cfg_max_dca2legs < eCpair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eXi/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);

          auto posLeg = cascade.template posTrack_as<TTracks>();
          auto negLeg = cascade.template negTrack_as<TTracks>();
          auto bachelor = cascade.template bachelor_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
          int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -2212, -3122, mcParticles);
          if (mcLambdaId > 0) { // true Lambda
            const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
            int mcXiId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, -3122, 211, -3312, mcParticles);
            if (mcXiId > 0) { // true xiMinus
              const auto& mcXi = mcParticles.rawIteratorAt(mcXiId);
              int mcXic0Id = FindCommonMotherFrom2Prongs(mcele, mcXi, 11, -3312, -4132, mcParticles);
              if (mcXic0Id > 0) { // true Xic0
                // const auto& mcXic0 = mcParticles.rawIteratorAt(mcXic0Id);
                fRegistry.fill(HIST("MC/eXifromXic0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);
              }
            }
          }
        } // end of Xi+ loop

        for (const auto& cascadeId : omegaPlusIds) {
          const auto& cascade = cascades.rawIteratorAt(cascadeId);
          float lxyz_cascade = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
          const auto& eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(fitter_eCascade, collision, ele, cascade, o2::track::PID::Electron, o2::track::PID::OmegaMinus);

          if (!eCpair.isOK) {
            continue;
          }
          if (!(lCPairCut.cfg_min_mass < eCpair.mass && eCpair.mass < lCPairCut.cfg_max_mass) || eCpair.cospa < lCPairCut.cfg_min_cospa || lCPairCut.cfg_max_lxyz < eCpair.lxyz / eCpair.lxyzErr || lCPairCut.cfg_max_dca2legs < eCpair.dca2legs) {
            continue;
          }
          fRegistry.fill(HIST("Data/eOmega/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);

          auto posLeg = cascade.template posTrack_as<TTracks>();
          auto negLeg = cascade.template negTrack_as<TTracks>();
          auto bachelor = cascade.template bachelor_as<TTracks>();
          const auto& mcposLeg = posLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcnegLeg = negLeg.template mcParticle_as<aod::McParticles>();
          const auto& mcbachelor = bachelor.template mcParticle_as<aod::McParticles>();
          int mcLambdaId = FindCommonMotherFrom2Prongs(mcposLeg, mcnegLeg, 211, -2212, -3122, mcParticles);
          if (mcLambdaId > 0) { // true Lambda
            const auto& mcLambda = mcParticles.rawIteratorAt(mcLambdaId);
            int mcOmegaId = FindCommonMotherFrom2Prongs(mcLambda, mcbachelor, -3122, 321, -3334, mcParticles);
            if (mcOmegaId > 0) { // true omegaMinus
              const auto& mcOmega = mcParticles.rawIteratorAt(mcOmegaId);
              int mcOmegac0Id = FindCommonMotherFrom2Prongs(mcele, mcOmega, 11, -3334, -4332, mcParticles);
              if (mcOmegac0Id > 0) { // true Omegac0
                // const auto& mcOmegac0 = mcParticles.rawIteratorAt(mcOmegac0Id);
                fRegistry.fill(HIST("MC/eOmegafromOmegac0/hs"), leptonParCov.getPt(), dca3DinSigma_lepton, cascade.pt(), lxyz_cascade, eCpair.mass, eCpair.lxyz / eCpair.lxyzErr, eCpair.cospa, eCpair.dca2legs);
              }
            }
          }
        } // end of Omega+ loop

      } // end of main electron sample

      electronIds.clear();
      electronIds.shrink_to_fit();
      positronIds.clear();
      positronIds.shrink_to_fit();

      kaonPlusIds.clear();
      kaonPlusIds.shrink_to_fit();
      kaonMinusIds.clear();
      kaonMinusIds.shrink_to_fit();

      k0Ids.clear();
      k0Ids.shrink_to_fit();

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

  template <typename TMCCollisions, typename TMCParticles>
  void runGen(TMCCollisions const& mcCollisions, TMCParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcCollisions) {

      auto mcDpms_per_mccollision = mcDpms.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (const auto& mcParticle : mcDpms_per_mccollision) {
        // for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
        //   auto daughter = mcParticles.rawIteratorAt(d);
        //   LOGF(info, "daughter.pdgCode() = %d", daughter.pdgCode());
        // }

        if (isSemiLeptonic(mcParticle, mcParticles, -cfgPdgLepton, cfgPdgLepton + 1, -321) || isSemiLeptonic(mcParticle, mcParticles, cfgPdgLepton, -cfgPdgLepton - 1, 321)) { // D+ -> l+ nul K- pi+ or D- -> l- anti-nul K+ pi-
          // LOGF(info, "semileptonic decay is found.");
          float ptLepton = 0, ptHadron = 0, etaLepton = 999.f, etaHadron = 999.f;
          for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
            auto daughter = mcParticles.rawIteratorAt(d);
            if (std::abs(daughter.pdgCode()) == cfgPdgLepton) {
              ptLepton = daughter.pt();
              etaLepton = daughter.eta();
            } else if (std::abs(daughter.pdgCode()) == 321) {
              ptHadron = daughter.pt();
              etaHadron = daughter.eta();
            }
          }
          fRegistry.fill(HIST("Generated/Dpm/hsAcc"), ptLepton, ptHadron, etaLepton, etaHadron);
        }

      } // end of Dpm loop per mcCollision

      auto mcD0s_per_mccollision = mcD0s.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (const auto& mcParticle : mcD0s_per_mccollision) {
        // for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
        //   auto daughter = mcParticles.rawIteratorAt(d);
        //   LOGF(info, "daughter.pdgCode() = %d", daughter.pdgCode());
        // }

        if (isSemiLeptonic(mcParticle, mcParticles, -cfgPdgLepton, cfgPdgLepton + 1, -321) || isSemiLeptonic(mcParticle, mcParticles, cfgPdgLepton, -cfgPdgLepton - 1, 321)) { // D0 -> l+ nul K- or D0bar -> l- anti-nul K+
          // LOGF(info, "semileptonic decay is found.");
          float ptLepton = 0, ptHadron = 0, etaLepton = 999.f, etaHadron = 999.f;
          for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
            auto daughter = mcParticles.rawIteratorAt(d);
            if (std::abs(daughter.pdgCode()) == cfgPdgLepton) {
              ptLepton = daughter.pt();
              etaLepton = daughter.eta();
            } else if (std::abs(daughter.pdgCode()) == 321) {
              ptHadron = daughter.pt();
              etaHadron = daughter.eta();
            }
          }
          fRegistry.fill(HIST("Generated/D0/hsAcc"), ptLepton, ptHadron, etaLepton, etaHadron);
        }

      } // end of D0 loop per mcCollision

      auto mcLcs_per_mccollision = mcLcs.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (const auto& mcParticle : mcLcs_per_mccollision) {
        // for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
        //   auto daughter = mcParticles.rawIteratorAt(d);
        //   LOGF(info, "daughter.pdgCode() = %d", daughter.pdgCode());
        // }

        if (isSemiLeptonic(mcParticle, mcParticles, -std::abs(cfgPdgLepton), std::abs(cfgPdgLepton) + 1, 3122) || isSemiLeptonic(mcParticle, mcParticles, std::abs(cfgPdgLepton), -std::abs(cfgPdgLepton) - 1, -3122)) { // Lc+ -> l+ nul L or Lc- -> l- anti-nul anti-L
          // LOGF(info, "semileptonic decay is found.");
          float ptLepton = 0, ptHadron = 0, etaLepton = 999.f, etaHadron = 999.f;
          for (int d = mcParticle.daughtersIds()[0]; d <= mcParticle.daughtersIds()[1]; ++d) {
            auto daughter = mcParticles.rawIteratorAt(d);
            if (std::abs(daughter.pdgCode()) == cfgPdgLepton) {
              ptLepton = daughter.pt();
              etaLepton = daughter.eta();
            } else if (std::abs(daughter.pdgCode()) == 3122) {
              ptHadron = daughter.pt();
              etaHadron = daughter.eta();
            }
          }
          fRegistry.fill(HIST("Generated/Lc/hsAcc"), ptLepton, ptHadron, etaLepton, etaHadron);
        }

      } // end of D0 loop per mcCollision
    }
  }

  Partition<aod::McParticles> mcDpms = nabs(o2::aod::mcparticle::pdgCode) == 411;
  Partition<aod::McParticles> mcD0s = nabs(o2::aod::mcparticle::pdgCode) == 421;
  Partition<aod::McParticles> mcDspms = nabs(o2::aod::mcparticle::pdgCode) == 431;
  Partition<aod::McParticles> mcLcs = nabs(o2::aod::mcparticle::pdgCode) == 4122;
  Partition<aod::McParticles> mcXic0s = nabs(o2::aod::mcparticle::pdgCode) == 4232;
  Partition<aod::McParticles> mcOmegac0s = nabs(o2::aod::mcparticle::pdgCode) == 4332;

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;
  Preslice<aod::V0Datas> perCol_v0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> perCol_casc = o2::aod::cascdata::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && (eventCut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventCut.cfgZvtxMax);
  Filter collisionFilter_centrality = (eventCut.cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < eventCut.cfgCentMax) || (eventCut.cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < eventCut.cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;
  using FilteredMyCollisionsWithMCLabel = soa::Filtered<MyCollisionsWithMCLabel>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // std::vector<std::pair<int, int>> stored_trackIds;

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cut), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > v0Cut.cfg_min_cospa&& o2::aod::v0data::dcaV0daughters<v0Cut.cfg_max_dca2legs && nabs(o2::aod::v0data::dcanegtopv)> v0Cut.cfg_min_dcaxy&& nabs(o2::aod::v0data::dcanegtopv) > v0Cut.cfg_min_dcaxy;
  using filteredV0s = soa::Filtered<MyV0s>;

  Filter cascadeFilter = nabs(o2::aod::cascdata::dcanegtopv) > cascadeCut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcanegtopv) > cascadeCut.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcabachtopv) > cascadeCut.cfg_min_dcaxy_bachelor;
  Filter cascadeFilter_dca = o2::aod::cascdata::dcacascdaughters < cascadeCut.cfg_max_dcadau && o2::aod::cascdata::dcaV0daughters < cascadeCut.cfg_max_dcadau_v0;
  using filteredMyCascades = soa::Filtered<MyCascades>;

  std::vector<int> electronIds;
  std::vector<int> positronIds;

  std::vector<int> kaonPlusIds;
  std::vector<int> kaonMinusIds;

  std::vector<int> k0Ids;

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
    runGen(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(taggingHFE, processMC, "process with TTCA", true);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(taggingHFE, processDummy, "process dummy", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taggingHFE>(cfgc, TaskName{"tagging-hfe"})};
}
