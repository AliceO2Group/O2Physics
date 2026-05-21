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

#ifndef PWGEM_DILEPTON_UTILS_ELECTRONMODULE_H_
#define PWGEM_DILEPTON_UTILS_ELECTRONMODULE_H_

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"
#include "PWGEM/Dilepton/Utils/MlResponsePID.h"
#include "PWGEM/Dilepton/Utils/MlResponseSCT.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/Dilepton/Utils/SemiCharmTag.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
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
#include <DCAFitter/DCAFitterN.h>
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
#include <TMath.h>

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

namespace o2
{
namespace pwgem::dilepton::utils
{

struct ElectronProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<aod::EMPrimaryElectrons> electronTable;
  o2::framework::Produces<aod::EMPrimaryElectronsCov> electronCovTable;
  o2::framework::Produces<aod::EMPrimaryElectronsPrefilterBit> electronPFTable;
  o2::framework::Produces<aod::EMPrimaryElectronsBDTSCT> sctTable;
};

struct electronCut : o2::framework::ConfigurableGroup {
  std::string prefix = "electronCut";
  o2::framework::Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  o2::framework::Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  o2::framework::Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  o2::framework::Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  o2::framework::Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  o2::framework::Configurable<float> minchi2tpc{"minchi2tpc", 0.0, "min. chi2/NclsTPC"};
  o2::framework::Configurable<float> minchi2its{"minchi2its", 0.0, "min. chi2/NclsITS"};
  o2::framework::Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  o2::framework::Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  o2::framework::Configurable<float> minpt{"minpt", 0.05, "min pt"};
  o2::framework::Configurable<float> maxeta{"maxeta", 0.9, "max eta"};
  o2::framework::Configurable<float> dca_xy_max{"dca_xy_max", 1.0, "max DCAxy in cm"};
  o2::framework::Configurable<float> dca_z_max{"dca_z_max", 1.0, "max DCAz in cm"};
  o2::framework::Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.5, "min. TPC n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.5, "max. TPC n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", 3.5, "max. TOF n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 2.5, "max. TPC n sigma for pion exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", -1e+10, "min. TPC n sigma for pion exclusion"}; // set to -2 for lowB, -1e+10 for nominalB
  o2::framework::Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 2.5, "max. TPC n sigma for kaon exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaKa{"minTPCNsigmaKa", -2.5, "min. TPC n sigma for kaon exclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 2.5, "max. TPC n sigma for proton exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaPr{"minTPCNsigmaPr", -2.5, "min. TPC n sigma for proton exclusion"};
  o2::framework::Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  o2::framework::Configurable<float> min_pin_for_pion_rejection{"min_pin_for_pion_rejection", 0.0, "pion rejection is applied above this pin"}; // this is used only in TOFreq
  o2::framework::Configurable<float> max_pin_for_pion_rejection{"max_pin_for_pion_rejection", 0.5, "pion rejection is applied below this pin"};
  o2::framework::Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  o2::framework::Configurable<float> maxMeanITSClusterSize{"maxMeanITSClusterSize", 16, "max <ITS cluster size> x cos(lambda)"};
  o2::framework::Configurable<bool> storeOnlyTrueElectronMC{"storeOnlyTrueElectronMC", false, "Flag to store only true electron in MC"};
  o2::framework::Configurable<int> minNelectron{"minNelectron", 0, "min number of electron candidates per collision"};
  o2::framework::Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to include ITSsa tracks only for MC. switch ON only if needed."};
  o2::framework::Configurable<bool> useTOFNSigmaDeltaBC{"useTOFNSigmaDeltaBC", false, "Flag to shift delta BC for TOF n sigma (only with TTCA)"};

  // configuration for PID ML
  o2::framework::Configurable<bool> usePIDML{"usePIDML", false, "Flag to use PID ML"};
  o2::framework::Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
  o2::framework::Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
  o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.15, 0.2, 0.25, 0.4, 0.8, 1.6, 2.0, 20}, "Bin limits for ML application"};
  o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.95, 0.95, 0.7, 0.7, 0.8, 0.8, 0.7, 0.7}, "ML cuts per bin"};
  o2::framework::Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"tpcInnerParam", "tpcNClsFound", "tpcChi2NCl", "tpcNSigmaEl", "tofNSigmaEl", "meanClusterSizeITSobCosTgl"}, "Names of ML model input features"};
  o2::framework::Configurable<std::string> nameBinningFeature{"nameBinningFeature", "tpcInnerParam", "Names of ML model binning feature"};
  o2::framework::Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  o2::framework::Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
};

struct electronPFCut : o2::framework::ConfigurableGroup {
  std::string prefix = "electronPFCut";
  o2::framework::Configurable<int> min_ncluster_tpc{"min_ncluster_tpc", 0, "min ncluster tpc"};
  o2::framework::Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  o2::framework::Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  o2::framework::Configurable<int> min_ncluster_its{"min_ncluster_its", 4, "min ncluster its"};
  o2::framework::Configurable<int> min_ncluster_itsib{"min_ncluster_itsib", 1, "min ncluster itsib"};
  o2::framework::Configurable<float> minchi2tpc{"minchi2tpc", 0.0, "min. chi2/NclsTPC"};
  o2::framework::Configurable<float> minchi2its{"minchi2its", 0.0, "min. chi2/NclsITS"};
  o2::framework::Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  o2::framework::Configurable<float> maxchi2its{"maxchi2its", 36.0, "max. chi2/NclsITS"};
  o2::framework::Configurable<float> minpt{"minpt", 0.05, "min pt"};
  o2::framework::Configurable<float> maxeta{"maxeta", 0.9, "max eta"};
  o2::framework::Configurable<float> dca_xy_max{"dca_xy_max", 1.0, "max DCAxy in cm"};
  o2::framework::Configurable<float> dca_z_max{"dca_z_max", 1.0, "max DCAz in cm"};
  o2::framework::Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 3.0, "max. TPC n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTOFNsigmaEl{"maxTOFNsigmaEl", 1e+10, "max. TOF n sigma for electron inclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaPi{"maxTPCNsigmaPi", 0.0, "max. TPC n sigma for pion exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaPi{"minTPCNsigmaPi", 0.0, "min. TPC n sigma for pion exclusion"}; // set to -2 for lowB, -1e+10 for nominalB
  o2::framework::Configurable<float> maxTPCNsigmaKa{"maxTPCNsigmaKa", 0.0, "max. TPC n sigma for kaon exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaKa{"minTPCNsigmaKa", 0.0, "min. TPC n sigma for kaon exclusion"};
  o2::framework::Configurable<float> maxTPCNsigmaPr{"maxTPCNsigmaPr", 0.0, "max. TPC n sigma for proton exclusion"};
  o2::framework::Configurable<float> minTPCNsigmaPr{"minTPCNsigmaPr", 0.0, "min. TPC n sigma for proton exclusion"};
  // o2::framework::Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  o2::framework::Configurable<float> min_pin_for_pion_rejection{"min_pin_for_pion_rejection", 0.0, "pion rejection is applied above this pin"}; // this is used only in TOFreq
  o2::framework::Configurable<float> max_pin_for_pion_rejection{"max_pin_for_pion_rejection", 0.5, "pion rejection is applied below this pin"};
  o2::framework::Configurable<float> max_frac_shared_clusters_tpc{"max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  o2::framework::Configurable<float> maxMeanITSClusterSize{"maxMeanITSClusterSize", 16, "max <ITS cluster size> x cos(lambda)"};

  // configuration for PID ML
  o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.15, 0.2, 0.25, 0.4, 0.8, 1.6, 2.0, 20}, "Bin limits for ML application"};
  o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.9, 0.9, 0.65, 0.65, 0.75, 0.75, 0.65, 0.65}, "ML cuts per bin"};

  // for pair
  o2::framework::Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  o2::framework::Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};
  o2::framework::Configurable<bool> doPF{"doPF", false, "flag to set pion prefilter"};
};

struct hadronCut : o2::framework::ConfigurableGroup {
  std::string prefix = "hadronCut";
  o2::framework::Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.01, "min pT for single track"};
  o2::framework::Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
  o2::framework::Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
  o2::framework::Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
  o2::framework::Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.f, "min. TPC Ncr/Nf ratio"};
  o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  o2::framework::Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
  o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
  o2::framework::Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 4, "min ncluster its"};
  o2::framework::Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
  o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
  o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
  o2::framework::Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
  o2::framework::Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};
  o2::framework::Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3, "min n sigma pi in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3, "max n sigma pi in TPC"};
  o2::framework::Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -3, "min n sigma pi in TOF"};
  o2::framework::Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +3, "max n sigma pi in TOF"};
  o2::framework::Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3, "min n sigma ka in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3, "max n sigma ka in TPC"};
  o2::framework::Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -3, "min n sigma ka in TOF"};
  o2::framework::Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +3, "max n sigma ka in TOF"};
  o2::framework::Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3, "min n sigma pr in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3, "max n sigma pr in TPC"};
  o2::framework::Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -3, "min n sigma pr in TOF"};
  o2::framework::Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +3, "max n sigma pr in TOF"};
  o2::framework::Configurable<bool> requirePiKaPr{"requirePiKaPr", true, "require hadron to be pion or kaon or proton"};
};

struct v0Cut : o2::framework::ConfigurableGroup {
  std::string prefix = "v0Cut";
  o2::framework::Configurable<float> cfg_min_mass_k0s{"cfg_min_mass_k0s", 0.48, "min mass for K0S"};
  o2::framework::Configurable<float> cfg_max_mass_k0s{"cfg_max_mass_k0s", 0.51, "max mass for K0S"};
  o2::framework::Configurable<float> cfg_min_mass_k0s_veto{"cfg_min_mass_k0s_veto", 0.48, "min mass for K0S veto for Lambda"};
  o2::framework::Configurable<float> cfg_max_mass_k0s_veto{"cfg_max_mass_k0s_veto", 0.51, "max mass for K0S veto for Lambda"};
  o2::framework::Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for Lambda"};
  o2::framework::Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for Lambda"};
  o2::framework::Configurable<float> cfg_min_mass_lambda_veto{"cfg_min_mass_lambda_veto", 1.11, "min mass for Lambda veto for K0S"};
  o2::framework::Configurable<float> cfg_max_mass_lambda_veto{"cfg_max_mass_lambda_veto", 1.12, "max mass for Lambda veto for K0S"};
  o2::framework::Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.95, "min cospa for v0"};
  o2::framework::Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 0.1, "max distance between 2 legs for v0"};
  o2::framework::Configurable<float> cfg_min_radius{"cfg_min_radius", 0.1, "min rxy for v"};
  o2::framework::Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.f, "min. TPC Ncr/Nf ratio"};
  o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
  o2::framework::Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
  o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
  o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
  o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
  o2::framework::Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
  o2::framework::Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
  o2::framework::Configurable<float> cfg_min_dcaxy{"cfg_min_dcaxy", 0.1, "min dca XY for v0 legs in cm"};

  o2::framework::Configurable<float> cfg_max_alpha_veto{"cfg_max_alpha_veto", 0.95, "max alpha for photon conversion rejection"};
  o2::framework::Configurable<float> cfg_max_qt_veto{"cfg_max_qt_veto", 0.01, "max qT for photon conversion rejection"};

  // for both v0 and cascade
  o2::framework::Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3, "min n sigma pi in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3, "max n sigma pi in TPC"};
  o2::framework::Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3, "min n sigma ka in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3, "max n sigma ka in TPC"};
  o2::framework::Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3, "min n sigma pr in TPC"};
  o2::framework::Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3, "max n sigma pr in TPC"};
  // o2::framework::Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -3, "min n sigma pi in TOF"};
  // o2::framework::Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +3, "max n sigma pi in TOF"};
  // o2::framework::Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -3, "min n sigma ka in TOF"};
  // o2::framework::Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +3, "max n sigma ka in TOF"};
  // o2::framework::Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -3, "min n sigma pr in TOF"};
  // o2::framework::Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +3, "max n sigma pr in TOF"};
  // o2::framework::Configurable<bool> applyTOFif{"applyTOFif", false, "apply TOFif for hadron identification"};
};

struct cascadeCut : o2::framework::ConfigurableGroup {
  std::string prefix = "cascadeCut";
  o2::framework::Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for lambda in cascade"};
  o2::framework::Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for lambda in cascade"};
  o2::framework::Configurable<float> cfg_min_mass_Xi{"cfg_min_mass_Xi", 1.314, "min mass for Xi"};
  o2::framework::Configurable<float> cfg_max_mass_Xi{"cfg_max_mass_Xi", 1.328, "max mass for Xi"};
  o2::framework::Configurable<float> cfg_min_mass_Xi_veto{"cfg_min_mass_Xi_veto", 1.31, "min mass for Xi veto"};
  o2::framework::Configurable<float> cfg_max_mass_Xi_veto{"cfg_max_mass_Xi_veto", 1.33, "max mass for Xi veto"};
  o2::framework::Configurable<float> cfg_min_mass_Omega{"cfg_min_mass_Omega", 1.668, "min mass for Omega"};
  o2::framework::Configurable<float> cfg_max_mass_Omega{"cfg_max_mass_Omega", 1.678, "max mass for Omega"};
  o2::framework::Configurable<float> cfg_min_mass_Omega_veto{"cfg_min_mass_Omega_veto", 1.665, "min mass for Omega veto"};
  o2::framework::Configurable<float> cfg_max_mass_Omega_veto{"cfg_max_mass_Omega_veto", 1.680, "max mass for Omega veto"};
  o2::framework::Configurable<float> cfg_min_cospa_v0{"cfg_min_cospa_v0", 0.95, "minimum V0 CosPA in cascade"};
  o2::framework::Configurable<float> cfg_max_dcadau_v0{"cfg_max_dcadau_v0", 0.1, "max distance between V0 Daughters in cascade"};
  o2::framework::Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.95, "minimum cascade CosPA"};
  o2::framework::Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.1, "max distance between bachelor and V0"};
  o2::framework::Configurable<float> cfg_min_rxy_v0{"cfg_min_rxy_v0", 0.1, "minimum V0 rxy in cascade"};
  o2::framework::Configurable<float> cfg_min_rxy{"cfg_min_rxy", 0.1, "minimum V0 rxy in cascade"};
  o2::framework::Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};
  o2::framework::Configurable<float> cfg_min_dcaxy_bachelor{"cfg_min_dcaxy_bachelor", 0.05, "min dca XY for bachelor in cm"};
  o2::framework::Configurable<float> cfg_min_dcaxy_v0{"cfg_min_dcaxy_v0", 0.0, "min dca XY for V0 in cm"};
};

struct cfgDFeT : o2::framework::ConfigurableGroup {
  std::string prefix = "cfgDFeT";
  o2::framework::Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  o2::framework::Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  o2::framework::Configurable<float> maxDCA2legs{"maxDCA2legs", 1.0, "max distance between 2 legs in cm"};
  // configuration for PID ML
  o2::framework::Configurable<bool> useML{"useML", false, "Flag to use PID ML"};
  o2::framework::Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
  o2::framework::Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
  o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4, 20}, "Bin limits for ML application"};
  // o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "ML cuts per bin"};
  o2::framework::Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"ptH", "impPar3DHinSigma", "tpcNSigmaKa", "signedMassLH", "cpa", "cpaXY", "dcaLH", "impPar3DinSigma", "decayLength3DinSigma", "decayLengthXYinSigma"}, "Names of ML model input features"};
  o2::framework::Configurable<std::string> nameBinningFeature{"nameBinningFeature", "ptL", "Names of ML model binning feature"};
  o2::framework::Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  o2::framework::Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
};

struct cfgDFeV0 : o2::framework::ConfigurableGroup {
  std::string prefix = "cfgDFeV0";
  o2::framework::Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  o2::framework::Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  o2::framework::Configurable<float> maxDCA2legs{"maxDCA2legs", 1.0, "max distance between 2 legs in cm"};
  // configuration for PID ML
  o2::framework::Configurable<bool> useML{"useML", false, "Flag to use PID ML"};
  o2::framework::Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
  o2::framework::Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
  o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4, 20}, "Bin limits for ML application"};
  // o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "ML cuts per bin"};
  o2::framework::Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"ptH", "impPar3DHinSigma", "massLH", "cpa", "cpaXY", "dcaLH", "impPar3DinSigma", "decayLength3DinSigma", "decayLengthXYinSigma"}, "Names of ML model input features"};
  o2::framework::Configurable<std::string> nameBinningFeature{"nameBinningFeature", "ptL", "Names of ML model binning feature"};
  o2::framework::Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  o2::framework::Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
};

struct cfgDFeC : o2::framework::ConfigurableGroup {
  std::string prefix = "cfgDFeC";
  o2::framework::Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  o2::framework::Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  o2::framework::Configurable<float> maxDCA2legs{"maxDCA2legs", 1.0, "max distance between 2 legs in cm"};
  // configuration for PID ML
  o2::framework::Configurable<bool> useML{"useML", false, "Flag to use PID ML"};
  o2::framework::Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"filename"}, "ONNX file names for each bin (if not from CCDB full path)"};
  o2::framework::Configurable<std::vector<std::string>> onnxPathsCCDB{"onnxPathsCCDB", std::vector<std::string>{"path"}, "Paths of models on CCDB"};
  o2::framework::Configurable<std::vector<double>> binsMl{"binsMl", std::vector<double>{0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4, 20}, "Bin limits for ML application"};
  // o2::framework::Configurable<std::vector<double>> cutsMl{"cutsMl", std::vector<double>{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}, "ML cuts per bin"};
  o2::framework::Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"ptH", "impPar3DHinSigma", "massLH", "cpa", "cpaXY", "dcaLH", "impPar3DinSigma", "decayLength3DinSigma", "decayLengthXYinSigma"}, "Names of ML model input features"};
  o2::framework::Configurable<std::string> nameBinningFeature{"nameBinningFeature", "ptL", "Names of ML model binning feature"};
  o2::framework::Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};
  o2::framework::Configurable<bool> enableOptimizations{"enableOptimizations", false, "Enables the ONNX extended model-optimization: sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED)"};
};

class ElectronModule
{
 public:
  ElectronModule()
  {
    // constructor
  }
  ~ElectronModule()
  {
    // destructor
  }

  struct looseElectron {
    float pt{1e+10};
    float eta{1e+10};
    float phi{1e+10};
  };

  template <typename TElectronCut, typename TElectronPFCut, typename THadronCut, typename TV0Cut, typename TCascadeCut, typename TDFConfigET, typename TDFConfigEV0, typename TDFConfigEC, typename TInitContext, typename TCCDB, typename TTOFResponse>
  void init(TElectronCut const& eCut, TElectronPFCut const& ePFCut, THadronCut const& hCut, TV0Cut const& v0Cut, TCascadeCut const& cascadeCut, TDFConfigET const& cfgET, TDFConfigEV0 const& cfgEV0, TDFConfigEC const& cfgEC, TInitContext& initContext, TCCDB& ccdb, TTOFResponse const& tofResponse, std::string const& ccdburl)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdbApi.init(ccdburl);

    LOGF(info, "intializing configurations");
    fElectronCut = eCut;
    fElectronPFCut = ePFCut;
    fHadronCut = hCut;
    fV0Cut = v0Cut;
    fCascadeCut = cascadeCut;
    fConfigDFeT = cfgET;
    fConfigDFeV0 = cfgEV0;
    fConfigDFeC = cfgEC;

    LOGF(info, "intializing TOFResponse");
    mTOFResponse = tofResponse;
    // mTOFResponse->initSetup(&ccdb->instance(), initContext);
    mTOFResponse->initSetup(ccdb, initContext);

    dfeT.setPropagateToPCA(true);
    dfeT.setMaxR(200.f);
    dfeT.setMinParamChange(1e-3);
    dfeT.setMinRelChi2Change(0.9);
    dfeT.setMaxDZIni(1e9);
    dfeT.setMaxChi2(1e9);
    dfeT.setUseAbsDCA(fConfigDFeT.useAbsDCA);
    dfeT.setWeightedFinalPCA(fConfigDFeT.useWeightedFinalPCA);
    dfeT.setMatCorrType(matCorr);

    dfeV0.setPropagateToPCA(true);
    dfeV0.setMaxR(200.f);
    dfeV0.setMinParamChange(1e-3);
    dfeV0.setMinRelChi2Change(0.9);
    dfeV0.setMaxDZIni(1e9);
    dfeV0.setMaxChi2(1e9);
    dfeV0.setUseAbsDCA(fConfigDFeV0.useAbsDCA);
    dfeV0.setWeightedFinalPCA(fConfigDFeV0.useWeightedFinalPCA);
    dfeV0.setMatCorrType(matCorr);

    dfeC.setPropagateToPCA(true);
    dfeC.setMaxR(200.f);
    dfeC.setMinParamChange(1e-3);
    dfeC.setMinRelChi2Change(0.9);
    dfeC.setMaxDZIni(1e9);
    dfeC.setMaxChi2(1e9);
    dfeC.setUseAbsDCA(fConfigDFeC.useAbsDCA);
    dfeC.setWeightedFinalPCA(fConfigDFeC.useWeightedFinalPCA);
    dfeC.setMatCorrType(matCorr);
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Configuring for timestamp " << bc.timestamp() << " with magnetic field of " << d_bz << " kG";

    // initialize MLResponse
    if (fElectronCut.usePIDML) {
      static constexpr int nClassesMl = 2;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      const std::vector<std::string> labelsClasses = {"Background", "Signal"};
      const uint32_t nBinsMl = fElectronCut.binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = 0.0;
        cutsMlArr[i][1] = fElectronCut.cutsMl.value[i];
      }
      o2::framework::LabeledArray<double> cutsMltmp = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponsePID.configure(fElectronCut.binsMl.value, cutsMltmp, cutDirMl, nClassesMl);
      if (fElectronCut.loadModelsFromCCDB) {
        mlResponsePID.setModelPathsCCDB(fElectronCut.onnxFileNames, ccdbApi, fElectronCut.onnxPathsCCDB, bc.timestamp());
      } else {
        mlResponsePID.setModelPathsLocal(fElectronCut.onnxFileNames);
      }
      mlResponsePID.cacheInputFeaturesIndices(fElectronCut.namesInputFeatures);
      mlResponsePID.cacheBinningIndex(fElectronCut.nameBinningFeature);
      mlResponsePID.init(fElectronCut.enableOptimizations);
    } // end of ML PID

    if (fConfigDFeT.useML) {
      static constexpr int nClassesMl = 4;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      const std::vector<std::string> labelsClasses = {"prompt", "prompthc", "nonprompthc", "hb"};
      const uint32_t nBinsMl = fConfigDFeT.binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = 0.0;
        cutsMlArr[i][1] = 0.0;
        cutsMlArr[i][2] = 0.0;
        cutsMlArr[i][3] = 0.0;
      }
      o2::framework::LabeledArray<double> cutsMltmp = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponseSCTeT.configure(fConfigDFeT.binsMl.value, cutsMltmp, cutDirMl, nClassesMl);
      if (fConfigDFeT.loadModelsFromCCDB) {
        mlResponseSCTeT.setModelPathsCCDB(fConfigDFeT.onnxFileNames, ccdbApi, fConfigDFeT.onnxPathsCCDB, bc.timestamp());
      } else {
        mlResponseSCTeT.setModelPathsLocal(fConfigDFeT.onnxFileNames);
      }
      mlResponseSCTeT.cacheInputFeaturesIndices(fConfigDFeT.namesInputFeatures);
      mlResponseSCTeT.cacheBinningIndex(fConfigDFeT.nameBinningFeature);
      mlResponseSCTeT.init(fConfigDFeT.enableOptimizations);
    } // end of ML SCTeT

    if (fConfigDFeV0.useML) {
      static constexpr int nClassesMl = 4;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      const std::vector<std::string> labelsClasses = {"prompt", "prompthc", "nonprompthc", "hb"};
      const uint32_t nBinsMl = fConfigDFeV0.binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = 0.0;
        cutsMlArr[i][1] = 0.0;
        cutsMlArr[i][2] = 0.0;
        cutsMlArr[i][3] = 0.0;
      }
      o2::framework::LabeledArray<double> cutsMltmp = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponseSCTeV0.configure(fConfigDFeV0.binsMl.value, cutsMltmp, cutDirMl, nClassesMl);
      if (fConfigDFeV0.loadModelsFromCCDB) {
        mlResponseSCTeV0.setModelPathsCCDB(fConfigDFeV0.onnxFileNames, ccdbApi, fConfigDFeV0.onnxPathsCCDB, bc.timestamp());
      } else {
        mlResponseSCTeV0.setModelPathsLocal(fConfigDFeV0.onnxFileNames);
      }
      mlResponseSCTeV0.cacheInputFeaturesIndices(fConfigDFeV0.namesInputFeatures);
      mlResponseSCTeV0.cacheBinningIndex(fConfigDFeV0.nameBinningFeature);
      mlResponseSCTeV0.init(fConfigDFeV0.enableOptimizations);
    } // end of ML SCTeV0

    if (fConfigDFeC.useML) {
      static constexpr int nClassesMl = 4;
      const std::vector<int> cutDirMl = {o2::cuts_ml::CutNot, o2::cuts_ml::CutSmaller};
      const std::vector<std::string> labelsClasses = {"prompt", "prompthc", "nonprompthc", "hb"};
      const uint32_t nBinsMl = fConfigDFeC.binsMl.value.size() - 1;
      const std::vector<std::string> labelsBins(nBinsMl, "bin");
      double cutsMlArr[nBinsMl][nClassesMl];
      for (uint32_t i = 0; i < nBinsMl; i++) {
        cutsMlArr[i][0] = 0.0;
        cutsMlArr[i][1] = 0.0;
        cutsMlArr[i][2] = 0.0;
        cutsMlArr[i][3] = 0.0;
      }
      o2::framework::LabeledArray<double> cutsMltmp = {cutsMlArr[0], nBinsMl, nClassesMl, labelsBins, labelsClasses};

      mlResponseSCTeC.configure(fConfigDFeC.binsMl.value, cutsMltmp, cutDirMl, nClassesMl);
      if (fConfigDFeC.loadModelsFromCCDB) {
        mlResponseSCTeC.setModelPathsCCDB(fConfigDFeC.onnxFileNames, ccdbApi, fConfigDFeC.onnxPathsCCDB, bc.timestamp());
      } else {
        mlResponseSCTeC.setModelPathsLocal(fConfigDFeC.onnxFileNames);
      }
      mlResponseSCTeC.cacheInputFeaturesIndices(fConfigDFeC.namesInputFeatures);
      mlResponseSCTeC.cacheBinningIndex(fConfigDFeC.nameBinningFeature);
      mlResponseSCTeC.init(fConfigDFeC.enableOptimizations);
    } // end of ML SCTeC

    mRunNumber = bc.runNumber();
    mTOFResponse->processSetup(bc);
  }

  template <typename THistoregistry>
  void addHistograms(THistoregistry& registry)
  {
    registry.add("Track/hPt", "pT;p_{T} (GeV/c)", o2::framework::HistType::kTH1F, {{1000, 0.0f, 10}}, false);
    registry.add("Track/hQoverPt", "q/pT;q/p_{T} (GeV/c)^{-1}", o2::framework::HistType::kTH1F, {{4000, -20, 20}}, false);
    registry.add("Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", o2::framework::HistType::kTH2F, {{180, 0, 2 * M_PI}, {20, -1.0f, 1.0f}}, false);
    registry.add("Track/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", o2::framework::HistType::kTH2F, {{200, -1.0f, 1.0f}, {200, -1.0f, 1.0f}}, false);
    registry.add("Track/hDCAxyzSigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", o2::framework::HistType::kTH2F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}}, false);
    registry.add("Track/hDCAxyRes_Pt", "DCA_{xy} resolution vs. pT;p_{T} (GeV/c);DCA_{xy} resolution (#mum)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    registry.add("Track/hDCAzRes_Pt", "DCA_{z} resolution vs. pT;p_{T} (GeV/c);DCA_{z} resolution (#mum)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {500, 0., 500}}, false);
    registry.add("Track/hNclsTPC", "number of TPC clusters", o2::framework::HistType::kTH1F, {{161, -0.5, 160.5}}, false);
    registry.add("Track/hNcrTPC", "number of TPC crossed rows", o2::framework::HistType::kTH1F, {{161, -0.5, 160.5}}, false);
    registry.add("Track/hChi2TPC", "chi2/number of TPC clusters", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("Track/hChi2TOF", "chi2 of TOF", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("Track/hTPCNcr2Nf", "TPC Ncr/Nfindable", o2::framework::HistType::kTH1F, {{200, 0, 2}}, false);
    registry.add("Track/hTPCNcls2Nf", "TPC Ncls/Nfindable", o2::framework::HistType::kTH1F, {{200, 0, 2}}, false);
    registry.add("Track/hTPCNclsShared", "TPC Ncls shared/Ncls;p_{T} (GeV/c);N_{cls}^{shared}/N_{cls} in TPC", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    registry.add("Track/hNclsITS", "number of ITS clusters", o2::framework::HistType::kTH1F, {{8, -0.5, 7.5}}, false);
    registry.add("Track/hChi2ITS", "chi2/number of ITS clusters", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("Track/hITSClusterMap", "ITS cluster map", o2::framework::HistType::kTH1F, {{128, -0.5, 127.5}}, false);
    registry.add("Track/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    registry.add("Track/hTPCdEdxMC", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    registry.add("Track/hTPCNsigmaEl", "TPC n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TPC}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTPCNsigmaPi", "TPC n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTPCNsigmaKa", "TPC n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TPC}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTPCNsigmaPr", "TPC n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TPC}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTOFbeta", "TOF beta;p_{pv} (GeV/c);#beta", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);
    registry.add("Track/hTOFNsigmaEl", "TOF n sigma el;p_{in} (GeV/c);n #sigma_{e}^{TOF}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTOFNsigmaPi", "TOF n sigma pi;p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTOFNsigmaKa", "TOF n sigma ka;p_{in} (GeV/c);n #sigma_{K}^{TOF}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hTOFNsigmaPr", "TOF n sigma pr;p_{in} (GeV/c);n #sigma_{p}^{TOF}", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, -5, +5}}, false);
    registry.add("Track/hMeanClusterSizeITS", "mean cluster size ITS;p_{pv} (GeV/c);<ITS cluster size> #times cos(#lambda)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
    registry.add("Track/hMeanClusterSizeITSib", "mean cluster size ITSib;p_{pv} (GeV/c);<ITSib cluster size> #times cos(#lambda)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
    registry.add("Track/hMeanClusterSizeITSob", "mean cluster size ITSob;p_{pv} (GeV/c);<ITSob cluster size> #times cos(#lambda)", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {150, 0, 15}}, false);
    registry.add("Track/hProbElBDT", "probability to be e from BDT;p_{in} (GeV/c);BDT score;", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {100, 0, 1}}, false);
    registry.add("Track/hNe", "electron counts;N_{e} per collision", o2::framework::HistType::kTH1F, {{51, -0.5, 50.5}}, false);

    registry.add("Prefilter/Track/hPt", "pT;p_{T} (GeV/c)", o2::framework::HistType::kTH1F, {{1000, 0.0f, 10}}, false);
    registry.add("Prefilter/Track/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", o2::framework::HistType::kTH2F, {{180, 0, 2 * M_PI}, {40, -2, 2}}, false);
    registry.add("Prefilter/Track/hTPCdEdx", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    registry.add("Prefilter/Track/hTOFbeta", "TOF #beta vs. p;p_{pv} (GeV/c);TOF #beta", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);
    registry.add("Prefilter/Pair/hMvsPhiV", "mass vs. phiv;#varphi_{V} (rad.);m_{ee} (GeV/c^{2})", o2::framework::HistType::kTH2F, {{90, 0.f, M_PI}, {100, 0, 0.1}});

    // for charged hadron
    registry.add("SCT/Track/hs", "hs;p_{T} (GeV/c);#eta;#varphi (rad.)", o2::framework::HistType::kTHnSparseF, {{100, 0, 10}, {80, -2, 2}, {36, 0, 2 * M_PI}}, false);
    registry.add("SCT/Track/hDCA", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", o2::framework::HistType::kTH2F, {{200, -1, 1}, {200, -1, 1}}, false);
    registry.add("SCT/Track/hTPCdEdx", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    registry.add("SCT/Track/hTOFbeta", "TOF #beta vs. p;p_{pv} (GeV/c);TOF #beta", o2::framework::HistType::kTH2F, {{1000, 0, 10}, {600, 0, 1.2}}, false);

    // for V0
    registry.add("SCT/V0/hPt", "pT of V0;p_{T} (GeV/c)", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("SCT/V0/hYPhi_K0S", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{K0S}", o2::framework::HistType::kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    registry.add("SCT/V0/hYPhi_Lambda", "rapidity vs. #varphi of V0;#varphi (rad.);rapidity_{#Lambda}", o2::framework::HistType::kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    registry.add("SCT/V0/hAP", "AP plot;#alpha;q_{T} (GeV/c)", o2::framework::HistType::kTH2F, {{200, -1, 1}, {250, 0, 0.25}}, false);
    registry.add("SCT/V0/hLxy", "decay length from PV;L_{xy} (cm)", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("SCT/V0/hCosPA", "cosPA;cosine of pointing angle", o2::framework::HistType::kTH1F, {{100, 0.9, 1}}, false);
    registry.add("SCT/V0/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", o2::framework::HistType::kTH1F, {{100, 0, 1}}, false);
    registry.add("SCT/V0/hMassK0S", "K0S mass;m_{#pi#pi} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 0.45, 0.55}}, false);
    registry.add("SCT/V0/hMassLambda", "Lambda mass;m_{p#pi^{#minus}} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.08, 1.18}}, false);
    registry.add("SCT/V0/hMassAntiLambda", "Anti-Lambda mass;m_{#bar{p}#pi^{+}} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.08, 1.18}}, false);
    registry.add("SCT/V0/hMassGamma_misid", "#gamma mass;m_{ee} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 0, 0.1}}, false);
    registry.add("SCT/V0/hMassK0S_misid", "K0S mass;m_{#pi#pi} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 0.45, 0.55}}, false);
    registry.add("SCT/V0/hMassLambda_misid", "Lambda mass;m_{p#pi^{#minus}} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.08, 1.18}}, false);
    registry.add("SCT/V0/hMassAntiLambda_misid", "Anti-Lambda mass;m_{#bar{p}#pi^{+}} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.08, 1.18}}, false);

    // for cascadeSCT/
    registry.add("SCT/Cascade/hPt", "pT of cascade;p_{T} (GeV/c)", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.add("SCT/Cascade/hYPhi_Xi", "rapidity vs. #varphi of cascade;#varphi (rad.);rapidity_{#Xi}", o2::framework::HistType::kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    registry.add("SCT/Cascade/hYPhi_Omega", "rapidity vs. #varphi of cascade;#varphi (rad.);rapidity_{#Omega}", o2::framework::HistType::kTH2F, {{90, 0, 2 * M_PI}, {80, -2, +2}}, false);
    registry.add("SCT/Cascade/hCosPA", "cosPA;cosine of pointing angle", o2::framework::HistType::kTH1F, {{100, 0.9, 1}}, false);
    registry.add("SCT/Cascade/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm)", o2::framework::HistType::kTH1F, {{100, 0, 1}}, false);
    registry.add("SCT/Cascade/hV0CosPA", "cosPA of V0 in cascade;cosine of pointing angle", o2::framework::HistType::kTH1F, {{100, 0.9, 1}}, false);
    registry.add("SCT/Cascade/hV0DCA2Legs", "distance between 2 legs at PCA of V0 in cascade;distance between 2 legs (cm)", o2::framework::HistType::kTH1F, {{100, 0, 1}}, false);
    registry.add("SCT/Cascade/hMassLambda", "Lambda mass;m_{p#pi^{-}} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.08, 1.18}}, false);
    registry.add("SCT/Cascade/hMassXi", "#Xi mass;m_{#Lambda#pi} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.27, 1.37}}, false);
    registry.add("SCT/Cascade/hMassOmega", "#Omega mass;m_{#LambdaK} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.62, 1.72}}, false);
    registry.add("SCT/Cascade/hMassXi_misid", "#Xi mass;m_{#Lambda#pi} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.27, 1.37}}, false);
    registry.add("SCT/Cascade/hMassOmega_misid", "#Omega mass;m_{#LambdaK} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 1.62, 1.72}}, false);

    registry.add("SCT/eT/hImpactParameter", "eH impact parameter;IP in XY (cm);IP in Z (cm);", o2::framework::HistType::kTH2F, {{200, -1, 1}, {200, -1, 1}}, false);
    registry.add("SCT/eT/hDecayLength", "eH decay length;L_{xy} (cm);L_{z} (cm);", o2::framework::HistType::kTH2F, {{100, 0, 1}, {200, -1, 1}}, false);
    registry.add("SCT/eT/hCosPA", "eH cosPA;cosine of pointing angle", o2::framework::HistType::kTH1F, {{200, -1, 1}}, false);
    registry.add("SCT/eT/hDCA2legs", "dca between e and hadron;dca between 2legs (cm)", o2::framework::HistType::kTH1F, {{100, 0, 1}}, false);
    registry.add("SCT/eT/hMass", "invariant mass of e-h;m_{eh} (GeV/c^{2})", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    registry.addClone("SCT/eT/", "SCT/eV0/");
    registry.addClone("SCT/eT/", "SCT/eC/");
  }

  template <bool withTTCA, typename TCollisions, typename TBCs, typename TTracks, typename TTrackAssoc, typename TSliceCache, typename TPresliceTrack, typename TPresliceTrackAssoc>
  void calculateTOFNSigmaWithReassociation(TCollisions const& collisions, TBCs const&, TTracks const& tracks, TTrackAssoc const& trackIndices, TSliceCache const&, TPresliceTrack const& perColTrack, TPresliceTrackAssoc const& trackIndicesPerCollision)
  {
    if (fElectronCut.useTOFNSigmaDeltaBC) {
      if constexpr (withTTCA) {
        std::unordered_map<int, double> mapCollisionTime;
        std::unordered_map<int, double> mapCollisionTimeError;
        for (const auto& track : tracks) {
          if (mapCollisionTime.find(track.collisionId()) == mapCollisionTime.end()) {
            mapCollisionTime[track.collisionId()] = track.tofEvTime();
            mapCollisionTimeError[track.collisionId()] = track.tofEvTimeErr();
          }
        }

        for (const auto& collision : collisions) {
          if (mapCollisionTime.find(collision.globalIndex()) == mapCollisionTime.end()) {
            continue;
          }
          auto bcCollision = collision.template bc_as<TBCs>();
          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
          for (const auto& trackId : trackIdsThisCollision) {
            auto track = trackId.template track_as<TTracks>();
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }

            if (track.hasTOF() && track.has_collision()) { // TTCA may use orphan tracks.
              auto bcTrack = track.template collision_as<TCollisions>().template bc_as<TBCs>();
              fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = mTOFResponse->nSigma<o2::track::PID::Electron>(track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()), track.tofExpMom(), track.length(), track.p(), track.eta(), mapCollisionTime[collision.globalIndex()], mapCollisionTimeError[collision.globalIndex()]);
              fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = mTOFResponse->nSigma<o2::track::PID::Pion>(track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()), track.tofExpMom(), track.length(), track.p(), track.eta(), mapCollisionTime[collision.globalIndex()], mapCollisionTimeError[collision.globalIndex()]);
              fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = mTOFResponse->nSigma<o2::track::PID::Kaon>(track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()), track.tofExpMom(), track.length(), track.p(), track.eta(), mapCollisionTime[collision.globalIndex()], mapCollisionTimeError[collision.globalIndex()]);
              fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = mTOFResponse->nSigma<o2::track::PID::Proton>(track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()), track.tofExpMom(), track.length(), track.p(), track.eta(), mapCollisionTime[collision.globalIndex()], mapCollisionTimeError[collision.globalIndex()]);
              fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.length() / (track.tofSignalInAnotherBC(bcTrack.globalBC(), bcCollision.globalBC()) - mapCollisionTime[collision.globalIndex()]) / (TMath::C() * 1e+2 * 1e-12);
              ;
            } else {
              fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
              fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPi();
              fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaKa();
              fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPr();
              fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
            }
          } // end of track loop
        } // end of collision loop
        mapCollisionTime.clear();
        mapCollisionTimeError.clear();
      } else { // without TTCA
        for (const auto& collision : collisions) {
          auto tracks_per_coll = tracks.sliceBy(perColTrack, collision.globalIndex());
          for (const auto& track : tracks_per_coll) {
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPi();
            fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaKa();
            fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPr();
            fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
          }
        } // end of track loop
      } // end of collision loop
    } else {
      if constexpr (withTTCA) {
        for (const auto& collision : collisions) {
          auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
          for (const auto& trackId : trackIdsThisCollision) {
            auto track = trackId.template track_as<TTracks>();
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPi();
            fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaKa();
            fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPr();
            fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
          } // end of track loop
        } // end of collision loop
      } else {
        for (const auto& collision : collisions) {
          auto tracks_per_coll = tracks.sliceBy(perColTrack, collision.globalIndex());
          for (const auto& track : tracks_per_coll) {
            if (!track.hasITS() || !track.hasTPC()) { // apply only minimal cut
              continue;
            }
            fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaEl();
            fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPi();
            fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaKa();
            fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.tofNSigmaPr();
            fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())] = track.beta();
          }
        } // end of track loop
      } // end of collision loop
    }
  }

  template <typename TCollision, typename TTrack>
  void fillMapMLPID(TCollision const& collision, TTrack const& track)
  {
    if (fElectronCut.usePIDML) {
      o2::dataformats::DCA mDcaInfoCov;
      mDcaInfoCov.set(999, 999, 999, 999, 999);
      auto trackParCov = getTrackParCov(track);
      trackParCov.setPID(o2::track::PID::Electron);
      bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
      if (!isPropOK) {
        return;
      }

      if (!isElectron_TOFif(track, collision, fElectronCut)) { // minimal n sigma cut is taken from the main electron cut.
        return;
      }

      o2::analysis::pwgem::dilepton::mlpid::candidate candidate;
      candidate.tpcInnerParam = track.tpcInnerParam();
      candidate.tpcNClsFound = track.tpcNClsFound();
      candidate.tpcChi2NCl = track.tpcChi2NCl();
      candidate.tpcNSigmaEl = track.tpcNSigmaEl();
      candidate.tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
      int total_cluster_size_ob = 0, nl_ob = 0;
      for (unsigned int layer = 3; layer < 7; layer++) {
        int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
        if (cluster_size_per_layer > 0) {
          nl_ob++;
        }
        total_cluster_size_ob += cluster_size_per_layer;
      }
      candidate.meanClusterSizeITSobCosTgl = static_cast<float>(total_cluster_size_ob) / static_cast<float>(nl_ob) * std::cos(std::atan(trackParCov.getTgl()));

      std::vector<float> inputFeatures = mlResponsePID.getInputFeatures(candidate);
      float binningFeature = mlResponsePID.getBinningFeature(candidate);

      int pbin = lower_bound(fElectronCut.binsMl.value.begin(), fElectronCut.binsMl.value.end(), binningFeature) - fElectronCut.binsMl.value.begin() - 1;
      if (pbin < 0) {
        pbin = 0;
      } else if (static_cast<int>(fElectronCut.binsMl.value.size()) - 2 < pbin) {
        pbin = static_cast<int>(fElectronCut.binsMl.value.size()) - 2;
      }

      float probaEl = mlResponsePID.getModelOutput(inputFeatures, pbin)[1]; // 0: hadron, 1:electron
      fMapProbaEl[std::make_pair(collision.globalIndex(), track.globalIndex())] = probaEl;
    } else {
      fMapProbaEl[std::make_pair(collision.globalIndex(), track.globalIndex())] = 1.0;
    }
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrack(TCollision const& collision, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
      if (fElectronCut.storeOnlyTrueElectronMC) {
        const auto& mcParticle = track.template mcParticle_as<aod::McParticles>();
        if (std::abs(mcParticle.pdgCode()) != 11) {
          return false;
        }
      }
    }

    float tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    if (fElectronCut.requireTOF && !(track.hasTOF() && std::fabs(tofNSigmaEl) < fElectronCut.maxTOFNsigmaEl)) {
      return false;
    }

    if (!track.hasITS()) {
      return false;
    }

    if (track.itsChi2NCl() < fElectronCut.minchi2its || fElectronCut.maxchi2its < track.itsChi2NCl()) {
      return false;
    }
    if (track.itsNCls() < fElectronCut.min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < fElectronCut.min_ncluster_itsib) {
      return false;
    }

    if (!fElectronCut.includeITSsa && (!track.hasITS() || !track.hasTPC())) {
      return false;
    }

    if (track.hasTPC()) {
      if (track.tpcChi2NCl() < fElectronCut.minchi2tpc || fElectronCut.maxchi2tpc < track.tpcChi2NCl()) {
        return false;
      }

      if (track.tpcNClsFound() < fElectronCut.min_ncluster_tpc) {
        return false;
      }

      if (track.tpcNClsCrossedRows() < fElectronCut.mincrossedrows) {
        return false;
      }

      if (track.tpcCrossedRowsOverFindableCls() < fElectronCut.min_tpc_cr_findable_ratio) {
        return false;
      }

      if (track.tpcFractionSharedCls() > fElectronCut.max_frac_shared_clusters_tpc) {
        return false;
      }
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > fElectronCut.dca_xy_max || std::fabs(dcaZ) > fElectronCut.dca_z_max) {
      return false;
    }

    if (trackParCov.getPt() < fElectronCut.minpt || std::fabs(trackParCov.getEta()) > fElectronCut.maxeta) {
      return false;
    }

    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 0; layer < 7; layer++) {
      int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }

    if (fElectronCut.maxMeanITSClusterSize < static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(trackParCov.getTgl()))) {
      return false;
    }

    return true;
  }

  template <bool isMC, typename TCollision, typename TTrack>
  bool checkTrackPF(TCollision const&, TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsChi2NCl() < fElectronPFCut.minchi2its || fElectronPFCut.maxchi2its < track.itsChi2NCl()) {
      return false;
    }
    if (track.itsNCls() < fElectronPFCut.min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < fElectronPFCut.min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < fElectronPFCut.minchi2tpc || fElectronPFCut.maxchi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < fElectronPFCut.min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < fElectronPFCut.mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < fElectronPFCut.min_tpc_cr_findable_ratio) {
      return false;
    }

    if (track.tpcFractionSharedCls() > fElectronPFCut.max_frac_shared_clusters_tpc) {
      return false;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov = getTrackParCov(track);
    trackParCov.setPID(o2::track::PID::Electron);
    bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
    if (!isPropOK) {
      return false;
    }
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(dcaXY) > fElectronPFCut.dca_xy_max || std::fabs(dcaZ) > fElectronPFCut.dca_z_max) {
      return false;
    }

    if (trackParCov.getPt() < fElectronPFCut.minpt || std::fabs(trackParCov.getEta()) > fElectronPFCut.maxeta) {
      return false;
    }

    int total_cluster_size = 0, nl = 0;
    for (unsigned int layer = 0; layer < 7; layer++) {
      int cluster_size_per_layer = track.itsClsSizeInLayer(layer);
      if (cluster_size_per_layer > 0) {
        nl++;
      }
      total_cluster_size += cluster_size_per_layer;
    }

    if (fElectronPFCut.maxMeanITSClusterSize < static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(trackParCov.getTgl()))) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TTrack, typename TConfig>
  bool isElectron(TCollision const& collision, TTrack const& track, TConfig const& eCut)
  {
    if (fElectronCut.usePIDML) { // keep fElectronCut here intentinoally
      if (isElectron_TOFif(track, collision, eCut)) {
        o2::analysis::pwgem::dilepton::mlpid::candidate candidate;
        candidate.tpcInnerParam = track.tpcInnerParam();
        float binningFeature = mlResponsePID.getBinningFeature(candidate);
        int pbin = lower_bound(eCut.binsMl.value.begin(), eCut.binsMl.value.end(), binningFeature) - eCut.binsMl.value.begin() - 1;
        if (pbin < 0) {
          pbin = 0;
        } else if (static_cast<int>(eCut.binsMl.value.size()) - 2 < pbin) {
          pbin = static_cast<int>(eCut.binsMl.value.size()) - 2;
        }
        return fMapProbaEl[std::make_pair(collision.globalIndex(), track.globalIndex())] > eCut.cutsMl.value[pbin];
      } else {
        return false;
      }
    } else {
      return isElectron_TPChadrej(track, collision, eCut) || isElectron_TOFreq(track, collision, eCut);
    }
  }

  template <typename TTrack, typename TCollision, typename TConfig>
  bool isElectron_TOFif(TTrack const& track, TCollision const& collision, TConfig const& eCut)
  {
    float tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    bool is_EL_TPC = eCut.minTPCNsigmaEl.value < track.tpcNSigmaEl() && track.tpcNSigmaEl() < eCut.maxTPCNsigmaEl.value;
    bool is_EL_TOF = track.hasTOF() ? (std::fabs(tofNSigmaEl) < eCut.maxTOFNsigmaEl.value) : true; // TOFif
    return is_EL_TPC && is_EL_TOF;
  }

  template <typename TTrack, typename TCollision, typename TConfig>
  bool isElectron_TPChadrej(TTrack const& track, TCollision const& collision, TConfig const& eCut)
  {
    float tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    if (track.tpcNSigmaEl() < eCut.minTPCNsigmaEl.value || eCut.maxTPCNsigmaEl.value < track.tpcNSigmaEl()) {
      return false;
    }
    if (eCut.minTPCNsigmaPi.value < track.tpcNSigmaPi() && track.tpcNSigmaPi() < eCut.maxTPCNsigmaPi.value && track.tpcInnerParam() < eCut.max_pin_for_pion_rejection.value) {
      return false;
    }
    if (eCut.minTPCNsigmaKa.value < track.tpcNSigmaKa() && track.tpcNSigmaKa() < eCut.maxTPCNsigmaKa.value) {
      return false;
    }
    if (eCut.minTPCNsigmaPr.value < track.tpcNSigmaPr() && track.tpcNSigmaPr() < eCut.maxTPCNsigmaPr.value) {
      return false;
    }
    if (track.hasTOF() && (eCut.maxTOFNsigmaEl.value < std::fabs(tofNSigmaEl))) {
      return false;
    }
    return true;
  }

  template <typename TTrack, typename TCollision, typename TConfig>
  bool isElectron_TOFreq(TTrack const& track, TCollision const& collision, TConfig const& eCut)
  {
    float tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    if (eCut.minTPCNsigmaPi.value < track.tpcNSigmaPi() && track.tpcNSigmaPi() < eCut.maxTPCNsigmaPi.value && (eCut.min_pin_for_pion_rejection.value < track.tpcInnerParam() && track.tpcInnerParam() < eCut.max_pin_for_pion_rejection.value)) {
      return false;
    }
    return eCut.minTPCNsigmaEl.value < track.tpcNSigmaEl() && track.tpcNSigmaEl() < eCut.maxTPCNsigmaEl.value && std::fabs(tofNSigmaEl) < eCut.maxTOFNsigmaEl.value;
  }

  template <typename TCollision, typename TTrack, typename TTrackParCov>
  bool isSelectedHadron(TCollision const&, TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < fHadronCut.cfg_min_pt_track || fHadronCut.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (trackParCov.getEta() < fHadronCut.cfg_min_eta_track || fHadronCut.cfg_max_eta_track < trackParCov.getEta()) {
      return false;
    }

    if (std::fabs(dcaXY) > fHadronCut.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > fHadronCut.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() < 0.f || fHadronCut.cfg_max_chi2its < track.itsChi2NCl()) {
      return false;
    }

    if (track.itsNCls() < fHadronCut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < fHadronCut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() < 0.f || fHadronCut.cfg_max_chi2tpc < track.tpcChi2NCl()) {
      return false;
    }

    if (track.tpcNClsFound() < fHadronCut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < fHadronCut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < fHadronCut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > fHadronCut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    if (fHadronCut.requirePiKaPr && !isPiKaPr(track)) {
      return false;
    }

    return true;
  }

  template <bool isMC, typename TCollision, typename TTrack, typename TTrackParCov, typename TProducts, typename THistoregistry>
  void fillElectronTable(TCollision const& collision, TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ, TProducts& products, THistoregistry& registry)
  {
    float pt = trackParCov.getPt();
    float eta = trackParCov.getEta();
    float phi = RecoDecay::constrainAngle(trackParCov.getPhi(), 0, 1U);

    float tofNSigmaEl = fMapTOFNsigmaElReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    float tofNSigmaPi = fMapTOFNsigmaPiReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    float tofNSigmaKa = fMapTOFNsigmaKaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    float tofNSigmaPr = fMapTOFNsigmaPrReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    float beta = fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())];
    float probaEl = fMapProbaEl[std::make_pair(collision.globalIndex(), track.globalIndex())];

    bool isAssociatedToMPC = collision.globalIndex() == track.collisionId();
    float mcTunedTPCSignal = 0.f;
    if constexpr (isMC) {
      mcTunedTPCSignal = track.mcTunedTPCSignal();
    }

    products.electronTable(collision.globalIndex(), track.globalIndex(), track.sign(),
                           pt, eta, phi,
                           dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(),
                           track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusPID(), track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(),
                           track.tpcChi2NCl(), track.tpcInnerParam(),
                           track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                           beta, tofNSigmaEl,
                           track.itsClusterSizes(),
                           track.itsChi2NCl(), track.tofChi2(), track.detectorMap(),
                           isAssociatedToMPC, false, probaEl, track.flags(), mcTunedTPCSignal);

    products.electronCovTable(
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

    registry.fill(HIST("Track/hPt"), pt);
    registry.fill(HIST("Track/hQoverPt"), track.sign() / pt);
    registry.fill(HIST("Track/hEtaPhi"), phi, eta);
    registry.fill(HIST("Track/hDCAxyz"), dcaXY, dcaZ);
    registry.fill(HIST("Track/hDCAxyzSigma"), dcaXY / std::sqrt(trackParCov.getSigmaY2()), dcaZ / std::sqrt(trackParCov.getSigmaZ2()));
    registry.fill(HIST("Track/hDCAxyRes_Pt"), pt, std::sqrt(trackParCov.getSigmaY2()) * 1e+4); // convert cm to um
    registry.fill(HIST("Track/hDCAzRes_Pt"), pt, std::sqrt(trackParCov.getSigmaZ2()) * 1e+4);  // convert cm to um
    registry.fill(HIST("Track/hNclsITS"), track.itsNCls());
    registry.fill(HIST("Track/hNclsTPC"), track.tpcNClsFound());
    registry.fill(HIST("Track/hNcrTPC"), track.tpcNClsCrossedRows());
    registry.fill(HIST("Track/hTPCNcr2Nf"), track.tpcCrossedRowsOverFindableCls());
    registry.fill(HIST("Track/hTPCNcls2Nf"), track.tpcFoundOverFindableCls());
    registry.fill(HIST("Track/hTPCNclsShared"), track.pt(), track.tpcFractionSharedCls());
    registry.fill(HIST("Track/hChi2TPC"), track.tpcChi2NCl());
    registry.fill(HIST("Track/hChi2ITS"), track.itsChi2NCl());
    registry.fill(HIST("Track/hChi2TOF"), track.tofChi2());
    registry.fill(HIST("Track/hITSClusterMap"), track.itsClusterMap());
    registry.fill(HIST("Track/hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    registry.fill(HIST("Track/hTPCdEdxMC"), track.tpcInnerParam(), mcTunedTPCSignal);
    registry.fill(HIST("Track/hTPCNsigmaEl"), track.tpcInnerParam(), track.tpcNSigmaEl());
    registry.fill(HIST("Track/hTPCNsigmaPi"), track.tpcInnerParam(), track.tpcNSigmaPi());
    registry.fill(HIST("Track/hTPCNsigmaKa"), track.tpcInnerParam(), track.tpcNSigmaKa());
    registry.fill(HIST("Track/hTPCNsigmaPr"), track.tpcInnerParam(), track.tpcNSigmaPr());
    registry.fill(HIST("Track/hTOFbeta"), trackParCov.getP(), beta);
    registry.fill(HIST("Track/hTOFNsigmaEl"), track.tpcInnerParam(), tofNSigmaEl);
    registry.fill(HIST("Track/hTOFNsigmaPi"), track.tpcInnerParam(), tofNSigmaPi);
    registry.fill(HIST("Track/hTOFNsigmaKa"), track.tpcInnerParam(), tofNSigmaKa);
    registry.fill(HIST("Track/hTOFNsigmaPr"), track.tpcInnerParam(), tofNSigmaPr);
    registry.fill(HIST("Track/hMeanClusterSizeITS"), trackParCov.getP(), static_cast<float>(total_cluster_size) / static_cast<float>(nl) * std::cos(std::atan(trackParCov.getTgl())));
    registry.fill(HIST("Track/hMeanClusterSizeITSib"), trackParCov.getP(), static_cast<float>(total_cluster_size_ib) / static_cast<float>(nl_ib) * std::cos(std::atan(trackParCov.getTgl())));
    registry.fill(HIST("Track/hMeanClusterSizeITSob"), trackParCov.getP(), static_cast<float>(total_cluster_size_ob) / static_cast<float>(nl_ob) * std::cos(std::atan(trackParCov.getTgl())));
    registry.fill(HIST("Track/hProbElBDT"), track.tpcInnerParam(), probaEl);
  }

  template <typename TTrack>
  bool isPiKaPr(TTrack const& track)
  {
    bool is_pi_included_TPC = fHadronCut.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < fHadronCut.cfg_max_TPCNsigmaPi;
    bool is_ka_included_TPC = fHadronCut.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < fHadronCut.cfg_max_TPCNsigmaKa;
    bool is_pr_included_TPC = fHadronCut.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < fHadronCut.cfg_max_TPCNsigmaPr;
    return is_pi_included_TPC || is_ka_included_TPC || is_pr_included_TPC;
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

    if (track.itsChi2NCl() > fV0Cut.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < fV0Cut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < fV0Cut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > fV0Cut.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < fV0Cut.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < fV0Cut.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < fV0Cut.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > fV0Cut.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    return fV0Cut.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < fV0Cut.cfg_max_TPCNsigmaPi;
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    return fV0Cut.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < fV0Cut.cfg_max_TPCNsigmaKa;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    return fV0Cut.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < fV0Cut.cfg_max_TPCNsigmaPr;
  }

  template <typename TV0>
  bool isK0S(TV0 const& v0)
  {
    return (fV0Cut.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < fV0Cut.cfg_max_mass_k0s) && (v0.mLambda() < fV0Cut.cfg_min_mass_lambda_veto || fV0Cut.cfg_max_mass_lambda_veto < v0.mLambda()) && (v0.mAntiLambda() < fV0Cut.cfg_min_mass_lambda_veto || fV0Cut.cfg_max_mass_lambda_veto < v0.mAntiLambda());
  }

  template <typename TV0>
  bool isLambda(TV0 const& v0)
  {
    return (fV0Cut.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < fV0Cut.cfg_max_mass_lambda) && (v0.mK0Short() < fV0Cut.cfg_min_mass_k0s_veto || fV0Cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TV0>
  bool isAntiLambda(TV0 const& v0)
  {
    return (fV0Cut.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < fV0Cut.cfg_max_mass_lambda) && (v0.mK0Short() < fV0Cut.cfg_min_mass_k0s_veto || fV0Cut.cfg_max_mass_k0s_veto < v0.mK0Short());
  }

  template <typename TCascade>
  bool isXi(TCascade const& cascade)
  {
    return (fCascadeCut.cfg_min_mass_Xi < cascade.mXi() && cascade.mXi() < fCascadeCut.cfg_max_mass_Xi) && (cascade.mOmega() < fCascadeCut.cfg_min_mass_Omega_veto || fCascadeCut.cfg_max_mass_Omega_veto < cascade.mOmega());
  }

  template <typename TCascade>
  bool isOmega(TCascade const& cascade)
  {
    return (fCascadeCut.cfg_min_mass_Omega < cascade.mOmega() && cascade.mOmega() < fCascadeCut.cfg_max_mass_Omega) && (cascade.mXi() < fCascadeCut.cfg_min_mass_Xi_veto || fCascadeCut.cfg_max_mass_Xi_veto < cascade.mXi());
  }

  template <typename TTracks, typename TCollision, typename TV0, typename THistoregistry>
  void fillV0Histograms(TCollision const&, TV0 const& v0, THistoregistry& registry)
  {
    auto pos = v0.template posTrack_as<TTracks>();
    auto neg = v0.template negTrack_as<TTracks>();
    registry.fill(HIST("SCT/V0/hPt"), v0.pt());
    registry.fill(HIST("SCT/V0/hAP"), v0.alpha(), v0.qtarm());
    registry.fill(HIST("SCT/V0/hCosPA"), v0.v0cosPA());
    registry.fill(HIST("SCT/V0/hLxy"), v0.v0radius());
    registry.fill(HIST("SCT/V0/hDCA2Legs"), v0.dcaV0daughters());

    if (isPion(pos) && isPion(neg)) {
      if ((v0.mLambda() < fV0Cut.cfg_min_mass_lambda_veto || fV0Cut.cfg_max_mass_lambda_veto < v0.mLambda()) && (v0.mAntiLambda() < fV0Cut.cfg_min_mass_lambda_veto || fV0Cut.cfg_max_mass_lambda_veto < v0.mAntiLambda())) {
        registry.fill(HIST("SCT/V0/hMassK0S"), v0.mK0Short());
        registry.fill(HIST("SCT/V0/hYPhi_K0S"), v0.phi(), v0.yK0Short());
      }
      registry.fill(HIST("SCT/V0/hMassGamma_misid"), v0.mGamma());
      registry.fill(HIST("SCT/V0/hMassLambda_misid"), v0.mLambda());
      registry.fill(HIST("SCT/V0/hMassAntiLambda_misid"), v0.mAntiLambda());
    }

    if (isProton(pos) && isPion(neg)) {
      if (v0.mK0Short() < fV0Cut.cfg_min_mass_k0s_veto || fV0Cut.cfg_max_mass_k0s_veto < v0.mK0Short()) {
        registry.fill(HIST("SCT/V0/hMassLambda"), v0.mLambda());
        registry.fill(HIST("SCT/V0/hYPhi_Lambda"), v0.phi(), v0.yLambda());
      }
      registry.fill(HIST("SCT/V0/hMassGamma_misid"), v0.mGamma());
      registry.fill(HIST("SCT/V0/hMassK0S_misid"), v0.mK0Short());
    }

    if (isProton(neg) && isPion(pos)) {
      if (v0.mK0Short() < fV0Cut.cfg_min_mass_k0s_veto || fV0Cut.cfg_max_mass_k0s_veto < v0.mK0Short()) {
        registry.fill(HIST("SCT/V0/hMassAntiLambda"), v0.mAntiLambda());
        registry.fill(HIST("SCT/V0/hYPhi_Lambda"), v0.phi(), v0.yLambda());
      }
      registry.fill(HIST("SCT/V0/hMassGamma_misid"), v0.mGamma());
      registry.fill(HIST("SCT/V0/hMassK0S_misid"), v0.mK0Short());
    }
  }

  template <typename TTracks, typename TCollision, typename TCascade, typename THistoregistry>
  void fillCascadeHistograms(TCollision const& collision, TCascade const& cascade, THistoregistry& registry)
  {
    auto pos = cascade.template posTrack_as<TTracks>();
    auto neg = cascade.template negTrack_as<TTracks>();
    auto bachelor = cascade.template bachelor_as<TTracks>();

    registry.fill(HIST("SCT/Cascade/hPt"), cascade.pt());
    registry.fill(HIST("SCT/Cascade/hMassLambda"), cascade.mLambda());
    registry.fill(HIST("SCT/Cascade/hCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
    registry.fill(HIST("SCT/Cascade/hDCA2Legs"), cascade.dcacascdaughters());
    registry.fill(HIST("SCT/Cascade/hV0CosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
    registry.fill(HIST("SCT/Cascade/hV0DCA2Legs"), cascade.dcaV0daughters());

    if (cascade.sign() < 0) { // Xi- or Omega-
      if (isPion(bachelor) && isProton(pos) && isPion(neg)) {
        if (cascade.mOmega() < fCascadeCut.cfg_min_mass_Omega_veto || fCascadeCut.cfg_max_mass_Omega_veto < cascade.mOmega()) {
          registry.fill(HIST("SCT/Cascade/hMassXi"), cascade.mXi());
          registry.fill(HIST("SCT/Cascade/hYPhi_Xi"), cascade.phi(), cascade.yXi());
        }
        registry.fill(HIST("SCT/Cascade/hMassOmega_misid"), cascade.mOmega());
      }
      if (isKaon(bachelor) && isProton(pos) && isPion(neg)) {
        if (cascade.mXi() < fCascadeCut.cfg_min_mass_Xi_veto || fCascadeCut.cfg_max_mass_Xi_veto < cascade.mXi()) {
          registry.fill(HIST("SCT/Cascade/hMassOmega"), cascade.mOmega());
          registry.fill(HIST("SCT/Cascade/hYPhi_Omega"), cascade.phi(), cascade.yOmega());
        }
        registry.fill(HIST("SCT/Cascade/hMassXi_misid"), cascade.mXi());
      }
    } else { // Xi+ or Omega+
      if (isPion(bachelor) && isProton(neg) && isPion(pos)) {
        if (cascade.mOmega() < fCascadeCut.cfg_min_mass_Omega_veto || fCascadeCut.cfg_max_mass_Omega_veto < cascade.mOmega()) {
          registry.fill(HIST("SCT/Cascade/hMassXi"), cascade.mXi());
          registry.fill(HIST("SCT/Cascade/hYPhi_Xi"), cascade.phi(), cascade.yXi());
        }
        registry.fill(HIST("SCT/Cascade/hMassOmega_misid"), cascade.mOmega());
      }
      if (isKaon(bachelor) && isProton(neg) && isPion(pos)) {
        if (cascade.mXi() < fCascadeCut.cfg_min_mass_Xi_veto || fCascadeCut.cfg_max_mass_Xi_veto < cascade.mXi()) {
          registry.fill(HIST("SCT/Cascade/hMassOmega"), cascade.mOmega());
          registry.fill(HIST("SCT/Cascade/hYPhi_Omega"), cascade.phi(), cascade.yOmega());
        }
        registry.fill(HIST("SCT/Cascade/hMassXi_misid"), cascade.mXi());
      }
    }
  }

  template <bool isMC, bool isTriggerAnalysis, typename TBCs, typename TCollisions, typename TTracks, typename TV0s, typename TCascades, typename TTrackAssoc, typename TMCParticles, typename TProducts, typename THistoregistry, typename TSliceCache, typename TPresliceTrack, typename TPresliceTrackAssoc, typename TPresliceV0, typename TPresliceCascade>
  void processWithTTCA(TBCs const& bcs, TCollisions const& collisions, TTracks const& tracks, TV0s const& v0s, TCascades const& cascades, TTrackAssoc const& trackIndices, TMCParticles const&, TProducts& products, THistoregistry& registry, TSliceCache& cache, TPresliceTrack const& perColTrack, TPresliceTrackAssoc const& trackIndicesPerCollision, TPresliceV0 const& perColV0, TPresliceCascade const& perColCasc)
  {
    initCCDB(bcs.begin());

    calculateTOFNSigmaWithReassociation<true>(collisions, bcs, tracks, trackIndices, cache, perColTrack, trackIndicesPerCollision);

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<TBCs>();
      initCCDB(bc);

      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      if (!collision.isSelected()) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      std::unordered_multimap<int, int> multiMapTracksPerCollision; // collisionId -> trackIds

      auto trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      std::vector<int> looseElectronIds; // for pion prefilter
      if (fElectronPFCut.doPF) {
        looseElectronIds.reserve(trackIdsThisCollision.size());
      }

      for (const auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<TTracks>();
        if (!checkTrack<isMC>(collision, track) && !checkTrackPF<isMC>(collision, track)) {
          continue;
        }

        fillMapMLPID(collision, track); // call this function for all track also for loose electron sample.

        if (checkTrack<isMC>(collision, track) && isElectron(collision, track, fElectronCut)) {
          multiMapTracksPerCollision.insert(std::make_pair(collision.globalIndex(), track.globalIndex()));
        }

        if (fElectronPFCut.doPF && checkTrackPF<isMC>(collision, track) && isElectron(collision, track, fElectronPFCut)) {
          looseElectronIds.emplace_back(track.globalIndex());
        }
      } // end of track loop

      int count_electrons = multiMapTracksPerCollision.count(collision.globalIndex());
      registry.fill(HIST("Track/hNe"), count_electrons);
      if (count_electrons < fElectronCut.minNelectron) {
        continue;
      }

      std::vector<int> hadronIds;
      std::vector<int> k0sIds;
      std::vector<int> lambdaIds;
      std::vector<int> antilambdaIds;
      std::vector<int> xiMinusIds;
      std::vector<int> xiPlusIds;
      std::vector<int> omegaMinusIds;
      std::vector<int> omegaPlusIds;

      if (fDoSCTwithTracks) {
        hadronIds.reserve(trackIdsThisCollision.size());
      }

      if (fDoSCTwithV0s) {
        auto v0s_per_collision = v0s.sliceBy(perColV0, collision.globalIndex());
        k0sIds.reserve(v0s_per_collision.size());
        lambdaIds.reserve(v0s_per_collision.size());
        antilambdaIds.reserve(v0s_per_collision.size());

        for (const auto& v0 : v0s_per_collision) {
          auto pos = v0.template posTrack_as<TTracks>();
          auto neg = v0.template negTrack_as<TTracks>();

          if (pos.sign() * neg.sign() > 0) {
            continue;
          }
          if (!isSelectedV0Leg<isMC>(pos) || !isSelectedV0Leg<isMC>(neg)) {
            continue;
          }

          if (v0.dcaV0daughters() > fV0Cut.cfg_max_dca2legs) {
            continue;
          }

          if (v0.v0radius() < fV0Cut.cfg_min_radius) {
            continue;
          }

          if (v0.v0cosPA() < fV0Cut.cfg_min_cospa) {
            continue;
          }

          if (std::sqrt(std::pow(v0.alpha() / fV0Cut.cfg_max_alpha_veto, 2) + std::pow(v0.qtarm() / fV0Cut.cfg_max_qt_veto, 2)) < 1.f) { // photon conversion rejection at small qT
            continue;
          }

          fillV0Histograms<TTracks>(collision, v0, registry);

          if (isK0S(v0) && isPion(pos) && isPion(neg)) {
            k0sIds.emplace_back(v0.globalIndex());
          }

          if (isLambda(v0) && isProton(pos) && isPion(neg)) {
            lambdaIds.emplace_back(v0.globalIndex());
          } else if (isAntiLambda(v0) && isProton(neg) && isPion(pos)) {
            antilambdaIds.emplace_back(v0.globalIndex());
          }

        } // end of v0 loop
      }

      if (fDoSCTwithCascades) {
        auto cascades_per_collision = cascades.sliceBy(perColCasc, collision.globalIndex());
        xiMinusIds.reserve(cascades_per_collision.size());
        xiPlusIds.reserve(cascades_per_collision.size());
        omegaMinusIds.reserve(cascades_per_collision.size());
        omegaPlusIds.reserve(cascades_per_collision.size());

        for (const auto& cascade : cascades_per_collision) {
          auto pos = cascade.template posTrack_as<TTracks>();
          auto neg = cascade.template negTrack_as<TTracks>();
          auto bachelor = cascade.template bachelor_as<TTracks>();
          if (pos.sign() * neg.sign() > 0) {
            continue;
          }
          if (cascade.mLambda() < fCascadeCut.cfg_min_mass_lambda || fCascadeCut.cfg_max_mass_lambda < cascade.mLambda()) {
            continue;
          }

          if (!isSelectedV0Leg<isMC>(pos) || !isSelectedV0Leg<isMC>(neg) || !isSelectedV0Leg<isMC>(bachelor)) {
            continue;
          }

          if (cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < fCascadeCut.cfg_min_cospa) {
            continue;
          }
          if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < fCascadeCut.cfg_min_cospa_v0) {
            continue;
          }

          if (cascade.cascradius() > cascade.v0radius()) {
            continue;
          }

          if (cascade.dcaV0daughters() > fCascadeCut.cfg_max_dcadau_v0) {
            continue;
          }
          if (cascade.v0radius() < fCascadeCut.cfg_min_rxy_v0) {
            continue;
          }
          if (cascade.cascradius() < fCascadeCut.cfg_min_rxy) {
            continue;
          }

          if (cascade.dcacascdaughters() > fCascadeCut.cfg_max_dcadau) {
            continue;
          }

          if (std::fabs(cascade.dcav0topv(collision.posX(), collision.posY(), collision.posZ())) < fCascadeCut.cfg_min_dcaxy_v0) {
            continue;
          }

          fillCascadeHistograms<TTracks>(collision, cascade, registry);

          if (cascade.sign() < 0) { // Xi- or Omega-
            if (isXi(cascade) && isPion(bachelor) && isProton(pos) && isPion(neg)) {
              xiMinusIds.emplace_back(cascade.globalIndex());
            }
            if (isOmega(cascade) && isKaon(bachelor) && isProton(pos) && isPion(neg)) {
              omegaMinusIds.emplace_back(cascade.globalIndex());
            }
          } else { // Xi+ or Omega+
            if (isXi(cascade) && isPion(bachelor) && isProton(neg) && isPion(pos)) {
              xiPlusIds.emplace_back(cascade.globalIndex());
            }
            if (isOmega(cascade) && isKaon(bachelor) && isProton(neg) && isPion(pos)) {
              omegaPlusIds.emplace_back(cascade.globalIndex());
            }
          }
        } // end of cascade loop
      }

      for (const auto& trackId : trackIdsThisCollision) {
        auto track = trackId.template track_as<TTracks>();
        auto trackParCov = getTrackParCov(track);
        o2::dataformats::DCA mDcaInfoCov;
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        trackParCov.setPID(track.pidForTracking());
        bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        if (!isPropOK) {
          continue;
        }
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();
        if (isSelectedHadron(collision, track, trackParCov, dcaXY, dcaZ)) {
          float tpcSignal = track.tpcSignal();
          if constexpr (isMC) {
            tpcSignal = track.mcTunedTPCSignal();
          }
          registry.fill(HIST("SCT/Track/hs"), trackParCov.getPt(), trackParCov.getEta(), RecoDecay::constrainAngle(trackParCov.getPhi(), 0, 1U));
          registry.fill(HIST("SCT/Track/hDCA"), dcaXY, dcaZ);
          registry.fill(HIST("SCT/Track/hTPCdEdx"), track.tpcInnerParam(), tpcSignal);
          registry.fill(HIST("SCT/Track/hTOFbeta"), trackParCov.getP(), fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), track.globalIndex())]);
          hadronIds.emplace_back(track.globalIndex());
        }
      } // end of track loop

      auto range_electrons = multiMapTracksPerCollision.equal_range(collision.globalIndex());
      for (auto it = range_electrons.first; it != range_electrons.second; it++) {
        auto electron = tracks.rawIteratorAt(it->second);
        o2::dataformats::DCA mDcaInfoCov;
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto trackParCov = getTrackParCov(electron);
        trackParCov.setPID(o2::track::PID::Electron);
        bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        if (!isPropOK) {
          continue;
        }
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();
        fillElectronTable<isMC>(collision, electron, trackParCov, dcaXY, dcaZ, products, registry);
        ROOT::Math::PtEtaPhiMVector v1(trackParCov.getPt(), trackParCov.getEta(), RecoDecay::constrainAngle(trackParCov.getPhi(), 0, 1U), o2::constants::physics::MassElectron); // main electron

        // apply pion prefilter
        uint8_t pfb = 0;
        for (const auto& looseElectronId : looseElectronIds) {
          auto looseElectron = tracks.rawIteratorAt(looseElectronId);
          if (looseElectron.globalIndex() == electron.globalIndex()) {
            continue;
          }
          if (looseElectron.sign() * electron.sign() > 0) { // reject LS, because pi0dalitz and photon conversion must be ULS.
            continue;
          }

          auto trackParCov2 = getTrackParCov(looseElectron);
          trackParCov2.setPID(o2::track::PID::Electron);
          bool isPropOK = o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov2, 2.f, matCorr, &mDcaInfoCov);
          if (!isPropOK) {
            continue;
          }

          float tpcSignal = looseElectron.tpcSignal();
          if constexpr (isMC) {
            tpcSignal = looseElectron.mcTunedTPCSignal();
          }
          registry.fill(HIST("Prefilter/Track/hPt"), trackParCov2.getPt());
          registry.fill(HIST("Prefilter/Track/hEtaPhi"), RecoDecay::constrainAngle(trackParCov2.getPhi(), 0, 1U), trackParCov2.getEta());
          registry.fill(HIST("Prefilter/Track/hTPCdEdx"), looseElectron.tpcInnerParam(), tpcSignal);
          registry.fill(HIST("Prefilter/Track/hTOFbeta"), trackParCov2.getP(), fMapTOFBetaReassociated[std::make_pair(collision.globalIndex(), looseElectron.globalIndex())]);

          ROOT::Math::PtEtaPhiMVector v2(trackParCov2.getPt(), trackParCov2.getEta(), RecoDecay::constrainAngle(trackParCov2.getPhi(), 0, 1U), o2::constants::physics::MassElectron); // loose electron
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(v1.Px(), v1.Py(), v1.Pz(), v2.Px(), v2.Py(), v2.Pz(), electron.sign(), looseElectron.sign(), d_bz);
          registry.fill(HIST("Prefilter/Pair/hMvsPhiV"), phiv, v12.M());

          if (v12.M() < fElectronPFCut.slope * phiv + fElectronPFCut.intercept) {
            pfb |= (uint8_t(1) << static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPC));
          }

          for (int i = 0; i < static_cast<int>(max_mee_vec.size()); i++) {
            if (v12.M() < max_mee_vec.at(i)) {
              pfb |= (uint8_t(1) << (static_cast<int>(o2::aod::pwgem::dilepton::utils::pairutil::DileptonPrefilterBit::kElFromPi0_20MeV) + i));
            }
          } // end of mee loop

        } // end of loose electron loop
        products.electronPFTable(pfb);

        // perform SCT from here
        std::vector<float> bdtScorePrompt;
        std::vector<float> bdtScorePromptHc;
        std::vector<float> bdtScoreNonpromptHc;
        std::vector<float> bdtScoreHb;
        std::vector<uint8_t> hadronType;

        bdtScorePrompt.reserve(hadronIds.size() + k0sIds.size() + lambdaIds.size() + antilambdaIds.size());
        bdtScorePromptHc.reserve(hadronIds.size() + k0sIds.size() + lambdaIds.size() + antilambdaIds.size());
        bdtScoreNonpromptHc.reserve(hadronIds.size() + k0sIds.size() + lambdaIds.size() + antilambdaIds.size());
        bdtScoreHb.reserve(hadronIds.size() + k0sIds.size() + lambdaIds.size() + antilambdaIds.size());
        hadronType.reserve(hadronIds.size() + k0sIds.size() + lambdaIds.size() + antilambdaIds.size());

        // eTrack pair
        for (const auto& hadronId : hadronIds) {
          auto hadron = tracks.rawIteratorAt(hadronId);
          if (hadron.globalIndex() == electron.globalIndex()) {
            continue;
          }

          auto hadronParCov = getTrackParCov(hadron);
          o2::dataformats::DCA mDcaInfoCov;
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          hadronParCov.setPID(hadron.pidForTracking());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, hadronParCov, 2.f, matCorr, &mDcaInfoCov);

          auto eTpair = o2::aod::pwgem::dilepton::utils::makePairLeptonTrack(dfeT, collision, electron, hadron, o2::track::PID::Electron, hadron.pidForTracking());
          registry.fill(HIST("SCT/eT/hImpactParameter"), eTpair.impParXY, eTpair.impParZ);
          registry.fill(HIST("SCT/eT/hDecayLength"), eTpair.lxy, eTpair.lz);
          registry.fill(HIST("SCT/eT/hCosPA"), eTpair.cospa);
          registry.fill(HIST("SCT/eT/hDCA2legs"), eTpair.dca2legs);
          registry.fill(HIST("SCT/eT/hMass"), eTpair.mass);
          if (eTpair.isOK && fConfigDFeT.useML) {
            o2::analysis::pwgem::dilepton::sct::candidate candidate;
            fillCandidate(candidate, eTpair, trackParCov, mDcaInfoCov);
            candidate.ptL = trackParCov.getPt();
            candidate.signLH = electron.sign() * hadron.sign();
            candidate.signedMassLH = electron.sign() * hadron.sign() * eTpair.mass;
            candidate.tpcNSigmaKa = hadron.tpcNSigmaKa();

            std::vector<float> inputFeatures = mlResponseSCTeT.getInputFeatures(candidate);
            float binningFeature = mlResponseSCTeT.getBinningFeature(candidate);

            int pbin = lower_bound(fConfigDFeT.binsMl.value.begin(), fConfigDFeT.binsMl.value.end(), binningFeature) - fConfigDFeT.binsMl.value.begin() - 1;
            if (pbin < 0) {
              pbin = 0;
            } else if (static_cast<int>(fConfigDFeT.binsMl.value.size()) - 2 < pbin) {
              pbin = static_cast<int>(fConfigDFeT.binsMl.value.size()) - 2;
            }

            float probaPrompt = mlResponseSCTeT.getModelOutput(inputFeatures, pbin)[0];
            float probaPromptHc = mlResponseSCTeT.getModelOutput(inputFeatures, pbin)[1];
            float probaNonpromptHc = mlResponseSCTeT.getModelOutput(inputFeatures, pbin)[2];
            float probaHb = mlResponseSCTeT.getModelOutput(inputFeatures, pbin)[3];

            bdtScorePrompt.emplace_back(probaPrompt);
            bdtScorePromptHc.emplace_back(probaPromptHc);
            bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
            bdtScoreHb.emplace_back(probaHb);
            hadronType.emplace_back(0);
          }
        } // end of charged track loop

        // eK0S pair
        for (const auto& k0sId : k0sIds) {
          auto v0 = v0s.rawIteratorAt(k0sId);
          if (v0.posTrackId() == electron.globalIndex() || v0.negTrackId() == electron.globalIndex()) {
            continue;
          }
          const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
          const std::array<float, 3> momV0 = {v0.px(), v0.py(), v0.pz()};
          std::array<float, 21> covV0 = {0.};
          constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
          for (int i = 0; i < 6; i++) {
            covV0[MomInd[i]] = v0.momentumCovMat()[i];
            covV0[i] = v0.positionCovMat()[i];
          }
          auto v0ParCov = o2::track::TrackParCov(vertexV0, momV0, covV0, 0, true);
          v0ParCov.setAbsCharge(0);
          v0ParCov.setPID(o2::track::PID::K0);
          o2::dataformats::DCA impactParameterV0;
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, v0ParCov, 2.f, matCorr, &impactParameterV0); // v0ParCov is TrackParCov object

          auto eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(dfeV0, collision, electron, v0, o2::track::PID::Electron, o2::track::PID::K0);
          registry.fill(HIST("SCT/eV0/hImpactParameter"), eV0pair.impParXY, eV0pair.impParZ);
          registry.fill(HIST("SCT/eV0/hDecayLength"), eV0pair.lxy, eV0pair.lz);
          registry.fill(HIST("SCT/eV0/hCosPA"), eV0pair.cospa);
          registry.fill(HIST("SCT/eV0/hDCA2legs"), eV0pair.dca2legs);
          registry.fill(HIST("SCT/eV0/hMass"), eV0pair.mass);
          if (eV0pair.isOK && fConfigDFeV0.useML) {
            o2::analysis::pwgem::dilepton::sct::candidate candidate;
            fillCandidate(candidate, eV0pair, v0ParCov, impactParameterV0);
            candidate.ptL = trackParCov.getPt();

            std::vector<float> inputFeatures = mlResponseSCTeV0.getInputFeatures(candidate);
            float binningFeature = mlResponseSCTeV0.getBinningFeature(candidate);

            int pbin = lower_bound(fConfigDFeV0.binsMl.value.begin(), fConfigDFeV0.binsMl.value.end(), binningFeature) - fConfigDFeV0.binsMl.value.begin() - 1;
            if (pbin < 0) {
              pbin = 0;
            } else if (static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2 < pbin) {
              pbin = static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2;
            }

            float probaPrompt = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[0];
            float probaPromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[1];
            float probaNonpromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[2];
            float probaHb = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[3];
            bdtScorePrompt.emplace_back(probaPrompt);
            bdtScorePromptHc.emplace_back(probaPromptHc);
            bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
            bdtScoreHb.emplace_back(probaHb);
            hadronType.emplace_back(1);
          }
        } // end of k0s loop

        if (electron.sign() > 0) { // positron
          // eL pair // sign is restricted in baryon decay: Lc+ -> e+ nu_e Lambda
          for (const auto& lambdaId : lambdaIds) {
            auto v0 = v0s.rawIteratorAt(lambdaId);
            if (v0.posTrackId() == electron.globalIndex() || v0.negTrackId() == electron.globalIndex()) {
              continue;
            }
            const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
            const std::array<float, 3> momV0 = {v0.px(), v0.py(), v0.pz()};
            std::array<float, 21> covV0 = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covV0[MomInd[i]] = v0.momentumCovMat()[i];
              covV0[i] = v0.positionCovMat()[i];
            }
            auto v0ParCov = o2::track::TrackParCov(vertexV0, momV0, covV0, 0, true);
            v0ParCov.setAbsCharge(0);
            v0ParCov.setPID(o2::track::PID::Lambda);
            o2::dataformats::DCA impactParameterV0;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, v0ParCov, 2.f, matCorr, &impactParameterV0); // v0ParCov is TrackParCov object

            auto eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(dfeV0, collision, electron, v0, o2::track::PID::Electron, o2::track::PID::Lambda);
            registry.fill(HIST("SCT/eV0/hImpactParameter"), eV0pair.impParXY, eV0pair.impParZ);
            registry.fill(HIST("SCT/eV0/hDecayLength"), eV0pair.lxy, eV0pair.lz);
            registry.fill(HIST("SCT/eV0/hCosPA"), eV0pair.cospa);
            registry.fill(HIST("SCT/eV0/hDCA2legs"), eV0pair.dca2legs);
            registry.fill(HIST("SCT/eV0/hMass"), eV0pair.mass);
            if (eV0pair.isOK && fConfigDFeV0.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eV0pair, v0ParCov, impactParameterV0);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeV0.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeV0.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeV0.binsMl.value.begin(), fConfigDFeV0.binsMl.value.end(), binningFeature) - fConfigDFeV0.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(2);
            }
          } // end of lambda loop

          // eXi pair // sign is restricted in baryon decay: Xic0 -> e+ nu_e Xi-
          for (const auto& xiId : xiMinusIds) {
            auto cascade = cascades.rawIteratorAt(xiId);
            if (cascade.posTrackId() == electron.globalIndex() || cascade.negTrackId() == electron.globalIndex() || cascade.bachelorId() == electron.globalIndex()) {
              continue;
            }

            const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
            const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};
            std::array<float, 21> covCasc = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
              covCasc[i] = cascade.positionCovMat()[i];
            }
            auto cascadeParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
            cascadeParCov.setAbsCharge(1);
            cascadeParCov.setPID(o2::track::PID::XiMinus);
            o2::dataformats::DCA impactParameterCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, cascadeParCov, 2.f, matCorr, &impactParameterCasc); // cascadeParCov is TrackParCov object

            auto eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(dfeC, collision, electron, cascade, o2::track::PID::Electron, o2::track::PID::XiMinus);
            registry.fill(HIST("SCT/eC/hImpactParameter"), eCpair.impParXY, eCpair.impParZ);
            registry.fill(HIST("SCT/eC/hDecayLength"), eCpair.lxy, eCpair.lz);
            registry.fill(HIST("SCT/eC/hCosPA"), eCpair.cospa);
            registry.fill(HIST("SCT/eC/hDCA2legs"), eCpair.dca2legs);
            registry.fill(HIST("SCT/eC/hMass"), eCpair.mass);
            if (eCpair.isOK && fConfigDFeC.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eCpair, cascadeParCov, impactParameterCasc);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeC.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeC.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeC.binsMl.value.begin(), fConfigDFeC.binsMl.value.end(), binningFeature) - fConfigDFeC.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(4);
            }
          } // end of Xi- loop

          // eOmega pair // sign is restricted in baryon decay: Omegac0 -> e+ nu_e Omega-
          for (const auto& omegaId : omegaMinusIds) {
            auto cascade = cascades.rawIteratorAt(omegaId);
            if (cascade.posTrackId() == electron.globalIndex() || cascade.negTrackId() == electron.globalIndex() || cascade.bachelorId() == electron.globalIndex()) {
              continue;
            }

            const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
            const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};
            std::array<float, 21> covCasc = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
              covCasc[i] = cascade.positionCovMat()[i];
            }
            auto cascadeParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
            cascadeParCov.setAbsCharge(1);
            cascadeParCov.setPID(o2::track::PID::OmegaMinus);
            o2::dataformats::DCA impactParameterCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, cascadeParCov, 2.f, matCorr, &impactParameterCasc); // cascadeParCov is TrackParCov object

            auto eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(dfeC, collision, electron, cascade, o2::track::PID::Electron, o2::track::PID::OmegaMinus);
            registry.fill(HIST("SCT/eC/hImpactParameter"), eCpair.impParXY, eCpair.impParZ);
            registry.fill(HIST("SCT/eC/hDecayLength"), eCpair.lxy, eCpair.lz);
            registry.fill(HIST("SCT/eC/hCosPA"), eCpair.cospa);
            registry.fill(HIST("SCT/eC/hDCA2legs"), eCpair.dca2legs);
            registry.fill(HIST("SCT/eC/hMass"), eCpair.mass);
            if (eCpair.isOK && fConfigDFeC.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eCpair, cascadeParCov, impactParameterCasc);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeC.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeC.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeC.binsMl.value.begin(), fConfigDFeC.binsMl.value.end(), binningFeature) - fConfigDFeC.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(6);
            }
          } // end of Omega- loop

        } else { // electron
          // eL pair // sign is restricted in baryon decay: Lc- -> e- anti_nu_e antiLambda
          for (const auto& antilambdaId : antilambdaIds) {
            auto v0 = v0s.rawIteratorAt(antilambdaId);
            if (v0.posTrackId() == electron.globalIndex() || v0.negTrackId() == electron.globalIndex()) {
              continue;
            }
            const std::array<float, 3> vertexV0 = {v0.x(), v0.y(), v0.z()};
            const std::array<float, 3> momV0 = {v0.px(), v0.py(), v0.pz()};
            std::array<float, 21> covV0 = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covV0[MomInd[i]] = v0.momentumCovMat()[i];
              covV0[i] = v0.positionCovMat()[i];
            }
            auto v0ParCov = o2::track::TrackParCov(vertexV0, momV0, covV0, 0, true);
            v0ParCov.setAbsCharge(0);
            v0ParCov.setPID(o2::track::PID::Lambda);
            o2::dataformats::DCA impactParameterV0;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, v0ParCov, 2.f, matCorr, &impactParameterV0); // v0ParCov is TrackParCov object

            auto eV0pair = o2::aod::pwgem::dilepton::utils::makePairLeptonV0(dfeV0, collision, electron, v0, o2::track::PID::Electron, o2::track::PID::Lambda);
            registry.fill(HIST("SCT/eV0/hImpactParameter"), eV0pair.impParXY, eV0pair.impParZ);
            registry.fill(HIST("SCT/eV0/hDecayLength"), eV0pair.lxy, eV0pair.lz);
            registry.fill(HIST("SCT/eV0/hCosPA"), eV0pair.cospa);
            registry.fill(HIST("SCT/eV0/hDCA2legs"), eV0pair.dca2legs);
            registry.fill(HIST("SCT/eV0/hMass"), eV0pair.mass);
            if (eV0pair.isOK && fConfigDFeV0.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eV0pair, v0ParCov, impactParameterV0);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeV0.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeV0.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeV0.binsMl.value.begin(), fConfigDFeV0.binsMl.value.end(), binningFeature) - fConfigDFeV0.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeV0.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeV0.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(3);
            }
          } // end of antilambda loop

          // eXi pair // sign is restricted in baryon decay: Xic0bar -> e- anti_nu_e Xi+
          for (const auto& xiId : xiPlusIds) {
            auto cascade = cascades.rawIteratorAt(xiId);
            if (cascade.posTrackId() == electron.globalIndex() || cascade.negTrackId() == electron.globalIndex() || cascade.bachelorId() == electron.globalIndex()) {
              continue;
            }

            const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
            const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};
            std::array<float, 21> covCasc = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
              covCasc[i] = cascade.positionCovMat()[i];
            }
            auto cascadeParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
            cascadeParCov.setAbsCharge(1);
            cascadeParCov.setPID(o2::track::PID::XiMinus);
            o2::dataformats::DCA impactParameterCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, cascadeParCov, 2.f, matCorr, &impactParameterCasc); // cascadeParCov is TrackParCov object

            auto eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(dfeC, collision, electron, cascade, o2::track::PID::Electron, o2::track::PID::XiMinus);
            registry.fill(HIST("SCT/eC/hImpactParameter"), eCpair.impParXY, eCpair.impParZ);
            registry.fill(HIST("SCT/eC/hDecayLength"), eCpair.lxy, eCpair.lz);
            registry.fill(HIST("SCT/eC/hCosPA"), eCpair.cospa);
            registry.fill(HIST("SCT/eC/hDCA2legs"), eCpair.dca2legs);
            registry.fill(HIST("SCT/eC/hMass"), eCpair.mass);
            if (eCpair.isOK && fConfigDFeC.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eCpair, cascadeParCov, impactParameterCasc);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeC.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeC.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeC.binsMl.value.begin(), fConfigDFeC.binsMl.value.end(), binningFeature) - fConfigDFeC.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(5);
            }
          } // end of Xi- loop

          // eOmega pair // sign is restricted in baryon decay: Omegac0bar -> e- anti_nu_e Omega+
          for (const auto& omegaId : omegaPlusIds) {
            auto cascade = cascades.rawIteratorAt(omegaId);
            if (cascade.posTrackId() == electron.globalIndex() || cascade.negTrackId() == electron.globalIndex() || cascade.bachelorId() == electron.globalIndex()) {
              continue;
            }

            const std::array<float, 3> vertexCasc = {cascade.x(), cascade.y(), cascade.z()};
            const std::array<float, 3> momCasc = {cascade.px(), cascade.py(), cascade.pz()};
            std::array<float, 21> covCasc = {0.};
            constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
            for (int i = 0; i < 6; i++) {
              covCasc[MomInd[i]] = cascade.momentumCovMat()[i];
              covCasc[i] = cascade.positionCovMat()[i];
            }
            auto cascadeParCov = o2::track::TrackParCov(vertexCasc, momCasc, covCasc, cascade.sign(), true);
            cascadeParCov.setAbsCharge(1);
            cascadeParCov.setPID(o2::track::PID::OmegaMinus);
            o2::dataformats::DCA impactParameterCasc;
            o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, cascadeParCov, 2.f, matCorr, &impactParameterCasc); // cascadeParCov is TrackParCov object

            auto eCpair = o2::aod::pwgem::dilepton::utils::makePairLeptonCascade(dfeC, collision, electron, cascade, o2::track::PID::Electron, o2::track::PID::OmegaMinus);
            registry.fill(HIST("SCT/eC/hImpactParameter"), eCpair.impParXY, eCpair.impParZ);
            registry.fill(HIST("SCT/eC/hDecayLength"), eCpair.lxy, eCpair.lz);
            registry.fill(HIST("SCT/eC/hCosPA"), eCpair.cospa);
            registry.fill(HIST("SCT/eC/hDCA2legs"), eCpair.dca2legs);
            registry.fill(HIST("SCT/eC/hMass"), eCpair.mass);
            if (eCpair.isOK && fConfigDFeC.useML) {
              o2::analysis::pwgem::dilepton::sct::candidate candidate;
              fillCandidate(candidate, eCpair, cascadeParCov, impactParameterCasc);
              candidate.ptL = trackParCov.getPt();

              std::vector<float> inputFeatures = mlResponseSCTeC.getInputFeatures(candidate);
              float binningFeature = mlResponseSCTeC.getBinningFeature(candidate);

              int pbin = lower_bound(fConfigDFeC.binsMl.value.begin(), fConfigDFeC.binsMl.value.end(), binningFeature) - fConfigDFeC.binsMl.value.begin() - 1;
              if (pbin < 0) {
                pbin = 0;
              } else if (static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2 < pbin) {
                pbin = static_cast<int>(fConfigDFeC.binsMl.value.size()) - 2;
              }

              float probaPrompt = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[0];
              float probaPromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[1];
              float probaNonpromptHc = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[2];
              float probaHb = mlResponseSCTeC.getModelOutput(inputFeatures, pbin)[3];
              bdtScorePrompt.emplace_back(probaPrompt);
              bdtScorePromptHc.emplace_back(probaPromptHc);
              bdtScoreNonpromptHc.emplace_back(probaNonpromptHc);
              bdtScoreHb.emplace_back(probaHb);
              hadronType.emplace_back(7);
            }
          } // end of Omega- loop
        }

        products.sctTable(bdtScorePrompt, bdtScorePromptHc, bdtScoreNonpromptHc, bdtScoreHb, hadronType);

        bdtScorePrompt.clear();
        bdtScorePromptHc.clear();
        bdtScoreNonpromptHc.clear();
        bdtScoreHb.clear();
        hadronType.clear();
        bdtScorePrompt.shrink_to_fit();
        bdtScorePromptHc.shrink_to_fit();
        bdtScoreNonpromptHc.shrink_to_fit();
        bdtScoreHb.shrink_to_fit();
        hadronType.shrink_to_fit();
      } // end of main electron loop

      hadronIds.clear();
      hadronIds.shrink_to_fit();

      k0sIds.clear();
      k0sIds.shrink_to_fit();

      lambdaIds.clear();
      lambdaIds.shrink_to_fit();
      antilambdaIds.clear();
      antilambdaIds.shrink_to_fit();

      xiMinusIds.clear();
      xiMinusIds.shrink_to_fit();
      xiPlusIds.clear();
      xiPlusIds.shrink_to_fit();

      omegaMinusIds.clear();
      omegaMinusIds.shrink_to_fit();
      omegaPlusIds.clear();
      omegaPlusIds.shrink_to_fit();

      looseElectronIds.clear();
      looseElectronIds.shrink_to_fit();
      multiMapTracksPerCollision.clear();
    } // end of collision loop
    clear();
  }

  template <bool isMC, bool isTriggerAnalysis, typename TBCs, typename TCollisions, typename TTracks, typename TV0s, typename TCascades, typename TMCParticles, typename TProducts, typename THistoregistry>
  void processWithoutTTCA(TBCs const& bcs, TCollisions const& collisions, TTracks const&, TV0s const&, TCascades const&, TMCParticles const&, TProducts&, THistoregistry&)
  {
    LOGF(info, "processWithoutTTCA is not supported. Bye.");
    return;

    initCCDB(bcs.begin());
    // calculateTOFNSigmaWithReassociation<false>(collisions, bcs, tracks, nullptr);

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<TBCs>();
      initCCDB(bc);

      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
      }

      if (!collision.isSelected()) {
        continue;
      }

      if constexpr (isTriggerAnalysis) {
        if (collision.swtaliastmp_raw() == 0) {
          continue;
        }
      }

    } // end of collision loop

    clear();
  }

  void clear()
  {
    // this should be called at the end of DF.
    fMapProbaEl.clear();
    fMapTOFNsigmaElReassociated.clear();
    fMapTOFNsigmaPiReassociated.clear();
    fMapTOFNsigmaKaReassociated.clear();
    fMapTOFNsigmaPrReassociated.clear();
    fMapTOFBetaReassociated.clear();
  }

  void doSCTwithTracks(const bool flag) { fDoSCTwithTracks = flag; }
  void doSCTwithV0s(const bool flag) { fDoSCTwithV0s = flag; }
  void doSCTwithCascades(const bool flag) { fDoSCTwithCascades = flag; }

  float dca3DinSigmaOTF(const float dcaXY, const float dcaZ, const float cYY, const float cZZ, const float cZY)
  {
    float det = cYY * cZZ - cZY * cZY; // determinant
    if (det < 0) {
      return 999.f;
    } else {
      return std::sqrt(std::fabs((dcaXY * dcaXY * cZZ + dcaZ * dcaZ * cYY - 2. * dcaXY * dcaZ * cZY) / det / 2.)); // dca 3d in sigma
    }
  }

  template <typename TCandidate, typename TPair, typename TTrackParCov, typename TDCA>
  void fillCandidate(TCandidate& candidate, TPair const& pair, TTrackParCov const& trackParCov, TDCA const& dcaInfo)
  {
    candidate.ptH = trackParCov.getPt();
    candidate.tpcNSigmaKa = 0;
    candidate.impParXYH = dcaInfo.getY();
    candidate.impParZH = dcaInfo.getZ();
    candidate.impPar3DH = std::sqrt(std::pow(candidate.impParXYH, 2) + std::pow(candidate.impParZH, 2));
    candidate.impParXYHinSigma = candidate.impParXYH / std::sqrt(trackParCov.getSigmaY2());
    candidate.impParZHinSigma = candidate.impParZH / std::sqrt(trackParCov.getSigmaZ2());
    candidate.impPar3DHinSigma = dca3DinSigmaOTF(candidate.impParXYH, candidate.impParZH, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());
    candidate.signLH = 0;
    candidate.dcaLH = pair.dca2legs;
    candidate.massLH = pair.mass;
    candidate.ptLH = pair.pt;
    candidate.signedMassLH = pair.mass;
    candidate.cpa = pair.cospa;
    candidate.cpaXY = pair.cospaXY;
    candidate.impParXY = pair.impParXY;
    candidate.impParZ = pair.impParZ;
    candidate.impPar3D = std::sqrt(std::pow(candidate.impParXY, 2) + std::pow(candidate.impParZ, 2));
    candidate.impParXYinSigma = candidate.impParXY / std::sqrt(pair.impParCYY);
    candidate.impParZinSigma = candidate.impParZ / std::sqrt(pair.impParCZZ);
    candidate.impPar3DinSigma = dca3DinSigmaOTF(candidate.impParXY, candidate.impParZ, pair.impParCYY, pair.impParCZY, pair.impParCZZ);
    candidate.decayLengthXY = pair.lxy;
    candidate.decayLengthZ = pair.lz;
    candidate.decayLength3D = pair.lxyz;
    candidate.decayLengthXYinSigma = pair.lxy / pair.lxyErr;
    candidate.decayLengthZinSigma = pair.lz / pair.lzErr;
    candidate.decayLength3DinSigma = pair.lxyz / pair.lxyzErr;
  }

 protected:
  electronCut fElectronCut;
  electronPFCut fElectronPFCut;
  hadronCut fHadronCut;
  v0Cut fV0Cut;
  cascadeCut fCascadeCut;
  cfgDFeT fConfigDFeT;
  cfgDFeV0 fConfigDFeV0;
  cfgDFeC fConfigDFeC;

  std::map<std::pair<int, int>, float> fMapProbaEl;                 // map pair(collisionId, trackId) -> probaEl
  std::map<std::pair<int, int>, float> fMapTOFNsigmaElReassociated; // map pair(collisionId, trackId) -> tof n sigma el
  std::map<std::pair<int, int>, float> fMapTOFNsigmaPiReassociated; // map pair(collisionId, trackId) -> tof n sigma pi
  std::map<std::pair<int, int>, float> fMapTOFNsigmaKaReassociated; // map pair(collisionId, trackId) -> tof n sigma ka
  std::map<std::pair<int, int>, float> fMapTOFNsigmaPrReassociated; // map pair(collisionId, trackId) -> tof n sigma pr
  std::map<std::pair<int, int>, float> fMapTOFBetaReassociated;     // map pair(collisionId, trackId) -> tof beta

  int mRunNumber{0};
  float d_bz{0};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  o2::framework::Service<o2::pid::tof::TOFResponse> mTOFResponse;
  o2::analysis::MlResponsePID<float> mlResponsePID;
  o2::analysis::MlResponseSCT<float> mlResponseSCTeT;
  o2::analysis::MlResponseSCT<float> mlResponseSCTeV0;
  o2::analysis::MlResponseSCT<float> mlResponseSCTeC;

  bool fDoSCTwithTracks{false};
  bool fDoSCTwithV0s{false};
  bool fDoSCTwithCascades{false};
  const std::vector<float> max_mee_vec{0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14};
  o2::vertexing::DCAFitterN<2> dfeT;
  o2::vertexing::DCAFitterN<2> dfeV0;
  o2::vertexing::DCAFitterN<2> dfeC;
  o2::ccdb::CcdbApi ccdbApi;

}; // end ElectronModule

} // namespace pwgem::dilepton::utils
} // namespace o2

#endif // PWGEM_DILEPTON_UTILS_ELECTRONMODULE_H_
