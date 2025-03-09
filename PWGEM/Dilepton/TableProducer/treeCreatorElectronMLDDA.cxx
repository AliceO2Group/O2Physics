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
//
// ========================
//
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include <vector>
#include <array>
#include <random>
#include <string>
#include "Math/Vector4D.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;

struct TreeCreatorElectronMLDDA {
  enum class EM_V0_Label : int { // Reconstructed V0
    kUndef = -1,
    kGamma = 0,
    kK0S = 1,
    kLambda = 2,
    kAntiLambda = 3,
  };

  SliceCache cache;
  Produces<o2::aod::EMPrimaryTracks> emprimarytracks; // flat table containing collision + track information

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"V0/hAP", "Armenteros Podolanski", {HistType::kTH2F, {{200, -1.f, +1.f}, {250, 0, 0.25}}}},
      {"V0/hXY_Gamma", "photon conversion point in XY;X (cm);Y (cm)", {HistType::kTH2F, {{400, -100, +100}, {400, -100, +100}}}},
      {"V0/hMassGamma_Rxy", "V0 mass gamma", {HistType::kTH2F, {{200, 0, 100}, {100, 0, 0.1}}}},
      {"V0/hCosPA", "V0 cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"V0/hPCA", "V0 distance between 2 legs", {HistType::kTH1F, {{200, 0.f, 2.f}}}},
      {"V0/hMassGamma", "V0 mass gamma", {HistType::kTH1F, {{100, 0, 0.1}}}},
      {"V0/hMassK0Short", "V0 mass K0S", {HistType::kTH1F, {{200, 0.4, 0.6}}}},
      {"V0/hMassLambda", "V0 mass Lambda", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"V0/hMassAntiLambda", "V0 mass AntiLambda", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"hMvsPhiV", "mee vs. phiv", {HistType::kTH2F, {{72, 0, M_PI}, {100, 0, 0.1}}}},

      {"V0/hTPCdEdx_P_El", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Mu", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Pi", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Ka", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      {"V0/hTPCdEdx_P_Pr", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},

      {"V0/hTOFbeta_P_El", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Mu", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Pi", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Ka", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      {"V0/hTOFbeta_P_Pr", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},

      {"V0/hITSClusterSize_P_El", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},
      {"V0/hITSClusterSize_P_Mu", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},
      {"V0/hITSClusterSize_P_Pi", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},
      {"V0/hITSClusterSize_P_Ka", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},
      {"V0/hITSClusterSize_P_Pr", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},

      // {"PrimaryTrack/hTPCdEdx_P", "TPC dEdx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0, 5}, {200, 0, 200}}}},
      // {"PrimaryTrack/hTOFbeta_P", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0, 5}, {220, 0, 1.1}}}},
      // {"PrimaryTrack/hITSClusterSize_P", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {150, 0.0, 15}}}},
      // {"PrimaryTrack/hTPCNsigmaEl_P", "TPC n#sigma_{e} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTPCNsigmaMu_P", "TPC n#sigma_{#mu} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{#mu}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTPCNsigmaPi_P", "TPC n#sigma_{#pi} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTPCNsigmaKa_P", "TPC n#sigma_{K} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{K}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTPCNsigmaPr_P", "TPC n#sigma_{p} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{p}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTOFNsigmaEl_P", "TOF n#sigma_{e} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTOFNsigmaMu_P", "TOF n#sigma_{#mu} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{#mu}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTOFNsigmaPi_P", "TOF n#sigma_{#pi} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTOFNsigmaKa_P", "TOF n#sigma_{K} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{K}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      // {"PrimaryTrack/hTOFNsigmaPr_P", "TOF n#sigma_{p} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{p}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},

      {"Cascade/hRxy_Xi", "R_{xy} of cascade vs. mass;m_{#Lambda#pi};R_{xy} (cm)", {HistType::kTH2F, {{200, 1.2, 1.4}, {200, 0, 20.f}}}},
      {"Cascade/hRxy_Omega", "R_{xy} of cascade vs. mass;m_{#LambdaK};R_{xy} (cm)", {HistType::kTH2F, {{200, 1.6, 1.8}, {200, 0, 20.f}}}},
      {"Cascade/hCTau_Xi", "c#tau vs. mass;m_{#Lambda#pi};c#tau (cm)", {HistType::kTH2F, {{200, 1.2, 1.4}, {200, 0, 20.f}}}},
      {"Cascade/hCTau_Omega", "c#tau vs. mass;m_{#LambdaK};c#tau (cm)", {HistType::kTH2F, {{200, 1.6, 1.8}, {200, 0, 20.f}}}},
      {"Cascade/hV0CosPA", "V0 cosine of pointing angle", {HistType::kTH1F, {{50, 0.95, 1.f}}}},
      {"Cascade/hV0PCA", "V0 distance between 2 legs", {HistType::kTH1F, {{200, 0.f, 2.f}}}},
      {"Cascade/hCosPA", "cascade cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"Cascade/hPCA", "cascade distance between 2 legs", {HistType::kTH1F, {{200, 0.f, 2.f}}}},
      {"Cascade/hMassLambda", "V0 mass Lambda in cascade", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"Cascade/hMassAntiLambda", "V0 mass AntiLambda in cascade", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"Cascade/hMassXi", "cascade mass #Xi", {HistType::kTH1F, {{200, 1.2, 1.4}}}},
      {"Cascade/hMassOmega", "cascade mass #Omega", {HistType::kTH1F, {{200, 1.6, 1.8}}}},
      {"Cascade/hMassPt_Xi", "cascade mass #Xi^{#pm};m_{#Lambda#pi} (GeV/c^{2});p_{T,#Lambda#pi} (GeV/c)", {HistType::kTH2F, {{200, 1.2, 1.4}, {100, 0, 10}}}},
      {"Cascade/hMassPt_Xi_bachelor", "cascade mass #Xi^{#pm};m_{#Lambda#pi} (GeV/c^{2});p_{T,#pi} (GeV/c)", {HistType::kTH2F, {{200, 1.2, 1.4}, {100, 0, 10}}}},
      {"Cascade/hMassPt_Omega", "cascade mass #Omega^{#pm};m_{#LambdaK} (GeV/c^{2});p_{T,#LambdaK} (GeV/c)", {HistType::kTH2F, {{200, 1.6, 1.8}, {100, 0, 10}}}},
      {"Cascade/hMassPt_Omega_bachelor", "cascade mass #Omega^{#pm};m_{#LambdaK} (GeV/c^{2});p_{T,K} (GeV/c)", {HistType::kTH2F, {{200, 1.6, 1.8}, {100, 0, 10}}}},
    },
  };

  // Configurables

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz_input", -999, "bz field, -999 is automatic"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  Configurable<float> downscaling_electron{"downscaling_electron", 0.005, "down scaling factor to store electron"};
  Configurable<float> downscaling_pion{"downscaling_pion", 0.001, "down scaling factor to store pion"};
  Configurable<float> downscaling_kaon{"downscaling_kaon", 1.1, "down scaling factor to store kaon"};
  Configurable<float> downscaling_proton{"downscaling_proton", 0.01, "down scaling factor to store proton"};
  Configurable<bool> store_v0photons{"store_v0photons", true, "create training data from v0 photons"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0380, "intercept for m vs. phiv"};

  struct : ConfigurableGroup {
    std::string prefix = "trackcut_group";
    Configurable<float> cfg_min_pt{"cfg_min_pt", 0.05, "min pt for v0 legs"};
    Configurable<float> cfg_max_eta{"cfg_max_eta", 0.9, "max. eta for v0 legs"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 5.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 6.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.3, "max dca XY in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.3, "max dca Z in cm"};
  } trackcuts;

  struct : ConfigurableGroup {
    std::string prefix = "v0cut_group";
    Configurable<float> cfg_min_pt{"cfg_min_pt", 0.05, "min pt for v0 legs"};
    Configurable<float> cfg_max_eta{"cfg_max_eta", 0.9, "max. eta for v0 legs"};
    Configurable<float> cfg_min_mass_photon{"cfg_min_mass_photon", 0.00, "min mass for photon conversion"};
    Configurable<float> cfg_max_mass_photon{"cfg_max_mass_photon", 0.02, "max mass for photon conversion"};
    Configurable<float> cfg_min_mass_k0s{"cfg_min_mass_k0s", 0.490, "min mass for K0S"};
    Configurable<float> cfg_max_mass_k0s{"cfg_max_mass_k0s", 0.505, "max mass for K0S"};
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for Lambda rejection"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for Lambda rejection"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min cospa for v0"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.5, "max distance between 2 legs for v0"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 2, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 0, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 5.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 6.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};

    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -4, "min n sigma e in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +4, "max n sigma e in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -4, "min n sigma pi in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +4, "max n sigma pi in TPC"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -4, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +4, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -4, "min n sigma pr in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +4, "max n sigma pr in TPC"};

    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -4, "min n sigma e in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +4, "max n sigma e in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -4, "min n sigma pi in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +4, "max n sigma pi in TOF"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -4, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +4, "max n sigma ka in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -4, "min n sigma pr in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +4, "max n sigma pr in TOF"};

    Configurable<float> cfg_min_TPCNsigmaPi_tight{"cfg_min_TPCNsigmaPi_tight", -2, "min n sigma pi in TPC for Lambda and cascade"};
    Configurable<float> cfg_max_TPCNsigmaPi_tight{"cfg_max_TPCNsigmaPi_tight", +2, "max n sigma pi in TPC for Lambda and cascade"};
    Configurable<float> cfg_min_TPCNsigmaPr_tight{"cfg_min_TPCNsigmaPr_tight", -2, "min n sigma pr in TPC for cascade"};
    Configurable<float> cfg_max_TPCNsigmaPr_tight{"cfg_max_TPCNsigmaPr_tight", +2, "max n sigma pr in TPC for cascade"};

    Configurable<float> cfg_min_TOFNsigmaPi_tight{"cfg_min_TOFNsigmaPi_tight", -2, "min n sigma pi in TOF for Lambda and cascade"};
    Configurable<float> cfg_max_TOFNsigmaPi_tight{"cfg_max_TOFNsigmaPi_tight", +2, "max n sigma pi in TOF for Lambda and cascade"};
    Configurable<float> cfg_min_TOFNsigmaPr_tight{"cfg_min_TOFNsigmaPr_tight", -2, "min n sigma pr in TOF for cascade"};
    Configurable<float> cfg_max_TOFNsigmaPr_tight{"cfg_max_TOFNsigmaPr_tight", +2, "max n sigma pr in TOF for cascade"};
  } v0cuts;

  struct : ConfigurableGroup {
    std::string prefix = "cascadecut_group";
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for lambda in cascade"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for lambda in cascade"};
    Configurable<float> cfg_min_mass_Xi{"cfg_min_mass_Xi", 1.317, "min mass for Xi"};
    Configurable<float> cfg_max_mass_Xi{"cfg_max_mass_Xi", 1.327, "max mass for Xi"};
    Configurable<float> cfg_min_mass_Omega{"cfg_min_mass_Omega", 1.667, "min mass for Omega"};
    Configurable<float> cfg_max_mass_Omega{"cfg_max_mass_Omega", 1.677, "max mass for Omega"};
    Configurable<float> cfg_min_cospa_v0{"cfg_min_cospa_v0", 0.97, "minimum V0 CosPA in cascade"};
    Configurable<float> cfg_max_dcadau_v0{"cfg_max_dcadau_v0", 0.2, "max distance between V0 Daughters in cascade"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "minimum cascade CosPA"};
    Configurable<float> cfg_max_dcadau{"cfg_max_dcadau", 0.4, "max distance between bachelor and V0"};
    Configurable<float> cfg_min_rxy_v0{"cfg_min_rxy_v0", 1.2, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_rxy{"cfg_min_rxy", 0.5, "minimum V0 rxy in cascade"};
    Configurable<float> cfg_min_dcaxy_v0leg{"cfg_min_dcaxy_v0leg", 0.1, "min dca XY for v0 legs in cm"};
    Configurable<float> cfg_min_dcaxy_bachelor{"cfg_min_dcaxy_bachelor", 0.1, "min dca XY for bachelor in cm"};
  } cascadecuts;

  int mRunNumber;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::dataformats::DCA mDcaInfoCov;

  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
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
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
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

  template <int il0, int il1, typename TTrack>
  float meanClusterSizeITS(TTrack const& track)
  {
    int total_size = 0;
    int nl = 0;
    for (int il = il0; il < il1; il++) {
      if (track.itsClsSizeInLayer(il) > 0) {
        total_size += track.itsClsSizeInLayer(il);
        nl++;
      }
    }
    if (nl > 0) {
      return static_cast<float>(total_size) / static_cast<float>(nl);
    } else {
      return 0.f;
    }
  }

  template <typename TCollision, typename TTrack>
  bool isSelectedTrack(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < trackcuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < trackcuts.cfg_min_ncluster_itsib) {
      return false;
    }
    if (track.itsChi2NCl() > trackcuts.cfg_max_chi2its) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < trackcuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }
    if (track.tpcNClsFound() < trackcuts.cfg_min_ncluster_tpc) {
      return false;
    }
    if (track.tpcChi2NCl() > trackcuts.cfg_max_chi2tpc) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < trackcuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }
    if (track.tpcFractionSharedCls() > trackcuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(track_par_cov_recalc.getEta()) > trackcuts.cfg_max_eta || track_par_cov_recalc.getPt() < trackcuts.cfg_min_pt) {
      return false;
    }

    if (std::fabs(dcaXY) > trackcuts.cfg_max_dcaxy) {
      return false;
    }
    if (std::fabs(dcaZ) > trackcuts.cfg_max_dcaz) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TTrack>
  bool isSelectedV0Leg(TCollision const& collision, TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < v0cuts.cfg_min_ncluster_its) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < v0cuts.cfg_min_ncluster_itsib) {
      return false;
    }
    if (track.itsChi2NCl() > v0cuts.cfg_max_chi2its) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < v0cuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }
    if (track.tpcNClsFound() < v0cuts.cfg_min_ncluster_tpc) {
      return false;
    }
    if (track.tpcChi2NCl() > v0cuts.cfg_max_chi2tpc) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < v0cuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }
    if (track.tpcFractionSharedCls() > v0cuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto track_par_cov_recalc = getTrackParCov(track);
    track_par_cov_recalc.setPID(track.pidForTracking());
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    // float dcaZ = mDcaInfoCov.getZ();

    if (std::fabs(track_par_cov_recalc.getEta()) > v0cuts.cfg_max_eta || track_par_cov_recalc.getPt() < v0cuts.cfg_min_pt) {
      return false;
    }

    if (std::fabs(dcaXY) < v0cuts.cfg_min_dcaxy_v0leg) { // this is applied in filter.
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isElectron(TTrack const& track)
  {
    bool is_El_TPC = v0cuts.cfg_min_TPCNsigmaEl < track.tpcNSigmaEl() && track.tpcNSigmaEl() < v0cuts.cfg_max_TPCNsigmaEl;
    bool is_El_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaEl < track.tofNSigmaEl() && track.tofNSigmaEl() < v0cuts.cfg_max_TOFNsigmaEl : true; // TOFif
    return is_El_TPC && is_El_TOF;
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    bool is_Pi_TPC = v0cuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0cuts.cfg_max_TPCNsigmaPi;
    bool is_Pi_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < v0cuts.cfg_max_TOFNsigmaPi : true; // TOFif
    return is_Pi_TPC && is_Pi_TOF;
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    bool is_Ka_TPC = v0cuts.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < v0cuts.cfg_max_TPCNsigmaKa;
    bool is_Ka_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < v0cuts.cfg_max_TOFNsigmaKa : true; // TOFif
    return is_Ka_TPC && is_Ka_TOF;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    bool is_Pr_TPC = v0cuts.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0cuts.cfg_max_TPCNsigmaPr;
    bool is_Pr_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < v0cuts.cfg_max_TOFNsigmaPr : true; // TOFif
    return is_Pr_TPC && is_Pr_TOF;
  }

  template <typename TTrack>
  bool isPionTight(TTrack const& track)
  {
    bool is_Pi_TPC = v0cuts.cfg_min_TPCNsigmaPi_tight < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0cuts.cfg_max_TPCNsigmaPi_tight;
    bool is_Pi_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPi_tight < track.tofNSigmaPi() && track.tofNSigmaPi() < v0cuts.cfg_max_TOFNsigmaPi_tight : true; // TOFif
    return is_Pi_TPC && is_Pi_TOF;
  }

  template <typename TTrack>
  bool isProtonTight(TTrack const& track)
  {
    bool is_Pr_TPC = v0cuts.cfg_min_TPCNsigmaPr_tight < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0cuts.cfg_max_TPCNsigmaPr_tight;
    bool is_Pr_TOF = track.hasTOF() ? v0cuts.cfg_min_TOFNsigmaPr_tight < track.tofNSigmaPr() && track.tofNSigmaPr() < v0cuts.cfg_max_TOFNsigmaPr_tight : true; // TOFif
    return is_Pr_TPC && is_Pr_TOF;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track, const int pidlabel, const int tracktype)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), track.globalIndex()) == stored_trackIds.end()) {
      mDcaInfoCov.set(999, 999, 999, 999, 999);
      auto track_par_cov_recalc = getTrackParCov(track);
      track_par_cov_recalc.setPID(track.pidForTracking());
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
      o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
      float dcaXY = mDcaInfoCov.getY();
      float dcaZ = mDcaInfoCov.getZ();

      emprimarytracks(collision.globalIndex(), collision.posZ(), collision.numContrib(),
                      track.pt(), track.eta(), track.phi(), track.tgl(), track.signed1Pt(),
                      dcaXY, dcaZ, track_par_cov_recalc.getSigmaY2(), track_par_cov_recalc.getSigmaZ2(), track_par_cov_recalc.getSigmaZY(),
                      track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                      track.tpcChi2NCl(), track.tpcInnerParam(),
                      track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                      track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                      track.itsClusterSizes(), track.itsChi2NCl(), track.tofChi2(), track.detectorMap(), pidlabel, tracktype);
      stored_trackIds.emplace_back(track.globalIndex());
    }
  }

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > v0cuts.cfg_min_cospa.value&& o2::aod::v0data::dcaV0daughters<v0cuts.cfg_max_dcadau.value && nabs(o2::aod::v0data::dcanegtopv)> v0cuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::v0data::dcanegtopv) > v0cuts.cfg_min_dcaxy_v0leg;
  using filteredV0s = soa::Filtered<aod::V0Datas>;

  Filter cascadeFilter = o2::aod::cascdata::dcacascdaughters < cascadecuts.cfg_max_dcadau.value && nabs(o2::aod::cascdata::dcanegtopv) > cascadecuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcanegtopv) > cascadecuts.cfg_min_dcaxy_v0leg&& nabs(o2::aod::cascdata::dcabachtopv) > cascadecuts.cfg_min_dcaxy_bachelor;
  using filteredCascades = soa::Filtered<aod::CascDatas>;

  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::evsel::sel8 == true;
  using filteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<aod::V0Datas> perCollision_v0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDatas> perCollision_cascade = o2::aod::cascdata::collisionId;
  Preslice<MyTracks> perCollision_track = o2::aod::track::collisionId;

  Partition<MyTracks> posTracks = o2::aod::track::signed1Pt > 0.f && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  Partition<MyTracks> negTracks = o2::aod::track::signed1Pt < 0.f && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  std::vector<uint64_t> stored_trackIds;

  void processPID(filteredMyCollisions const& collisions, aod::BCsWithTimestamps const&, filteredV0s const& v0s, filteredCascades const& cascades, MyTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());
    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1.0); // all

      auto bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      auto v0s_coll = v0s.sliceBy(perCollision_v0, collision.globalIndex());
      for (const auto& v0 : v0s_coll) {
        auto pos = v0.template posTrack_as<MyTracks>();
        auto neg = v0.template negTrack_as<MyTracks>();
        if (!isSelectedV0Leg(collision, pos) || !isSelectedV0Leg(collision, neg)) {
          continue;
        }
        if (pos.sign() * neg.sign() > 0) {
          continue;
        }

        registry.fill(HIST("V0/hPCA"), v0.dcaV0daughters());
        registry.fill(HIST("V0/hCosPA"), v0.v0cosPA());
        registry.fill(HIST("V0/hAP"), v0.alpha(), v0.qtarm());

        if (v0.dcaV0daughters() > v0cuts.cfg_max_dcadau) {
          continue;
        }
        if (v0.v0cosPA() < v0cuts.cfg_min_cospa) {
          continue;
        }

        if (isPion(pos) && isPion(neg)) {
          registry.fill(HIST("V0/hMassK0Short"), v0.mK0Short());
          if (v0cuts.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < v0cuts.cfg_max_mass_k0s) {
            registry.fill(HIST("V0/hTPCdEdx_P_Pi"), neg.p(), neg.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Pi"), neg.p(), neg.beta());
            registry.fill(HIST("V0/hTPCdEdx_P_Pi"), pos.p(), pos.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Pi"), pos.p(), pos.beta());
            if (dist01(engine) < downscaling_pion) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
            if (dist01(engine) < downscaling_pion) {
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
          }
        }
        if (isProton(pos) && isPionTight(neg)) {
          registry.fill(HIST("V0/hMassLambda"), v0.mLambda());
          if (v0cuts.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < v0cuts.cfg_max_mass_lambda) {
            if (dist01(engine) < downscaling_proton) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
            registry.fill(HIST("V0/hTPCdEdx_P_Pr"), pos.p(), pos.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Pr"), pos.p(), pos.beta());
          }
        }
        if (isPionTight(pos) && isProton(neg)) {
          registry.fill(HIST("V0/hMassAntiLambda"), v0.mAntiLambda());
          if (v0cuts.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < v0cuts.cfg_max_mass_lambda) {
            if (dist01(engine) < downscaling_proton) {
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
            registry.fill(HIST("V0/hTPCdEdx_P_Pr"), neg.p(), neg.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Pr"), neg.p(), neg.beta());
          }
        }
        if (isElectron(pos) && isElectron(neg) && store_v0photons) {
          registry.fill(HIST("V0/hMassGamma"), v0.mGamma());
          registry.fill(HIST("V0/hXY_Gamma"), v0.x(), v0.y());
          registry.fill(HIST("V0/hMassGamma_Rxy"), v0.v0radius(), v0.mGamma());
          if ((v0cuts.cfg_min_mass_photon < v0.mGamma() && v0.mGamma() < v0cuts.cfg_max_mass_photon)) {
            if (dist01(engine) < downscaling_electron) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
            if (dist01(engine) < downscaling_electron) {
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kSecondary));
            }
            registry.fill(HIST("V0/hTPCdEdx_P_El"), neg.p(), neg.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), neg.p(), neg.beta());
            registry.fill(HIST("V0/hTPCdEdx_P_El"), pos.p(), pos.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), pos.p(), pos.beta());
          }
        }
      } // end of v0 loop

      if (!store_v0photons) {
        auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
        auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
        for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
          if (!isSelectedTrack(collision, pos) || !isSelectedTrack(collision, neg)) {
            continue;
          }

          if (!isElectron(pos) || !isElectron(neg)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v1(neg.pt(), neg.eta(), neg.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
          float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(pos.px(), pos.py(), pos.pz(), neg.px(), neg.py(), neg.pz(), pos.sign(), neg.sign(), d_bz);
          registry.fill(HIST("hMvsPhiV"), phiv, v12.M());

          if (v12.M() < slope * phiv + intercept) { // photon conversion is found.
            registry.fill(HIST("V0/hTPCdEdx_P_El"), neg.p(), neg.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), neg.p(), neg.beta());
            registry.fill(HIST("V0/hTPCdEdx_P_El"), pos.p(), pos.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_El"), pos.p(), pos.beta());

            if (dist01(engine) < downscaling_electron) {
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kPrimary)); // secondary in primary electron candidates
            }
            if (dist01(engine) < downscaling_electron) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kPrimary)); // secondary in primary electron candidates
            }
          }
        } // end of ULS pair loop
      }

      auto cascades_coll = cascades.sliceBy(perCollision_cascade, collision.globalIndex());
      for (const auto& cascade : cascades_coll) {
        // Track casting
        auto bachelor = cascade.template bachelor_as<MyTracks>();
        auto pos = cascade.template posTrack_as<MyTracks>();
        auto neg = cascade.template negTrack_as<MyTracks>();
        if (!isSelectedV0Leg(collision, pos) || !isSelectedV0Leg(collision, neg) || !isSelectedV0Leg(collision, bachelor)) {
          continue;
        }

        if (pos.sign() * neg.sign() > 0) {
          continue;
        }

        if (bachelor.sign() < 0) { // Omega- -> L + K- -> p + pi- + K-
          if (!isProtonTight(pos) || !isPionTight(neg)) {
            continue;
          }
        } else { // Omegabar+ -> Lbar + K+ -> pbar + pi+ + K+
          if (!isProtonTight(neg) || !isPionTight(pos)) {
            continue;
          }
        }

        registry.fill(HIST("Cascade/hMassLambda"), cascade.mLambda());
        if (!(v0cuts.cfg_min_mass_lambda < cascade.mLambda() && cascade.mLambda() < v0cuts.cfg_max_mass_lambda)) {
          continue;
        }

        if (cascade.cascradius() > cascade.v0radius()) {
          continue;
        }

        registry.fill(HIST("Cascade/hV0PCA"), cascade.dcaV0daughters());
        registry.fill(HIST("Cascade/hV0CosPA"), cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        registry.fill(HIST("Cascade/hPCA"), cascade.dcacascdaughters()); // distance between bachelor and V0.
        registry.fill(HIST("Cascade/hCosPA"), cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()));

        if (cascade.dcaV0daughters() > cascadecuts.cfg_max_dcadau_v0) {
          continue;
        }
        if (cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecuts.cfg_min_cospa_v0) {
          continue;
        }
        if (cascade.v0radius() < cascadecuts.cfg_min_rxy_v0) {
          continue;
        }

        if (cascade.cascradius() < cascadecuts.cfg_min_rxy) {
          continue;
        }

        if (cascade.dcacascdaughters() > cascadecuts.cfg_max_dcadau) {
          continue;
        }
        if (cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascadecuts.cfg_min_cospa) {
          continue;
        }

        float length = std::sqrt(std::pow(cascade.x() - collision.posX(), 2) + std::pow(cascade.y() - collision.posY(), 2) + std::pow(cascade.z() - collision.posZ(), 2));
        float mom = cascade.p();
        float ctauXi = length / mom * o2::constants::physics::MassXiMinus;       // 4.91 cm in PDG
        float ctauOmega = length / mom * o2::constants::physics::MassOmegaMinus; // 2.46 cm in PDG

        if (isPion(bachelor)) {
          registry.fill(HIST("Cascade/hMassXi"), cascade.mXi());
          registry.fill(HIST("Cascade/hMassPt_Xi"), cascade.mXi(), cascade.pt());
          registry.fill(HIST("Cascade/hMassPt_Xi_bachelor"), cascade.mXi(), bachelor.p());
          registry.fill(HIST("Cascade/hRxy_Xi"), cascade.mXi(), cascade.cascradius());
          registry.fill(HIST("Cascade/hCTau_Xi"), cascade.mXi(), ctauXi);
        }
        if (!(cascadecuts.cfg_min_mass_Xi < cascade.mXi() && cascade.mXi() < cascadecuts.cfg_max_mass_Xi) && isKaon(bachelor)) { // reject Xi candidates
          registry.fill(HIST("Cascade/hMassOmega"), cascade.mOmega());
          registry.fill(HIST("Cascade/hMassPt_Omega"), cascade.mOmega(), cascade.pt());
          registry.fill(HIST("Cascade/hMassPt_Omega_bachelor"), cascade.mOmega(), bachelor.p());
          registry.fill(HIST("Cascade/hRxy_Omega"), cascade.mOmega(), cascade.cascradius());
          registry.fill(HIST("Cascade/hCTau_Omega"), cascade.mOmega(), ctauOmega);
          if (cascadecuts.cfg_min_mass_Omega < cascade.mOmega() && cascade.mOmega() < cascadecuts.cfg_max_mass_Omega) { // select Omega candidates
            registry.fill(HIST("V0/hTPCdEdx_P_Ka"), bachelor.p(), bachelor.tpcSignal());
            registry.fill(HIST("V0/hTOFbeta_P_Ka"), bachelor.p(), bachelor.beta());
            if (dist01(engine) < downscaling_kaon) {
              fillTrackTable(collision, bachelor, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kKaon), static_cast<int>(o2::aod::pwgem::dilepton::Track_Type::kPrimary));
            }
          }
        }
      } // end of cascade loop

      // auto tracks_coll = tracks.sliceBy(perCollision_track, collision.globalIndex());
      // for (const auto& track : tracks_coll) {
      //   if (!isSelectedTrack(collision, track)) {
      //     continue;
      //   }
      //   registry.fill(HIST("PrimaryTrack/hTPCdEdx_P"), track.p(), track.tpcSignal());
      //   registry.fill(HIST("PrimaryTrack/hTOFbeta_P"), track.p(), track.beta());
      //   registry.fill(HIST("PrimaryTrack/hITSClusterSize_P"), track.p(), meanClusterSizeITS<0, 7>(track) * std::cos(std::atan(track.tgl())));
      //   registry.fill(HIST("PrimaryTrack/hTPCNsigmaEl_P"), track.p(), track.tpcNSigmaEl());
      //   registry.fill(HIST("PrimaryTrack/hTOFNsigmaEl_P"), track.p(), track.tofNSigmaEl());
      //   registry.fill(HIST("PrimaryTrack/hTPCNsigmaMu_P"), track.p(), track.tpcNSigmaMu());
      //   registry.fill(HIST("PrimaryTrack/hTOFNsigmaMu_P"), track.p(), track.tofNSigmaMu());
      //   registry.fill(HIST("PrimaryTrack/hTPCNsigmaPi_P"), track.p(), track.tpcNSigmaPi());
      //   registry.fill(HIST("PrimaryTrack/hTOFNsigmaPi_P"), track.p(), track.tofNSigmaPi());
      //   registry.fill(HIST("PrimaryTrack/hTPCNsigmaKa_P"), track.p(), track.tpcNSigmaKa());
      //   registry.fill(HIST("PrimaryTrack/hTOFNsigmaKa_P"), track.p(), track.tofNSigmaKa());
      //   registry.fill(HIST("PrimaryTrack/hTPCNsigmaPr_P"), track.p(), track.tpcNSigmaPr());
      //   registry.fill(HIST("PrimaryTrack/hTOFNsigmaPr_P"), track.p(), track.tofNSigmaPr());
      // } // end of track loop

    } // end of collision loop
    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  } // end of process

  // please choose only 1 process function.
  void processDummy(filteredMyCollisions const&) {}

  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processPID, "produce ML input for single track level", false); // this is for eID with ITSsa. e/pi/k/p are selected by TOF, and these can be used for ITS-TPC PID.
  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processDummy, "process dummy", true);
};

struct MLTrackQC {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hTPCdEdx_P_All", "TPC dE/dx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0.f, 5.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_All", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0.f, 5.f}, {220, 0.0, 1.1}}}},
      {"hITSClusterSize_P_All", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {64, 0.0, 16}}}},
      {"hTPCdEdx_P_Electron", "TPC dE/dx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0.f, 5.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Electron", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0.f, 5.f}, {220, 0.0, 1.1}}}},
      {"hITSClusterSize_P_Electron", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {64, 0.0, 16}}}},
      {"hTPCdEdx_P_Pion", "TPC dE/dx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0.f, 5.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Pion", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0.f, 5.f}, {220, 0.0, 1.1}}}},
      {"hITSClusterSize_P_Pion", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {64, 0.0, 16}}}},
      {"hTPCdEdx_P_Kaon", "TPC dE/dx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0.f, 5.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Kaon", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0.f, 5.f}, {220, 0.0, 1.1}}}},
      {"hITSClusterSize_P_Kaon", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {64, 0.0, 16}}}},
      {"hTPCdEdx_P_Proton", "TPC dE/dx vs. p;p^{ITS-TPC} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{500, 0.f, 5.f}, {200, 0, 200}}}},
      {"hTOFbeta_P_Proton", "TOF beta vs. p;p^{ITS-TPC} (GeV/c);TOF #beta", {HistType::kTH2F, {{500, 0.f, 5.f}, {220, 0.0, 1.1}}}},
      {"hITSClusterSize_P_Proton", "mean ITS cluster size vs. p;p^{ITS-TPC} (GeV/c);<ITS cluster size> #times cos(#lambda)", {HistType::kTH2F, {{500, 0.f, 5.f}, {64, 0.0, 16}}}},

      {"hTPCNsigmaEl_P", "TPC n#sigma_{e} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTPCNsigmaPi_P", "TPC n#sigma_{#pi} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTPCNsigmaKa_P", "TPC n#sigma_{K} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{K}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTPCNsigmaPr_P", "TPC n#sigma_{p} vs. p;p^{ITS-TPC} (GeV/c);n #sigma_{p}^{TPC}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTOFNsigmaEl_P", "TOF n#sigma_{e} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTOFNsigmaPi_P", "TOF n#sigma_{#pi} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTOFNsigmaKa_P", "TOF n#sigma_{K} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{K}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
      {"hTOFNsigmaPr_P", "TOF n#sigma_{p} vs. p;p^{ITS-TOF} (GeV/c);n #sigma_{p}^{TOF}", {HistType::kTH2F, {{500, 0.f, 5.f}, {100, -5, +5}}}},
    },
  };

  void processQC(aod::EMPrimaryTracks const& tracks)
  {
    for (const auto& track : tracks) {
      registry.fill(HIST("hTPCdEdx_P_All"), track.p(), track.tpcSignal());
      registry.fill(HIST("hTOFbeta_P_All"), track.p(), track.beta());
      registry.fill(HIST("hITSClusterSize_P_All"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
      if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron)) {
        registry.fill(HIST("hTPCdEdx_P_Electron"), track.p(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Electron"), track.p(), track.beta());
        registry.fill(HIST("hITSClusterSize_P_Electron"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaEl_P"), track.p(), track.tpcNSigmaEl());
        registry.fill(HIST("hTOFNsigmaEl_P"), track.p(), track.tofNSigmaEl());
      } else if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion)) {
        registry.fill(HIST("hTPCdEdx_P_Pion"), track.p(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Pion"), track.p(), track.beta());
        registry.fill(HIST("hITSClusterSize_P_Pion"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaPi_P"), track.p(), track.tpcNSigmaPi());
        registry.fill(HIST("hTOFNsigmaPi_P"), track.p(), track.tofNSigmaPi());
      } else if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kKaon)) {
        registry.fill(HIST("hTPCdEdx_P_Kaon"), track.p(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Kaon"), track.p(), track.beta());
        registry.fill(HIST("hITSClusterSize_P_Kaon"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaKa_P"), track.p(), track.tpcNSigmaKa());
        registry.fill(HIST("hTOFNsigmaKa_P"), track.p(), track.tofNSigmaKa());
      } else if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton)) {
        registry.fill(HIST("hTPCdEdx_P_Proton"), track.p(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_P_Proton"), track.p(), track.beta());
        registry.fill(HIST("hITSClusterSize_P_Proton"), track.p(), track.meanClusterSizeITS() * std::cos(std::atan(track.tgl())));
        registry.fill(HIST("hTPCNsigmaPr_P"), track.p(), track.tpcNSigmaPr());
        registry.fill(HIST("hTOFNsigmaPr_P"), track.p(), track.tofNSigmaPr());
      }
    } // end of track loop
  }
  PROCESS_SWITCH(MLTrackQC, processQC, "process QC for single track level", false);

  void processDummy(aod::EMPrimaryTracks const&) {}
  PROCESS_SWITCH(MLTrackQC, processDummy, "process dummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronMLDDA>(cfgc, TaskName{"tree-creator-ele-ml-dda"}),
    adaptAnalysisTask<MLTrackQC>(cfgc, TaskName{"ml-track-qc"})};
}
