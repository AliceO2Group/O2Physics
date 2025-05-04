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

#include <vector>
#include <string>
#include <array>

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/RecoDecay.h"
#include "DCAFitter/DCAFitterN.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct taggingHFE {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;

  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                             aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                             aod::McTrackLabels>;

  using MyV0s = soa::Join<aod::V0Datas, aod::V0Covs>;

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
    std::string prefix = "trackcut_group";
    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.05, "min pT for single track"};
    Configurable<float> cfg_max_pt_track{"cfg_max_pt_track", 1e+10, "max pT for single track"};
    Configurable<float> cfg_min_eta_track{"cfg_min_eta_track", -0.9, "min eta for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", +0.9, "max eta for single track"};
    Configurable<float> cfg_min_eta_track_wide{"cfg_min_eta_track_wide", -1.2, "min eta for single track (tagger)"};
    Configurable<float> cfg_max_eta_track_wide{"cfg_max_eta_track_wide", +1.2, "max eta for single track (tagger)"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 0.7, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 70, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2/NclsTOF"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.3, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.3, "max dca Z for single track in cm"};

    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3, "min n sigma e in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3, "max n sigma e in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3, "min n sigma pi in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3, "max n sigma pi in TPC"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3, "min n sigma pr in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3, "max n sigma pr in TPC"};

    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3, "min n sigma el in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3, "max n sigma el in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -3, "min n sigma pi in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +3, "max n sigma pi in TOF"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -3, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +3, "max n sigma ka in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -3, "min n sigma pr in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +3, "max n sigma pr in TOF"};
  } trackcuts;

  struct : ConfigurableGroup {
    std::string prefix = "svcut_group";
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.8, "min cospa"};
    Configurable<float> cfg_min_cospaXY{"cfg_min_cospaXY", 0.0, "min cospaXY"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 1.0, "max distance between 2 legs"};
    Configurable<float> cfg_min_lxy{"cfg_min_lxy", -1, "min lxy for charm hadron candidate"};
    Configurable<float> cfg_max_mass_eK{"cfg_max_mass_eK", 2.0, "max mass for eK pair"};
    Configurable<float> cfg_max_mass_eL{"cfg_max_mass_eL", 2.3, "max mass for eL pair"};
  } svcuts;

  struct : ConfigurableGroup {
    std::string prefix = "v0cut_group";
    Configurable<float> cfg_min_mass_k0s{"cfg_min_mass_k0s", 0.49, "min mass for K0S"};
    Configurable<float> cfg_max_mass_k0s{"cfg_max_mass_k0s", 0.51, "max mass for K0S"};
    Configurable<float> cfg_min_mass_lambda{"cfg_min_mass_lambda", 1.11, "min mass for Lambda rejection"};
    Configurable<float> cfg_max_mass_lambda{"cfg_max_mass_lambda", 1.12, "max mass for Lambda rejection"};
    Configurable<float> cfg_min_cospa_v0hadron{"cfg_min_cospa_v0hadron", 0.95, "min cospa for v0hadron"};
    Configurable<float> cfg_max_pca_v0hadron{"cfg_max_pca_v0hadron", 0.5, "max distance between 2 legs for v0hadron"};
    // Configurable<float> cfg_min_radius_v0hadron{"cfg_min_radius_v0hadron", 0.1, "min rxy for v0hadron"};
    Configurable<float> cfg_min_cr2findable_ratio_tpc{"cfg_min_cr2findable_ratio_tpc", 0.8, "min. TPC Ncr/Nf ratio"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 40, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 4, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 1, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tof{"cfg_max_chi2tof", 1e+10, "max chi2 for TOF"};
    Configurable<float> cfg_min_dcaxy{"cfg_min_dcaxy", -1, "min dca XY for v0 legs in cm"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3, "min n sigma e in TPC"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3, "max n sigma e in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -3, "min n sigma pi in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3, "max n sigma pi in TPC"};
    Configurable<float> cfg_min_TPCNsigmaKa{"cfg_min_TPCNsigmaKa", -3, "min n sigma ka in TPC"};
    Configurable<float> cfg_max_TPCNsigmaKa{"cfg_max_TPCNsigmaKa", +3, "max n sigma ka in TPC"};
    Configurable<float> cfg_min_TPCNsigmaPr{"cfg_min_TPCNsigmaPr", -3, "min n sigma pr in TPC"};
    Configurable<float> cfg_max_TPCNsigmaPr{"cfg_max_TPCNsigmaPr", +3, "max n sigma pr in TPC"};

    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3, "min n sigma el in TOF"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3, "max n sigma el in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPi{"cfg_min_TOFNsigmaPi", -3, "min n sigma pi in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPi{"cfg_max_TOFNsigmaPi", +3, "max n sigma pi in TOF"};
    Configurable<float> cfg_min_TOFNsigmaKa{"cfg_min_TOFNsigmaKa", -3, "min n sigma ka in TOF"};
    Configurable<float> cfg_max_TOFNsigmaKa{"cfg_max_TOFNsigmaKa", +3, "max n sigma ka in TOF"};
    Configurable<float> cfg_min_TOFNsigmaPr{"cfg_min_TOFNsigmaPr", -3, "min n sigma pr in TOF"};
    Configurable<float> cfg_max_TOFNsigmaPr{"cfg_max_TOFNsigmaPr", +3, "max n sigma pr in TOF"};
  } v0cuts;

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};

  HistogramRegistry fRegistry{"fRegistry"};
  static constexpr std::string_view hadron_names[7] = {"LF/", "promptJpsi/", "nonpromptJpsi/", "D0/", "Dpm/", "Ds/", "Lc/"};
  static constexpr std::string_view pair_names[3] = {"e_Kpm/", "e_K0S/", "e_Lambda/"};
  static constexpr std::string_view hTypes[4] = {"findable/", "correct/", "fake/", "miss/"};

  void init(o2::framework::InitContext&)
  {
    if (doprocessSA && doprocessTTCA) {
      LOGF(fatal, "Cannot enable doprocessWithoutFTTCA and doprocessWithFTTCA at the same time. Please choose one.");
    }

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
  o2::dataformats::VertexBase mVtx;
  const o2::dataformats::MeanVertexObject* mMeanVtx = nullptr;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<2> fitter;

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

    // for charm hadrons
    fRegistry.add("e_Kpm/all/hLxy", "decay length XY from PV;L_{xy} (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("e_Kpm/all/hLz", "decay length Z from PV;L_{z} (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("e_Kpm/all/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("e_Kpm/all/hCosPAXY", "cosPA in XY;cosine of pointing angle in XY", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("e_Kpm/all/hDCA2Legs", "distance between 2 legs;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.add("e_Kpm/all/hMass", "mass;mass (GeV/c^{2});p_{T} (GeV/c)", kTH2F, {{40, 0.5, 2.5}, {100, 0, 10}}, false);

    fRegistry.addClone("e_Kpm/all/", "e_Kpm/D0/");
    fRegistry.addClone("e_Kpm/all/", "e_Kpm/Dpm/");
    fRegistry.addClone("e_Kpm/all/", "e_Kpm/Ds/");
    fRegistry.addClone("e_Kpm/all/", "e_Kpm/fake/");

    fRegistry.addClone("e_Kpm/all/", "e_K0S/all/");
    fRegistry.addClone("e_Kpm/all/", "e_K0S/D0/");
    fRegistry.addClone("e_Kpm/all/", "e_K0S/Dpm/");
    fRegistry.addClone("e_Kpm/all/", "e_K0S/Ds/");
    fRegistry.addClone("e_Kpm/all/", "e_K0S/fake/");

    fRegistry.addClone("e_Kpm/all/", "e_Lambda/all/");
    fRegistry.addClone("e_Kpm/all/", "e_Lambda/Lc/");
    fRegistry.addClone("e_Kpm/all/", "e_Lambda/fake/");

    // for V0s
    fRegistry.add("V0/K0S/hPt", "pT of V0;p_{T} (GeV/c)", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("V0/K0S/hYPhi", "Y vs. #varphi of V0;#varphi (rad.);rapidity", kTH2F, {{36, 0, 2 * M_PI}, {80, -2, +2}}, false);
    fRegistry.add("V0/K0S/hAP", "Ap plot;#alpha;q_{T} (GeV/c)", kTH2F, {{200, -1, 1}, {250, 0, 0.25}}, false);
    fRegistry.add("V0/K0S/hLxy", "decay length from PV;L_{xy} (cm)", kTH1F, {{200, 0, 2}}, false);
    fRegistry.add("V0/K0S/hCosPA", "cosPA;cosine of pointing angle", kTH1F, {{100, 0.9, 1}}, false);
    fRegistry.add("V0/K0S/hDCA2Legs", "distance between 2 legs;distance between 2 legs (cm)", kTH1F, {{100, 0, 1}}, false);
    fRegistry.addClone("V0/K0S/", "V0/Lambda/");
    fRegistry.addClone("V0/K0S/", "V0/AntiLambda/");
    fRegistry.add("V0/K0S/hMassK0S", "K0S mass;m_{#pi#pi} (GeV/c^{2})", kTH1F, {{200, 0.4, 0.6}}, false);
    fRegistry.add("V0/Lambda/hMassLambda", "Lambda mass;m_{p#pi} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);
    fRegistry.add("V0/AntiLambda/hMassAntiLambda", "Anti-Lambda mass;m_{p#pi} (GeV/c^{2})", kTH1F, {{100, 1.08, 1.18}}, false);

    const AxisSpec axis_pt{{0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "p_{T,e} (GeV/c)"};
    const AxisSpec axis_dca_sigma{{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10}, "DCA_{e}^{3D} (#sigma)"};

    // for tracks
    fRegistry.add("LF/electron/findable/hs", "electron;p_{T,e} (GeV/c);#eta_{e};#varphi_{e} (rad.);DCA_{e}^{3D} (#sigma)", kTHnSparseF, {{axis_pt}, {80, -2, +2}, {36, 0, 2 * M_PI}, {axis_dca_sigma}}, false);
    fRegistry.add("LF/electron/findable/hTPCdEdx", "TPC dE/dx;p_{in} (GeV/c);TPC dE/dx (a.u.)", kTH2F, {{1000, 0, 10}, {200, 0, 200}}, false);
    fRegistry.add("LF/electron/findable/hTOFbeta", "TOF #beta;p_{pv} (GeV/c);#beta", kTH2F, {{1000, 0, 10}, {240, 0, 1.2}}, false);

    fRegistry.addClone("LF/electron/findable/", "LF/electron/correct/");
    fRegistry.addClone("LF/electron/findable/", "LF/electron/fake/");
    fRegistry.addClone("LF/electron/", "promptJpsi/electron/");
    fRegistry.addClone("LF/electron/", "nonpromptJpsi/electron/");

    fRegistry.addClone("LF/electron/", "D0/electron/");  // D0 -> K- e+ nu, Br = 0.03549 | D0 -> K- e+ pi0 nu, Br = 0.016 | D0 -> K*(892)- e+ nu, Br = 0.0215 // D0 -> anti-K0S e+ pi- nu, Br = 0.0144
    fRegistry.addClone("LF/electron/", "Dpm/electron/"); // D+ -> K- pi+ e+ nu, Br = 0.0402 | D+ -> anti-K*(892)0 e+ nu, Br = 0.0540 // D+ -> anti-K0S e+ nu, Br = 0.0872
    fRegistry.addClone("LF/electron/", "Ds/electron/");  // Ds+ -> K0S e+ nu, Br = 0.0034 // Ds+ -> phi e+ nu, Br = 0.0239
    fRegistry.addClone("LF/electron/", "Lc/electron/");  // Lc+ -> L e+ nu, Br = 0.0356
  }

  template <typename TTrack>
  bool isKaon(TTrack const& track)
  {
    // TOFif
    bool is_ka_included_TPC = trackcuts.cfg_min_TPCNsigmaKa < track.tpcNSigmaKa() && track.tpcNSigmaKa() < trackcuts.cfg_max_TPCNsigmaKa;
    bool is_ka_included_TOF = track.hasTOF() ? (trackcuts.cfg_min_TOFNsigmaKa < track.tofNSigmaKa() && track.tofNSigmaKa() < trackcuts.cfg_max_TOFNsigmaKa && track.tofChi2() < trackcuts.cfg_max_chi2tof) : true;
    return is_ka_included_TPC && is_ka_included_TOF;
  }

  template <typename TTrack>
  bool isPion(TTrack const& track)
  {
    // TOFif
    bool is_pi_included_TPC = v0cuts.cfg_min_TPCNsigmaPi < track.tpcNSigmaPi() && track.tpcNSigmaPi() < v0cuts.cfg_max_TPCNsigmaPi;
    bool is_pi_included_TOF = track.hasTOF() ? (v0cuts.cfg_min_TOFNsigmaPi < track.tofNSigmaPi() && track.tofNSigmaPi() < v0cuts.cfg_max_TOFNsigmaPi && track.tofChi2() < v0cuts.cfg_max_chi2tof) : true;
    return is_pi_included_TPC && is_pi_included_TOF;
  }

  template <typename TTrack>
  bool isProton(TTrack const& track)
  {
    // TOFif
    bool is_pr_included_TPC = v0cuts.cfg_min_TPCNsigmaPr < track.tpcNSigmaPr() && track.tpcNSigmaPr() < v0cuts.cfg_max_TPCNsigmaPr;
    bool is_pr_included_TOF = track.hasTOF() ? (v0cuts.cfg_min_TOFNsigmaPr < track.tofNSigmaPr() && track.tofNSigmaPr() < v0cuts.cfg_max_TOFNsigmaPr && track.tofChi2() < v0cuts.cfg_max_chi2tof) : true;
    return is_pr_included_TPC && is_pr_included_TOF;
  }

  template <typename TTrack, typename TTrackParCov>
  bool isSelectedTrack(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (trackParCov.getPt() < trackcuts.cfg_min_pt_track || trackcuts.cfg_max_pt_track < trackParCov.getPt()) {
      return false;
    }

    if (std::fabs(dcaXY) > trackcuts.cfg_max_dcaxy) {
      return false;
    }

    if (std::fabs(dcaZ) > trackcuts.cfg_max_dcaz) {
      return false;
    }

    if (track.itsChi2NCl() > trackcuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < trackcuts.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < trackcuts.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > trackcuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < trackcuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < trackcuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < trackcuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > trackcuts.cfg_max_frac_shared_clusters_tpc) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isSelectedV0Leg(TTrack const& track, const float dcaXY)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (std::fabs(dcaXY) < v0cuts.cfg_min_dcaxy) {
      return false;
    }

    if (track.itsChi2NCl() > v0cuts.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < v0cuts.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < v0cuts.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > v0cuts.cfg_max_chi2tpc) {
      return false;
    }

    if (track.tpcNClsFound() < v0cuts.cfg_min_ncluster_tpc) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < v0cuts.cfg_min_ncrossedrows_tpc) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < v0cuts.cfg_min_cr2findable_ratio_tpc) {
      return false;
    }

    if (track.tpcFractionSharedCls() > v0cuts.cfg_max_frac_shared_clusters_tpc) {
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

  template <int charmHadronId, int findId, typename TTrack, typename TTrackParCov>
  void fillElectronHistograms(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    float dca3DinSigma = dca3DinSigmaOTF(dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());
    fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(hTypes[findId]) + HIST("hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
    fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(hTypes[findId]) + HIST("hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
    fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(hTypes[findId]) + HIST("hTOFbeta"), trackParCov.getP(), track.beta());
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

  template <int pairId, typename TCollision, typename TElectron, typename TTrackIds, typename TTracks, typename TMCParticles, typename TMCCollisions>
  void runPairEandTrack(TCollision const& collision, TElectron const& ele, TTrackIds const& trackIds, TTracks const& tracks, TMCParticles const& mcParticles, TMCCollisions const&)
  {
    std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    const auto& mcele = ele.template mcParticle_as<aod::McParticles>();
    const auto& mcCollision1 = mcele.template mcCollision_as<aod::McCollisions>();

    if (cfgEventGeneratorType >= 0 && mcCollision1.getSubGeneratorId() != cfgEventGeneratorType) {
      return;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto eleParCov = getTrackParCov(ele);
    eleParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, eleParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    for (const auto& trackId : trackIds) {
      const auto& track = tracks.rawIteratorAt(trackId);
      const auto& mctrack = track.template mcParticle_as<aod::McParticles>();
      const auto& mcCollision2 = mctrack.template mcCollision_as<aod::McCollisions>();

      if (cfgEventGeneratorType >= 0 && mcCollision2.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }
      if (mcCollision1.globalIndex() != mcCollision2.globalIndex()) {
        continue;
      }

      auto trackParCov = getTrackParCov(track);
      std::array<float, 3> svpos = {0.}; // secondary vertex position
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

      int nCand = 0;
      try {
        nCand = fitter.process(eleParCov, trackParCov);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      fitter.propagateTracksToVertex(); // propagate e and K to D vertex
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        svpos[i] = vtx[i];
      }
      fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
      fitter.getTrack(1).getPxPyPzGlo(pvec1); // strange hadron
      std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

      float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
      float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
      float lz = std::fabs(svpos[2] - collision.posZ());
      float mEK = RecoDecay::m(std::array{pvec0, pvec1}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassKaonCharged});
      float cpa = RecoDecay::cpa(pVtx, svpos, pvecSum);
      float cpaXY = RecoDecay::cpaXY(pVtx, svpos, pvecSum);
      float ptEK = RecoDecay::sqrtSumOfSquares(pvec0[0] + pvec1[0], pvec0[1] + pvec1[1]);

      if (cpa < svcuts.cfg_min_cospa || cpaXY < svcuts.cfg_min_cospaXY || svcuts.cfg_max_mass_eK < mEK || lxy < svcuts.cfg_min_lxy || svcuts.cfg_max_dca2legs < dca2legs) {
        continue;
      }

      fRegistry.fill(HIST("e_Kpm/all/hDCA2Legs"), dca2legs);
      fRegistry.fill(HIST("e_Kpm/all/hLxy"), lxy);
      fRegistry.fill(HIST("e_Kpm/all/hLz"), lz);
      fRegistry.fill(HIST("e_Kpm/all/hCosPAXY"), cpaXY);
      fRegistry.fill(HIST("e_Kpm/all/hCosPA"), cpa);
      fRegistry.fill(HIST("e_Kpm/all/hMass"), mEK, ptEK);

      int commonMotherId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2ProngsWithoutPDG(mcele, mctrack); // e and K+/-
      if (commonMotherId < 0 && mctrack.has_mothers()) {
        const auto& mctrack_mother = mctrack.template mothers_first_as<aod::McParticles>(); // mother particle of Kaon. For example K*(892)+ -> K+ pi0 or K*(892)0 -> K+ pi- and CC, or phi->K+K-
        if (std::abs(mctrack_mother.pdgCode()) == 313 || std::abs(mctrack_mother.pdgCode()) == 323 || std::abs(mctrack_mother.pdgCode()) == 333) {
          commonMotherId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2ProngsWithoutPDG(mcele, mctrack_mother); // e and K*(892)0 or K*(892)+/- or phi(1019)
        }
      }
      if (commonMotherId >= 0) { // common mother is correctly found by DCAFitterN.
        const auto& cmp = mcParticles.rawIteratorAt(commonMotherId);
        if (std::abs(cmp.pdgCode()) == 421) { // D0
          fillElectronHistograms<3, 1>(ele, eleParCov, dcaXY, dcaZ);
          fRegistry.fill(HIST("e_Kpm/D0/hDCA2Legs"), dca2legs);
          fRegistry.fill(HIST("e_Kpm/D0/hLxy"), lxy);
          fRegistry.fill(HIST("e_Kpm/D0/hLz"), lz);
          fRegistry.fill(HIST("e_Kpm/D0/hCosPAXY"), cpaXY);
          fRegistry.fill(HIST("e_Kpm/D0/hCosPA"), cpa);
          fRegistry.fill(HIST("e_Kpm/D0/hMass"), mEK, ptEK);
        } else if (std::abs(cmp.pdgCode()) == 411) { // Dpm
          fillElectronHistograms<4, 1>(ele, eleParCov, dcaXY, dcaZ);
          fRegistry.fill(HIST("e_Kpm/Dpm/hDCA2Legs"), dca2legs);
          fRegistry.fill(HIST("e_Kpm/Dpm/hLxy"), lxy);
          fRegistry.fill(HIST("e_Kpm/Dpm/hLz"), lz);
          fRegistry.fill(HIST("e_Kpm/Dpm/hCosPAXY"), cpaXY);
          fRegistry.fill(HIST("e_Kpm/Dpm/hCosPA"), cpa);
          fRegistry.fill(HIST("e_Kpm/Dpm/hMass"), mEK, ptEK);
        } else if (std::abs(cmp.pdgCode()) == 431) { // Ds
          fillElectronHistograms<5, 1>(ele, eleParCov, dcaXY, dcaZ);
          fRegistry.fill(HIST("e_Kpm/Ds/hDCA2Legs"), dca2legs);
          fRegistry.fill(HIST("e_Kpm/Ds/hLxy"), lxy);
          fRegistry.fill(HIST("e_Kpm/Ds/hLz"), lz);
          fRegistry.fill(HIST("e_Kpm/Ds/hCosPAXY"), cpaXY);
          fRegistry.fill(HIST("e_Kpm/Ds/hCosPA"), cpa);
          fRegistry.fill(HIST("e_Kpm/Ds/hMass"), mEK, ptEK);
        }
      } else { // common mother does not exist, but DCAFitterN found something. i.e. fake
        const auto& mp = mcele.template mothers_first_as<aod::McParticles>();
        if (std::abs(mp.pdgCode()) == 111 || std::abs(mp.pdgCode()) == 221 || std::abs(mp.pdgCode()) == 331 || std::abs(mp.pdgCode()) == 113 || std::abs(mp.pdgCode()) == 223 || std::abs(mp.pdgCode()) == 333) { // LF
          fillElectronHistograms<0, 2>(ele, eleParCov, dcaXY, dcaZ);
        }
        fRegistry.fill(HIST("e_Kpm/fake/hDCA2Legs"), dca2legs);
        fRegistry.fill(HIST("e_Kpm/fake/hLxy"), lxy);
        fRegistry.fill(HIST("e_Kpm/fake/hLz"), lz);
        fRegistry.fill(HIST("e_Kpm/fake/hCosPAXY"), cpaXY);
        fRegistry.fill(HIST("e_Kpm/fake/hCosPA"), cpa);
        fRegistry.fill(HIST("e_Kpm/fake/hMass"), mEK, ptEK);
      }
    } // end of kaon loop
  }

  template <int pairId, typename TCollision, typename TElectron, typename TV0Ids, typename TV0s, typename TMCParticles, typename TMCCollisions>
  void runPairEandV0(TCollision const& collision, TElectron const& ele, TV0Ids const& v0Ids, TV0s const& v0s, TMCParticles const& mcParticles, TMCCollisions const&)
  {
    std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    const auto& mcele = ele.template mcParticle_as<aod::McParticles>();
    const auto& mcCollision1 = mcele.template mcCollision_as<aod::McCollisions>();
    if (cfgEventGeneratorType >= 0 && mcCollision1.getSubGeneratorId() != cfgEventGeneratorType) {
      return;
    }

    o2::dataformats::DCA mDcaInfoCov;
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto eleParCov = getTrackParCov(ele);
    eleParCov.setPID(o2::track::PID::Electron);
    mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
    mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, eleParCov, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY = mDcaInfoCov.getY();
    float dcaZ = mDcaInfoCov.getZ();

    for (const auto& v0Id : v0Ids) {
      const auto& v0 = v0s.rawIteratorAt(v0Id);
      auto pos = v0.template posTrack_as<MyTracks>();
      // auto neg = v0.template negTrack_as<MyTracks>();

      const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
      // const auto& mcneg = neg.template mcParticle_as<aod::McParticles>();
      const auto& mcv0 = mcpos.template mothers_first_as<aod::McParticles>(); // check mother of K0S. namely, K0 [311 or -311].
      const auto& mcCollision2 = mcv0.template mcCollision_as<aod::McCollisions>();

      if (cfgEventGeneratorType >= 0 && mcCollision2.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }
      if (mcCollision1.globalIndex() != mcCollision2.globalIndex()) {
        continue;
      }

      const auto& mcv0_mother = mcv0.template mothers_first_as<aod::McParticles>(); // mother particle of K0S.

      const std::array<float, 3> vertex = {v0.x(), v0.y(), v0.z()};
      const std::array<float, 3> momentum = {v0.px(), v0.py(), v0.pz()};
      std::array<float, 21> covV0 = {0.f};

      constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
      for (int i = 0; i < 6; i++) {
        covV0[MomInd[i]] = v0.momentumCovMat()[i];
        covV0[i] = v0.positionCovMat()[i];
      }

      auto tV0 = o2::track::TrackParCov(vertex, momentum, covV0, 0, true);
      tV0.setAbsCharge(0);
      tV0.setPID(o2::track::PID::K0);

      std::array<float, 3> svpos = {0.}; // secondary vertex position
      std::array<float, 3> pvec0 = {0.};
      std::array<float, 3> pvec1 = {0.};

      int nCand = 0;
      try {
        nCand = fitter.process(eleParCov, tV0);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter process call!";
        continue;
      }
      if (nCand == 0) {
        continue;
      }

      fitter.propagateTracksToVertex(); // propagate e and K to D vertex
      const auto& vtx = fitter.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        svpos[i] = vtx[i];
      }
      fitter.getTrack(0).getPxPyPzGlo(pvec0); // electron
      fitter.getTrack(1).getPxPyPzGlo(pvec1); // v0
      std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

      float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
      float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2));
      float lz = std::fabs(svpos[2] - collision.posZ());
      float mEK = 0;
      float ptEK = RecoDecay::sqrtSumOfSquares(pvec0[0] + pvec1[0], pvec0[1] + pvec1[1]);
      if constexpr (pairId == 1) {
        mEK = RecoDecay::m(std::array{pvec0, pvec1}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassK0Short});
        if (svcuts.cfg_max_mass_eK < mEK) {
          continue;
        }
      } else if constexpr (pairId == 2) {
        mEK = RecoDecay::m(std::array{pvec0, pvec1}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassLambda});
        if (svcuts.cfg_max_mass_eL < mEK) {
          continue;
        }
      }
      float cpa = RecoDecay::cpa(pVtx, svpos, pvecSum);
      float cpaXY = RecoDecay::cpaXY(pVtx, svpos, pvecSum);

      if (cpa < svcuts.cfg_min_cospa || cpaXY < svcuts.cfg_min_cospaXY || lxy < svcuts.cfg_min_lxy || svcuts.cfg_max_dca2legs < dca2legs) {
        continue;
      }

      int commonMotherId = -1;
      if constexpr (pairId == 1) {
        commonMotherId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2ProngsWithoutPDG(mcele, mcv0_mother); // K0, not K0S
      } else if constexpr (pairId == 2) {
        commonMotherId = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2ProngsWithoutPDG(mcele, mcv0); // lambda
      }

      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hDCA2Legs"), dca2legs);
      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hLxy"), lxy);
      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hLz"), lz);
      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hCosPAXY"), cpaXY);
      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hCosPA"), cpa);
      fRegistry.fill(HIST(pair_names[pairId]) + HIST("all/hMass"), mEK, ptEK);

      if (commonMotherId >= 0) { // common mother is correctly found by DCAFitterN.
        const auto& cmp = mcParticles.rawIteratorAt(commonMotherId);
        if constexpr (pairId == 1) {
          if (std::abs(cmp.pdgCode()) == 421) { // D0
            fillElectronHistograms<3, 1>(ele, eleParCov, dcaXY, dcaZ);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hDCA2Legs"), dca2legs);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hLxy"), lxy);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hLz"), lz);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hCosPAXY"), cpaXY);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hCosPA"), cpa);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("D0/hMass"), mEK, ptEK);
          } else if (std::abs(cmp.pdgCode()) == 411) { // Dpm
            fillElectronHistograms<4, 1>(ele, eleParCov, dcaXY, dcaZ);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hDCA2Legs"), dca2legs);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hLxy"), lxy);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hLz"), lz);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hCosPAXY"), cpaXY);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hCosPA"), cpa);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Dpm/hMass"), mEK, ptEK);
          } else if (std::abs(cmp.pdgCode()) == 431) { // Ds
            fillElectronHistograms<5, 1>(ele, eleParCov, dcaXY, dcaZ);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hDCA2Legs"), dca2legs);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hLxy"), lxy);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hLz"), lz);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hCosPAXY"), cpaXY);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hCosPA"), cpa);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Ds/hMass"), mEK, ptEK);
          }
        } else if constexpr (pairId == 2) {
          if (std::abs(cmp.pdgCode()) == 4122) { // Lc
            fillElectronHistograms<6, 1>(ele, eleParCov, dcaXY, dcaZ);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hDCA2Legs"), dca2legs);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hLxy"), lxy);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hLz"), lz);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hCosPAXY"), cpaXY);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hCosPA"), cpa);
            fRegistry.fill(HIST(pair_names[pairId]) + HIST("Lc/hMass"), mEK, ptEK);
          }
        }
      } else { // common mother does not exist, but DCAFitterN found something. i.e. fake
        const auto& mp = mcele.template mothers_first_as<aod::McParticles>();
        if (std::abs(mp.pdgCode()) == 111 || std::abs(mp.pdgCode()) == 221 || std::abs(mp.pdgCode()) == 331 || std::abs(mp.pdgCode()) == 113 || std::abs(mp.pdgCode()) == 223 || std::abs(mp.pdgCode()) == 333) { // LF
          fillElectronHistograms<0, 2>(ele, eleParCov, dcaXY, dcaZ);
        }
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hDCA2Legs"), dca2legs);
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hLxy"), lxy);
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hLz"), lz);
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hCosPAXY"), cpaXY);
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hCosPA"), cpa);
        fRegistry.fill(HIST(pair_names[pairId]) + HIST("fake/hMass"), mEK, ptEK);
      }

    } // end of v0 loop
  }

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Preslice<aod::V0Datas> perCol_v0 = o2::aod::v0data::collisionId;

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && (cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < cfgZvtxMax);
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  std::vector<std::pair<int, int>> stored_trackIds;
  Filter trackFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  using MyFilteredTracks = soa::Filtered<MyTracks>;

  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1.
  Filter v0Filter = o2::aod::v0data::v0Type == uint8_t(1) && o2::aod::v0data::v0cosPA > v0cuts.cfg_min_cospa_v0hadron.value&& o2::aod::v0data::dcaV0daughters < v0cuts.cfg_max_pca_v0hadron.value;
  using filteredV0s = soa::Filtered<MyV0s>;

  std::vector<int> electronIds;
  std::vector<int> positronIds;
  std::vector<int> negKaonIds;
  std::vector<int> posKaonIds;

  std::vector<int> k0sIds;
  std::vector<int> antik0sIds;
  std::vector<int> lambdaIds;
  std::vector<int> antilambdaIds;

  void processSA(FilteredMyCollisions const&, aod::BCsWithTimestamps const&, MyTracks const&, filteredV0s const&, aod::McParticles const&, aod::McCollisions const&) {}
  PROCESS_SWITCH(taggingHFE, processSA, "process without TTCA", false);

  void processTTCA(FilteredMyCollisions const& collisions, aod::BCsWithTimestamps const&, MyTracks const& tracks, aod::TrackAssoc const& trackIndices, filteredV0s const& v0s, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    for (const auto& collision : collisions) {
      const auto& bc = collision.template foundBC_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (!collision.has_mcCollision()) {
        continue;
      }
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      const auto& mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (cfgEventGeneratorType < 0 || mcCollision.getSubGeneratorId() == cfgEventGeneratorType) {
        fillEventHistograms(collision);
      }
      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      electronIds.reserve(trackIdsThisCollision.size());
      positronIds.reserve(trackIdsThisCollision.size());
      negKaonIds.reserve(trackIdsThisCollision.size());
      posKaonIds.reserve(trackIdsThisCollision.size());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = trackId.template track_as<MyTracks>();
        if (!track.hasITS() || !track.hasTPC()) {
          continue;
        }
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto& mctrack = track.template mcParticle_as<aod::McParticles>();
        const auto& mp = mctrack.template mothers_first_as<aod::McParticles>(); // mother particle of electron
        if (!mctrack.has_mothers() || !(mctrack.isPhysicalPrimary() || mctrack.producedByGenerator())) {
          continue;
        }
        if (std::abs(mctrack.pdgCode()) != 11 && std::abs(mctrack.pdgCode()) != 321) {
          continue;
        }

        const auto& mcCollision1 = mctrack.template mcCollision_as<aod::McCollisions>();
        if (cfgEventGeneratorType >= 0 && mcCollision1.getSubGeneratorId() != cfgEventGeneratorType) {
          continue;
        }

        o2::dataformats::DCA mDcaInfoCov;
        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto track_par_cov_recalc = getTrackParCov(track);
        track_par_cov_recalc.setPID(o2::track::PID::Electron);
        mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
        mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();

        if (!isSelectedTrack(track, track_par_cov_recalc, dcaXY, dcaZ)) {
          continue;
        }

        if (std::abs(mctrack.pdgCode()) == 11) {
          if (track_par_cov_recalc.getEta() < trackcuts.cfg_min_eta_track || trackcuts.cfg_max_eta_track < track_par_cov_recalc.getEta()) {
            continue;
          }
          if (std::abs(mp.pdgCode()) == 111 || std::abs(mp.pdgCode()) == 221 || std::abs(mp.pdgCode()) == 331 || std::abs(mp.pdgCode()) == 113 || std::abs(mp.pdgCode()) == 223 || std::abs(mp.pdgCode()) == 333) { // LF
            fillElectronHistograms<0, 0>(track, track_par_cov_recalc, dcaXY, dcaZ);
          } else if (std::abs(mp.pdgCode()) == 443) { // Jpsi
            if (o2::aod::pwgem::dilepton::utils::mcutil::IsFromBeauty(mp, mcParticles) < 0) {
              fillElectronHistograms<1, 0>(track, track_par_cov_recalc, dcaXY, dcaZ); // prompt Jpsi
            } else {
              fillElectronHistograms<2, 0>(track, track_par_cov_recalc, dcaXY, dcaZ); // nonprompt Jpsi
            }
          } else if (std::abs(mp.pdgCode()) == 421) { // D0
            fillElectronHistograms<3, 0>(track, track_par_cov_recalc, dcaXY, dcaZ);
          } else if (std::abs(mp.pdgCode()) == 411) { // Dpm
            fillElectronHistograms<4, 0>(track, track_par_cov_recalc, dcaXY, dcaZ);
          } else if (std::abs(mp.pdgCode()) == 431) { // Ds
            fillElectronHistograms<5, 0>(track, track_par_cov_recalc, dcaXY, dcaZ);
          } else if (std::abs(mp.pdgCode()) == 4122) { // Lc
            fillElectronHistograms<6, 0>(track, track_par_cov_recalc, dcaXY, dcaZ);
          }

          if (track.sign() > 0) { // positron
            positronIds.emplace_back(trackId.trackId());
          } else { // electron
            electronIds.emplace_back(trackId.trackId());
          }
        } else if (std::abs(mctrack.pdgCode()) == 321 && isKaon(track)) {
          if (track_par_cov_recalc.getEta() < trackcuts.cfg_min_eta_track_wide || trackcuts.cfg_max_eta_track_wide < track_par_cov_recalc.getEta()) {
            continue;
          }
          if (track.sign() > 0) { // positive kaon
            posKaonIds.emplace_back(trackId.trackId());
          } else { // negative kaon
            negKaonIds.emplace_back(trackId.trackId());
          }
        }
      } // end of track loop for electron and kaon selection

      const auto& v0s_per_coll = v0s.sliceBy(perCol_v0, collision.globalIndex());
      k0sIds.reserve(v0s_per_coll.size());
      antik0sIds.reserve(v0s_per_coll.size());
      lambdaIds.reserve(v0s_per_coll.size());
      antilambdaIds.reserve(v0s_per_coll.size());

      for (const auto& v0 : v0s_per_coll) {
        if (v0cuts.cfg_min_mass_k0s < v0.mK0Short() && v0.mK0Short() < v0cuts.cfg_max_mass_k0s) {
          auto pos = v0.template posTrack_as<MyTracks>();
          auto neg = v0.template negTrack_as<MyTracks>();

          const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
          const auto& mcneg = neg.template mcParticle_as<aod::McParticles>();
          if (!mcpos.has_mothers() || !mcneg.has_mothers()) {
            continue;
          }

          if (o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(mcpos, mcneg, 211, -211, 310, mcParticles) < 0) {
            continue;
          }

          const auto& mcK0S = mcpos.template mothers_first_as<aod::McParticles>(); // mother particle of pion
          if (!mcK0S.has_mothers()) {
            continue;
          }
          const auto& mcK0 = mcK0S.template mothers_first_as<aod::McParticles>(); // mother particle of K0S
          const auto& mcCollision1 = mcK0S.template mcCollision_as<aod::McCollisions>();
          if (cfgEventGeneratorType >= 0 && mcCollision1.getSubGeneratorId() != cfgEventGeneratorType) {
            continue;
          }

          if (!isPion(pos) || !isPion(neg)) {
            continue;
          }

          float dcaXY = 999.f;
          o2::dataformats::DCA mDcaInfoCov;
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_pos = getTrackParCov(pos);
          track_par_cov_recalc_pos.setPID(o2::track::PID::Pion);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_pos, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(pos, dcaXY)) {
            continue;
          }

          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_neg = getTrackParCov(neg);
          track_par_cov_recalc_neg.setPID(o2::track::PID::Pion);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_neg, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(neg, dcaXY)) {
            continue;
          }

          fRegistry.fill(HIST("V0/K0S/hPt"), v0.pt());
          fRegistry.fill(HIST("V0/K0S/hYPhi"), v0.phi(), v0.yK0Short());
          fRegistry.fill(HIST("V0/K0S/hCosPA"), v0.v0cosPA());
          fRegistry.fill(HIST("V0/K0S/hLxy"), v0.v0radius());
          fRegistry.fill(HIST("V0/K0S/hDCA2Legs"), v0.dcaV0daughters());
          fRegistry.fill(HIST("V0/K0S/hAP"), v0.alpha(), v0.qtarm());
          fRegistry.fill(HIST("V0/K0S/hMassK0S"), v0.mK0Short());
          if (mcK0.pdgCode() == 311) {
            k0sIds.emplace_back(v0.globalIndex());
          } else if (mcK0.pdgCode() == -311) {
            antik0sIds.emplace_back(v0.globalIndex());
          }
        } else if (v0cuts.cfg_min_mass_lambda < v0.mLambda() && v0.mLambda() < v0cuts.cfg_max_mass_lambda) {
          auto pos = v0.template posTrack_as<MyTracks>();
          auto neg = v0.template negTrack_as<MyTracks>();

          const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
          const auto& mcneg = neg.template mcParticle_as<aod::McParticles>();
          if (!mcpos.has_mothers() || !mcneg.has_mothers()) {
            continue;
          }
          if (o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(mcpos, mcneg, 2212, -211, 3122, mcParticles) < 0) {
            continue;
          }

          if (!isProton(pos) || !isPion(neg)) {
            continue;
          }

          float dcaXY = 999.f;
          o2::dataformats::DCA mDcaInfoCov;
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_pos = getTrackParCov(pos);
          track_par_cov_recalc_pos.setPID(o2::track::PID::Proton);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_pos, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(pos, dcaXY)) {
            continue;
          }

          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_neg = getTrackParCov(neg);
          track_par_cov_recalc_neg.setPID(o2::track::PID::Pion);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_neg, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(neg, dcaXY)) {
            continue;
          }

          fRegistry.fill(HIST("V0/Lambda/hPt"), v0.pt());
          fRegistry.fill(HIST("V0/Lambda/hYPhi"), v0.phi(), v0.yLambda());
          fRegistry.fill(HIST("V0/Lambda/hCosPA"), v0.v0cosPA());
          fRegistry.fill(HIST("V0/Lambda/hLxy"), v0.v0radius());
          fRegistry.fill(HIST("V0/Lambda/hDCA2Legs"), v0.dcaV0daughters());
          fRegistry.fill(HIST("V0/Lambda/hAP"), v0.alpha(), v0.qtarm());
          fRegistry.fill(HIST("V0/Lambda/hMassLambda"), v0.mLambda());
          lambdaIds.emplace_back(v0.globalIndex());
        } else if (v0cuts.cfg_min_mass_lambda < v0.mAntiLambda() && v0.mAntiLambda() < v0cuts.cfg_max_mass_lambda) {
          auto pos = v0.template posTrack_as<MyTracks>();
          auto neg = v0.template negTrack_as<MyTracks>();

          const auto& mcpos = pos.template mcParticle_as<aod::McParticles>();
          const auto& mcneg = neg.template mcParticle_as<aod::McParticles>();
          if (!mcpos.has_mothers() || !mcneg.has_mothers()) {
            continue;
          }
          if (o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(mcpos, mcneg, 211, -2212, -3122, mcParticles) < 0) {
            continue;
          }

          if (!isPion(pos) || !isProton(neg)) {
            continue;
          }
          float dcaXY = 999.f;
          o2::dataformats::DCA mDcaInfoCov;
          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_pos = getTrackParCov(pos);
          track_par_cov_recalc_pos.setPID(o2::track::PID::Pion);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_pos, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(pos, dcaXY)) {
            continue;
          }

          mDcaInfoCov.set(999, 999, 999, 999, 999);
          auto track_par_cov_recalc_neg = getTrackParCov(neg);
          track_par_cov_recalc_neg.setPID(o2::track::PID::Proton);
          mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
          mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
          o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, track_par_cov_recalc_neg, 2.f, matCorr, &mDcaInfoCov);
          dcaXY = mDcaInfoCov.getY();
          if (!isSelectedV0Leg(neg, dcaXY)) {
            continue;
          }
          fRegistry.fill(HIST("V0/AntiLambda/hPt"), v0.pt());
          fRegistry.fill(HIST("V0/AntiLambda/hYPhi"), v0.phi(), v0.yLambda());
          fRegistry.fill(HIST("V0/AntiLambda/hCosPA"), v0.v0cosPA());
          fRegistry.fill(HIST("V0/AntiLambda/hLxy"), v0.v0radius());
          fRegistry.fill(HIST("V0/AntiLambda/hDCA2Legs"), v0.dcaV0daughters());
          fRegistry.fill(HIST("V0/AntiLambda/hAP"), v0.alpha(), v0.qtarm());
          fRegistry.fill(HIST("V0/AntiLambda/hMassAntiLambda"), v0.mAntiLambda());
          antilambdaIds.emplace_back(v0.globalIndex());
        }

      } // end of v0 loop

      for (const auto& trackId : electronIds) {
        const auto& ele = tracks.rawIteratorAt(trackId);
        runPairEandTrack<0>(collision, ele, posKaonIds, tracks, mcParticles, mcCollisions);
        runPairEandTrack<0>(collision, ele, negKaonIds, tracks, mcParticles, mcCollisions); // only for Ds
        runPairEandV0<1>(collision, ele, k0sIds, v0s, mcParticles, mcCollisions);
        runPairEandV0<1>(collision, ele, antik0sIds, v0s, mcParticles, mcCollisions); // only for Ds
        runPairEandV0<2>(collision, ele, antilambdaIds, v0s, mcParticles, mcCollisions);
      } // end of electron loop

      for (const auto& trackId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(trackId);
        runPairEandTrack<0>(collision, pos, negKaonIds, tracks, mcParticles, mcCollisions);
        runPairEandTrack<0>(collision, pos, posKaonIds, tracks, mcParticles, mcCollisions); // only for Ds
        runPairEandV0<1>(collision, pos, antik0sIds, v0s, mcParticles, mcCollisions);
        runPairEandV0<1>(collision, pos, k0sIds, v0s, mcParticles, mcCollisions); // only for Ds
        runPairEandV0<2>(collision, pos, lambdaIds, v0s, mcParticles, mcCollisions);
      } // end of positron loop

      electronIds.clear();
      electronIds.shrink_to_fit();
      positronIds.clear();
      positronIds.shrink_to_fit();
      negKaonIds.clear();
      negKaonIds.shrink_to_fit();
      posKaonIds.clear();
      posKaonIds.shrink_to_fit();

      k0sIds.clear();
      k0sIds.shrink_to_fit();
      antik0sIds.clear();
      antik0sIds.shrink_to_fit();
      lambdaIds.clear();
      lambdaIds.shrink_to_fit();
      antilambdaIds.clear();
      antilambdaIds.shrink_to_fit();

    } // end of collision loop
  }
  PROCESS_SWITCH(taggingHFE, processTTCA, "process with TTCA", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<taggingHFE>(cfgc, TaskName{"tagging-hfe"})};
}
