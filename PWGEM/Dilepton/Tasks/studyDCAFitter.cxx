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

/// \file studyDCAFitter.cxx
/// \brief a task to study tagging e from charm hadron decays in MC
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

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
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

struct studyDCAFitter {
  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::McCollisionLabels>;

  using MyTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta,
                             aod::McTrackLabels>;

  struct DielectronAtSV { // ee pair at SV
    bool isfound{false};
    float mass{-999.f};
    float pt{-999.f};
    float dca2legs{-999.f};
    float cospa{-999.f};
    float lxy{-999.f};
    float lz{-999.f};
    float lxyz = std::sqrt(std::pow(lxy, 2) + std::pow(lz, 2));
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
    Configurable<int> cfg_min_ncrossedrows_tpc{"cfg_min_ncrossedrows_tpc", 80, "min ncrossed rows"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncluster_itsib{"cfg_min_ncluster_itsib", 3, "min ncluster itsib"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 5.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 1.0, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 1.0, "max dca Z for single track in cm"};

    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2, "min TPC n sigma el inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3, "max TPC n sigma el inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -1e+10, "min TPC n sigma pi exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +3, "max TPC n sigma pi exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3, "min TOF n sigma el inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3, "max TOF n sigma el inclusion"};
  } electroncut;

  struct : ConfigurableGroup {
    std::string prefix = "svcut";
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min cospa"};
    Configurable<float> cfg_max_dca2legs{"cfg_max_dca2legs", 0.1, "max distance between 2 legs"};
  } svcut;

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", -1.f, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};
  Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
  Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};

  HistogramRegistry fRegistry{"fRegistry"};
  static constexpr std::string_view hadron_names[6] = {"LF/", "Jpsi/", "D0/", "Dpm/", "Ds/", "Lc/"};
  static constexpr std::string_view pair_names[3] = {"e_Kpm/", "e_K0S/", "e_Lambda/"};
  static constexpr std::string_view hTypes[4] = {"findable/", "correct/", "fake/", "miss/"};
  static constexpr std::string_view promptTypes[2] = {"prompt/", "nonprompt/"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(5.f);
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

  template <typename TBC>
  void initCCDB(TBC const& bc)
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

    // for pairs
    fRegistry.add("Pair/PV/Data/uls/hs", "hs;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);DCA_{ee}^{3D} (#sigma);", kTHnSparseF, {{500, 0, 5}, {100, 0, 10}, {100, 0, 10}}, true);
    fRegistry.addClone("Pair/PV/Data/uls/", "Pair/PV/Data/lspp/");
    fRegistry.addClone("Pair/PV/Data/uls/", "Pair/PV/Data/lsmm/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/PromptPhi/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/NonPromptPhi/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/PromptOmega/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/NonPromptOmega/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/PromptJpsi/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/NonPromptJpsi/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/c2e_c2e/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/b2e_b2e/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/b2c2e_b2c2e/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/b2c2e_b2e_sameb/");
    fRegistry.addClone("Pair/PV/Data/", "Pair/PV/MC/b2c2e_b2e_diffb/");

    fRegistry.add("Pair/SV/Data/uls/hs", "hs;m_{ee} (GeV/c^{2});p_{T,ee} (GeV/c);L_{xy} m_{ee}/p_{T,ee} (mm);", kTHnSparseF, {{500, 0, 5}, {100, 0, 10}, {200, -10, 10}}, true);
    fRegistry.add("Pair/SV/Data/uls/hCosPA", "cosPA;cosPA;", kTH1F, {{200, -1, 1}}, false);
    fRegistry.add("Pair/SV/Data/uls/hDCA2Legs", "distance between 2 legs at PCA;distance between 2 legs (cm);", kTH1F, {{100, 0, 0.1}}, false);
    fRegistry.addClone("Pair/SV/Data/uls/", "Pair/SV/Data/lspp/");
    fRegistry.addClone("Pair/SV/Data/uls/", "Pair/SV/Data/lsmm/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/PromptPhi/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/NonPromptPhi/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/PromptOmega/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/NonPromptOmega/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/PromptJpsi/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/NonPromptJpsi/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/c2e_c2e/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/b2e_b2e/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/b2c2e_b2c2e/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/b2c2e_b2e_sameb/");
    fRegistry.addClone("Pair/SV/Data/", "Pair/SV/MC/b2c2e_b2e_diffb/");
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
  bool isSelectedTrack(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
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

    if (track.itsChi2NCl() > electroncut.cfg_max_chi2its) {
      return false;
    }

    if (track.itsNCls() < electroncut.cfg_min_ncluster_its) {
      return false;
    }

    if (track.itsNClsInnerBarrel() < electroncut.cfg_min_ncluster_itsib) {
      return false;
    }

    if (track.tpcChi2NCl() > electroncut.cfg_max_chi2tpc) {
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

  template <int charmHadronId, int findId, int promptId, typename TTrack, typename TTrackParCov>
  void fillElectronHistograms(TTrack const& track, TTrackParCov const& trackParCov, const float dcaXY, const float dcaZ)
  {
    if (std::find(used_electronIds.begin(), used_electronIds.end(), std::make_pair(findId, track.globalIndex())) == used_electronIds.end()) {
      float dca3DinSigma = dca3DinSigmaOTF(dcaXY, dcaZ, trackParCov.getSigmaY2(), trackParCov.getSigmaZ2(), trackParCov.getSigmaZY());
      fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(promptTypes[promptId]) + HIST(hTypes[findId]) + HIST("hs"), trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), dca3DinSigma);
      fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(promptTypes[promptId]) + HIST(hTypes[findId]) + HIST("hTPCdEdx"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST(hadron_names[charmHadronId]) + HIST("electron/") + HIST(promptTypes[promptId]) + HIST(hTypes[findId]) + HIST("hTOFbeta"), trackParCov.getP(), track.beta());
      used_electronIds.emplace_back(std::make_pair(findId, track.globalIndex()));
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

  template <typename TTrack, typename TMCParticles>
  int FindLF(TTrack const& posmc, TTrack const& negmc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 22, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 111, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 100443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 100553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 200553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -11, 11, 300553, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <bool isMC, uint8_t signType, typename TTrack, typename TMCParticles>
  void runPairingAtPV(TTrack const& t1, TTrack const& t2, TMCParticles const& mcParticles)
  {
    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov1 = getTrackParCov(t1);
    trackParCov1.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov1, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY1 = mDcaInfoCov.getY();
    float dcaZ1 = mDcaInfoCov.getZ();

    mDcaInfoCov.set(999, 999, 999, 999, 999);
    auto trackParCov2 = getTrackParCov(t2);
    trackParCov2.setPID(o2::track::PID::Electron);
    o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov2, 2.f, matCorr, &mDcaInfoCov);
    float dcaXY2 = mDcaInfoCov.getY();
    float dcaZ2 = mDcaInfoCov.getZ();

    ROOT::Math::PtEtaPhiMVector v1(trackParCov1.getPt(), trackParCov1.getEta(), trackParCov1.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v2(trackParCov2.getPt(), trackParCov2.getEta(), trackParCov2.getPhi(), o2::constants::physics::MassElectron);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    float dca3DinSigma1 = dca3DinSigmaOTF(dcaXY1, dcaZ1, trackParCov1.getSigmaY2(), trackParCov1.getSigmaZ2(), trackParCov1.getSigmaZY());
    float dca3DinSigma2 = dca3DinSigmaOTF(dcaXY2, dcaZ2, trackParCov2.getSigmaY2(), trackParCov2.getSigmaZ2(), trackParCov2.getSigmaZY());
    float pair_dca = pairDCAQuadSum(dca3DinSigma1, dca3DinSigma2);

    if constexpr (signType == 0) { // ULS
      fRegistry.fill(HIST("Pair/PV/Data/uls/hs"), v12.M(), v12.Pt(), pair_dca);
    } else if constexpr (signType == 1) { // LS++
      fRegistry.fill(HIST("Pair/PV/Data/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
    } else if constexpr (signType == 2) { // LS--
      fRegistry.fill(HIST("Pair/PV/Data/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
    }

    if constexpr (isMC) {
      const auto& t1mc = t1.template mcParticle_as<aod::McParticles>();
      const auto& t2mc = t2.template mcParticle_as<aod::McParticles>();
      if (std::abs(t1mc.pdgCode()) != 11) {
        return;
      }
      if (!(t1mc.isPhysicalPrimary() || t1mc.producedByGenerator())) {
        return;
      }
      if (std::abs(t2mc.pdgCode()) != 11) {
        return;
      }
      if (!(t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
        return;
      }
      // const auto& mp1 = t1mc.template mothers_first_as<aod::McParticles>(); // mother particle of t1
      // const auto& mp2 = t2mc.template mothers_first_as<aod::McParticles>(); // mother particle of t2
      int mcCommonMotherid = FindLF(t1mc, t2mc, mcParticles);
      int hfee_type = IsHF(t1mc, t2mc, mcParticles);

      if (mcCommonMotherid > -1) {
        const auto cmp = mcParticles.rawIteratorAt(mcCommonMotherid);
        switch (std::abs(cmp.pdgCode())) {
          case 223:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/PV/MC/PromptOmega/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/PromptOmega/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/PromptOmega/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptOmega/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptOmega/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptOmega/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            }
            break;
          case 333:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/PV/MC/PromptPhi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/PromptPhi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/PromptPhi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptPhi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptPhi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptPhi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            }
            break;
          case 443:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/PV/MC/PromptJpsi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/PromptJpsi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/PromptJpsi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptJpsi/uls/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptJpsi/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/PV/MC/NonPromptJpsi/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
              }
            }
            break;
          default:
            break;
        } // end of switch for LF
      } else if (hfee_type > -1) {
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/MC/c2e_c2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/MC/c2e_c2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/MC/c2e_c2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBe_Be):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/MC/b2e_b2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/MC/b2e_b2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/MC/b2e_b2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_BCe):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2c2e/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2c2e/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2c2e/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_sameb/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_sameb/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_sameb/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_diffb/uls/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 1) { // LS+diff
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_diffb/lspp/hs"), v12.M(), v12.Pt(), pair_dca);
            } else if constexpr (signType == 2) { // LS-diff
              fRegistry.fill(HIST("Pair/PV/MC/b2c2e_b2e_diffb/lsmm/hs"), v12.M(), v12.Pt(), pair_dca);
            }
            break;

          default:
            break;
        } // end of switch for HFee
      } // end of HFee
    } // end of isMC
  }

  template <bool isMC, uint8_t signType, typename TCollision, typename TTrack, typename TMCParticles>
  void runSVFinder(TCollision const& collision, TTrack const& t1, TTrack const& t2, TMCParticles const& mcParticles)
  {
    DielectronAtSV eeatsv;

    auto trackParCov1 = getTrackParCov(t1);
    trackParCov1.setPID(o2::track::PID::Electron);
    auto trackParCov2 = getTrackParCov(t2);
    trackParCov2.setPID(o2::track::PID::Electron);

    std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
    std::array<float, 3> svpos = {0.}; // secondary vertex position
    std::array<float, 3> pvec0 = {0.};
    std::array<float, 3> pvec1 = {0.};

    int nCand = 0;
    try {
      nCand = fitter.process(trackParCov1, trackParCov2);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return;
    }
    if (nCand == 0) {
      return;
    }

    fitter.propagateTracksToVertex(); // propagate e and K to D vertex
    if (!fitter.isPropagateTracksToVertexDone()) {
      return;
    }
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      svpos[i] = vtx[i];
    }
    fitter.getTrack(0).getPxPyPzGlo(pvec0);
    fitter.getTrack(1).getPxPyPzGlo(pvec1);
    std::array<float, 3> pvecSum = {pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

    float cpa = RecoDecay::cpa(pVtx, svpos, pvecSum);
    float dca2legs = std::sqrt(fitter.getChi2AtPCACandidate());
    // float lxy = std::sqrt(std::pow(svpos[0] - collision.posX(), 2) + std::pow(svpos[1] - collision.posY(), 2)); // in cm
    float lz = std::fabs(svpos[2] - collision.posZ()); // in cm

    float meeAtSV = RecoDecay::m(std::array{pvec0, pvec1}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
    float pteeAtSV = RecoDecay::sqrtSumOfSquares(pvecSum[0], pvecSum[1]);
    float lxy = RecoDecay::dotProd(std::array{pvecSum[0], pvecSum[1]}, std::array{svpos[0] - collision.posX(), svpos[1] - collision.posY()}) / pteeAtSV;
    float ppdl = lxy * 1e-2 * meeAtSV / pteeAtSV * 1e+3; // pseudo-proper decay length in mm

    if (cpa < svcut.cfg_min_cospa || svcut.cfg_max_dca2legs < dca2legs) {
      return;
    }

    if constexpr (signType == 0) { // ULS
      fRegistry.fill(HIST("Pair/SV/Data/uls/hs"), meeAtSV, pteeAtSV, ppdl);
      fRegistry.fill(HIST("Pair/SV/Data/uls/hCosPA"), cpa);
      fRegistry.fill(HIST("Pair/SV/Data/uls/hDCA2Legs"), dca2legs);
    } else if constexpr (signType == 1) { // LS++
      fRegistry.fill(HIST("Pair/SV/Data/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
      fRegistry.fill(HIST("Pair/SV/Data/lspp/hCosPA"), cpa);
      fRegistry.fill(HIST("Pair/SV/Data/lspp/hDCA2Legs"), dca2legs);
    } else if constexpr (signType == 2) { // LS--
      fRegistry.fill(HIST("Pair/SV/Data/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
      fRegistry.fill(HIST("Pair/SV/Data/lsmm/hCosPA"), cpa);
      fRegistry.fill(HIST("Pair/SV/Data/lsmm/hDCA2Legs"), dca2legs);
    }

    if constexpr (isMC) {
      const auto& t1mc = t1.template mcParticle_as<aod::McParticles>();
      const auto& t2mc = t2.template mcParticle_as<aod::McParticles>();
      if (std::abs(t1mc.pdgCode()) != 11) {
        return;
      }
      if (!(t1mc.isPhysicalPrimary() || t1mc.producedByGenerator())) {
        return;
      }
      if (std::abs(t2mc.pdgCode()) != 11) {
        return;
      }
      if (!(t2mc.isPhysicalPrimary() || t2mc.producedByGenerator())) {
        return;
      }
      // const auto& mp1 = t1mc.template mothers_first_as<aod::McParticles>(); // mother particle of t1
      // const auto& mp2 = t2mc.template mothers_first_as<aod::McParticles>(); // mother particle of t2
      int mcCommonMotherid = FindLF(t1mc, t2mc, mcParticles);
      int hfee_type = IsHF(t1mc, t2mc, mcParticles);

      if (mcCommonMotherid > -1) {
        const auto cmp = mcParticles.rawIteratorAt(mcCommonMotherid);
        switch (cmp.pdgCode()) {
          case 223:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptOmega/lsmm/hDCA2Legs"), dca2legs);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptOmega/lsmm/hDCA2Legs"), dca2legs);
              }
            }
            break;
          case 333:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptPhi/lsmm/hDCA2Legs"), dca2legs);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptPhi/lsmm/hDCA2Legs"), dca2legs);
              }
            }
            break;
          case 443:
            if (IsFromCharm(cmp, mcParticles) < 0 && IsFromBeauty(cmp, mcParticles) < 0) { // prompt
              if constexpr (signType == 0) {                                               // ULS
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/PromptJpsi/lsmm/hDCA2Legs"), dca2legs);
              }
            } else {                         // nonprompt
              if constexpr (signType == 0) { // ULS
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/uls/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/uls/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/uls/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 1) { // LS++
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lspp/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lspp/hDCA2Legs"), dca2legs);
              } else if constexpr (signType == 2) { // LS--
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lsmm/hCosPA"), cpa);
                fRegistry.fill(HIST("Pair/SV/MC/NonPromptJpsi/lsmm/hDCA2Legs"), dca2legs);
              }
            }
            break;
          default:
            break;
        } // end of switch for LF
      } else if (hfee_type > -1) {
        switch (hfee_type) {
          case static_cast<int>(EM_HFeeType::kCe_Ce):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/uls/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/uls/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lspp/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/c2e_c2e/lsmm/hDCA2Legs"), dca2legs);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBe_Be):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/uls/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/uls/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lspp/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2e_b2e/lsmm/hDCA2Legs"), dca2legs);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_BCe):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/uls/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/uls/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lspp/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2c2e/lsmm/hDCA2Legs"), dca2legs);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/uls/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/uls/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 1) { // LS++
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lspp/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 2) { // LS--
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_sameb/lsmm/hDCA2Legs"), dca2legs);
            }
            break;
          case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
            if constexpr (signType == 0) { // ULS
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/uls/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/uls/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/uls/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 1) { // LS+diff
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lspp/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lspp/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lspp/hDCA2Legs"), dca2legs);
            } else if constexpr (signType == 2) { // LS-diff
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lsmm/hs"), meeAtSV, pteeAtSV, ppdl);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lsmm/hCosPA"), cpa);
              fRegistry.fill(HIST("Pair/SV/MC/b2c2e_b2e_diffb/lsmm/hDCA2Legs"), dca2legs);
            }
            break;

          default:
            break;
        } // end of switch for HFee
      } // end of HFee
    } // end of isMC

    eeatsv.isfound = true;
    eeatsv.mass = meeAtSV;
    eeatsv.pt = pteeAtSV;
    eeatsv.cospa = cpa;
    eeatsv.dca2legs = dca2legs;
    eeatsv.lxy = lxy;
    eeatsv.lz = lz;
    return;
  }

  template <bool isMC, typename TBCs, typename TCollisions, typename TTracks, typename TTrackAssoc, typename TMCCollisions, typename TMCParticles>
  void run(TBCs const&, TCollisions const& collisions, TTracks const& tracks, TTrackAssoc const& trackIndices, TMCCollisions const&, TMCParticles const& mcParticles)
  {
    used_electronIds.reserve(tracks.size());

    for (const auto& collision : collisions) {
      const auto& bc = collision.template foundBC_as<TBCs>();
      initCCDB(bc);
      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);

      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);
      fillEventHistograms(collision);

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());
      electronIds.reserve(trackIdsThisCollision.size());
      positronIds.reserve(trackIdsThisCollision.size());
      mVtx.setPos({collision.posX(), collision.posY(), collision.posZ()});
      mVtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = trackId.template track_as<MyTracks>();
        if (!track.hasITS() || !track.hasTPC()) {
          continue;
        }

        if constexpr (isMC) {
          if (!track.has_mcParticle()) {
            continue;
          }
          const auto& mctrack = track.template mcParticle_as<aod::McParticles>();
          const auto& mcCollision = mctrack.template mcCollision_as<aod::McCollisions>();
          if (cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != cfgEventGeneratorType) {
            continue;
          }
          if (!mctrack.has_mothers()) {
            continue;
          }
        }

        mDcaInfoCov.set(999, 999, 999, 999, 999);
        auto trackParCov = getTrackParCov(track);
        trackParCov.setPID(o2::track::PID::Electron);
        o2::base::Propagator::Instance()->propagateToDCABxByBz(mVtx, trackParCov, 2.f, matCorr, &mDcaInfoCov);
        float dcaXY = mDcaInfoCov.getY();
        float dcaZ = mDcaInfoCov.getZ();

        if (isSelectedTrack(track, trackParCov, dcaXY, dcaZ) && isElectron(track)) {
          if (track.sign() > 0) { // positron
            positronIds.emplace_back(trackId.trackId());
          } else { // electron
            electronIds.emplace_back(trackId.trackId());
          }
        }
      } // end of track loop for electron selection

      for (const auto& posId : positronIds) {
        const auto& pos = tracks.rawIteratorAt(posId);
        for (const auto& eleId : electronIds) {
          const auto& ele = tracks.rawIteratorAt(eleId);
          runSVFinder<isMC, 0>(collision, pos, ele, mcParticles);
          runPairingAtPV<isMC, 0>(pos, ele, mcParticles);
        } // end of electron loop
      } // end of positron loop

      for (const auto& posId1 : positronIds) {
        const auto& pos1 = tracks.rawIteratorAt(posId1);
        for (const auto& posId2 : positronIds) {
          const auto& pos2 = tracks.rawIteratorAt(posId2);
          if (pos1.globalIndex() == pos2.globalIndex()) {
            continue;
          }
          runSVFinder<isMC, 1>(collision, pos1, pos2, mcParticles);
          runPairingAtPV<isMC, 1>(pos1, pos2, mcParticles);
        } // end of positron loop
      } // end of positron loop

      for (const auto& eleId1 : electronIds) {
        const auto& ele1 = tracks.rawIteratorAt(eleId1);
        for (const auto& eleId2 : electronIds) {
          const auto& ele2 = tracks.rawIteratorAt(eleId2);
          if (ele1.globalIndex() == ele2.globalIndex()) {
            continue;
          }
          runSVFinder<isMC, 2>(collision, ele1, ele2, mcParticles);
          runPairingAtPV<isMC, 2>(ele1, ele2, mcParticles);
        } // end of electron loop
      } // end of electron loop

      electronIds.clear();
      electronIds.shrink_to_fit();
      positronIds.clear();
      positronIds.shrink_to_fit();
    } // end of collision loop

    used_electronIds.clear();
    used_electronIds.shrink_to_fit();
  }

  SliceCache cache;
  Preslice<aod::TracksIU> perCol = o2::aod::track::collisionId;

  Filter collisionFilter_evsel = o2::aod::evsel::sel8 == true && (cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < cfgZvtxMax);
  Filter collisionFilter_centrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  // Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  // Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;

  std::vector<int> electronIds;
  std::vector<int> positronIds;
  std::vector<std::pair<int, int>> used_electronIds; // pair of hTypeId and electronId

  void processMC(FilteredMyCollisions const& collisions, aod::BCsWithTimestamps const& bcs, MyTracks const& tracks, aod::TrackAssoc const& trackIndices, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    run<true>(bcs, collisions, tracks, trackIndices, mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(studyDCAFitter, processMC, "processMC", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<studyDCAFitter>(cfgc, TaskName{"study-dcafitter"})};
}
