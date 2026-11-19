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
/// \file Pi0EtaMesonDalitzFormFactor.cxx
/// \brief This code builds PCM-Dalitz (gamma + e+e-) candidates and fills a THnSparse of (pT, m(ee gamma), m(ee)) to measure the Dalitz form factor of pi0/eta.
/// \author Joshua Koenig (joshua.konig@cern.ch)

#include "PWGEM/PhotonMeson/Core/DalitzEECut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/EventTables.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/DalitzBinning.h"
#include "PWGEM/PhotonMeson/Utils/EMPhoton.h"
// Dilepton headers
#include "PWGEM/Dilepton/Utils/EventMixingHandler.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

enum class PairKind {
  kSame,
  kMix,
  kTruePi0, // reconstructed, all three decay products from the same pi0
  kTrueEta  // reconstructed, all three decay products from the same eta
};

struct Pi0EtaMesonDalitzFormFactor {
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<float> d_bz_input{"d_bz_input", -999, "bz field in kG, -999 is automatic"};
  Configurable<float> maxY{"maxY", 0.8, "maximum rapidity for reconstructed particles"};
  Configurable<bool> cfgRequireTrueAssociation{"cfgRequireTrueAssociation", false, "flag to require true mc collision association for validated mesons"};
  Configurable<bool> cfgDoMix{"cfgDoMix", true, "flag for event mixing"};
  Configurable<int> ndepth{"ndepth", 10, "depth for event mixing"};
  Configurable<uint64_t> ndiff_bc_mix{"ndiff_bc_mix", 594, "difference in global BC required in mixed events"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.f, 0.5f, 1.f, 1.5f, 2.f, 2.5f, 3.f, 4.f, 5.f, 6.f, 8.f, 10.f, 20.f}, "pT of the meson candidate (GeV/c)"};
  ConfigurableAxis axisMmeson{"axisMmeson", {400, 0.f, 0.8f}, "m(ee#gamma) (GeV/c^{2})"};
  ConfigurableAxis axisMee{"axisMee", {VARIABLE_WIDTH, 0.f, 0.005f, 0.01f, 0.02f, 0.03f, 0.04f, 0.05f, 0.06f, 0.08f, 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f}, "m(ee) of the virtual photon (GeV/c^{2})"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", +10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
  } eventcuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min V0 CosPA"};
    Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to V0"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
  } pcmcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfg_min_mass{"cfg_min_mass", 0.0, "min mass"};
    Configurable<float> cfg_max_mass{"cfg_max_mass", 0.6, "max mass"};
    Configurable<bool> cfg_apply_phiv{"cfg_apply_phiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfg_require_itsib_any{"cfg_require_itsib_any", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfg_require_itsib_1st{"cfg_require_itsib_1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfg_phiv_slope{"cfg_phiv_slope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfg_phiv_intercept{"cfg_phiv_intercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfg_min_pt_track{"cfg_min_pt_track", 0.1, "min pT for single track"};
    Configurable<float> cfg_max_eta_track{"cfg_max_eta_track", 0.8, "max eta for single track"};
    Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    Configurable<int> cfg_min_ncluster_its{"cfg_min_ncluster_its", 5, "min ncluster its"};
    Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfg_max_dcaxy{"cfg_max_dcaxy", 0.05, "max dca XY for single track in cm"};
    Configurable<float> cfg_max_dcaz{"cfg_max_dcaz", 0.05, "max dca Z for single track in cm"};
    Configurable<float> cfg_max_dca3dsigma_track{"cfg_max_dca3dsigma_track", 1.5, "max DCA 3D in sigma"};
    Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to electron"};

    Configurable<int> cfg_pid_scheme{"cfg_pid_scheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfg_min_TPCNsigmaPi{"cfg_min_TPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_max_TPCNsigmaPi{"cfg_max_TPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfg_min_TOFNsigmaEl{"cfg_min_TOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfg_max_TOFNsigmaEl{"cfg_max_TOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  using MyCollisions = Join<aod::PMEvents, aod::EMEventsAlias>;
  using MyV0Photons = Filtered<Join<aod::V0PhotonsKF, aod::V0KFEMEventIds, aod::V0PhotonsKFPrefilterBitDerived>>;
  using MyPrimaryElectrons = Filtered<Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronDaEMEventIds, aod::EMPrimaryElectronsPrefilterBitDerived>>;

  using MyMCCollisions = Join<aod::PMEvents, aod::EMEventsAlias, aod::EMMCEventLabels>;
  using MyMCV0Legs = Join<aod::V0Legs, aod::V0LegMCLabels>;
  using MyMCPrimaryElectrons = Filtered<Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronDaEMEventIds, aod::EMPrimaryElectronsPrefilterBitDerived, aod::EMPrimaryElectronMCLabels>>;

  SliceCache cache;
  Preslice<MyV0Photons> perCollision_pcm = aod::v0photonkf::pmeventId;
  Preslice<MyPrimaryElectrons> perCollision_electron = aod::emprimaryelectronda::pmeventId;
  Filter prefilter_pcm = ifnode(pcmcuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::v0photonkf::pfbderived == static_cast<uint16_t>(0), true);
  Filter prefilter_primaryelectron = ifnode(dileptoncuts.cfg_apply_cuts_from_prefilter_derived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > static_cast<int8_t>(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt&& nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < static_cast<int8_t>(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;

  Partition<MyMCPrimaryElectrons> positronsMC = o2::aod::emprimaryelectron::sign > static_cast<int8_t>(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt&& nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;
  Partition<MyMCPrimaryElectrons> electronsMC = o2::aod::emprimaryelectron::sign < static_cast<int8_t>(0) && dileptoncuts.cfg_min_pt_track < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfg_max_eta_track;

  using EventKey = std::pair<int, int64_t>; // <ndf, global collision index>
  o2::aod::pwgem::dilepton::utils::EventMixingHandler<int, EventKey, o2::aod::pwgem::photonmeson::utils::EMPhoton> emhPhoton;
  o2::aod::pwgem::dilepton::utils::EventMixingHandler<int, EventKey, o2::aod::pwgem::photonmeson::utils::EMPhoton> emhDilepton;
  std::map<EventKey, uint64_t> mapEventToGlobalBC;
  std::vector<float> zvtx_bin_edges;
  int ndf = 0;
  std::unique_ptr<TH3D> hGenBinning; // (pT, |y|, m_ee) binning of the generated data, see DalitzBinning.h
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  int mRunNumber = -1;
  float d_bz = 0.f;

  void init(InitContext&)
  {
    const AxisSpec axPt{axisPt, "#it{p}_{T,ee#gamma} (GeV/#it{c})"};
    const AxisSpec axMmeson{axisMmeson, "#it{m}_{ee#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec axMee{axisMee, "#it{m}_{ee} (GeV/#it{c}^{2})"};
    fRegistry.add("Pair/same/hFormFactor", "ee#gamma candidates;#it{p}_{T,ee#gamma};#it{m}_{ee#gamma};#it{m}_{ee}", HistType::kTHnSparseD, {axMmeson, axPt, axMee}, true);
    if (doprocessAnalysis && cfgDoMix) {
      fRegistry.addClone("Pair/same/", "Pair/mix/");
      fRegistry.add("Pair/mix/hDiffBC", "diff. global BC in mixed event;|BC_{current} - BC_{mixed}|", HistType::kTH1D, {{10001, -0.5, 10000.5}}, true);
    }
    if (doprocessAnalysisMC) {
      fRegistry.addClone("Pair/same/", "Pair/true_pi0/");
      fRegistry.addClone("Pair/same/", "Pair/true_eta/");
      fRegistry.addClone("Pair/same/", "Gen/pi0/");
      fRegistry.addClone("Pair/same/", "Gen/eta/");

      namespace dalitzbinning = o2::aod::pwgem::photonmeson::utils::dalitzbinning;
      const auto ptBins = dalitzbinning::ptBins();
      const auto yBins = dalitzbinning::yBins();
      const auto meeBins = dalitzbinning::meeBins();
      hGenBinning = std::make_unique<TH3D>("hGenBinning", "binning of the generated Dalitz decays", ptBins.size() - 1, ptBins.data(), yBins.size() - 1, yBins.data(), meeBins.size() - 1, meeBins.data());
    }

    zvtx_bin_edges = std::vector<float>(ConfVtxBins.value.begin(), ConfVtxBins.value.end());
    zvtx_bin_edges.erase(zvtx_bin_edges.begin());
    emhPhoton.SetNdepth(ndepth);
    emhDilepton.SetNdepth(ndepth);

    DefineEMEventCut();
    DefinePCMCut();
    DefineDileptonCut();

    mRunNumber = 0;
    d_bz = 0.f;
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  template <typename TCollision>
  void initCCDB(TCollision const& collision)
  {
    if (mRunNumber == collision.runNumber()) {
      return;
    }
    mRunNumber = collision.runNumber();

    if (d_bz_input > -990) { // o2-linter: disable=magic-number (override value)
      d_bz = d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(d_bz) > 1e-5) {                   // o2-linter: disable=magic-number (override value)
        grpmag.setL3Current(30000.f / (d_bz / 5.0f)); // o2-linter: disable=magic-number (override value)
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
    } else {
      auto* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", collision.runNumber());
      if (grpmag == nullptr) {
        LOGF(fatal, "Could not retrieve GRPMagField for run %d", collision.runNumber());
      }
      d_bz = grpmag->getNominalL3Field();
    }
    fV0PhotonCut.SetD_Bz(d_bz);
  }

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(eventcuts.cfgZvtxMin, eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
  }

  void DefinePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    fV0PhotonCut.SetV0PtRange(pcmcuts.cfg_min_pt_v0, pcmcuts.cfg_max_pt_v0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfg_min_eta_v0, pcmcuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfg_min_v0radius, pcmcuts.cfg_max_v0radius, -1.f, -1.f);
    fV0PhotonCut.SetAPRange(pcmcuts.cfg_max_alpha_ap, pcmcuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(pcmcuts.cfg_reject_v0_on_itsib);

    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfg_min_TPCNsigmaEl, pcmcuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfg_max_chi2its);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfg_disable_tpconly_track);
  }

  void DefineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    fDileptonCut.SetMeeRange(dileptoncuts.cfg_min_mass, dileptoncuts.cfg_max_mass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfg_phiv_intercept) / dileptoncuts.cfg_phiv_slope; });
    fDileptonCut.ApplyPhiV(dileptoncuts.cfg_apply_phiv);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfg_require_itsib_any);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfg_require_itsib_1st);

    fDileptonCut.SetTrackPtRange(dileptoncuts.cfg_min_pt_track, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfg_max_eta_track, +dileptoncuts.cfg_max_eta_track);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfg_min_ncluster_tpc);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfg_min_ncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetMaxFracSharedClustersTPC(dileptoncuts.cfg_max_frac_shared_clusters_tpc);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfg_max_chi2tpc);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfg_max_chi2its);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfg_min_ncluster_its, 7);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfg_max_dcaxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfg_max_dcaz);
    fDileptonCut.SetTrackDca3DRange(0.f, dileptoncuts.cfg_max_dca3dsigma_track); // in sigma

    fDileptonCut.SetPIDScheme(dileptoncuts.cfg_pid_scheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfg_min_TPCNsigmaEl, dileptoncuts.cfg_max_TPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfg_min_TPCNsigmaPi, dileptoncuts.cfg_max_TPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfg_min_TOFNsigmaEl, dileptoncuts.cfg_max_TOFNsigmaEl);
  }

  /// \brief common pairing step (data, mixed, MC true): applies the rapidity cut and fills the sparse
  template <PairKind kind>
  bool fillPair(ROOT::Math::PtEtaPhiMVector const& v_gamma, ROOT::Math::PtEtaPhiMVector const& v_ee)
  {
    ROOT::Math::PtEtaPhiMVector v_eeg = v_gamma + v_ee;
    if (std::fabs(v_eeg.Rapidity()) > maxY) {
      return false;
    }
    if constexpr (kind == PairKind::kSame) {
      fRegistry.fill(HIST("Pair/same/hFormFactor"), v_eeg.M(), v_eeg.Pt(), v_ee.M());
    } else if constexpr (kind == PairKind::kMix) {
      fRegistry.fill(HIST("Pair/mix/hFormFactor"), v_eeg.M(), v_eeg.Pt(), v_ee.M());
    } else if constexpr (kind == PairKind::kTruePi0) {
      fRegistry.fill(HIST("Pair/true_pi0/hFormFactor"), v_eeg.M(), v_eeg.Pt(), v_ee.M());
    } else if constexpr (kind == PairKind::kTrueEta) {
      fRegistry.fill(HIST("Pair/true_eta/hFormFactor"), v_eeg.M(), v_eeg.Pt(), v_ee.M());
    }
    return true;
  }

  static ROOT::Math::PtEtaPhiMVector toVector(o2::aod::pwgem::photonmeson::utils::EMPhoton const& p)
  {
    return {p.pt(), p.eta(), p.phi(), p.mass()};
  }

  /// \brief pairs photons of the current event with dileptons from pooled events, and vice versa
  template <typename TCollision>
  void runMixing(TCollision const& collision, int zbin, EventKey const& key_df_collision)
  {
    const auto& photons_this = emhPhoton.GetTracksPerCollision(key_df_collision);
    const auto& dileptons_this = emhDilepton.GetTracksPerCollision(key_df_collision);

    for (const auto& key_mix : emhPhoton.GetCollisionIdsFromEventPool(zbin)) {
      if (key_mix == key_df_collision) {
        continue;
      }
      const uint64_t globalBC_mix = mapEventToGlobalBC[key_mix];
      const uint64_t diffBC = std::max(collision.globalBC(), globalBC_mix) - std::min(collision.globalBC(), globalBC_mix);
      fRegistry.fill(HIST("Pair/mix/hDiffBC"), diffBC);
      if (diffBC < ndiff_bc_mix) {
        continue;
      }

      for (const auto& g : photons_this) {
        for (const auto& ee : emhDilepton.GetTracksPerCollision(key_mix)) {
          fillPair<PairKind::kMix>(toVector(g), toVector(ee));
        }
      }
      for (const auto& ee : dileptons_this) {
        for (const auto& g : emhPhoton.GetTracksPerCollision(key_mix)) {
          fillPair<PairKind::kMix>(toVector(g), toVector(ee));
        }
      }
    }
  }

  void processAnalysis(MyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const&, MyPrimaryElectrons const&)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      int zbin = std::lower_bound(zvtx_bin_edges.begin(), zvtx_bin_edges.end(), collision.posZ()) - zvtx_bin_edges.begin() - 1;
      zbin = std::clamp(zbin, 0, static_cast<int>(zvtx_bin_edges.size()) - 2);
      const EventKey key_df_collision{ndf, collision.globalIndex()};
      std::vector<int64_t> used_photonIds;
      std::vector<std::pair<int, int>> used_dileptonIds;
      int npair = 0;

      auto photons_per_collision = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
      auto positrons_per_collision = positrons->sliceByCached(aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
      auto electrons_per_collision = electrons->sliceByCached(aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);

      for (const auto& g1 : photons_per_collision) {
        if (!fV0PhotonCut.IsSelected<decltype(g1), aod::V0Legs>(g1)) {
          continue;
        }
        auto pos1 = g1.template posTrack_as<aod::V0Legs>();
        auto ele1 = g1.template negTrack_as<aod::V0Legs>();
        ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

        for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
          if (pos2.trackId() == ele2.trackId()) { // protection against pairing identical tracks
            continue;
          }
          if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) { // electrons already used by the PCM photon
            continue;
          }
          if (!fDileptonCut.IsSelected(pos2, ele2, d_bz)) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;
          if (!fillPair<PairKind::kSame>(v_gamma, v_ee)) {
            continue;
          }
          npair++;

          if (!cfgDoMix) {
            continue;
          }
          if (std::find(used_photonIds.begin(), used_photonIds.end(), g1.globalIndex()) == used_photonIds.end()) {
            emhPhoton.AddTrackToEventPool(key_df_collision, o2::aod::pwgem::photonmeson::utils::EMPhoton(g1.pt(), g1.eta(), g1.phi(), 0.f));
            used_photonIds.emplace_back(g1.globalIndex());
          }
          const std::pair<int, int> ids_ee{pos2.trackId(), ele2.trackId()};
          if (std::find(used_dileptonIds.begin(), used_dileptonIds.end(), ids_ee) == used_dileptonIds.end()) {
            emhDilepton.AddTrackToEventPool(key_df_collision, o2::aod::pwgem::photonmeson::utils::EMPhoton(v_ee.Pt(), v_ee.Eta(), v_ee.Phi(), v_ee.M()));
            used_dileptonIds.emplace_back(ids_ee);
          }
        } // end of dielectron loop
      } // end of photon loop

      if (!cfgDoMix || npair == 0) {
        continue;
      }
      runMixing(collision, zbin, key_df_collision);

      emhPhoton.AddCollisionIdAtLast(zbin, key_df_collision);
      if (auto evicted = emhDilepton.AddCollisionIdAtLast(zbin, key_df_collision)) {
        mapEventToGlobalBC.erase(*evicted);
      }
      mapEventToGlobalBC[key_df_collision] = collision.globalBC();
    } // end of collision loop
    ndf++;
  }
  PROCESS_SWITCH(Pi0EtaMesonDalitzFormFactor, processAnalysis, "process ee-gamma pair analysis", true);

  /// \brief fills the generated pi0/eta -> e+ e- gamma distributions from the binned data of the MC event
  template <bool isPi0, typename TBinnedData>
  void fillBinnedGen(TBinnedData const& binned_data, const float mass)
  {
    int xbin = 0, ybin = 0, zbin = 0;
    for (int ibin = 0; ibin < hGenBinning->GetNcells(); ibin++) {
      const int nentry = binned_data[ibin];
      if (nentry == 0) {
        continue;
      }
      if (hGenBinning->IsBinUnderflow(ibin) || hGenBinning->IsBinOverflow(ibin)) {
        continue;
      }
      hGenBinning->GetBinXYZ(ibin, xbin, ybin, zbin);
      if (hGenBinning->GetYaxis()->GetBinCenter(ybin) > maxY) {
        continue;
      }
      const double pt = hGenBinning->GetXaxis()->GetBinCenter(xbin);
      const double mee = hGenBinning->GetZaxis()->GetBinCenter(zbin);
      if constexpr (isPi0) {
        fRegistry.fill(HIST("Gen/pi0/hFormFactor"), mass, pt, mee, static_cast<double>(nentry));
      } else {
        fRegistry.fill(HIST("Gen/eta/hFormFactor"), mass, pt, mee, static_cast<double>(nentry));
      }
    }
  }

  void processAnalysisMC(MyMCCollisions const& collisions, MyV0Photons const& v0photons, MyMCV0Legs const&, MyMCPrimaryElectrons const&, Join<aod::EMMCEvents, aod::BinnedGenDalitzs> const&, aod::EMMCParticles const& mcparticles)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      auto mccollision = collision.template emmcevent_as<Join<aod::EMMCEvents, aod::BinnedGenDalitzs>>();
      fillBinnedGen<true>(mccollision.generatedPi0Dalitz(), o2::constants::physics::MassPionNeutral);
      fillBinnedGen<false>(mccollision.generatedEtaDalitz(), o2::constants::physics::MassEta);

      auto photons_per_collision = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
      auto positrons_per_collision = positronsMC->sliceByCached(aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);
      auto electrons_per_collision = electronsMC->sliceByCached(aod::emprimaryelectronda::pmeventId, collision.globalIndex(), cache);

      for (const auto& g1 : photons_per_collision) {
        if (!fV0PhotonCut.IsSelected<decltype(g1), MyMCV0Legs>(g1)) {
          continue;
        }
        auto pos1 = g1.template posTrack_as<MyMCV0Legs>();
        auto ele1 = g1.template negTrack_as<MyMCV0Legs>();
        auto pos1mc = pos1.template emmcparticle_as<aod::EMMCParticles>();
        auto ele1mc = ele1.template emmcparticle_as<aod::EMMCParticles>();
        const int photonid = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom2Prongs(pos1mc, ele1mc, -11, 11, 22, mcparticles);
        if (photonid < 0) {
          continue;
        }
        auto g1mc = mcparticles.iteratorAt(photonid);
        ROOT::Math::PtEtaPhiMVector v_gamma(g1.pt(), g1.eta(), g1.phi(), 0.);

        for (const auto& [pos2, ele2] : combinations(CombinationsFullIndexPolicy(positrons_per_collision, electrons_per_collision))) {
          if (pos2.trackId() == ele2.trackId()) { // protection against pairing identical tracks
            continue;
          }
          if (pos1.trackId() == pos2.trackId() || ele1.trackId() == ele2.trackId()) { // electrons already used by the PCM photon
            continue;
          }
          if (!fDileptonCut.IsSelected(pos2, ele2, d_bz)) {
            continue;
          }

          auto pos2mc = mcparticles.iteratorAt(pos2.emmcparticleId());
          auto ele2mc = mcparticles.iteratorAt(ele2.emmcparticleId());
          const int pi0id = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 111, mcparticles);
          const int etaid = o2::aod::pwgem::dilepton::utils::mcutil::FindCommonMotherFrom3Prongs(g1mc, pos2mc, ele2mc, 22, -11, 11, 221, mcparticles);
          if (pi0id < 0 && etaid < 0) {
            continue;
          }

          ROOT::Math::PtEtaPhiMVector v_pos(pos2.pt(), pos2.eta(), pos2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ele(ele2.pt(), ele2.eta(), ele2.phi(), o2::constants::physics::MassElectron);
          ROOT::Math::PtEtaPhiMVector v_ee = v_pos + v_ele;

          if (pi0id > 0) {
            if (cfgRequireTrueAssociation && mcparticles.iteratorAt(pi0id).emmceventId() != collision.emmceventId()) {
              continue;
            }
            fillPair<PairKind::kTruePi0>(v_gamma, v_ee);
          } else {
            if (cfgRequireTrueAssociation && mcparticles.iteratorAt(etaid).emmceventId() != collision.emmceventId()) {
              continue;
            }
            fillPair<PairKind::kTrueEta>(v_gamma, v_ee);
          }
        } // end of dielectron loop
      } // end of photon loop
    } // end of collision loop
  }
  PROCESS_SWITCH(Pi0EtaMesonDalitzFormFactor, processAnalysisMC, "process ee-gamma pairs with MC truth (true and generated mesons)", false);

  void processDummy(aod::PMEvents const&) {}
  PROCESS_SWITCH(Pi0EtaMesonDalitzFormFactor, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{
    adaptAnalysisTask<Pi0EtaMesonDalitzFormFactor>(context)};
}
