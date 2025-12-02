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
// TableProducer to generate Trees for pure protons ans pions (from Lambda), Kaons (from Omegas), deuterons (identified with TOF) and He3 (identified with TPC).
// The output trees contain the ITS cluster size information.
//
// Author: Giorgio Alberto Lucia

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGLF/DataModel/LFClusterStudiesTable.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"
#include "TPDGCode.h"

#include <cstdint>
#include <cstring>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace ::o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU, aod::pidTPCEl, aod::pidTOFEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFDe, aod::TOFSignal, aod::TOFEvTime>;
using TracksFullIUMc = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TracksCovIU, aod::pidTPCEl, aod::pidTOFEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFDe, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;

using CollisionsCustom = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionsCustomMc = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

namespace BetheBloch
{

constexpr double defaultParams[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> parNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

} // namespace BetheBloch

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda,
  Photon, // deprecated, electrons are now selected from pi0 photons
  V0TypeAll
};

enum CascadeType : uint8_t {
  XiMinus = 0,
  OmegaMinus
};

enum Selections {
  kNoCut = 0,
  kSel8,
  kVtxZ,
  kAll
};

enum V0Selections {
  kV0NoCut = 0,
  kV0DaughterQuality,
  kV0Topology,
  kV0PID,
  kV0All
};

enum CascSelections {
  kCascNoCut = 0,
  kCascTopology,
  kRejectedXi,
  kAcceptedOmega,
  kNSigmaTPC,
  kCascAll
};

enum NucleiSelections {
  kNucleiNoCut = 0,
  kNucleiNClsIts,
  kNucleiPIDtpc,
  kNucleiPIDtof,
  kNucleiAll
};

enum ESelections {
  kENoCut = 0,
  kETrackQuality,
  kEPrimary,
  kEPid,
  kEPi0,
  kEAll
};

enum PartID {
  none = 0,
  el,
  pi,
  ka,
  pr,
  de,
  he,
  all
};

static constexpr std::string_view cNames[] = {"none", "electron", "pion", "kaon", "proton", "deuteron", "He3"};

struct Candidate {
  float p = -999.f; // momentum * charge
  float eta = -999.f;
  float phi = -999.f;
  uint32_t itsClusterSize = 0;
  uint8_t partID = PartID::none;
  float pTPC = -999.f;
  uint32_t pidInTrk = 0; // PID in tracking
  float nsigmaTPC = -999.f;
  float nsigmaTOF = -999.f;
  float tofMass = -999.f;
  float cosPAMother = -999.f; // Cosine of the pointing angle of the mother
  float massMother = -999.f;  // Invariant mass of the mother
  int pdgCode = 0;
};

struct LfTreeCreatorClusterStudies {

  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  SliceCache m_cache;
  int m_runNumber;
  int m_collisionCounter = 0;
  float m_d_bz;
  uint32_t m_randomSeed = 0.;

  Configurable<bool> setting_fillV0{"fillV0", true, "Fill the V0 tree"};
  Configurable<bool> setting_fillExtraTable{"fillExtraTable", false, "Fill the extra table"};
  Configurable<bool> setting_fillCollTable{"fillCollTable", false, "Fill the collision table"};

  Configurable<int> setting_materialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  Configurable<float> setting_zVtxMax{"zVtxMax", 10.f, "Maximum z vertex position"};

  Configurable<float> setting_downscaleFactor{"downscaleFactor", 1.f, "Downscale factor for the V0 candidates"};
  Configurable<bool> setting_applyAdditionalEvSel{"applyAdditionalEvSel", false, "Apply additional event selection"};

  Configurable<float> track_nClsItsMin{"track_NclsItsMin", 0.f, "Minimum number of ITS clusters for the V0 daughters"};
  Configurable<float> track_nClsTpcMin{"track_NclsTpcMin", 100.f, "Minimum number of TPC clusters for the V0 daughters"};
  Configurable<float> track_nClsTpcMaxShared{"track_NclsTpcMaxShared", 5.f, "Maximum number of shared TPC clusters for the V0 daughters"};
  Configurable<float> track_etaMax{"etaMax", 0.8f, "Maximum eta"};
  Configurable<float> track_tpcChi2Min{"track_tpcChi2Min", 0.5f, "Minimum TPC chi2 per cluster"};

  // Configurable<float> v0setting_etaMaxV0{"etaMaxV0", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0setting_dcaV0daughters{"v0setting_dcaV0daughters", 1.f, "DCA between the V0 daughters"};
  Configurable<float> v0setting_dcaMinV0DaughterToPv{"v0setting_dcaMinV0DaughterToPv", 0.06f, "DCA of the daughters to the primary vertex"};
  Configurable<float> v0setting_radiusV0{"v0setting_radiusV0", 0.5f, "Maximum radius of the V0 accepted"};
  Configurable<float> v0setting_cosPA{"v0setting_cosPA", 0.98f, "Cosine of the pointing angle of the V0"};
  Configurable<float> v0setting_massWindowLambda{"v0setting_massWindowLambda", 0.02f, "Mass window for the Lambda"};
  Configurable<float> v0setting_massWindowK0s{"v0setting_massWindowK0s", 0.02f, "Mass window for the K0s"};
  Configurable<float> v0setting_nsigmatpcPi{"v0setting_nsigmaTPCPi", 2.f, "Number of sigmas for the TPC PID for pions"};
  Configurable<float> v0setting_nsigmatpcPr{"v0setting_nsigmaTPCPr", 2.f, "Number of sigmas for the TPC PID for protons"};
  Configurable<float> lambdasetting_qtAPcut{"lambdasetting_qtAPcut", 0.02f, "Cut on the qt for the Armenteros-Podolanski plot for photon rejection"};
  Configurable<float> lambdasetting_pmin{"lambdasetting_pmin", 0.0f, "Minimum momentum for the V0 daughters"};

  Configurable<float> cascsetting_dcaMinV0DaughterToPv{"cascsetting_dcaMinV0DaughterToPv", 0.03f, "DCA of one of the daugthers of Lambda to pv"};
  Configurable<float> cascsetting_dcaMinProtonToPv{"cascsetting_dcaMinProtonToPv", 0.03f, "DCA of the proton coming from Lambda to pv"};
  Configurable<float> cascsetting_dcaMinBachelorToPv{"cascsetting_dcaMinBachelorToPv", 0.04f, "DCA of the bachelor to pv"};
  Configurable<float> cascsetting_dcaMinV0ToPv{"cascsetting_dcaMinV0ToPv", 0.04f, "DCA of the V0 to pv"};
  Configurable<float> cascsetting_dcaV0Daughters{"cascsetting_dcaV0daughters", 0.4f, "DCA between the V0 daughters"};
  Configurable<float> cascsetting_cascCosPA{"cascsetting_cascCosPA", 0.99f, "Minimum cCosine of the pointing angle of the cascade"};
  Configurable<float> cascsetting_v0cosPA{"cascsetting_v0cosPA", 0.97f, "Minimum cCosine of the pointing angle of the v0"};
  Configurable<float> cascsetting_dcaMaxBachelorToV0{"cascsetting_dcaMaxBachelorToV0", 0.8f, "DCA of the bachelor to V0"};
  // Configurable<float> cascsetting_dcaMinBachelorToProton{"cascsetting_dcaMinBachelorToProton", 0.015f, "DCA of the bachelor to proton"};
  Configurable<float> cascsetting_radiusV0{"cascsetting_radiusV0", 1.2f, "Minimum radius of the V0 accepted"};
  Configurable<float> cascsetting_radiusCasc{"cascsetting_radiusCasc", 0.5f, "Minimum radius of the cascade accepted"};

  Configurable<float> cascsetting_massWindowOmega{"cascsetting_massWindowOmega", 0.01f, "Mass window for the Omega"};
  Configurable<float> cascsetting_massWindowXi{"cascsetting_massWindowXi", 0.01f, "Mass window for the Xi"};
  Configurable<float> cascsetting_nsigmatpc{"cascsetting_nsigmaTPC", 3.f, "Number of sigmas for the TPC PID"};

  Configurable<float> electronsetting_conversion_rmin{"electron_conversion_rmin", 1.76f, "Minimum radius for the photon conversion (cm)"};
  Configurable<float> electronsetting_conversion_rmax{"electron_conversion_rmax", 19.77f, "Maximum radius for the photon conversion (cm)"};
  Configurable<float> electronsetting_maxDcaxy{"electronsetting_maxDcaxy", 0.1f, "Maximum value for the DCAxy"};
  Configurable<float> electronsetting_maxDcaz{"electronsetting_maxDcaz", 0.5f, "Maximum value for the DCAz"};
  Configurable<float> electronsetting_minNsigmatpcEl{"electronsetting_minNsigmaTPCEl", -2.5f, "Minimum value for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsetting_maxNsigmatpcEl{"electronsetting_maxNsigmaTPCEl", 3.5f, "Maximum number for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsetting_maxNsigmatpcPi{"electronsetting_maxNsigmaTPCPi", 2.f, "Maximum number for the number of sigmas for pi rejection for the TPC PID for electrons"};
  Configurable<float> electronsetting_maxNsigmatpcKa{"electronsetting_maxNsigmaTPCKa", 2.f, "Maximum number for the number of sigmas for K rejection for the TPC PID for electrons"};
  Configurable<float> electronsetting_maxNsigmatpcPr{"electronsetting_maxNsigmaTPCPr", 2.f, "Maximum number for the number of sigmas for p rejection for the TPC PID for electrons"};
  Configurable<float> electronsetting_maxNsigmatofEl{"electronsetting_maxNsigmaTOFEl", 4.f, "Minimum value for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsetting_minPt{"electronsetting_minPt", 0.f, "Minimum pT accepted for electrons"};

  Configurable<int> desetting_nClsIts{"desetting_nClsIts", 6, "Minimum number of ITS clusters"};
  Configurable<float> desetting_nsigmatpc{"desetting_nsigmaCutTPC", 2.f, "Number of sigmas for the TPC PID"};
  Configurable<float> desetting_nsigmatof{"desetting_nsigmaCutTOF", 2.f, "Number of sigmas for the TOF PID"};
  Configurable<int> he3setting_nClsIts{"he3setting_nClsIts", 6, "Minimum number of ITS clusters"};
  Configurable<bool> he3setting_compensatePIDinTracking{"he3setting_compensatePIDinTracking", true, "Compensate PID in tracking"};
  Configurable<float> he3setting_nsigmatpc{"he3setting_nsigmaCutTPC", 2.f, "Number of sigmas for the TPC PID"};
  Configurable<float> he3setting_tofmasslow{"he3setting_tofmasslow", 1.8f, "Lower limit for the TOF mass"};
  Configurable<float> he3setting_tofmasshigh{"he3setting_tofmasshigh", 4.2f, "Upper limit for the TOF mass"};

  // Bethe Bloch parameters
  std::array<float, 6> m_BBparamsDe, m_BBparamsHe;
  Configurable<LabeledArray<double>> setting_BetheBlochParams{"setting_BetheBlochParams", {BetheBloch::defaultParams[0], 2, 6, {"De", "He3"}, BetheBloch::parNames}, "TPC Bethe-Bloch parameterisation for nuclei"};

  Preslice<aod::V0s> m_perCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::Cascades> m_perCollisionCascade = o2::aod::cascade::collisionId;
  Preslice<TracksFullIU> m_perCol = aod::track::collisionId;
  Preslice<TracksFullIUMc> m_perColMC = aod::track::collisionId;

  HistogramRegistry m_hAnalysis{
    "LFTreeCreator",
    {
      {"collision_selections", "Collision selection; selection; counts", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
      {"v0_selections", "V0 selection; selection; counts", {HistType::kTH1F, {{V0Selections::kV0All, -0.5, static_cast<double>(V0Selections::kV0All) - 0.5}}}},
      {"casc_selections", "Cascade selection; selection; counts", {HistType::kTH1F, {{CascSelections::kCascAll, -0.5, static_cast<double>(CascSelections::kCascAll) - 0.5}}}},
      {"e_selections", "e^{#pm} selection; selection; counts", {HistType::kTH1F, {{ESelections::kEAll, -0.5, static_cast<double>(ESelections::kEAll) - 0.5}}}},
      {"v0_type", "Selected V0; particle; counts", {HistType::kTH1F, {{V0Type::V0TypeAll, -0.5, static_cast<double>(V0Type::V0TypeAll) - 0.5}}}},
      {"radiusV0", "Decay radius (xy) V0; radius (cm); counts", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"massLambda", "#Lambda invariant mass; signed #it{p}_{T} (GeV/#it{c}); #it{m}_{#Lambda} (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {200, 1.08f, 1.18f}}}},
      {"massLambdaMc", "#Lambda invariant mass (MC); signed #it{p}_{T} (GeV/#it{c}); #it{m}_{#Lambda} (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {200, 1.08f, 1.18f}}}},
      {"Lambda_vs_K0s", "Mass #Lambda vs K^{0}_s; #it{m}_{K^{0}_{s}} (GeV/#it{c}^{2}); #it{m}_{#Lambda} (GeV/#it{c}^{2})", {HistType::kTH2F, {{50, 0.f, 1.f}, {70, 0.6f, 2.f}}}},
      {"armenteros_plot_before_selections", "Armenteros-Podolanski plot; #alpha; #it{q}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
      {"armenteros_plot", "Armenteros-Podolanski plot; #alpha; #it{q}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
      {"armenteros_plot_lambda", "Armenteros-Podolanski plot (#Lambda only); #alpha; #it{q}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
      {"armenteros_plot_gamma", "Armenteros-Podolanski plot (#gamma only); #alpha; #it{q}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
      {"photon_radiusV0", "Photon conversion radius (xy) V0; #it{r} (cm); counts", {HistType::kTH1F, {{100, 0., 100.}}}},
      {"photon_conversion_position", "Photon conversion position; #it{x} (cm); #it{y} (cm)", {HistType::kTH2F, {{250, -5.f, 5.f}, {250, -5.f, 5.f}}}},
      {"photon_conversion_position_layer", "Photon conversion position (ITS layers); #it{x} (cm); #it{y} (cm)", {HistType::kTH2F, {{100, -5.f, 5.f}, {100, -5.f, 5.f}}}},
      {"casc_dca_daughter_pairs", "DCA (xy) for cascade daughter pairs; DCA_{#it{xy}} (cm); counts", {HistType::kTH1F, {{100, -0.1, 0.1}}}},
      {"Xi_vs_Omega", "Mass Xi vs Omega; mass Omega (GeV/#it{c}^{2}); #it{m}_#Xi (GeV/#it{c}^{2})", {HistType::kTH2F, {{50, 1.f, 2.f}, {50, 1.f, 2.f}}}},
      {"massOmega", "Mass #Omega; signed #it{p}_{T} (GeV/#it{c}); #it{m}_{#Omega} (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {400, 1.62f, 1.72f}}}},
      {"massOmegaMc", "Mass #Omega (MC); signed #it{p}_{T} (GeV/#it{c}); #it{m}_{#Omega} (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {400, 1.62f, 1.72f}}}},
      {"massPi0", "Mass #pi^{0}; #it{m}_{#pi^{0}} (GeV/#it{c}^{2})", {HistType::kTH1F, {{100, 0.0f, 0.200f}}}},
      {"massPi0Mc", "Mass #pi^{0} (MC); #it{m}_{#pi^{0}} (GeV/#it{c}^{2})", {HistType::kTH1F, {{100, 0.0f, 0.200f}}}},
      {"massPi0WithBkg", "Mass #pi^{0} with Background; #it{m}_{#pi^{0}} (GeV/#it{c}^{2}); counts", {HistType::kTH1F, {{100, 0.0f, 0.200f}}}},
      {"zVtx", "Binning for the vertex z in cm; #it{z}_{vertex} (cm)", {HistType::kTH1F, {{100, -20.f, 20.f}}}},
      {"isPositive", "is the candidate positive?; isPositive; counts", {HistType::kTH1F, {{2, -0.5f, 1.5f}}}},

      {"electron/DCAxyBeforeSelection", "DCA (xy) for cascade daughter pairs; DCA_{#it{xy}} (cm); counts", {HistType::kTH1F, {{100, -0.1, 0.1}}}},
      {"electron/DCAzBeforeSelection", "DCA (z) for cascade daughter pairs; DCA_{#it{z}} (cm); counts", {HistType::kTH1F, {{200, -0.2, 0.2}}}},
    },
    OutputObjHandlingPolicy::AnalysisObject,
    false};

  Produces<o2::aod::ClStTable> m_ClusterStudiesTable;
  Produces<o2::aod::ClStTableExtra> m_ClusterStudiesTableExtra;
  Produces<o2::aod::ClStTableColl> m_ClusterStudiesTableCollision;
  Produces<o2::aod::ClStTableMc> m_ClusterStudiesTableMc;

  o2::aod::ITSResponse m_responseITS;

  template <typename Tcollision>
  bool collisionSelection(const Tcollision& collision)
  {
    m_hAnalysis.fill(HIST("collision_selections"), Selections::kNoCut);
    if (!collision.sel8()) {
      return false;
    }
    m_hAnalysis.fill(HIST("collision_selections"), Selections::kSel8);
    // if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
    //   return false;
    // }
    if (std::abs(collision.posZ()) > setting_zVtxMax) {
      return false;
    }
    m_hAnalysis.fill(HIST("collision_selections"), Selections::kVtxZ);
    return true;
  }

  // =========================================================================================================

  /**
   * Select the V0 daughters based on the quality cuts
   */
  template <typename Track>
  bool qualityTrackSelection(const Track& track)
  {
    if (std::abs(track.eta()) > track_etaMax ||
        track.itsNCls() < track_nClsItsMin ||
        track.tpcNClsFound() < track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > track_nClsTpcMaxShared ||
        track.tpcChi2NCl() < track_tpcChi2Min ||
        track.tpcChi2NCl() > 4.0f) {
      return false;
    }
    return true;
  }

  bool qualitySelectionV0(aod::V0Datas::iterator const& v0)
  {
    if (std::abs(v0.dcapostopv()) < v0setting_dcaMinV0DaughterToPv)
      return false;
    if (std::abs(v0.dcanegtopv()) < v0setting_dcaMinV0DaughterToPv)
      return false;
    if (std::abs(v0.dcaV0daughters()) > v0setting_dcaV0daughters)
      return false;
    if (v0.v0radius() < v0setting_radiusV0)
      return false;
    if (std::abs(v0.v0cosPA()) < v0setting_cosPA)
      return false;

    return true;
  }

  bool qualitySelectionCascade(aod::CascDatas::iterator const& cascade, const std::array<float, 3>& pv)
  {
    if (std::abs(cascade.dcapostopv()) < cascsetting_dcaMinV0DaughterToPv)
      return false;
    if (std::abs(cascade.dcanegtopv()) < cascsetting_dcaMinV0DaughterToPv)
      return false;
    if (std::abs(cascade.dcabachtopv()) < cascsetting_dcaMinBachelorToPv)
      return false;
    if (std::abs(cascade.dcav0topv(pv[0], pv[1], pv[2])) < cascsetting_dcaMinV0ToPv)
      return false;
    if (std::abs(cascade.dcaV0daughters()) > cascsetting_dcaV0Daughters)
      return false;
    if (std::abs(cascade.casccosPA(pv[0], pv[1], pv[2])) < cascsetting_cascCosPA)
      return false;
    if (std::abs(cascade.v0cosPA(pv[0], pv[1], pv[2])) < cascsetting_v0cosPA)
      return false;
    if (std::abs(cascade.v0radius()) < cascsetting_radiusV0)
      return false;
    if (std::abs(cascade.cascradius()) < cascsetting_radiusCasc)
      return false;

    m_hAnalysis.fill(HIST("casc_dca_daughter_pairs"), cascade.dcacascdaughters());
    if (std::abs(cascade.dcacascdaughters()) > cascsetting_dcaMaxBachelorToV0)
      return false;

    return true;
  }

  uint8_t selectV0MotherHypothesis(aod::V0Datas::iterator const& v0)
  {
    uint8_t v0Bitmask(0);
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0setting_massWindowK0s) {
      SETBIT(v0Bitmask, K0s);
    }
    if ((std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0setting_massWindowLambda) && (v0.alpha() > 0)) {
      SETBIT(v0Bitmask, Lambda);
    }
    if ((std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < v0setting_massWindowLambda) && (v0.alpha() < 0)) {
      SETBIT(v0Bitmask, AntiLambda);
    }
    return v0Bitmask;
  }

  template <typename Track>
  bool selectPidV0Daughters(Candidate& candidatePos, Candidate& candidateNeg,
                            aod::V0Datas::iterator const& v0, const Track& posTrack, const Track& negTrack, uint8_t v0Bitmask)
  {
    if (TESTBIT(v0Bitmask, Lambda)) {
      if (v0.qtarm() < lambdasetting_qtAPcut)
        return false;
      if (std::abs(posTrack.tpcNSigmaPr()) > v0setting_nsigmatpcPr || std::abs(negTrack.tpcNSigmaPi()) > v0setting_nsigmatpcPi)
        return false;
      if (v0.p() < lambdasetting_pmin)
        return false;
      candidatePos.partID = PartID::pr;
      candidateNeg.partID = PartID::pi;
      candidatePos.nsigmaTPC = posTrack.tpcNSigmaPr();
      candidateNeg.nsigmaTPC = negTrack.tpcNSigmaPi();
      m_hAnalysis.fill(HIST("v0_type"), V0Type::Lambda);

    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      if (v0.qtarm() < lambdasetting_qtAPcut)
        return false;
      if (std::abs(posTrack.tpcNSigmaPi()) > v0setting_nsigmatpcPi || std::abs(negTrack.tpcNSigmaPr()) > v0setting_nsigmatpcPr)
        return false;
      if (v0.p() < lambdasetting_pmin)
        return false;
      candidatePos.partID = PartID::pi;
      candidateNeg.partID = PartID::pr;
      candidatePos.nsigmaTPC = posTrack.tpcNSigmaPi();
      candidateNeg.nsigmaTPC = negTrack.tpcNSigmaPr();
      m_hAnalysis.fill(HIST("v0_type"), V0Type::AntiLambda);

    } else {
      return false;
    }

    return true;
  }

  /**
   * Fill the histograms for the V0 candidate and return the mass of the V0
   */
  template <typename Track>
  float fillHistogramsV0(aod::V0Datas::iterator const& v0, const Track& trackPos, const Track& trackNeg, const uint8_t v0Bitmask)
  {
    float massV0{0.f};
    if (TESTBIT(v0Bitmask, Lambda)) {
      massV0 = v0.mLambda();
      m_hAnalysis.fill(HIST("massLambda"), v0.pt(), v0.mLambda());
      fillHistogramsParticle<PartID::pr, /*isMC*/ false>(trackPos);
      fillHistogramsParticle<PartID::pi, /*isMC*/ false>(trackNeg);
    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      massV0 = v0.mAntiLambda();
      m_hAnalysis.fill(HIST("massLambda"), v0.pt() * -1.f, v0.mAntiLambda());
      fillHistogramsParticle<PartID::pi, /*isMC*/ false>(trackPos);
      fillHistogramsParticle<PartID::pr, /*isMC*/ false>(trackNeg);
    }

    m_hAnalysis.fill(HIST("radiusV0"), v0.v0radius());
    m_hAnalysis.fill(HIST("armenteros_plot_lambda"), v0.alpha(), v0.qtarm());
    m_hAnalysis.fill(HIST("armenteros_plot"), v0.alpha(), v0.qtarm());

    return massV0;
  }

  template <const int partID, const bool isMC, typename Track>
  void fillHistogramsParticle(const Track& track)
  {
    float nsigmaTpc = -999.f;
    switch (partID) {
      case PartID::el:
        nsigmaTpc = track.tpcNSigmaEl();
        break;
      case PartID::pi:
        nsigmaTpc = track.tpcNSigmaPi();
        break;
      case PartID::ka:
        nsigmaTpc = track.tpcNSigmaKa();
        break;
      case PartID::pr:
        nsigmaTpc = track.tpcNSigmaPr();
        break;
      case PartID::de:
        nsigmaTpc = track.tpcNSigmaDe();
        break;
      case PartID::he:
        nsigmaTpc = computeNSigmaTPCHe3(track);
        break;
      default:
        nsigmaTpc = -999.f;
        break;
    }

    float nsigmaTof = -999.f;
    switch (partID) {
      case PartID::el:
        nsigmaTof = track.tofNSigmaEl();
        break;
      case PartID::de:
        nsigmaTof = track.tofNSigmaDe();
        break;
      default:
        nsigmaTof = -999.f;
        break;
    }

    float massTof = -999.f;
    if (track.hasTOF()) {
      switch (partID) {
        case PartID::de:
          massTof = computeTOFmassDe<isMC>(track);
          break;
        case PartID::he:
          massTof = computeTOFmassHe3<isMC>(track);
          break;
        default:
          massTof = -999.f;
          break;
      }
    }

    float nsigmaIts = -999.f;
    switch (partID) {
      case PartID::el:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Electron>(track);
        break;
      case PartID::pi:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Pion>(track);
        break;
      case PartID::ka:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Kaon>(track);
        break;
      case PartID::pr:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Proton>(track);
        break;
      case PartID::de:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Deuteron>(track);
        break;
      case PartID::he:
        nsigmaIts = m_responseITS.nSigmaITS<o2::track::PID::Helium3>(track);
        break;
      default:
        nsigmaIts = -999.f;
        break;
    }

    float correctedTpcInnerParam = track.tpcInnerParam();
    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    correctedTpcInnerParam = (partID == PartID::he && he3setting_compensatePIDinTracking && heliumPID) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();

    m_hAnalysis.fill(HIST(cNames[partID]) + HIST("/nSigmaTPC"), track.p() * track.sign(), nsigmaTpc);
    m_hAnalysis.fill(HIST(cNames[partID]) + HIST("/nSigmaITS"), track.p() * track.sign(), nsigmaIts);
    m_hAnalysis.fill(HIST(cNames[partID]) + HIST("/nSigmaTOF"), track.p() * track.sign(), nsigmaTof);
    if (partID == static_cast<int>(PartID::de) || partID == static_cast<int>(PartID::he))
      m_hAnalysis.fill(HIST(cNames[partID]) + HIST("/TOFmass"), track.p() * track.sign(), massTof);
    m_hAnalysis.fill(HIST(cNames[partID]) + HIST("/pmatching"), correctedTpcInnerParam * track.sign(), (correctedTpcInnerParam - track.p()) / correctedTpcInnerParam);
  }

  template <typename McPart>
  void fillMcHistogramsV0(aod::V0Datas::iterator const& v0, const McPart& posDaughter, const McPart& negDaughter)
  {
    if ((std::abs(posDaughter.pdgCode()) != PDG_t::kProton && std::abs(posDaughter.pdgCode()) != PDG_t::kPiPlus) ||
        (std::abs(negDaughter.pdgCode()) != PDG_t::kProton && std::abs(negDaughter.pdgCode()) != PDG_t::kPiPlus)) {
      return;
    }

    int motherPdgCode = 0;
    for (const auto& posMother : posDaughter.template mothers_as<aod::McParticles>()) {
      for (const auto& negMother : negDaughter.template mothers_as<aod::McParticles>()) {
        if (negMother.globalIndex() == posMother.globalIndex()) {
          motherPdgCode = posMother.pdgCode();
        }
      }
    }

    if (motherPdgCode == PDG_t::kLambda0) {
      m_hAnalysis.fill(HIST("massLambdaMc"), v0.pt(), v0.mLambda());
    } else if (motherPdgCode == PDG_t::kLambda0Bar) {
      m_hAnalysis.fill(HIST("massLambdaMc"), v0.pt() * -1.f, v0.mAntiLambda());
    }
  }

  template <typename McPart>
  void fillMcHistogramsCascade(aod::CascDatas::iterator const& cascade,
                               const McPart& bachelorDaughter, const McPart& posV0Daughter)
  {
    McPart v0Daughter;
    for (const auto& iterV0Daughter : posV0Daughter.template mothers_as<aod::McParticles>()) {
      if (std::abs(iterV0Daughter.pdgCode()) != PDG_t::kLambda0) {
        continue;
      }
      v0Daughter = iterV0Daughter;
    }
    if (std::abs(bachelorDaughter.pdgCode()) != PDG_t::kKPlus) {
      return;
    }

    int motherPdgCode = 0;
    for (const auto& bachelorMother : bachelorDaughter.template mothers_as<aod::McParticles>()) {
      for (const auto& v0Mother : v0Daughter.template mothers_as<aod::McParticles>()) {
        if (v0Mother.globalIndex() == bachelorMother.globalIndex()) {
          motherPdgCode = bachelorMother.pdgCode();
        }
      }
    }

    if (motherPdgCode == PDG_t::kOmegaMinus) {
      m_hAnalysis.fill(HIST("massOmegaMc"), cascade.pt(), cascade.mOmega());
    } else if (motherPdgCode == -PDG_t::kOmegaMinus) {
      m_hAnalysis.fill(HIST("massOmegaMc"), cascade.pt() * -1.f, cascade.mOmega());
    }
  }

  template <bool isMC = false>
  void fillTable(const Candidate& candidate)
  {
    m_ClusterStudiesTable(
      candidate.p, candidate.eta, candidate.phi,
      candidate.itsClusterSize, static_cast<uint8_t>(candidate.partID));
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        candidate.pTPC, candidate.pidInTrk,
        candidate.nsigmaTPC, candidate.nsigmaTOF, candidate.tofMass,
        candidate.cosPAMother, candidate.massMother);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }

    if constexpr (isMC) {
      m_ClusterStudiesTableMc(
        candidate.pdgCode);
    }
  }

  // =========================================================================================================

  template <bool isMC = false, typename T>
  float computeTOFmassDe(const T& candidate)
  {
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    return candidate.tpcInnerParam() * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  template <typename Track>
  float computeNSigmaTPCHe3(const Track& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && he3setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), m_BBparamsHe[0], m_BBparamsHe[1], m_BBparamsHe[2], m_BBparamsHe[3], m_BBparamsHe[4]);
    double resoTPC{expTPCSignal * m_BBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <bool isMC = false, typename Track>
  float computeTOFmassHe3(const Track& candidate)
  {
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && he3setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    return correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  // =========================================================================================================

  template <typename Track>
  bool electronPrimarySelection(const Track& track)
  {
    m_hAnalysis.fill(HIST("electron/DCAxyBeforeSelection"), track.dcaXY());
    m_hAnalysis.fill(HIST("electron/DCAzBeforeSelection"), track.dcaZ());

    if (track.dcaXY() > electronsetting_maxDcaxy ||
        track.dcaZ() > electronsetting_maxDcaz) {
      return false;
    }
    return true;
  }

  template <typename Track>
  bool electronPidSelection(const Track& track)
  {
    if (track.tpcNSigmaEl() < electronsetting_minNsigmatpcEl ||
        electronsetting_maxNsigmatpcEl < track.tpcNSigmaEl() ||
        std::abs(track.tpcNSigmaPi()) < electronsetting_maxNsigmatpcPi ||
        std::abs(track.tpcNSigmaKa()) < electronsetting_maxNsigmatpcKa ||
        std::abs(track.tpcNSigmaPr()) < electronsetting_maxNsigmatpcPr)
      return false;

    if (electronsetting_maxNsigmatofEl != 0 &&
        std::abs(track.tofNSigmaEl()) < electronsetting_maxNsigmatofEl)
      return false;

    return true;
  }

  // =========================================================================================================

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (m_runNumber == bc.runNumber()) {
      return;
    }

    const auto& timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;

    auto grpmagPath{"GLO/Config/GRPMagField"};
    grpmag = m_ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    m_d_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << m_d_bz << " kG";
    m_runNumber = bc.runNumber();

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void init(o2::framework::InitContext&)
  {
    m_runNumber = 0;
    m_d_bz = 0;

    m_ccdb->setURL("http://alice-ccdb.cern.ch");
    m_ccdb->setCaching(true);
    m_ccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    for (int ipartid = 0; ipartid < static_cast<int>(PartID::all); ipartid++) {
      if (ipartid == 0)
        continue;

      m_hAnalysis.add(fmt::format("{}/nSigmaITS", cNames[ipartid]).c_str(), (fmt::format("nSigma ITS {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{ITS}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      m_hAnalysis.add(fmt::format("{}/nSigmaTPC", cNames[ipartid]).c_str(), (fmt::format("nSigma TPC {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{TPC}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}});
      m_hAnalysis.add(fmt::format("{}/nSigmaTOF", cNames[ipartid]).c_str(), (fmt::format("nSigma TOF {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{TOF}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}});
      if (ipartid == static_cast<int>(PartID::de) || ipartid == static_cast<int>(PartID::he)) {
        m_hAnalysis.add(fmt::format("{}/TOFmass", cNames[ipartid]).c_str(), (fmt::format("TOF mass {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); #it{m}_{TOF} (GeV/#it{c}^{2})")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, 1.0f, 5.0f}});
        m_hAnalysis.add(fmt::format("{}/trackSelections", cNames[ipartid]).c_str(), (fmt::format("track selections {};", cNames[ipartid]) + std::string("Selections; Counts")).c_str(), HistType::kTH1F, {{NucleiSelections::kNucleiAll, -0.5, static_cast<double>(NucleiSelections::kNucleiAll) - 0.5}});
      }
      m_hAnalysis.add(fmt::format("{}/pmatching", cNames[ipartid]).c_str(), (fmt::format("p matching {};", cNames[ipartid]) + std::string("signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}});
    }

    std::vector<std::string> trackSelectionLabels = {"All", "n clusters ITS", "TPC", "TOF"};
    for (int i = 0; i < NucleiSelections::kNucleiAll; i++) {
      m_hAnalysis.get<TH1>(HIST(cNames[static_cast<int>(PartID::de)]) + HIST("/trackSelections"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionLabels[i].c_str());
      m_hAnalysis.get<TH1>(HIST(cNames[static_cast<int>(PartID::he)]) + HIST("/trackSelections"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionLabels[i].c_str());
    }

    LOG(info) << "Bethe-Bloch parameters for He3:";
    for (int i = 0; i < 5; i++) {
      m_BBparamsHe[i] = setting_BetheBlochParams->get("He3", Form("p%i", i));
      LOG(info) << "p" << i << ": " << m_BBparamsHe[i];
    }
    m_BBparamsHe[5] = setting_BetheBlochParams->get("He3", "resolution");
    LOG(info) << "resolution: " << m_BBparamsHe[5];

    LOG(info) << "Bethe-Bloch parameters for De:";
    for (int i = 0; i < 5; i++) {
      m_BBparamsDe[i] = setting_BetheBlochParams->get("De", Form("p%i", i));
      LOG(info) << "p" << i << ": " << m_BBparamsDe[i];
    }
    m_BBparamsDe[5] = setting_BetheBlochParams->get("De", "resolution");
    LOG(info) << "resolution: " << m_BBparamsDe[5];

    std::vector<std::string> collision_selection_labels = {"All", "sel8", "z_{VTX} < 10 cm"};
    for (int i = 0; i < Selections::kAll; i++)
      m_hAnalysis.get<TH1>(HIST("collision_selections"))->GetXaxis()->SetBinLabel(i + 1, collision_selection_labels[i].c_str());

    std::vector<std::string> V0_selection_labels = {"All", "daughter track quality", "V0 topology", "V0 mass selection"};
    for (int i = 0; i < V0Selections::kV0All; i++)
      m_hAnalysis.get<TH1>(HIST("v0_selections"))->GetXaxis()->SetBinLabel(i + 1, V0_selection_labels[i].c_str());

    std::vector<std::string> Casc_selection_labels = {"All", "Topology", "Veto Xi", "Accepted Omega", "n#sigma_{TPC} K"};
    for (int i = 0; i < CascSelections::kCascAll; i++)
      m_hAnalysis.get<TH1>(HIST("casc_selections"))->GetXaxis()->SetBinLabel(i + 1, Casc_selection_labels[i].c_str());

    std::vector<std::string> E_selections_labels = {"All", "Track quality", "Primary", "Pid", "#pi^{0}"};
    for (int i = 0; i < ESelections::kEAll; i++)
      m_hAnalysis.get<TH1>(HIST("e_selections"))->GetXaxis()->SetBinLabel(i + 1, E_selections_labels[i].c_str());

    std::vector<std::string> V0Type_labels = {"K0s", "#Lambda", "#bar{#Lambda}", "Photon"};
    for (int i = 0; i < V0Type::V0TypeAll; i++)
      m_hAnalysis.get<TH1>(HIST("v0_type"))->GetXaxis()->SetBinLabel(i + 1, V0Type_labels[i].c_str());
  }

  template <bool isMC = false, typename Tracks>
  void fillV0Cand(const std::array<float, 3>& /*pv*/, const aod::V0Datas::iterator& v0, const Tracks&)
  {
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0NoCut);

    const auto& posTrack = v0.template posTrack_as<Tracks>();
    const auto& negTrack = v0.template negTrack_as<Tracks>();

    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack))
      return;
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0DaughterQuality);

    if (!qualitySelectionV0(v0))
      return;
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0Topology);

    std::array<float, 3> momPos{v0.pxpos(), v0.pypos(), v0.pzpos()},
      momNeg{v0.pxneg(), v0.pyneg(), v0.pzneg()};

    m_hAnalysis.fill(HIST("armenteros_plot_before_selections"), v0.alpha(), v0.qtarm());
    m_hAnalysis.fill(HIST("Lambda_vs_K0s"), v0.mK0Short(), v0.mAntiLambda());

    uint8_t v0Bitmask = selectV0MotherHypothesis(v0);
    if (v0Bitmask == 0 || (v0Bitmask & (v0Bitmask - 1)) != 0)
      return;
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0PID);

    Candidate candidatePos(std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(),
                           RecoDecay::eta(momPos), RecoDecay::phi(momPos), posTrack.itsClusterSizes(),
                           0, posTrack.tpcInnerParam() * posTrack.sign(), posTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);
    Candidate candidateNeg(std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(),
                           RecoDecay::eta(momNeg), RecoDecay::phi(momNeg), negTrack.itsClusterSizes(),
                           0, negTrack.tpcInnerParam() * negTrack.sign(), negTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);

    if (!selectPidV0Daughters(candidatePos, candidateNeg, v0, posTrack, negTrack, v0Bitmask))
      return;

    const float massV0 = fillHistogramsV0(v0, posTrack, negTrack, v0Bitmask);
    candidatePos.massMother = massV0;
    candidateNeg.massMother = massV0;

    if (!setting_fillV0)
      return;

    if constexpr (isMC) { // MC
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle())
        return;

      const auto& posMcParticle = posTrack.mcParticle();
      const auto& negMcParticle = negTrack.mcParticle();

      candidatePos.pdgCode = posMcParticle.pdgCode();
      candidateNeg.pdgCode = negMcParticle.pdgCode();

      fillMcHistogramsV0(v0, posMcParticle, negMcParticle);
    }

    fillTable<isMC>(candidatePos);
    fillTable<isMC>(candidateNeg);

    m_hAnalysis.fill(HIST("isPositive"), true);
    m_hAnalysis.fill(HIST("isPositive"), false);
  }

  template <bool isMC = false, typename Track>
  void fillKCand(const std::array<float, 3>& pv, aod::CascDatas::iterator const& cascade, const Track&)
  {
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kCascNoCut);
    const auto& bachelorTrack = cascade.template bachelor_as<Track>();

    std::array<float, 3> momBachelor{cascade.pxbach(), cascade.pybach(), cascade.pzbach()};

    if (!qualitySelectionCascade(cascade, pv)) {
      return;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kCascTopology);

    const float& massXi = cascade.mXi();
    const float& massOmega = cascade.mOmega();
    m_hAnalysis.fill(HIST("Xi_vs_Omega"), massOmega, massXi);

    if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsetting_massWindowXi) {
      return;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kRejectedXi);

    m_hAnalysis.fill(HIST("massOmega"), cascade.pt() * bachelorTrack.sign(), massOmega);
    if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsetting_massWindowOmega) {
      return;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kAcceptedOmega);

    if (std::abs(bachelorTrack.tpcNSigmaKa()) > cascsetting_nsigmatpc) {
      return;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kNSigmaTPC);
    fillHistogramsParticle<PartID::ka, isMC>(bachelorTrack);

    m_ClusterStudiesTable(
      std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(),
      RecoDecay::eta(momBachelor), RecoDecay::phi(momBachelor),
      bachelorTrack.itsClusterSizes(), PartID::ka);
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        bachelorTrack.tpcInnerParam() * bachelorTrack.sign(),
        bachelorTrack.pidForTracking(),
        bachelorTrack.tpcNSigmaKa(), /*TOF nsigma*/ -999.f, /*TOF mass*/ -999.f,
        cascade.casccosPA(pv[0], pv[1], pv[2]), massOmega);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }

    if constexpr (isMC) {
      if (!bachelorTrack.has_mcParticle()) {
        return;
      }
      const auto& mcParticle = bachelorTrack.mcParticle();

      m_ClusterStudiesTableMc(
        mcParticle.pdgCode());

      const auto& posV0Daughter = cascade.template posTrack_as<Track>();
      if (!posV0Daughter.has_mcParticle()) {
        return;
      }
      const auto& mcPosParticleV0 = posV0Daughter.mcParticle();
      fillMcHistogramsCascade(cascade, mcParticle, mcPosParticleV0);
    }

    m_hAnalysis.fill(HIST("isPositive"), bachelorTrack.sign() > 0);
  }

  template <int partID, bool isMC = false, typename Track>
  void fillNucleusTable(const Track& track)
  {
    constexpr int kPartID = partID;

    if (kPartID == static_cast<int>(PartID::de) && track.sign() > 0)
      return;
    m_hAnalysis.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiNoCut);

    if (track.itsNCls() < desetting_nClsIts)
      return;
    m_hAnalysis.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiNClsIts);

    const float tpcNsigma = kPartID == static_cast<int>(PartID::de) ? track.tpcNSigmaDe() : computeNSigmaTPCHe3(track);
    const float tpcNsigmaMax = kPartID == static_cast<int>(PartID::de) ? desetting_nsigmatpc : he3setting_nsigmatpc;
    if (std::abs(tpcNsigma) > tpcNsigmaMax)
      return;
    m_hAnalysis.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiPIDtpc);

    const float tofMass = track.hasTOF() ? (kPartID == static_cast<int>(PartID::de) ? computeTOFmassDe<isMC>(track) : computeTOFmassHe3<isMC>(track)) : -999.f;
    const float tofNsigma = kPartID == static_cast<int>(PartID::de) ? track.tofNSigmaDe() : -999.f;
    float correctedTPCinnerParam = track.tpcInnerParam();

    if (kPartID == static_cast<int>(PartID::de)) {
      if (!track.hasTOF() || std::abs(tofNsigma) > desetting_nsigmatof)
        return;
    } else {
      if (track.hasTOF() && (tofMass < he3setting_tofmasslow || tofMass > he3setting_tofmasshigh))
        return;

      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      correctedTPCinnerParam = (heliumPID && he3setting_compensatePIDinTracking) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();
    }
    m_hAnalysis.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiPIDtof);

    fillHistogramsParticle<kPartID, isMC>(track);

    m_ClusterStudiesTable(
      track.p() * track.sign(), track.eta(), track.phi(),
      track.itsClusterSizes(), kPartID);
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        correctedTPCinnerParam * track.sign(), track.pidForTracking(),
        tpcNsigma, tofNsigma, tofMass,
        /*cosPA*/ -999.f, /*mass mother*/ -999.f);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }

    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return;
      }

      const auto& mcParticle = track.mcParticle();

      m_ClusterStudiesTableMc(
        mcParticle.pdgCode());
    }

    m_hAnalysis.fill(HIST("isPositive"), track.sign() > 0);
  }

  void fillPKPiTable(const TracksFullIU::iterator& track)
  {
    uint8_t partID = 0;
    float tpcNSigma = 0.f;
    if (std::abs(track.tpcNSigmaPi()) < v0setting_nsigmatpcPi && std::abs(track.tpcNSigmaKa()) > 3) {
      partID = PartID::pi;
      tpcNSigma = track.tpcNSigmaPi();
      m_hAnalysis.fill(HIST("nSigmaTPCPi"), track.p() * track.sign(), tpcNSigma);
    } else if (std::abs(track.tpcNSigmaKa()) < cascsetting_nsigmatpc && (std::abs(track.tpcNSigmaPi()) > 3 /*&& std::abs(track.tpcNSigmaPr()) > 3*/)) {
      partID = PartID::ka;
      tpcNSigma = track.tpcNSigmaKa();
      m_hAnalysis.fill(HIST("nSigmaTPCKa"), track.p() * track.sign(), tpcNSigma);
    } else if (std::abs(track.tpcNSigmaPr()) < v0setting_nsigmatpcPr && std::abs(track.tpcNSigmaKa()) > 3) {
      partID = PartID::pr;
      tpcNSigma = track.tpcNSigmaPr();
      m_hAnalysis.fill(HIST("nSigmaTPCPr"), track.p() * track.sign(), tpcNSigma);
    } else {
      return;
    }

    m_ClusterStudiesTable(
      track.p() * track.sign(), track.eta(), track.phi(),
      track.itsClusterSizes(), partID);
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        track.tpcInnerParam() * track.sign(), track.pidForTracking(),
        tpcNSigma, /*TOF nsigma*/ -999.f, /*TOF mass*/ -999.f,
        /*cosPA*/ -999.f, /*mass mother*/ -999.f);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }
  }

  template <const bool isMC = false, typename Track>
  void fillElectronTable(const Track& posTrack, const Track& negTrack)
  {
    m_hAnalysis.fill(HIST("e_selections"), ESelections::kENoCut);
    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack)) {
      return;
    }
    m_hAnalysis.fill(HIST("e_selections"), ESelections::kETrackQuality);

    if (!electronPrimarySelection(posTrack) || !electronPrimarySelection(negTrack)) {
      return;
    }
    m_hAnalysis.fill(HIST("e_selections"), ESelections::kEPrimary);

    if (!electronPidSelection(posTrack) || !electronPidSelection(negTrack)) {
      return;
    }
    m_hAnalysis.fill(HIST("e_selections"), ESelections::kEPid);

    const float invariantMass = std::sqrt(RecoDecay::m2<2>(std::array<std::array<float, 3>, 2>{
                                                             std::array<float, 3>{posTrack.px(), posTrack.py(), posTrack.pz()},
                                                             std::array<float, 3>{negTrack.px(), negTrack.py(), negTrack.pz()}},
                                                           std::array<float, 2>{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}));

    m_hAnalysis.fill(HIST("massPi0WithBkg"), invariantMass);
    if (invariantMass > o2::constants::physics::MassPi0) {
      return;
    }
    m_hAnalysis.fill(HIST("e_selections"), ESelections::kEPi0);
    m_hAnalysis.fill(HIST("massPi0"), invariantMass);
    fillHistogramsParticle<PartID::el, isMC>(posTrack);
    fillHistogramsParticle<PartID::el, isMC>(negTrack);

    // float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz(), t1.sign(), t2.sign(), m_d_bz);
    // float opangle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

    m_ClusterStudiesTable(
      posTrack.p(), posTrack.eta(), posTrack.phi(),
      posTrack.itsClusterSizes(), PartID::el);
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        posTrack.tpcInnerParam(), posTrack.pidForTracking(),
        posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), -999.f /* TofMass */,
        -999.f /* cosPA */, invariantMass);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }

    m_ClusterStudiesTable(
      negTrack.p(), negTrack.eta(), negTrack.phi(),
      negTrack.itsClusterSizes(), PartID::el);
    if (setting_fillExtraTable) {
      m_ClusterStudiesTableExtra(
        negTrack.tpcInnerParam(), negTrack.pidForTracking(),
        negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), -999.f /* TofMass */,
        -999.f /* cosPA */, invariantMass);
    }
    if (setting_fillCollTable) {
      m_ClusterStudiesTableCollision(
        m_runNumber);
    }

    if constexpr (isMC) {
      const auto& posMcParticle = posTrack.mcParticle();
      const auto& negMcParticle = negTrack.mcParticle();

      m_ClusterStudiesTableMc(
        posMcParticle.pdgCode());
      m_ClusterStudiesTableMc(
        negMcParticle.pdgCode());

      if (!posMcParticle.has_mothers() || !negMcParticle.has_mothers())
        return;

      for (const auto& posMother : posMcParticle.template mothers_as<aod::McParticles>()) {
        for (const auto& negMother : negMcParticle.template mothers_as<aod::McParticles>()) {
          if (posMother.globalIndex() != negMother.globalIndex() || std::abs(posMother.pdgCode()) != PDG_t::kPi0)
            return;
          m_hAnalysis.fill(HIST("massPi0Mc"), std::sqrt((posMcParticle.e() + negMcParticle.e()) * (posMcParticle.e() + negMcParticle.e()) -
                                                        (posMcParticle.p() + negMcParticle.p()) * (posMcParticle.p() + posMcParticle.p())));
          break;
        }
      }
    }
  }

  // =========================================================================================================

  void processDataV0Casc(CollisionsCustom::iterator const& collision /*s*/, TracksFullIU const& tracks, aod::V0Datas const& v0s, aod::CascDatas const& cascades, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillV0Cand</*isMC*/ false>(PV, v0, tracks);
    }

    for (const auto& cascade : cascades) {
      fillKCand</*isMC*/ false>(PV, cascade, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataV0Casc, "process Data V0 and cascade", false);

  Partition<TracksFullIU> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<TracksFullIU> negTracks = o2::aod::track::signed1Pt < 0.f;
  void processDataElectrons(CollisionsCustom::iterator const& collision, TracksFullIU const& /*tracks*/, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

    const auto& posTracks_thisCollision = posTracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), m_cache);
    const auto& negTracks_thisCollision = negTracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), m_cache);

    for (const auto& [posTrack, negTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_thisCollision, negTracks_thisCollision))) {
      fillElectronTable(posTrack, negTrack);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataElectrons, "process Data Electrons", false);

  void processDataNuclei(CollisionsCustom::iterator const& collision, TracksFullIU const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

    for (const auto& track : tracks) {
      if (!qualityTrackSelection(track)) {
        return;
      }

      fillNucleusTable<PartID::de, /*isMC*/ false>(track);
      fillNucleusTable<PartID::he, /*isMC*/ false>(track);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataNuclei, "process Data Nuclei", false);

  /**
   * @brief Produce a dataset with high purity p, K, #pi
   */
  void processDataPKPi(CollisionsCustom::iterator const& collision, TracksFullIU const& tracks)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

    for (const auto& track : tracks) {
      if (!qualityTrackSelection(track)) {
        return;
      }

      fillPKPiTable(track);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataPKPi, "process Data p, K, pi", false);

  void processMcV0Casc(CollisionsCustomMc::iterator const& collision /*s*/, TracksFullIUMc const& tracks, aod::V0Datas const& v0s, aod::CascDatas const& cascades, aod::BCsWithTimestamps const&, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillV0Cand</*isMC*/ true>(PV, v0, tracks);
    }

    for (const auto& cascade : cascades) {
      fillKCand</*isMC*/ true>(PV, cascade, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcV0Casc, "process Mc V0 and cascade", false);

  Partition<TracksFullIUMc> posTracksMc = o2::aod::track::signed1Pt > 0.f;
  Partition<TracksFullIUMc> negTracksMc = o2::aod::track::signed1Pt < 0.f;
  void processMcElectrons(CollisionsCustomMc::iterator const& collision, TracksFullIUMc const& /*tracks*/, aod::BCsWithTimestamps const&, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

    const auto& posTracks_thisCollision = posTracksMc.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), m_cache);
    const auto& negTracks_thisCollision = negTracksMc.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), m_cache);

    for (const auto& [posTrack, negTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_thisCollision, negTracks_thisCollision))) {
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle())
        continue;

      fillElectronTable</*isMC*/ true>(posTrack, negTrack);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcElectrons, "process Mc Electrons", false);

  void processMcNuclei(CollisionsCustomMc::iterator const& collision, TracksFullIUMc const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collisionSelection(collision)) {
      return;
    }
    const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

    for (const auto& track : tracks) {
      if (!qualityTrackSelection(track)) {
        return;
      }

      fillNucleusTable<PartID::de, /*isMC*/ true>(track);
      fillNucleusTable<PartID::he, /*isMC*/ true>(track);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcNuclei, "process Mc Nuclei", false);

}; // LfTreeCreatorClusterStudies

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorClusterStudies>(cfgc)};
}
