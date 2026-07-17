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
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH1.h>
#include <TPDGCode.h>
#include <TString.h>

#include <fmt/format.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
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

constexpr std::array<std::array<double, 6>, 1> defaultParams = {{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}}};
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

static constexpr std::array<std::string_view, 8> cNames = {"none", "electron", "pion", "kaon", "proton", "deuteron", "He3"};

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

  Service<o2::ccdb::BasicCCDBManager> mccdb;
  SliceCache mcache;
  int mrunNumber = 0;
  int mcollisionCounter = 0;
  float md_bz = 0.f;
  uint32_t mrandomSeed = 0.;

  Configurable<bool> settingFillExtraTable{"fillExtraTable", false, "Fill the extra table"};
  Configurable<bool> settingFillCollTable{"fillCollTable", false, "Fill the collision table"};

  Configurable<float> settingZvtxMax{"zVtxMax", 10.f, "Maximum z vertex position"};

  // Configurable<float> settingDownscaleFactor{"downscaleFactor", 1.f, "Downscale factor for the V0 candidates"};

  Configurable<float> trackNclsItsMin{"trackNclsItsMin", 0.f, "Minimum number of ITS clusters for the V0 daughters"};
  Configurable<float> trackNclsTpcMin{"trackNclsTpcMin", 100.f, "Minimum number of TPC clusters for the V0 daughters"};
  Configurable<float> trackNclsTpcMaxShared{"trackNclsTpcMaxShared", 5.f, "Maximum number of shared TPC clusters for the V0 daughters"};
  Configurable<float> trackEtaMax{"etaMax", 0.8f, "Maximum eta"};
  Configurable<float> trackTpcChi2Min{"trackTpcChi2Min", 0.5f, "Minimum TPC chi2 per cluster"};

  // Configurable<float> v0settingEtaMaxV0{"etaMaxV0", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0settingDcaV0daughters{"v0settingDcaV0daughters", 1.f, "DCA between the V0 daughters"};
  Configurable<float> v0settingDcaMinV0DaughterToPv{"v0settingDcaMinV0DaughterToPv", 0.06f, "DCA of the daughters to the primary vertex"};
  Configurable<float> v0settingRadiusV0{"v0settingRadiusV0", 0.5f, "Maximum radius of the V0 accepted"};
  Configurable<float> v0settingCosPA{"v0settingCosPA", 0.98f, "Cosine of the pointing angle of the V0"};
  Configurable<float> v0settingMassWindowLambda{"v0settingMassWindowLambda", 0.02f, "Mass window for the Lambda"};
  Configurable<float> v0settingMassWindowK0s{"v0settingMassWindowK0s", 0.02f, "Mass window for the K0s"};
  Configurable<float> v0settingNsigmatpcPi{"v0settingNsigmaTPCPi", 2.f, "Number of sigmas for the TPC PID for pions"};
  Configurable<float> v0settingNsigmatpcPr{"v0settingNsigmaTPCPr", 2.f, "Number of sigmas for the TPC PID for protons"};
  Configurable<float> lambdasettingQtAPcut{"lambdasettingQtAPcut", 0.02f, "Cut on the qt for the Armenteros-Podolanski plot for photon rejection"};
  Configurable<float> lambdasettingPmin{"lambdasettingPmin", 0.0f, "Minimum momentum for the V0 daughters"};

  Configurable<float> cascsettingDcaMinV0DaughterToPv{"cascsettingDcaMinV0DaughterToPv", 0.03f, "DCA of one of the daugthers of Lambda to pv"};
  Configurable<float> cascsettingDcaMinProtonToPv{"cascsettingDcaMinProtonToPv", 0.03f, "DCA of the proton coming from Lambda to pv"};
  Configurable<float> cascsettingDcaMinBachelorToPv{"cascsettingDcaMinBachelorToPv", 0.04f, "DCA of the bachelor to pv"};
  Configurable<float> cascsettingDcaMinV0ToPv{"cascsettingDcaMinV0ToPv", 0.04f, "DCA of the V0 to pv"};
  Configurable<float> cascsettingDcaV0Daughters{"cascsettingDcaV0Daughters", 0.4f, "DCA between the V0 daughters"};
  Configurable<float> cascsettingCascCosPA{"cascsettingCascCosPA", 0.99f, "Minimum cCosine of the pointing angle of the cascade"};
  Configurable<float> cascsettingV0CosPA{"cascsettingV0CosPA", 0.97f, "Minimum cCosine of the pointing angle of the v0"};
  Configurable<float> cascsettingDcaMaxBachelorToV0{"cascsettingDcaMaxBachelorToV0", 0.8f, "DCA of the bachelor to V0"};
  // Configurable<float> cascsettingDcaMinBachelorToProton{"cascsettingDcaMinBachelorToProton", 0.015f, "DCA of the bachelor to proton"};
  Configurable<float> cascsettingRadiusV0{"cascsettingRadiusV0", 1.2f, "Minimum radius of the V0 accepted"};
  Configurable<float> cascsettingRadiusCasc{"cascsettingRadiusCasc", 0.5f, "Minimum radius of the cascade accepted"};

  Configurable<float> cascsettingMassWindowOmega{"cascsettingMassWindowOmega", 0.01f, "Mass window for the Omega"};
  Configurable<float> cascsettingMassWindowXi{"cascsettingMassWindowXi", 0.01f, "Mass window for the Xi"};
  Configurable<float> cascsettingNsigmatpc{"cascsettingNsigmaTPC", 3.f, "Number of sigmas for the TPC PID"};

  Configurable<bool> electronsettingFromPhotonConversion{"electronsettingFromPhotonConversion", false, "Flag to indicate if the electron is from photon conversion"};
  Configurable<float> electronsettingConversion_rmin{"electronsettingConversion_rmin", 1.76f, "Minimum radius for the photon conversion (cm)"};
  Configurable<float> electronsettingConversion_rmax{"electronsettingConversion_rmax", 19.77f, "Maximum radius for the photon conversion (cm)"};
  Configurable<float> electronsettingMaxDcaxy{"electronsettingMaxDcaxy", 0.1f, "Maximum value for the DCAxy"};
  Configurable<float> electronsettingMaxDcaz{"electronsettingMaxDcaz", 0.5f, "Maximum value for the DCAz"};
  Configurable<float> electronsettingMinNsigmatpcEl{"electronsettingMinNsigmaTPCEl", -2.5f, "Minimum value for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsettingMaxNsigmatpcEl{"electronsettingMaxNsigmaTPCEl", 3.5f, "Maximum number for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsettingMaxNsigmatpcPi{"electronsettingMaxNsigmaTPCPi", 2.f, "Maximum number for the number of sigmas for pi rejection for the TPC PID for electrons"};
  Configurable<float> electronsettingMaxNsigmatpcKa{"electronsettingMaxNsigmaTPCKa", 2.f, "Maximum number for the number of sigmas for K rejection for the TPC PID for electrons"};
  Configurable<float> electronsettingMaxNsigmatpcPr{"electronsettingMaxNsigmaTPCPr", 2.f, "Maximum number for the number of sigmas for p rejection for the TPC PID for electrons"};
  Configurable<float> electronsettingMaxNsigmatofEl{"electronsettingMaxNsigmaTOFEl", 4.f, "Minimum value for the number of sigmas for the TPC PID for electrons"};
  Configurable<float> electronsettingMinPt{"electronsettingMinPt", 0.f, "Minimum pT accepted for electrons"};
  Configurable<float> electronsettingDcaV0daughters{"electronsettingDcaV0daughters", 1.f, "DCA between the V0 daughters"};
  Configurable<float> electronsettingDcaMinV0DaughterToPv{"electronsettingDcaMinV0DaughterToPv", 0.06f, "DCA of the daughters to the primary vertex"};
  Configurable<float> electronsettingRadiusV0{"electronsettingRadiusV0", 0.5f, "Maximum radius of the V0 accepted"};
  Configurable<float> electronsettingCosPA{"electronsettingCosPA", 0.98f, "Cosine of the pointing angle of the V0"};
  Configurable<float> electronsettingMaxMassPi0{"electronsettingMaxMassPi0", o2::constants::physics::MassPi0, "Maximum mass for the Pi0"};
  Configurable<float> electronsettingMinMassPi0{"electronsettingMinMassPi0", 0.08f, "Minimum mass for the Pi0"};
  Configurable<float> electronsettingMaxPhiv{"electronsettingMaxPhiv", 1.5708f, "Maximum value for the Phiv"};
  Configurable<float> electronsettingMinPhiv{"electronsettingMinPhiv", 0.f, "Minimum value for the Phiv"};
  Configurable<float> electronsettingMaxMassPhoton{"electronsettingMaxMassPhoton", 0.f, "Maximum mass for the Photon"};
  Configurable<float> electronsettingMinMassPhoton{"electronsettingMinMassPhoton", 0.005f, "Minimum mass for the Photon"};

  Configurable<int> desettingNclsIts{"desettingNclsIts", 6, "Minimum number of ITS clusters"};
  Configurable<float> desettingNsigmatpc{"desettingNsigmaCutTPC", 2.f, "Number of sigmas for the TPC PID"};
  Configurable<float> desettingNsigmatof{"desettingNsigmaCutTOF", 2.f, "Number of sigmas for the TOF PID"};
  Configurable<int> he3settingNclsIts{"he3settingNclsIts", 6, "Minimum number of ITS clusters"};
  Configurable<bool> he3settingCompensatePIDinTracking{"he3settingCompensatePIDinTracking", true, "Compensate PID in tracking"};
  Configurable<float> he3settingNsigmatpc{"he3settingNsigmaCutTPC", 2.f, "Number of sigmas for the TPC PID"};
  Configurable<float> he3settingTofmasslow{"he3settingTofmasslow", 1.8f, "Lower limit for the TOF mass"};
  Configurable<float> he3settingTofmasshigh{"he3settingTofmasshigh", 4.2f, "Upper limit for the TOF mass"};

  // Bethe Bloch parameters
  std::array<float, 6> mBBparamsDe = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}}, mBBparamsHe = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
  Configurable<LabeledArray<double>> settingBetheBlochParams{"settingBetheBlochParams", {BetheBloch::defaultParams[0].data(), 2, 6, {"De", "He3"}, BetheBloch::parNames}, "TPC Bethe-Bloch parameterisation for nuclei"};

  Preslice<aod::V0s> mPerCollisionV0 = o2::aod::v0::collisionId;
  Preslice<aod::Cascades> mPerCollisionCascade = o2::aod::cascade::collisionId;
  Preslice<TracksFullIU> mPerCol = aod::track::collisionId;
  Preslice<TracksFullIUMc> mPerColMC = aod::track::collisionId;

  HistogramRegistry mHistograms{
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

  Produces<o2::aod::ClStTable> mClusterStudiesTable;
  Produces<o2::aod::ClStTableExtra> mClusterStudiesTableExtra;
  Produces<o2::aod::ClStTableColl> mClusterStudiesTableCollision;
  Produces<o2::aod::ClStTableMc> mClusterStudiesTableMc;

  o2::aod::ITSResponse mresponseITS;

  template <typename Tcollision>
  bool collisionSelection(const Tcollision& collision)
  {
    mHistograms.fill(HIST("collision_selections"), Selections::kNoCut);
    if (!collision.sel8()) {
      return false;
    }
    mHistograms.fill(HIST("collision_selections"), Selections::kSel8);
    // if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
    //   return false;
    // }
    if (std::abs(collision.posZ()) > settingZvtxMax) {
      return false;
    }
    mHistograms.fill(HIST("collision_selections"), Selections::kVtxZ);
    return true;
  }

  // =========================================================================================================

  /**
   * Select the V0 daughters based on the quality cuts
   */
  template <typename Track>
  bool qualityTrackSelection(const Track& track)
  {
    return (std::abs(track.eta()) < trackEtaMax &&
            track.itsNCls() >= trackNclsItsMin &&
            track.tpcNClsFound() >= trackNclsTpcMin &&
            track.tpcNClsCrossedRows() >= trackNclsTpcMin &&
            track.tpcNClsShared() <= trackNclsTpcMaxShared &&
            track.tpcChi2NCl() >= trackTpcChi2Min &&
            track.tpcChi2NCl() <= 4.0f);
  }

  bool qualitySelectionV0(aod::V0Datas::iterator const& v0)
  {
    if (std::abs(v0.dcapostopv()) < v0settingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(v0.dcanegtopv()) < v0settingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(v0.dcaV0daughters()) > v0settingDcaV0daughters) {
      return false;
    }
    if (v0.v0radius() < v0settingRadiusV0) {
      return false;
    }
    if (std::abs(v0.v0cosPA()) < v0settingCosPA) {
      return false;
    }

    return true;
  }

  bool qualitySelectionPhotons(aod::V0Datas::iterator const& v0)
  {
    if (std::abs(v0.dcapostopv()) < electronsettingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(v0.dcanegtopv()) < electronsettingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(v0.dcaV0daughters()) > electronsettingDcaV0daughters) {
      return false;
    }
    if (v0.v0radius() < electronsettingRadiusV0) {
      return false;
    }
    if (std::abs(v0.v0cosPA()) < electronsettingCosPA) {
      return false;
    }
    if (v0.qtarm() > lambdasettingQtAPcut) {
      return false;
    }

    return true;
  }

  bool qualitySelectionCascade(aod::CascDatas::iterator const& cascade, const std::array<float, 3>& pv)
  {
    if (std::abs(cascade.dcapostopv()) < cascsettingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(cascade.dcanegtopv()) < cascsettingDcaMinV0DaughterToPv) {
      return false;
    }
    if (std::abs(cascade.dcabachtopv()) < cascsettingDcaMinBachelorToPv) {
      return false;
    }
    if (std::abs(cascade.dcav0topv(pv[0], pv[1], pv[2])) < cascsettingDcaMinV0ToPv) {
      return false;
    }
    if (std::abs(cascade.dcaV0daughters()) > cascsettingDcaV0Daughters) {
      return false;
    }
    if (std::abs(cascade.casccosPA(pv[0], pv[1], pv[2])) < cascsettingCascCosPA) {
      return false;
    }
    if (std::abs(cascade.v0cosPA(pv[0], pv[1], pv[2])) < cascsettingV0CosPA) {
      return false;
    }
    if (std::abs(cascade.v0radius()) < cascsettingRadiusV0) {
      return false;
    }
    if (std::abs(cascade.cascradius()) < cascsettingRadiusCasc) {
      return false;
    }

    mHistograms.fill(HIST("casc_dca_daughter_pairs"), cascade.dcacascdaughters());
    if (std::abs(cascade.dcacascdaughters()) > cascsettingDcaMaxBachelorToV0) {
      return false;
    }

    return true;
  }

  uint8_t selectV0MotherHypothesis(aod::V0Datas::iterator const& v0)
  {
    uint8_t v0Bitmask(0);
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0settingMassWindowK0s) {
      SETBIT(v0Bitmask, K0s);
    }
    if ((std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0settingMassWindowLambda) && (v0.alpha() > 0)) {
      SETBIT(v0Bitmask, Lambda);
    }
    if ((std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < v0settingMassWindowLambda) && (v0.alpha() < 0)) {
      SETBIT(v0Bitmask, AntiLambda);
    }
    return v0Bitmask;
  }

  template <typename Track>
  bool selectPidV0Daughters(Candidate& candidatePos, Candidate& candidateNeg,
                            aod::V0Datas::iterator const& v0, const Track& posTrack, const Track& negTrack, uint8_t v0Bitmask)
  {
    if (TESTBIT(v0Bitmask, Lambda)) {
      if (v0.qtarm() < lambdasettingQtAPcut) {
        return false;
      }
      if (std::abs(posTrack.tpcNSigmaPr()) > v0settingNsigmatpcPr || std::abs(negTrack.tpcNSigmaPi()) > v0settingNsigmatpcPi) {
        return false;
      }
      if (v0.p() < lambdasettingPmin) {
        return false;
      }
      candidatePos.partID = PartID::pr;
      candidateNeg.partID = PartID::pi;
      candidatePos.nsigmaTPC = posTrack.tpcNSigmaPr();
      candidateNeg.nsigmaTPC = negTrack.tpcNSigmaPi();
      mHistograms.fill(HIST("v0_type"), V0Type::Lambda);

    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      if (v0.qtarm() < lambdasettingQtAPcut) {
        return false;
      }
      if (std::abs(posTrack.tpcNSigmaPi()) > v0settingNsigmatpcPi || std::abs(negTrack.tpcNSigmaPr()) > v0settingNsigmatpcPr) {
        return false;
      }
      if (v0.p() < lambdasettingPmin) {
        return false;
      }
      candidatePos.partID = PartID::pi;
      candidateNeg.partID = PartID::pr;
      candidatePos.nsigmaTPC = posTrack.tpcNSigmaPi();
      candidateNeg.nsigmaTPC = negTrack.tpcNSigmaPr();
      mHistograms.fill(HIST("v0_type"), V0Type::AntiLambda);

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
      mHistograms.fill(HIST("massLambda"), v0.pt(), v0.mLambda());
      fillHistogramsParticle<PartID::pr, /*isMC*/ false>(trackPos);
      fillHistogramsParticle<PartID::pi, /*isMC*/ false>(trackNeg);
    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      massV0 = v0.mAntiLambda();
      mHistograms.fill(HIST("massLambda"), v0.pt() * -1.f, v0.mAntiLambda());
      fillHistogramsParticle<PartID::pi, /*isMC*/ false>(trackPos);
      fillHistogramsParticle<PartID::pr, /*isMC*/ false>(trackNeg);
    }

    mHistograms.fill(HIST("radiusV0"), v0.v0radius());
    mHistograms.fill(HIST("armenteros_plot_lambda"), v0.alpha(), v0.qtarm());
    mHistograms.fill(HIST("armenteros_plot"), v0.alpha(), v0.qtarm());

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
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Electron>(track);
        break;
      case PartID::pi:
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Pion>(track);
        break;
      case PartID::ka:
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Kaon>(track);
        break;
      case PartID::pr:
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Proton>(track);
        break;
      case PartID::de:
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Deuteron>(track);
        break;
      case PartID::he:
        nsigmaIts = mresponseITS.nSigmaITS<o2::track::PID::Helium3>(track);
        break;
      default:
        nsigmaIts = -999.f;
        break;
    }

    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    float correctedTpcInnerParam = (partID == PartID::he && he3settingCompensatePIDinTracking && heliumPID) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();

    mHistograms.fill(HIST(cNames[partID]) + HIST("/nSigmaTPC"), track.p() * track.sign(), nsigmaTpc);
    mHistograms.fill(HIST(cNames[partID]) + HIST("/nSigmaITS"), track.p() * track.sign(), nsigmaIts);
    mHistograms.fill(HIST(cNames[partID]) + HIST("/nSigmaTOF"), track.p() * track.sign(), nsigmaTof);
    if (partID == static_cast<int>(PartID::de) || partID == static_cast<int>(PartID::he))
      mHistograms.fill(HIST(cNames[partID]) + HIST("/TOFmass"), track.p() * track.sign(), massTof);
    mHistograms.fill(HIST(cNames[partID]) + HIST("/pmatching"), correctedTpcInnerParam * track.sign(), (correctedTpcInnerParam - track.p()) / correctedTpcInnerParam);
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
      mHistograms.fill(HIST("massLambdaMc"), v0.pt(), v0.mLambda());
    } else if (motherPdgCode == PDG_t::kLambda0Bar) {
      mHistograms.fill(HIST("massLambdaMc"), v0.pt() * -1.f, v0.mAntiLambda());
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
      mHistograms.fill(HIST("massOmegaMc"), cascade.pt(), cascade.mOmega());
    } else if (motherPdgCode == -PDG_t::kOmegaMinus) {
      mHistograms.fill(HIST("massOmegaMc"), cascade.pt() * -1.f, cascade.mOmega());
    }
  }

  template <bool isMC = false>
  void fillTable(const Candidate& candidate)
  {
    mClusterStudiesTable(
      candidate.p, candidate.eta, candidate.phi,
      candidate.itsClusterSize, static_cast<uint8_t>(candidate.partID));
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        candidate.pTPC, candidate.pidInTrk,
        candidate.nsigmaTPC, candidate.nsigmaTOF, candidate.tofMass,
        candidate.cosPAMother, candidate.massMother);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }

    if constexpr (isMC) {
      mClusterStudiesTableMc(
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
    float correctedTPCinnerParam = (heliumPID && he3settingCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::common::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), mBBparamsHe[0], mBBparamsHe[1], mBBparamsHe[2], mBBparamsHe[3], mBBparamsHe[4]);
    double resoTPC{expTPCSignal * mBBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <bool isMC = false, typename Track>
  float computeTOFmassHe3(const Track& candidate)
  {
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && he3settingCompensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    return correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  // =========================================================================================================

  template <typename Track>
  bool electronPrimarySelection(const Track& track)
  {
    mHistograms.fill(HIST("electron/DCAxyBeforeSelection"), track.dcaXY());
    mHistograms.fill(HIST("electron/DCAzBeforeSelection"), track.dcaZ());

    if (track.dcaXY() > electronsettingMaxDcaxy ||
        track.dcaZ() > electronsettingMaxDcaz) {
      return false;
    }
    return true;
  }

  template <typename Track>
  bool electronPidSelection(const Track& track)
  {
    if (track.tpcNSigmaEl() < electronsettingMinNsigmatpcEl ||
        electronsettingMaxNsigmatpcEl < track.tpcNSigmaEl() ||
        std::abs(track.tpcNSigmaPi()) < electronsettingMaxNsigmatpcPi ||
        std::abs(track.tpcNSigmaKa()) < electronsettingMaxNsigmatpcKa ||
        std::abs(track.tpcNSigmaPr()) < electronsettingMaxNsigmatpcPr) {
      return false;
    }

    if (electronsettingMaxNsigmatofEl > 1e-7 &&
        std::abs(track.tofNSigmaEl()) < electronsettingMaxNsigmatofEl) {
      return false;
    }

    return true;
  }

  // =========================================================================================================

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mrunNumber == bc.runNumber()) {
      return;
    }

    const auto& timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;

    auto grpmagPath{"GLO/Config/GRPMagField"};
    grpmag = mccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    md_bz = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << md_bz << " kG";
    mrunNumber = bc.runNumber();

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void init(o2::framework::InitContext&)
  {
    mrunNumber = 0;
    md_bz = 0;

    mccdb->setURL("http://alice-ccdb.cern.ch");
    mccdb->setCaching(true);
    mccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    for (int ipartid = 0; ipartid < static_cast<int>(PartID::all); ipartid++) {
      if (ipartid == 0) {
        continue;
      }

      mHistograms.add(fmt::format("{}/nSigmaITS", cNames[ipartid]).c_str(), (fmt::format("nSigma ITS {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{ITS}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -5.0f, 5.0f}});
      mHistograms.add(fmt::format("{}/nSigmaTPC", cNames[ipartid]).c_str(), (fmt::format("nSigma TPC {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{TPC}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}});
      mHistograms.add(fmt::format("{}/nSigmaTOF", cNames[ipartid]).c_str(), (fmt::format("nSigma TOF {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); n#sigma_{TOF}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}});
      if (ipartid == static_cast<int>(PartID::de) || ipartid == static_cast<int>(PartID::he)) {
        mHistograms.add(fmt::format("{}/TOFmass", cNames[ipartid]).c_str(), (fmt::format("TOF mass {};", cNames[ipartid]) + std::string("signed #it{p} (GeV/#it{c}); #it{m}_{TOF} (GeV/#it{c}^{2})")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, 1.0f, 5.0f}});
        mHistograms.add(fmt::format("{}/trackSelections", cNames[ipartid]).c_str(), (fmt::format("track selections {};", cNames[ipartid]) + std::string("Selections; Counts")).c_str(), HistType::kTH1F, {{NucleiSelections::kNucleiAll, -0.5, static_cast<double>(NucleiSelections::kNucleiAll) - 0.5}});
      }
      mHistograms.add(fmt::format("{}/pmatching", cNames[ipartid]).c_str(), (fmt::format("p matching {};", cNames[ipartid]) + std::string("signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}")).c_str(), HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}});
    }

    std::vector<std::string> trackSelectionLabels = {"All", "n clusters ITS", "TPC", "TOF"};
    for (int i = 0; i < NucleiSelections::kNucleiAll; i++) {
      mHistograms.get<TH1>(HIST(cNames[static_cast<int>(PartID::de)]) + HIST("/trackSelections"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionLabels[i].c_str());
      mHistograms.get<TH1>(HIST(cNames[static_cast<int>(PartID::he)]) + HIST("/trackSelections"))->GetXaxis()->SetBinLabel(i + 1, trackSelectionLabels[i].c_str());
    }

    LOG(info) << "Bethe-Bloch parameters for He3:";
    for (int i = 0; i < 5; i++) {
      mBBparamsHe[i] = settingBetheBlochParams->get("He3", Form("p%i", i));
      LOG(info) << "p" << i << ": " << mBBparamsHe[i];
    }
    mBBparamsHe[5] = settingBetheBlochParams->get("He3", "resolution");
    LOG(info) << "resolution: " << mBBparamsHe[5];

    LOG(info) << "Bethe-Bloch parameters for De:";
    for (int i = 0; i < 5; i++) {
      mBBparamsDe[i] = settingBetheBlochParams->get("De", Form("p%i", i));
      LOG(info) << "p" << i << ": " << mBBparamsDe[i];
    }
    mBBparamsDe[5] = settingBetheBlochParams->get("De", "resolution");
    LOG(info) << "resolution: " << mBBparamsDe[5];

    std::vector<std::string> collisionSelectionLabels = {"All", "sel8", "z_{VTX} < 10 cm"};
    for (int i = 0; i < Selections::kAll; i++)
      mHistograms.get<TH1>(HIST("collision_selections"))->GetXaxis()->SetBinLabel(i + 1, collisionSelectionLabels[i].c_str());

    std::vector<std::string> V0selectionLabels = {"All", "daughter track quality", "V0 topology", "V0 mass selection"};
    for (int i = 0; i < V0Selections::kV0All; i++)
      mHistograms.get<TH1>(HIST("v0_selections"))->GetXaxis()->SetBinLabel(i + 1, V0selectionLabels[i].c_str());

    std::vector<std::string> CascSelectionLabels = {"All", "Topology", "Veto Xi", "Accepted Omega", "n#sigma_{TPC} K"};
    for (int i = 0; i < CascSelections::kCascAll; i++)
      mHistograms.get<TH1>(HIST("casc_selections"))->GetXaxis()->SetBinLabel(i + 1, CascSelectionLabels[i].c_str());

    std::vector<std::string> ESelectionLabels = {"All", "Track quality", "Primary", "Pid", "#pi^{0}"};
    for (int i = 0; i < ESelections::kEAll; i++)
      mHistograms.get<TH1>(HIST("e_selections"))->GetXaxis()->SetBinLabel(i + 1, ESelectionLabels[i].c_str());

    std::vector<std::string> V0TypeLabels = {"K0s", "#Lambda", "#bar{#Lambda}", "Photon"};
    for (int i = 0; i < V0Type::V0TypeAll; i++)
      mHistograms.get<TH1>(HIST("v0_type"))->GetXaxis()->SetBinLabel(i + 1, V0TypeLabels[i].c_str());
  }

  template <bool isMC = false, typename Tracks>
  void fillV0Cand(const std::array<float, 3>& /*pv*/, const aod::V0Datas::iterator& v0, const Tracks&)
  {
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0NoCut);

    const auto& posTrack = v0.template posTrack_as<Tracks>();
    const auto& negTrack = v0.template negTrack_as<Tracks>();

    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0DaughterQuality);

    if (!qualitySelectionV0(v0)) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0Topology);

    std::array<float, 3> momPos{v0.pxpos(), v0.pypos(), v0.pzpos()},
      momNeg{v0.pxneg(), v0.pyneg(), v0.pzneg()};

    mHistograms.fill(HIST("armenteros_plot_before_selections"), v0.alpha(), v0.qtarm());
    mHistograms.fill(HIST("Lambda_vs_K0s"), v0.mK0Short(), v0.mAntiLambda());

    uint8_t v0Bitmask = selectV0MotherHypothesis(v0);
    if (v0Bitmask == 0 || (v0Bitmask & (v0Bitmask - 1)) != 0) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0PID);

    Candidate candidatePos(std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(),
                           RecoDecay::eta(momPos), RecoDecay::phi(momPos), posTrack.itsClusterSizes(),
                           0, posTrack.tpcInnerParam() * posTrack.sign(), posTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);
    Candidate candidateNeg(std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(),
                           RecoDecay::eta(momNeg), RecoDecay::phi(momNeg), negTrack.itsClusterSizes(),
                           0, negTrack.tpcInnerParam() * negTrack.sign(), negTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);

    if (!selectPidV0Daughters(candidatePos, candidateNeg, v0, posTrack, negTrack, v0Bitmask)) {
      return;
    }

    const float massV0 = fillHistogramsV0(v0, posTrack, negTrack, v0Bitmask);
    candidatePos.massMother = massV0;
    candidateNeg.massMother = massV0;

    if constexpr (isMC) { // MC
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
        return;
      }

      const auto& posMcParticle = posTrack.mcParticle();
      const auto& negMcParticle = negTrack.mcParticle();

      candidatePos.pdgCode = posMcParticle.pdgCode();
      candidateNeg.pdgCode = negMcParticle.pdgCode();

      fillMcHistogramsV0(v0, posMcParticle, negMcParticle);
    }

    fillTable<isMC>(candidatePos);
    fillTable<isMC>(candidateNeg);

    mHistograms.fill(HIST("isPositive"), true);
    mHistograms.fill(HIST("isPositive"), false);
  }

  template <bool isMC = false, typename Track>
  void fillKCand(const std::array<float, 3>& pv, aod::CascDatas::iterator const& cascade, const Track&)
  {
    mHistograms.fill(HIST("casc_selections"), CascSelections::kCascNoCut);
    const auto& bachelorTrack = cascade.template bachelor_as<Track>();

    std::array<float, 3> momBachelor{cascade.pxbach(), cascade.pybach(), cascade.pzbach()};

    if (!qualitySelectionCascade(cascade, pv)) {
      return;
    }
    mHistograms.fill(HIST("casc_selections"), CascSelections::kCascTopology);

    const float& massXi = cascade.mXi();
    const float& massOmega = cascade.mOmega();
    mHistograms.fill(HIST("Xi_vs_Omega"), massOmega, massXi);

    if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsettingMassWindowXi) {
      return;
    }
    mHistograms.fill(HIST("casc_selections"), CascSelections::kRejectedXi);

    mHistograms.fill(HIST("massOmega"), cascade.pt() * bachelorTrack.sign(), massOmega);
    if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsettingMassWindowOmega) {
      return;
    }
    mHistograms.fill(HIST("casc_selections"), CascSelections::kAcceptedOmega);

    if (std::abs(bachelorTrack.tpcNSigmaKa()) > cascsettingNsigmatpc) {
      return;
    }
    mHistograms.fill(HIST("casc_selections"), CascSelections::kNSigmaTPC);
    fillHistogramsParticle<PartID::ka, isMC>(bachelorTrack);

    mClusterStudiesTable(
      std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(),
      RecoDecay::eta(momBachelor), RecoDecay::phi(momBachelor),
      bachelorTrack.itsClusterSizes(), PartID::ka);
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        bachelorTrack.tpcInnerParam() * bachelorTrack.sign(),
        bachelorTrack.pidForTracking(),
        bachelorTrack.tpcNSigmaKa(), /*TOF nsigma*/ -999.f, /*TOF mass*/ -999.f,
        cascade.casccosPA(pv[0], pv[1], pv[2]), massOmega);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }

    if constexpr (isMC) {
      if (!bachelorTrack.has_mcParticle()) {
        return;
      }
      const auto& mcParticle = bachelorTrack.mcParticle();

      mClusterStudiesTableMc(
        mcParticle.pdgCode());

      const auto& posV0Daughter = cascade.template posTrack_as<Track>();
      if (!posV0Daughter.has_mcParticle()) {
        return;
      }
      const auto& mcPosParticleV0 = posV0Daughter.mcParticle();
      fillMcHistogramsCascade(cascade, mcParticle, mcPosParticleV0);
    }

    mHistograms.fill(HIST("isPositive"), bachelorTrack.sign() > 0);
  }

  template <int partID, bool isMC = false, typename Track>
  void fillNucleusTable(const Track& track)
  {
    constexpr int kPartID = partID;

    if (kPartID == static_cast<int>(PartID::de) && track.sign() > 0) {
      return;
    }
    mHistograms.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiNoCut);

    if (track.itsNCls() < desettingNclsIts)
      return;
    mHistograms.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiNClsIts);

    const float tpcNsigma = kPartID == static_cast<int>(PartID::de) ? track.tpcNSigmaDe() : computeNSigmaTPCHe3(track);
    const float tpcNsigmaMax = kPartID == static_cast<int>(PartID::de) ? desettingNsigmatpc : he3settingNsigmatpc;
    if (std::abs(tpcNsigma) > tpcNsigmaMax) {
      return;
    }
    mHistograms.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiPIDtpc);

    const float tofMass = track.hasTOF() ? (kPartID == static_cast<int>(PartID::de) ? computeTOFmassDe<isMC>(track) : computeTOFmassHe3<isMC>(track)) : -999.f;
    const float tofNsigma = kPartID == static_cast<int>(PartID::de) ? track.tofNSigmaDe() : -999.f;
    float correctedTPCinnerParam = track.tpcInnerParam();

    if (kPartID == static_cast<int>(PartID::de)) {
      if (!track.hasTOF() || std::abs(tofNsigma) > desettingNsigmatof) {
        return;
      }
    } else {
      if (track.hasTOF() && (tofMass < he3settingTofmasslow || tofMass > he3settingTofmasshigh)) {
        return;
      }

      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      correctedTPCinnerParam = (heliumPID && he3settingCompensatePIDinTracking) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();
    }
    mHistograms.fill(HIST(cNames[kPartID]) + HIST("/trackSelections"), NucleiSelections::kNucleiPIDtof);

    fillHistogramsParticle<kPartID, isMC>(track);

    mClusterStudiesTable(
      track.p() * track.sign(), track.eta(), track.phi(),
      track.itsClusterSizes(), kPartID);
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        correctedTPCinnerParam * track.sign(), track.pidForTracking(),
        tpcNsigma, tofNsigma, tofMass,
        /*cosPA*/ -999.f, /*mass mother*/ -999.f);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }

    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return;
      }

      const auto& mcParticle = track.mcParticle();

      mClusterStudiesTableMc(
        mcParticle.pdgCode());
    }

    mHistograms.fill(HIST("isPositive"), track.sign() > 0);
  }

  void fillPKPiTable(const TracksFullIU::iterator& track)
  {
    uint8_t partID = 0;
    float tpcNSigma = 0.f;
    if (std::abs(track.tpcNSigmaPi()) < v0settingNsigmatpcPi && std::abs(track.tpcNSigmaKa()) > 3) {
      partID = PartID::pi;
      tpcNSigma = track.tpcNSigmaPi();
      mHistograms.fill(HIST("nSigmaTPCPi"), track.p() * track.sign(), tpcNSigma);
    } else if (std::abs(track.tpcNSigmaKa()) < cascsettingNsigmatpc && (std::abs(track.tpcNSigmaPi()) > 3 /*&& std::abs(track.tpcNSigmaPr()) > 3*/)) {
      partID = PartID::ka;
      tpcNSigma = track.tpcNSigmaKa();
      mHistograms.fill(HIST("nSigmaTPCKa"), track.p() * track.sign(), tpcNSigma);
    } else if (std::abs(track.tpcNSigmaPr()) < v0settingNsigmatpcPr && std::abs(track.tpcNSigmaKa()) > 3) {
      partID = PartID::pr;
      tpcNSigma = track.tpcNSigmaPr();
      mHistograms.fill(HIST("nSigmaTPCPr"), track.p() * track.sign(), tpcNSigma);
    } else {
      return;
    }

    mClusterStudiesTable(
      track.p() * track.sign(), track.eta(), track.phi(),
      track.itsClusterSizes(), partID);
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        track.tpcInnerParam() * track.sign(), track.pidForTracking(),
        tpcNSigma, /*TOF nsigma*/ -999.f, /*TOF mass*/ -999.f,
        /*cosPA*/ -999.f, /*mass mother*/ -999.f);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }
  }

  template <const bool isMC = false, typename Track>
  void fillElectronTableFromPi0Dalitz(const Track& posTrack, const Track& negTrack)
  {
    mHistograms.fill(HIST("e_selections"), ESelections::kENoCut);
    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("e_selections"), ESelections::kETrackQuality);

    if (!electronPrimarySelection(posTrack) || !electronPrimarySelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("e_selections"), ESelections::kEPrimary);

    if (!electronPidSelection(posTrack) || !electronPidSelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("e_selections"), ESelections::kEPid);

    const float invariantMass = std::sqrt(RecoDecay::m2<2>(std::array<std::array<float, 3>, 2>{
                                                             std::array<float, 3>{posTrack.px(), posTrack.py(), posTrack.pz()},
                                                             std::array<float, 3>{negTrack.px(), negTrack.py(), negTrack.pz()}},
                                                           std::array<float, 2>{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}));
    float phiv = o2::aod::pwgem::dilepton::utils::pairutil::getPhivPair(posTrack.px(), posTrack.py(), posTrack.pz(), negTrack.px(), negTrack.py(), negTrack.pz(), posTrack.sign(), negTrack.sign(), md_bz);

    mHistograms.fill(HIST("massPi0WithBkg"), invariantMass);
    if (invariantMass > electronsettingMaxMassPi0 || invariantMass < electronsettingMinMassPi0) {
      return;
    }
    if (phiv > electronsettingMaxPhiv || phiv < electronsettingMinPhiv) {
      return;
    }
    mHistograms.fill(HIST("e_selections"), ESelections::kEPi0);
    mHistograms.fill(HIST("massPi0"), invariantMass);
    fillHistogramsParticle<PartID::el, isMC>(posTrack);
    fillHistogramsParticle<PartID::el, isMC>(negTrack);

    // float opangle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(t1.px(), t1.py(), t1.pz(), t2.px(), t2.py(), t2.pz());

    mClusterStudiesTable(
      posTrack.p(), posTrack.eta(), posTrack.phi(),
      posTrack.itsClusterSizes(), PartID::el);
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        posTrack.tpcInnerParam(), posTrack.pidForTracking(),
        posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), -999.f /* TofMass */,
        -999.f /* cosPA */, invariantMass);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }

    mClusterStudiesTable(
      negTrack.p(), negTrack.eta(), negTrack.phi(),
      negTrack.itsClusterSizes(), PartID::el);
    if (settingFillExtraTable) {
      mClusterStudiesTableExtra(
        negTrack.tpcInnerParam(), negTrack.pidForTracking(),
        negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), -999.f /* TofMass */,
        -999.f /* cosPA */, invariantMass);
    }
    if (settingFillCollTable) {
      mClusterStudiesTableCollision(
        mrunNumber);
    }

    if constexpr (isMC) {
      const auto& posMcParticle = posTrack.mcParticle();
      const auto& negMcParticle = negTrack.mcParticle();

      mClusterStudiesTableMc(
        posMcParticle.pdgCode());
      mClusterStudiesTableMc(
        negMcParticle.pdgCode());

      if (!posMcParticle.has_mothers() || !negMcParticle.has_mothers()) {
        return;
      }

      for (const auto& posMother : posMcParticle.template mothers_as<aod::McParticles>()) {
        for (const auto& negMother : negMcParticle.template mothers_as<aod::McParticles>()) {
          if (posMother.globalIndex() != negMother.globalIndex() || std::abs(posMother.pdgCode()) != PDG_t::kPi0) {
            return;
          }
          mHistograms.fill(HIST("massPi0Mc"), std::sqrt((posMcParticle.e() + negMcParticle.e()) * (posMcParticle.e() + negMcParticle.e()) -
                                                        (posMcParticle.px() + negMcParticle.px()) * (posMcParticle.px() + negMcParticle.px()) -
                                                        (posMcParticle.py() + negMcParticle.py()) * (posMcParticle.py() + negMcParticle.py()) -
                                                        (posMcParticle.pz() + negMcParticle.pz()) * (posMcParticle.pz() + negMcParticle.pz())));
          break;
        }
      }
    }
  }

  template <const bool isMC = false, typename Tracks>
  void fillElectronTableFromPhotonConversion(const std::array<float, 3>& /*pv*/, const aod::V0Datas::iterator& v0, const Tracks&)
  {

    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0NoCut);

    const auto& posTrack = v0.template posTrack_as<Tracks>();
    const auto& negTrack = v0.template negTrack_as<Tracks>();

    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0DaughterQuality);

    if (!qualitySelectionPhotons(v0)) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0Topology);

    std::array<float, 3> momPos{v0.pxpos(), v0.pypos(), v0.pzpos()},
      momNeg{v0.pxneg(), v0.pyneg(), v0.pzneg()};

    mHistograms.fill(HIST("armenteros_plot_before_selections"), v0.alpha(), v0.qtarm());

    Candidate candidatePos(std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(),
                           RecoDecay::eta(momPos), RecoDecay::phi(momPos), posTrack.itsClusterSizes(),
                           0, posTrack.tpcInnerParam() * posTrack.sign(), posTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);
    Candidate candidateNeg(std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(),
                           RecoDecay::eta(momNeg), RecoDecay::phi(momNeg), negTrack.itsClusterSizes(),
                           0, negTrack.tpcInnerParam() * negTrack.sign(), negTrack.pidForTracking(),
                           -999.f, -999.f, -999.f, v0.v0cosPA(), -999.f, 0);

    float invariantMass = std::sqrt(RecoDecay::m2<2>(std::array<std::array<float, 3>, 2>{
                                                       std::array<float, 3>{posTrack.px(), posTrack.py(), posTrack.pz()},
                                                       std::array<float, 3>{negTrack.px(), negTrack.py(), negTrack.pz()}},
                                                     std::array<float, 2>{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron}));

    if (invariantMass > electronsettingMaxMassPhoton || invariantMass < electronsettingMinMassPhoton) {
      return;
    }

    if (!electronPidSelection(posTrack) || !electronPidSelection(negTrack)) {
      return;
    }
    mHistograms.fill(HIST("v0_selections"), V0Selections::kV0PID);

    mHistograms.fill(HIST("armenteros_plot_gamma"), v0.alpha(), v0.qtarm());
    mHistograms.fill(HIST("photon_radiusV0"), v0.v0radius());
    candidatePos.massMother = invariantMass;
    candidateNeg.massMother = invariantMass;

    if constexpr (isMC) { // MC
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
        return;
      }

      const auto& posMcParticle = posTrack.mcParticle();
      const auto& negMcParticle = negTrack.mcParticle();

      candidatePos.pdgCode = posMcParticle.pdgCode();
      candidateNeg.pdgCode = negMcParticle.pdgCode();

      if (candidatePos.pdgCode != PDG_t::kElectron || candidateNeg.pdgCode != PDG_t::kElectron) {
        return;
      }
      if (!posMcParticle.has_mothers() || !negMcParticle.has_mothers()) {
        return;
      }
    }

    fillHistogramsParticle<PartID::el, isMC>(posTrack);
    fillHistogramsParticle<PartID::el, isMC>(negTrack);

    fillTable<isMC>(candidatePos);
    fillTable<isMC>(candidateNeg);

    mHistograms.fill(HIST("isPositive"), true);
    mHistograms.fill(HIST("isPositive"), false);
  }

  // =========================================================================================================

  void processDataV0Casc(CollisionsCustom::iterator const& collision /*s*/, TracksFullIU const& tracks, aod::V0Datas const& v0s, aod::CascDatas const& cascades, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mHistograms.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillV0Cand</*isMC*/ false>(PV, v0, tracks);
    }

    for (const auto& cascade : cascades) {
      fillKCand</*isMC*/ false>(PV, cascade, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataV0Casc, "process Data V0 and cascade", false);

  void processDataElectronsFromPhotonConversion(CollisionsCustom::iterator const& collision /*s*/, TracksFullIU const& tracks, aod::V0Datas const& v0s, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mHistograms.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillElectronTableFromPhotonConversion</*isMC*/ false>(PV, v0, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataElectronsFromPhotonConversion, "process Data Electrons from Photon Conversion", false);

  Partition<TracksFullIU> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<TracksFullIU> negTracks = o2::aod::track::signed1Pt < 0.f;
  void processDataElectronsFromDalitz(CollisionsCustom::iterator const& collision, TracksFullIU const& /*tracks*/, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mHistograms.fill(HIST("zVtx"), collision.posZ());

    const auto& posTracks_thisCollision = posTracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), mcache);
    const auto& negTracks_thisCollision = negTracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), mcache);

    for (const auto& [posTrack, negTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_thisCollision, negTracks_thisCollision))) {
      fillElectronTableFromPi0Dalitz(posTrack, negTrack);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataElectronsFromDalitz, "process Data Electrons from Dalitz", false);

  void processDataNuclei(CollisionsCustom::iterator const& collision, TracksFullIU const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mHistograms.fill(HIST("zVtx"), collision.posZ());

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

    mHistograms.fill(HIST("zVtx"), collision.posZ());

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

    mHistograms.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillV0Cand</*isMC*/ true>(PV, v0, tracks);
    }

    for (const auto& cascade : cascades) {
      fillKCand</*isMC*/ true>(PV, cascade, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcV0Casc, "process Mc V0 and cascade", false);

  void processMcElectronsFromPhotonConversion(CollisionsCustomMc::iterator const& collision /*s*/, TracksFullIUMc const& tracks, aod::V0Datas const& v0s, aod::BCsWithTimestamps const&, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    mHistograms.fill(HIST("zVtx"), collision.posZ());
    std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

    for (const auto& v0 : v0s) {
      fillElectronTableFromPhotonConversion</*isMC*/ true>(PV, v0, tracks);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcElectronsFromPhotonConversion, "process Mc Electrons from Photon Conversion", false);

  Partition<TracksFullIUMc> posTracksMc = o2::aod::track::signed1Pt > 0.f;
  Partition<TracksFullIUMc> negTracksMc = o2::aod::track::signed1Pt < 0.f;
  void processMcElectronsFromDalitz(CollisionsCustomMc::iterator const& collision, TracksFullIUMc const& /*tracks*/, aod::BCsWithTimestamps const&, aod::McParticles const&, aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collisionSelection(collision)) {
      return;
    }

    const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    mHistograms.fill(HIST("zVtx"), collision.posZ());

    const auto& posTracks_thisCollision = posTracksMc.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), mcache);
    const auto& negTracks_thisCollision = negTracksMc.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), mcache);

    for (const auto& [posTrack, negTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posTracks_thisCollision, negTracks_thisCollision))) {
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle())
        continue;

      fillElectronTableFromPi0Dalitz</*isMC*/ true>(posTrack, negTrack);
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcElectronsFromDalitz, "process Mc Electrons from Dalitz", false);

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
    mHistograms.fill(HIST("zVtx"), collision.posZ());

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
