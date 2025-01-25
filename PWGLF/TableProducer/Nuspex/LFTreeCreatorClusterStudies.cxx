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

#include <vector>
#include <utility>
#include <random>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <numeric>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "DCAFitter/DCAFitterN.h"

#include "PWGLF/DataModel/LFClusterStudiesTable.h"

#include "TDatabasePDG.h"
#include "TPDGCode.h"

using namespace ::o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Track = o2::track::TrackParCov;
using TracksFullIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFDe, aod::TOFSignal, aod::TOFEvTime>;
using TracksFullIUMc = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFDe, aod::TOFSignal, aod::TOFEvTime, aod::McTrackLabels>;
using CollisionsCustom = soa::Join<aod::Collisions, aod::EvSels>;

namespace BetheBloch
{

constexpr double defaultParams[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> parNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

} // namespace BetheBloch

enum V0Type : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda,
  Photon,
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
  kV0DaughterDCA,
  // kV0DCA,
  kV0Radius,
  kV0CosPA,
  kV0PID,
  kV0DaughterDCAtoPV,
  kV0All
};

enum CascSelections {
  kCascNoCut = 0,
  kCascDCA,
  kCascCosPA,
  kAcceptedOmega,
  kRejectedXi,
  kNSigmaTPC,
  kCascAll
};

enum DeSelections {
  kDeNoCut = 0,
  kDeNClsIts,
  kDePIDtpc,
  kDePIDtof,
  kDeAll
};

enum He3Selections {
  kHe3NoCut = 0,
  kHe3NClsIts,
  kHe3PIDtpc,
  kHe3PIDtof,
  kHe3All
};

enum PartID {
  none = 0,
  el,
  pi,
  ka,
  pr,
  de,
  he
};

struct LfTreeCreatorClusterStudies {

  Service<o2::ccdb::BasicCCDBManager> m_ccdb;
  int m_runNumber;
  int m_collisionCounter = 0;
  float m_d_bz;
  uint32_t m_randomSeed = 0.;

  Configurable<bool> setting_fillV0{"fillV0", true, "Fill the V0 tree"};
  Configurable<bool> setting_fillK{"fillK", true, "Fill the K tree"};
  Configurable<bool> setting_fillDe{"fillDe", true, "Fill the De tree"};
  Configurable<bool> setting_fillHe3{"fillHe3", true, "Fill the He3 tree"};
  Configurable<bool> setting_fillPKPi{"fillPKPPi", true, "Fill the p, K, pi tree"};
  Configurable<bool> setting_smallTable{"smallTable", true, "Use a small table for testing"};

  Configurable<int> setting_materialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};

  Configurable<float> setting_zVtxMax{"zVtxMax", 10.f, "Maximum z vertex position"};

  Configurable<float> setting_downscaleFactor{"downscaleFactor", 1.f, "Downscale factor for the V0 candidates"};
  Configurable<bool> setting_applyAdditionalEvSel{"applyAdditionalEvSel", false, "Apply additional event selection"};

  Configurable<float> track_nClsItsMin{"track_NclsItsMin", 0.f, "Minimum number of ITS clusters for the V0 daughters"};
  Configurable<float> track_nClsTpcMin{"track_NclsTpcMin", 100.f, "Minimum number of TPC clusters for the V0 daughters"};
  Configurable<float> track_nClsTpcMaxShared{"track_NclsTpcMaxShared", 5.f, "Maximum number of shared TPC clusters for the V0 daughters"};

  // Configurable<float> v0setting_etaMaxV0{"etaMaxV0", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0setting_etaMaxV0dau{"etaMaxV0dau", 0.8f, "Maximum eta for the V0 daughters"};
  Configurable<float> v0setting_dcaV0daughters{"v0setting_dcaV0daughters", 0.5f, "DCA between the V0 daughters"};
  Configurable<float> v0setting_dcaV0toPV{"v0setting_dcaV0fromPV", 1.f, "DCA of the V0 to the primary vertex"};
  Configurable<float> v0setting_dcaDaughtersToPV{"v0setting_dcaDaughtersToPV", 1.f, "DCA of the daughters to the primary vertex"};
  Configurable<float> v0setting_radiusMax{"v0setting_radiusMax", 100.f, "Maximum radius of the V0 accepted"};
  Configurable<float> v0setting_radiusMin{"v0setting_radiusMin", 5.f, "Minimum radius of the V0 accepted"};
  Configurable<float> v0setting_cosPA{"v0setting_cosPA", 0.99f, "Cosine of the pointing angle of the V0"};
  Configurable<float> v0setting_nsigmatpc{"v0setting_nsigmaTPC", 4.f, "Number of sigmas for the TPC PID"};
  Configurable<float> v0setting_massWindowLambda{"v0setting_massWindowLambda", 0.02f, "Mass window for the Lambda"};
  Configurable<float> v0setting_massWindowK0s{"v0setting_massWindowK0s", 0.02f, "Mass window for the K0s"};
  Configurable<float> v0setting_nsigmatpcEl{"v0setting_nsigmaTPCEl", 1.f, "Number of sigmas for the TPC PID for electrons"};
  Configurable<float> v0setting_nsigmatpcPi{"v0setting_nsigmaTPCPi", 2.f, "Number of sigmas for the TPC PID for pions"};
  Configurable<float> v0setting_nsigmatpcPr{"v0setting_nsigmaTPCPr", 2.f, "Number of sigmas for the TPC PID for protons"};
  Configurable<float> lambdasetting_qtAPcut{"lambdasetting_qtAPcut", 0.02f, "Cut on the qt for the Armenteros-Podolanski plot for photon rejection"};
  Configurable<float> lambdasetting_pmin{"lambdasetting_pmin", 0.0f, "Minimum momentum for the V0 daughters"};

  Configurable<float> cascsetting_dcaCascDaughters{"casc_setting_dcaV0daughters", 0.1f, "DCA between the V0 daughters"};
  Configurable<float> cascsetting_cosPA{"casc_setting_cosPA", 0.99f, "Cosine of the pointing angle of the V0"};
  Configurable<float> cascsetting_massWindowOmega{"casc_setting_massWindowOmega", 0.01f, "Mass window for the Omega"};
  Configurable<float> cascsetting_massWindowXi{"casc_setting_massWindowXi", 0.01f, "Mass window for the Xi"};
  Configurable<float> cascsetting_nsigmatpc{"casc_setting_nsigmaTPC", 3.f, "Number of sigmas for the TPC PID"};

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
    {{"collision_selections", "Collision selection; selection; counts", {HistType::kTH1F, {{Selections::kAll, -0.5, static_cast<double>(Selections::kAll) - 0.5}}}},
     {"v0_selections", "V0 selection; selection; counts", {HistType::kTH1F, {{V0Selections::kV0All, -0.5, static_cast<double>(V0Selections::kV0All) - 0.5}}}},
     {"casc_selections", "Cascade selection; selection; counts", {HistType::kTH1F, {{CascSelections::kCascAll, -0.5, static_cast<double>(CascSelections::kCascAll) - 0.5}}}},
     {"de_selections", "Deuteron track selection; selection; counts", {HistType::kTH1F, {{DeSelections::kDeAll, -0.5, static_cast<double>(DeSelections::kDeAll) - 0.5}}}},
     {"he3_selections", "He3 track selection; selection; counts", {HistType::kTH1F, {{He3Selections::kHe3All, -0.5, static_cast<double>(He3Selections::kHe3All) - 0.5}}}},
     {"v0_type", "Selected V0; particle; counts", {HistType::kTH1F, {{V0Type::V0TypeAll, -0.5, static_cast<double>(V0Type::V0TypeAll) - 0.5}}}},
     {"radiusV0", "Decay radius (xy) V0; radius (cm); counts", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"massLambda", "#Lambda invariant mass; signed #it{p} (GeV/#it{c}); m (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {50, 1.08f, 1.18f}}}},
     {"Lambda_vs_K0s", "Mass #Lambda vs K^{0}_s; m_{K^{0}_{s}} (GeV/#it{c}^{2}); m_{#Lambda} (GeV/#it{c}^{2})", {HistType::kTH2F, {{50, 0.f, 1.f}, {70, 0.6f, 2.f}}}},
     {"armenteros_plot_before_selections", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
     {"armenteros_plot", "Armenteros-Podolanski plot; #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
     {"armenteros_plot_lambda", "Armenteros-Podolanski plot (#Lambda only); #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
     {"armenteros_plot_gamma", "Armenteros-Podolanski plot (#gamma only); #alpha; q_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1.f, 1.f}, {100, 0.f, 0.3f}}}},
     {"photon_radiusV0", "Photon conversion radius (xy) V0; radius (cm); counts", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"photon_conversion_position", "Photon conversion position; x (cm); y (cm)", {HistType::kTH2F, {{250, -5.f, 5.f}, {250, -5.f, 5.f}}}},
     {"photon_conversion_position_layer", "Photon conversion position (ITS layers); x (cm); y (cm)", {HistType::kTH2F, {{100, -5.f, 5.f}, {100, -5.f, 5.f}}}},
     {"Xi_vs_Omega", "Mass Xi vs Omega; mass Omega (GeV/#it{c}^{2}); mass Xi (GeV/#it{c}^{2})", {HistType::kTH2F, {{50, 1.f, 2.f}, {50, 1.f, 2.f}}}},
     {"massOmega", "Mass #Omega; signed #it{p}_{T} (GeV/#it{c}); mass (GeV/#it{c}^{2})", {HistType::kTH2F, {{100, -5.f, 5.f}, {100, 1.62f, 1.72f}}}},
     {"massOmegaWithBkg", "Mass Omega with Background; mass Omega (GeV/#it{c}^{2}); counts", {HistType::kTH1F, {{100, 1.62f, 1.72f}}}},
     {"nSigmaTPCEl", "nSigma TPC Electron; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} e", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -2.0f, 2.0f}}}},
     {"nSigmaTPCPi", "nSigma TPC Pion; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} #pi", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}}}},
     {"nSigmaTPCKa", "nSigma TPC Kaon; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} e", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -4.0f, 4.0f}}}},
     {"nSigmaTPCPr", "nSigma TPC Proton; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} p", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {60, -3.0f, 3.0f}}}},
     {"nSigmaTPCDe", "nSigma TPC Deuteron; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} d", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -3.0f, 3.0f}}}},
     {"nSigmaTPCHe", "nSigma TPC He3; signed #it{p} (GeV/#it{c}); n#sigma_{TPC} ^{3}He", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -3.0f, 3.0f}}}},
     {"nSigmaTOFDe", "nSigma TOF Deuteron; signed #it{p} (GeV/#it{c}); n#sigma_{TOF} d", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -3.0f, 3.0f}}}},
     {"TOFmassDe", "TOF mass De; signed #it{p}_{T} (GeV/#it{c}); mass_{TOF} ^{3}He (GeV/#it{c}^2)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, 1.0f, 5.0f}}}},
     {"TOFmassHe", "TOF mass He3; signed #it{p}_{T} (GeV/#it{c}); mass_{TOF} ^{3}He (GeV/#it{c}^2)", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, 1.0f, 5.0f}}}},
     {"nSigmaITSEl", "nSigma ITS Electron; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} e", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {60, -2.0f, 2.0f}}}},
     {"nSigmaITSPi", "nSigma ITS Pion; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} #pi", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {60, -3.0f, 3.0f}}}},
     {"nSigmaITSKa", "nSigma ITS Kaon; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} e", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {60, -4.0f, 4.0f}}}},
     {"nSigmaITSPr", "nSigma ITS Proton; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} p", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {60, -3.0f, 3.0f}}}},
     {"nSigmaITSDe", "nSigma ITS Deuteron; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} d", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {100, -3.0f, 3.0f}}}},
     {"nSigmaITSHe", "nSigma ITS He3; signed #it{p} (GeV/#it{c}); n#sigma_{ITS} ^{3}He", {HistType::kTH2F, {{50, -5.0f, 5.0f}, {100, -3.0f, 3.0f}}}},
     {"pmatchingEl", "#it{p} matching e; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"pmatchingPi", "#it{p} matching #pi; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"pmatchingKa", "#it{p} matching K; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"pmatchingPr", "#it{p} matching p; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"pmatchingDe", "#it{p} matching d; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"pmatchingHe", "#it{p} matching ^{3}He; signed #it{p}_{TPC} (GeV/#it{c}); #frac{#it{p}_{TPC} - #it{p}}{#it{p}_{TPC}}", {HistType::kTH2F, {{100, -5.0f, 5.0f}, {100, -0.5f, 0.5f}}}},
     {"zVtx", "Binning for the vertex z in cm", {HistType::kTH1F, {{100, -20.f, 20.f}}}},
     {"isPositive", "is the candidate positive?; isPositive; counts", {HistType::kTH1F, {{2, -0.5f, 1.5f}}}}},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true}; // check histograms

  Produces<o2::aod::ClStTable> m_ClusterStudiesTable;
  Produces<o2::aod::ClStTableExtra> m_ClusterStudiesTableExtra;
  Produces<o2::aod::ClStTableMc> m_ClusterStudiesTableMc;
  Produces<o2::aod::ClStTableMcExt> m_ClusterStudiesTableMcExtra;

  struct V0TrackParCov {
    int64_t globalIndex;
    Track trackParCov;
  };
  std::vector<V0TrackParCov> m_v0TrackParCovs;

  o2::vertexing::DCAFitterN<2> m_fitter;
  o2::aod::ITSResponse m_responseITS;

  template <typename T>
  bool initializeFitter(const T& trackParCovA, const T& trackParCovB)
  {
    int nCand = 0;
    try {
      nCand = m_fitter.process(trackParCovA, trackParCovB);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return false;
    }
    if (nCand == 0) {
      return false;
    }

    return true;
  }

  /**
   * Compute the momentum of the track using the fitter
   * @param itrack Index of the track in the fitter
   * @param mom Array to store the momentum
   */
  void computeTrackMomentum(const int itrack, std::array<float, 3>& mom)
  {
    auto fittedTrack = m_fitter.getTrack(itrack);
    fittedTrack.getPxPyPzGlo(mom);
  }

  void computeMotherMomentum(const std::array<float, 3>& momA, const std::array<float, 3>& momB, std::array<float, 3>& momMother)
  {
    momMother[0] = momA[0] + momB[0];
    momMother[1] = momA[1] + momB[1];
    momMother[2] = momA[2] + momB[2];
  }

  /**
   * Compute the alpha for the Armenteros-Podolanski plot
   */
  float computeAlphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momP, const std::array<float, 3>& momN)
  {
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momP.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momN.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  /**
   * Compute the qt for the Armenteros-Podolanski plot
   */
  float computeQtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momP)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momP.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momP.begin(), momP.end(), momP.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float dcaMotherToPV(const std::array<float, 3>& decayVtx, const std::array<float, 3>& PV, std::array<float, 3> momMother) const
  {
    std::array<float, 3> relPos = {decayVtx[0] - PV[0], decayVtx[1] - PV[1], decayVtx[2] - PV[2]};
    float lmomMotherl = std::hypot(momMother[0], momMother[1], momMother[2]);
    return std::sqrt((std::pow(relPos[1] * momMother[2] - relPos[2] * momMother[1], 2) + std::pow(relPos[2] * momMother[0] - relPos[0] * momMother[2], 2) + std::pow(relPos[0] * momMother[1] - relPos[1] * momMother[0], 2))) / lmomMotherl;
  }

  template <typename T>
  float dcaToPV(const std::array<float, 3>& PV, T& trackParCov, gpu::gpustd::array<float, 2>& dcaInfo)
  {
    o2::base::Propagator::Instance()->propagateToDCABxByBz({PV[0], PV[1], PV[2]}, trackParCov, 2.f, m_fitter.getMatCorrType(), &dcaInfo);
    return std::hypot(dcaInfo[0], dcaInfo[1]);
  }

  float computeMassMother(const float massA, const float massB, const std::array<float, 3>& momA, const std::array<float, 3>& momB, const std::array<float, 3>& momMother) const
  {
    float eA = std::hypot(massA, std::hypot(momA[0], momA[1], momA[2]));
    float eB = std::hypot(massB, std::hypot(momB[0], momB[1], momB[2]));
    float lmomMotherl = std::hypot(momMother[0], momMother[1], momMother[2]);
    float eMother = eA + eB;
    return std::sqrt(eMother * eMother - lmomMotherl * lmomMotherl);
  }

  bool collisionSelection(const CollisionsCustom::iterator& collision)
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
  template <typename T>
  bool qualityTrackSelection(const T& track)
  {
    if (std::abs(track.eta()) > v0setting_etaMaxV0dau) {
      return false;
    }
    if (track.itsNCls() < track_nClsItsMin ||
        track.tpcNClsFound() < track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < track_nClsTpcMin ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > track_nClsTpcMaxShared) {
      return false;
    }
    return true;
  }

  bool qualitySelectionV0(const double /*dcaV0toPV*/, const double dcaV0daughters, const double radiusV0, const double cosPA)
  {
    if (std::abs(dcaV0daughters) > v0setting_dcaV0daughters) {
      return false;
    }
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0DaughterDCA);
    if (radiusV0 > v0setting_radiusMax || radiusV0 < v0setting_radiusMin) {
      return false;
    }
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0Radius);
    if (std::abs(cosPA) < v0setting_cosPA) {
      return false;
    }
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0CosPA);
    return true;
  }

  bool qualitySelectionCascade(const double dcaCascDaughters, const double cosPA)
  {
    if (std::abs(dcaCascDaughters) > cascsetting_dcaCascDaughters) {
      return false;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kCascDCA);
    if (std::abs(cosPA) < cascsetting_cosPA) {
      return false;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kCascCosPA);
    return true;
  }

  // =========================================================================================================

  template <typename T>
  bool nucleiTrackSelection(const T& track)
  {
    if (track.tpcNClsFound() < 90) {
      return false;
    }
    return true;
  }

  // =========================================================================================================

  template <typename T>
  float computeNSigmaDe(const T& candidate)
  {
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() / constants::physics::MassDeuteron), m_BBparamsDe[0], m_BBparamsDe[1], m_BBparamsDe[2], m_BBparamsDe[3], m_BBparamsDe[4]);
    double resoTPC{expTPCSignal * m_BBparamsDe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDtpcDe(const T& candidate)
  {
    auto nSigmaDe = computeNSigmaDe(candidate);
    if (std::abs(nSigmaDe) < desetting_nsigmatpc) {
      return true;
    }
    return false;
  }

  template <bool isMC = false, typename T>
  float computeTOFmassDe(const T& candidate)
  {
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
    return candidate.tpcInnerParam() * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  // =========================================================================================================

  template <typename T>
  float computeNSigmaHe3(const T& candidate)
  {
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && he3setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(correctedTPCinnerParam * 2.f / constants::physics::MassHelium3), m_BBparamsHe[0], m_BBparamsHe[1], m_BBparamsHe[2], m_BBparamsHe[3], m_BBparamsHe[4]);
    double resoTPC{expTPCSignal * m_BBparamsHe[5]};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resoTPC);
  }

  template <typename T>
  bool selectionPIDtpcHe3(const T& candidate)
  {
    auto nSigmaHe3 = computeNSigmaHe3(candidate);
    if (std::abs(nSigmaHe3) < he3setting_nsigmatpc) {
      return true;
    }
    return false;
  }

  template <bool isMC = false, typename T>
  float computeTOFmassHe3(const T& candidate)
  {
    float beta = o2::pid::tof::Beta::GetBeta(candidate);
    beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
    bool heliumPID = candidate.pidForTracking() == o2::track::PID::Helium3 || candidate.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParamHe3 = (heliumPID && he3setting_compensatePIDinTracking) ? candidate.tpcInnerParam() / 2.f : candidate.tpcInnerParam();
    return correctedTPCinnerParamHe3 * 2.f * std::sqrt(1.f / (beta * beta) - 1.f);
  }

  // =========================================================================================================

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (m_runNumber == bc.runNumber()) {
      return;
    }

    auto timestamp = bc.timestamp();
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
    m_fitter.setBz(m_d_bz);

    // o2::base::Propagator::Instance()->setMatLUT(lut);
  }

  void init(o2::framework::InitContext&)
  {
    m_runNumber = 0;
    m_d_bz = 0;

    m_ccdb->setURL("http://alice-ccdb.cern.ch");
    m_ccdb->setCaching(true);
    m_ccdb->setLocalObjectValidityChecking();
    m_ccdb->setFatalWhenNull(false);
    // lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    m_fitter.setPropagateToPCA(true);
    m_fitter.setMaxR(200.);
    m_fitter.setMinParamChange(1e-3);
    m_fitter.setMinRelChi2Change(0.9);
    m_fitter.setMaxDZIni(4);
    m_fitter.setMaxDXYIni(4);
    m_fitter.setMaxChi2(1e9);
    m_fitter.setUseAbsDCA(true);
    m_fitter.setWeightedFinalPCA(false);
    int mat{static_cast<int>(setting_materialCorrection)};
    m_fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));

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

    std::vector<std::string> V0_selection_labels = {"All", "daughter track quality", "V0 daughters dca", "V0 radius", "V0 cosPA", "V0 mass selection", "V0 daughter DCA to PV"};
    for (int i = 0; i < V0Selections::kV0All; i++)
      m_hAnalysis.get<TH1>(HIST("v0_selections"))->GetXaxis()->SetBinLabel(i + 1, V0_selection_labels[i].c_str());

    std::vector<std::string> Casc_selection_labels = {"All", "Casc DCA", "Casc CosPA", "Accepted Omega", "Veto Xi", "n#sigma_{TPC} K"};
    for (int i = 0; i < CascSelections::kCascAll; i++)
      m_hAnalysis.get<TH1>(HIST("casc_selections"))->GetXaxis()->SetBinLabel(i + 1, Casc_selection_labels[i].c_str());

    std::vector<std::string> De_selection_labels = {"All", "n clusters ITS", "n#sigma_{TPC} d", "n#sigma_{TOF} d"};
    for (int i = 0; i < DeSelections::kDeAll; i++)
      m_hAnalysis.get<TH1>(HIST("de_selections"))->GetXaxis()->SetBinLabel(i + 1, De_selection_labels[i].c_str());

    std::vector<std::string> He3_selection_labels = {"All", "n clusters ITS", "n#sigma_{TPC} ^{3}He", "TOF mass ^{3}He"};
    for (int i = 0; i < He3Selections::kHe3All; i++)
      m_hAnalysis.get<TH1>(HIST("he3_selections"))->GetXaxis()->SetBinLabel(i + 1, He3_selection_labels[i].c_str());

    std::vector<std::string> V0Type_labels = {"K0s", "#Lambda", "#bar{#Lambda}", "Photon"};
    for (int i = 0; i < V0Type::V0TypeAll; i++)
      m_hAnalysis.get<TH1>(HIST("v0_type"))->GetXaxis()->SetBinLabel(i + 1, V0Type_labels[i].c_str());
  }

  template <bool isMC = false, typename Track>
  void fillV0Cand(const std::array<float, 3>& PV, const aod::V0s::iterator& v0, const Track&)
  {
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0NoCut);

    auto posTrack = v0.posTrack_as<Track>();
    auto negTrack = v0.negTrack_as<Track>();
    if (!qualityTrackSelection(posTrack) || !qualityTrackSelection(negTrack)) {
      return;
    }
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0DaughterQuality);

    auto daughterTrackCovarianceA = getTrackParCov(posTrack);
    auto daughterTrackCovarianceB = getTrackParCov(negTrack);
    if (!initializeFitter(daughterTrackCovarianceA, daughterTrackCovarianceB)) {
      return;
    }

    std::array<float, 3> momPos, momNeg, momMother;
    computeTrackMomentum(0, momPos);
    computeTrackMomentum(1, momNeg);
    computeMotherMomentum(momPos, momNeg, momMother);
    ROOT::Math::SVector<double, 3> vec_decayVtx = m_fitter.getPCACandidate();
    std::array<float, 3> decayVtx = {static_cast<float>(vec_decayVtx[0]), static_cast<float>(vec_decayVtx[1]), static_cast<float>(vec_decayVtx[2])};
    float alphaAP = computeAlphaAP(momMother, momPos, momNeg);
    float qtAP = computeQtAP(momMother, momPos);
    m_hAnalysis.fill(HIST("armenteros_plot_before_selections"), alphaAP, qtAP);

    gpu::gpustd::array<float, 2> dcaInfo;
    V0TrackParCov v0TrackParCov{v0.globalIndex(), m_fitter.createParentTrackParCov()};
    float dcaV0daughters = std::sqrt(std::abs(m_fitter.getChi2AtPCACandidate()));
    float radiusV0 = std::hypot(decayVtx[0], decayVtx[1]);
    float dcaV0toPV = dcaToPV(PV, v0TrackParCov.trackParCov, dcaInfo);
    float cosPA = RecoDecay::cpa(PV, decayVtx, momMother);
    if (!qualitySelectionV0(dcaV0toPV, dcaV0daughters, radiusV0, cosPA)) {
      return;
    }

    // mass hypothesis
    float massLambdaV0 = computeMassMother(o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);
    float massAntiLambdaV0 = computeMassMother(o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, momPos, momNeg, momMother);
    float massK0sV0 = computeMassMother(o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged, momPos, momNeg, momMother);
    m_hAnalysis.fill(HIST("Lambda_vs_K0s"), massK0sV0, massLambdaV0);
    // float massPhotonV0 = computeMassMother(o2::constants::physics::MassElectron, o2::constants::physics::MassElectron, momPos, momNeg, momMother);

    uint8_t v0Bitmask(0);
    if (v0.isPhotonV0()) {
      SETBIT(v0Bitmask, Photon);
    }
    if (std::abs(massK0sV0 - o2::constants::physics::MassK0Short) < v0setting_massWindowK0s) {
      SETBIT(v0Bitmask, K0s);
    }
    if ((std::abs(massLambdaV0 - o2::constants::physics::MassLambda0) < v0setting_massWindowLambda) && (alphaAP > 0)) {
      SETBIT(v0Bitmask, Lambda);
    }
    if ((std::abs(massAntiLambdaV0 - o2::constants::physics::MassLambda0) < v0setting_massWindowLambda) && (alphaAP < 0)) {
      SETBIT(v0Bitmask, AntiLambda);
    }
    if (v0Bitmask == 0 || (v0Bitmask & (v0Bitmask - 1)) != 0) {
      return;
    }
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0PID);

    uint8_t partID_pos{0}, partID_neg{0};
    if (TESTBIT(v0Bitmask, Lambda)) {
      if (qtAP < lambdasetting_qtAPcut)
        return;
      if (std::abs(posTrack.tpcNSigmaPr()) > v0setting_nsigmatpcPr || std::abs(negTrack.tpcNSigmaPi()) > v0setting_nsigmatpcPi)
        return;
      if (std::hypot(momMother[0], momMother[1], momMother[2]) < lambdasetting_pmin)
        return;
      partID_pos = PartID::pr;
      partID_neg = PartID::pi;
      m_hAnalysis.fill(HIST("v0_type"), V0Type::Lambda);
    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      if (qtAP < lambdasetting_qtAPcut)
        return;
      if (std::abs(posTrack.tpcNSigmaPi()) > v0setting_nsigmatpcPr || std::abs(negTrack.tpcNSigmaPr()) > v0setting_nsigmatpcPi)
        return;
      if (std::hypot(momMother[0], momMother[1], momMother[2]) < lambdasetting_pmin)
        return;
      partID_pos = PartID::pi;
      partID_neg = PartID::pr;
      m_hAnalysis.fill(HIST("v0_type"), V0Type::AntiLambda);
    } else if (TESTBIT(v0Bitmask, K0s)) {
      m_hAnalysis.fill(HIST("v0_type"), V0Type::K0s);
      return; // K0s not implemented
    } else if (TESTBIT(v0Bitmask, Photon)) {
      // require photon conversion to happen in one of the Inner Tracker layers (± 0.5 cm resolution)
      m_hAnalysis.fill(HIST("photon_conversion_position"), decayVtx[0], decayVtx[1]);
      m_hAnalysis.fill(HIST("photon_radiusV0"), radiusV0);
      if (!(radiusV0 > 1.76 && radiusV0 < 4.71))
        return;
      if (std::abs(posTrack.tpcNSigmaEl()) > v0setting_nsigmatpcEl || std::abs(negTrack.tpcNSigmaEl()) > v0setting_nsigmatpcEl)
        return;
      m_hAnalysis.fill(HIST("photon_conversion_position_layer"), decayVtx[0], decayVtx[1]);
      partID_pos = PartID::el;
      partID_neg = PartID::el;
      m_hAnalysis.fill(HIST("v0_type"), V0Type::Photon);
    } else {
      return;
    }

    float dcaToPVpos = dcaToPV(PV, daughterTrackCovarianceA, dcaInfo);
    if (std::abs(dcaToPVpos) < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
      return;
    }
    float dcaToPVneg = dcaToPV(PV, daughterTrackCovarianceB, dcaInfo);
    if (std::abs(dcaToPVneg) < v0setting_dcaDaughtersToPV /*&& std::abs(dcaInfo[0]) < v0setting_dcaDaughtersToPV*/) {
      return;
    }

    float massV0{0.f};
    m_hAnalysis.fill(HIST("v0_selections"), V0Selections::kV0DaughterDCAtoPV);
    if (TESTBIT(v0Bitmask, Lambda)) {
      massV0 = massLambdaV0;
      m_hAnalysis.fill(HIST("massLambda"), std::hypot(momMother[0], momMother[1], momMother[2]), massLambdaV0);
      m_hAnalysis.fill(HIST("armenteros_plot_lambda"), alphaAP, qtAP);
      m_hAnalysis.fill(HIST("nSigmaTPCPr"), std::hypot(momPos[0], momPos[1], momPos[2]), posTrack.tpcNSigmaPr());
      m_hAnalysis.fill(HIST("nSigmaITSPr"), std::hypot(momPos[0], momPos[1], momPos[2]), m_responseITS.nSigmaITS<o2::track::PID::Proton>(posTrack.itsClusterSizes(), posTrack.p(), posTrack.eta()));
      m_hAnalysis.fill(HIST("nSigmaTPCPi"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, negTrack.tpcNSigmaPi());
      m_hAnalysis.fill(HIST("nSigmaITSPi"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, m_responseITS.nSigmaITS<o2::track::PID::Pion>(negTrack.itsClusterSizes(), negTrack.p(), negTrack.eta()));
      m_hAnalysis.fill(HIST("pmatchingPr"), posTrack.tpcInnerParam(), (posTrack.tpcInnerParam() - posTrack.p()) / posTrack.tpcInnerParam());
      m_hAnalysis.fill(HIST("pmatchingPi"), -negTrack.tpcInnerParam(), (negTrack.tpcInnerParam() - negTrack.p()) / negTrack.tpcInnerParam());

    } else if (TESTBIT(v0Bitmask, AntiLambda)) {
      massV0 = massAntiLambdaV0;
      m_hAnalysis.fill(HIST("massLambda"), std::hypot(momMother[0], momMother[1], momMother[2]) * -1.f, massAntiLambdaV0);
      // "signed" pt for antimatter
      m_hAnalysis.fill(HIST("armenteros_plot_lambda"), alphaAP, qtAP);
      m_hAnalysis.fill(HIST("nSigmaTPCPi"), std::hypot(momPos[0], momPos[1], momPos[2]), posTrack.tpcNSigmaPi());
      m_hAnalysis.fill(HIST("nSigmaITSPi"), std::hypot(momPos[0], momPos[1], momPos[2]), m_responseITS.nSigmaITS<o2::track::PID::Pion>(posTrack.itsClusterSizes(), posTrack.p(), posTrack.eta()));
      m_hAnalysis.fill(HIST("nSigmaTPCPr"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, negTrack.tpcNSigmaPr());
      m_hAnalysis.fill(HIST("nSigmaITSPr"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, m_responseITS.nSigmaITS<o2::track::PID::Proton>(negTrack.itsClusterSizes(), negTrack.p(), negTrack.eta()));
      m_hAnalysis.fill(HIST("pmatchingPi"), posTrack.tpcInnerParam(), (posTrack.tpcInnerParam() - posTrack.p()) / posTrack.tpcInnerParam());
      m_hAnalysis.fill(HIST("pmatchingPr"), -negTrack.tpcInnerParam(), (negTrack.tpcInnerParam() - negTrack.p()) / negTrack.tpcInnerParam());

    } else if (TESTBIT(v0Bitmask, Photon)) {
      massV0 = 0.f;
      m_hAnalysis.fill(HIST("nSigmaTPCEl"), std::hypot(momPos[0], momPos[1], momPos[2]), posTrack.tpcNSigmaEl());
      m_hAnalysis.fill(HIST("nSigmaTPCEl"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, negTrack.tpcNSigmaEl());
      m_hAnalysis.fill(HIST("nSigmaITSEl"), std::hypot(momPos[0], momPos[1], momPos[2]), m_responseITS.nSigmaITS<o2::track::PID::Electron>(posTrack.itsClusterSizes(), posTrack.p(), posTrack.eta()));
      m_hAnalysis.fill(HIST("nSigmaITSEl"), std::hypot(momNeg[0], momNeg[1], momNeg[2]) * -1.f, m_responseITS.nSigmaITS<o2::track::PID::Electron>(negTrack.itsClusterSizes(), negTrack.p(), negTrack.eta()));
      m_hAnalysis.fill(HIST("armenteros_plot_gamma"), alphaAP, qtAP);
      m_hAnalysis.fill(HIST("pmatchingEl"), posTrack.tpcInnerParam(), (posTrack.tpcInnerParam() - posTrack.p()) / posTrack.tpcInnerParam());
      m_hAnalysis.fill(HIST("pmatchingEl"), -negTrack.tpcInnerParam(), (negTrack.tpcInnerParam() - negTrack.p()) / negTrack.tpcInnerParam());
    }
    m_hAnalysis.fill(HIST("radiusV0"), radiusV0);
    m_hAnalysis.fill(HIST("armenteros_plot"), alphaAP, qtAP);
    m_v0TrackParCovs.push_back(v0TrackParCov);

    if (!setting_fillV0) {
      return;
    }

    if constexpr (isMC) { // MC
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
        return;
      }

      auto posMcParticle = posTrack.mcParticle();
      auto negMcParticle = negTrack.mcParticle();

      if (setting_smallTable) {
        m_ClusterStudiesTableMc(
          std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(), // p_pos
          RecoDecay::eta(momPos),                                        // eta_pos
          RecoDecay::phi(momPos),                                        // phi_pos
          posTrack.itsClusterSizes(),                                    // itsClsize_pos
          partID_pos,                                                    // partID_pos
          posMcParticle.pdgCode());                                      // pdgCode_pos
        m_ClusterStudiesTableMc(
          std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(), // p_neg
          RecoDecay::eta(momNeg),                                        // eta_neg
          RecoDecay::phi(momNeg),                                        // phi_neg
          negTrack.itsClusterSizes(),                                    // itsClsize_neg
          partID_neg,                                                    // partID_neg
          negMcParticle.pdgCode());                                      // pdgCode_neg
      } else {
        m_ClusterStudiesTableMcExtra(
          std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(), // p_pos
          RecoDecay::eta(momPos),                                        // eta_pos
          RecoDecay::phi(momPos),                                        // phi_pos
          posTrack.itsClusterSizes(),                                    // itsClsize_pos
          partID_pos,                                                    // partID_pos
          posMcParticle.pdgCode(),                                       // pdgCode_pos
          posTrack.tpcInnerParam() * posTrack.sign(),                    // pTPC_pos
          posTrack.pidForTracking(),                                     // pidInTrk_pos
          -999.f,                                                        // TpcNSigma_pos
          -999.f,                                                        // TofNSigma_pos
          -999.f,                                                        // TofMass_pos
          cosPA,                                                         // cosPA
          massV0);                                                       // massV0
        m_ClusterStudiesTableMcExtra(
          std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(), // p_neg
          RecoDecay::eta(momNeg),                                        // eta_neg
          RecoDecay::phi(momNeg),                                        // phi_neg
          negTrack.itsClusterSizes(),                                    // itsClsize_neg
          partID_neg,                                                    // partID_neg
          negMcParticle.pdgCode(),                                       // pdgCode_pos
          negTrack.tpcInnerParam() * negTrack.sign(),                    // pTPC_neg
          negTrack.pidForTracking(),                                     // pidInTrk_neg
          -999.f,                                                        // TpcNSigma_neg
          -999.f,                                                        // TofNSigma_neg
          -999.f,                                                        // TofMass_neg
          cosPA,                                                         // cosPA
          massV0);                                                       // massV0
      }
    } else { // data
      if (setting_smallTable) {
        m_ClusterStudiesTable(
          std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(), // p_pos
          RecoDecay::eta(momPos),                                        // eta_pos
          RecoDecay::phi(momPos),                                        // phi_pos
          posTrack.itsClusterSizes(),                                    // itsClsize_pos
          partID_pos);                                                   // partID_pos
        m_ClusterStudiesTable(
          std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(), // p_neg
          RecoDecay::eta(momNeg),                                        // eta_neg
          RecoDecay::phi(momNeg),                                        // phi_neg
          negTrack.itsClusterSizes(),                                    // itsClsize_neg
          partID_neg);                                                   // partID_neg
      } else {
        m_ClusterStudiesTableExtra(
          std::hypot(momPos[0], momPos[1], momPos[2]) * posTrack.sign(), // p_pos
          RecoDecay::eta(momPos),                                        // eta_pos
          RecoDecay::phi(momPos),                                        // phi_pos
          posTrack.itsClusterSizes(),                                    // itsClsize_pos
          partID_pos,                                                    // partID_pos
          posTrack.tpcInnerParam() * posTrack.sign(),                    // pTPC_pos
          posTrack.pidForTracking(),                                     // pidInTrk_pos
          -999.f,                                                        // TpcNSigma_pos
          -999.f,                                                        // TofNSigma_pos
          -999.f,                                                        // TofMass_pos
          cosPA,                                                         // cosPA
          massV0);                                                       // massV0
        m_ClusterStudiesTableExtra(
          std::hypot(momNeg[0], momNeg[1], momNeg[2]) * negTrack.sign(), // p_neg
          RecoDecay::eta(momNeg),                                        // eta_neg
          RecoDecay::phi(momNeg),                                        // phi_neg
          negTrack.itsClusterSizes(),                                    // itsClsize_neg
          partID_neg,                                                    // partID_neg
          negTrack.tpcInnerParam() * negTrack.sign(),                    // pTPC_neg
          negTrack.pidForTracking(),                                     // pidInTrk_neg
          -999.f,                                                        // TpcNSigma_neg
          -999.f,                                                        // TofNSigma_neg
          -999.f,                                                        // TofMass_neg
          cosPA,                                                         // cosPA
          massV0);                                                       // massV0
      }
    }

    m_hAnalysis.fill(HIST("isPositive"), true);
    m_hAnalysis.fill(HIST("isPositive"), false);
  }

  template <bool isMC = false, typename Track>
  void fillKCand(const std::array<float, 3>& PV, const aod::Cascades::iterator& cascade, const Track&)
  {
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kCascNoCut);

    auto v0Track = cascade.template v0_as<aod::V0s>();
    auto bachelorTrack = cascade.template bachelor_as<Track>();

    auto itv0 = std::find_if(m_v0TrackParCovs.begin(), m_v0TrackParCovs.end(), [&](const V0TrackParCov& v0) { return v0.globalIndex == v0Track.globalIndex(); });
    if (itv0 == m_v0TrackParCovs.end()) {
      return;
    }

    auto v0TrackCovariance = itv0->trackParCov;
    auto bachelorTrackCovariance = getTrackParCov(bachelorTrack);
    if (!initializeFitter(v0TrackCovariance, bachelorTrackCovariance)) {
      return;
    }

    std::array<float, 3> momV0, momBachelor, momMother;
    computeTrackMomentum(0, momV0);
    computeTrackMomentum(1, momBachelor);
    computeMotherMomentum(momV0, momBachelor, momMother);

    ROOT::Math::SVector<double, 3> vec_decayVtx = m_fitter.getPCACandidate();
    std::array<float, 3> decayVtx = {static_cast<float>(vec_decayVtx[0]), static_cast<float>(vec_decayVtx[1]), static_cast<float>(vec_decayVtx[2])};

    float dcaV0daughters = std::sqrt(std::abs(m_fitter.getChi2AtPCACandidate()));
    float cosPA = RecoDecay::cpa(PV, decayVtx, momMother);

    if (!qualitySelectionCascade(dcaV0daughters, cosPA)) {
      return;
    }
    // gpu::gpustd::array<float, 2> dcaInfo;
    // float dcaToPVbachelor = dcaToPV(PV, bachelorTrackCovariance, dcaInfo);

    float massXi = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassPionCharged, momV0, momBachelor, momMother);
    float massOmega = computeMassMother(o2::constants::physics::MassLambda0, o2::constants::physics::MassKaonCharged, momV0, momBachelor, momMother);
    m_hAnalysis.fill(HIST("Xi_vs_Omega"), massOmega, massXi);
    if (std::abs(massOmega - o2::constants::physics::MassOmegaMinus) > cascsetting_massWindowOmega) {
      return;
    }
    m_hAnalysis.fill(HIST("massOmegaWithBkg"), massOmega);
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kAcceptedOmega);
    if (std::abs(massXi - o2::constants::physics::MassXiMinus) < cascsetting_massWindowXi) {
      return;
    } // enhance purity by rejecting Xi background
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kRejectedXi);
    if (std::abs(bachelorTrack.tpcNSigmaKa()) > cascsetting_nsigmatpc) {
      return;
    }
    m_hAnalysis.fill(HIST("casc_selections"), CascSelections::kNSigmaTPC);
    m_hAnalysis.fill(HIST("massOmega"), std::hypot(momMother[0], momMother[1]) * bachelorTrack.sign(), massOmega);
    m_hAnalysis.fill(HIST("pmatchingKa"), bachelorTrack.sign() * bachelorTrack.tpcInnerParam(), (bachelorTrack.tpcInnerParam() - bachelorTrack.p()) / bachelorTrack.tpcInnerParam());
    m_hAnalysis.fill(HIST("nSigmaTPCKa"), bachelorTrack.sign() * std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]), bachelorTrack.tpcNSigmaKa());

    uint8_t partID_bachelor = PartID::ka;

    if constexpr (isMC) {
      if (!bachelorTrack.has_mcParticle()) {
        return;
      }
      auto mcParticle = bachelorTrack.mcParticle();

      if (setting_smallTable) {
        m_ClusterStudiesTableMc(
          std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(), // p_K
          RecoDecay::eta(momBachelor),                                                       // eta_K
          RecoDecay::phi(momBachelor),                                                       // phi_K
          bachelorTrack.itsClusterSizes(),                                                   // itsClSize_K
          partID_bachelor,                                                                   // partID_K
          mcParticle.pdgCode());                                                             // pdgCode_K
      } else {
        m_ClusterStudiesTableMcExtra(
          std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(), // p_K
          RecoDecay::eta(momBachelor),                                                       // eta_K
          RecoDecay::phi(momBachelor),                                                       // phi_K
          bachelorTrack.itsClusterSizes(),                                                   // itsClSize_K
          partID_bachelor,                                                                   // partID_K
          mcParticle.pdgCode(),                                                              // pdgCode_K
          bachelorTrack.tpcInnerParam() * bachelorTrack.sign(),                              // pTPC_K
          bachelorTrack.pidForTracking(),                                                    // PIDinTrk_K
          -999.f,                                                                            // TpcNSigma_K
          -999.f,                                                                            // TofNSigma_K
          -999.f,                                                                            // TofMass_K
          cosPA,                                                                             // cosPA
          massOmega);                                                                        // massMother
      }
    } else {
      if (setting_smallTable) {
        m_ClusterStudiesTable(
          std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(), // p_K
          RecoDecay::eta(momBachelor),                                                       // eta_K
          RecoDecay::phi(momBachelor),                                                       // phi_K
          bachelorTrack.itsClusterSizes(),                                                   // itsClSize_K
          partID_bachelor);                                                                  // partID_K
      } else {
        m_ClusterStudiesTableExtra(
          std::hypot(momBachelor[0], momBachelor[1], momBachelor[2]) * bachelorTrack.sign(), // p_K
          RecoDecay::eta(momBachelor),                                                       // eta_K
          RecoDecay::phi(momBachelor),                                                       // phi_K
          bachelorTrack.itsClusterSizes(),                                                   // itsClSize_K
          partID_bachelor,                                                                   // partID_K
          bachelorTrack.tpcInnerParam() * bachelorTrack.sign(),                              // pTPC_K
          bachelorTrack.pidForTracking(),                                                    // PIDinTrk_K
          -999.f,                                                                            // TpcNSigma_K
          -999.f,                                                                            // TofNSigma_K
          -999.f,                                                                            // TofMass_K
          cosPA,                                                                             // cosPA
          massOmega);
      }
    }

    m_hAnalysis.fill(HIST("isPositive"), bachelorTrack.p() > 0);
  }

  template <bool isMC = false, typename Track>
  void fillDeTable(const Track& track)
  {
    if (track.sign() > 0) {
      return;
    }
    m_hAnalysis.fill(HIST("de_selections"), DeSelections::kDeNoCut);
    if (track.itsNCls() < desetting_nClsIts) {
      return;
    }
    m_hAnalysis.fill(HIST("de_selections"), DeSelections::kDeNClsIts);
    if (!selectionPIDtpcDe(track)) {
      return;
    }
    m_hAnalysis.fill(HIST("de_selections"), DeSelections::kDePIDtpc);
    if (!track.hasTOF() || std::abs(track.tofNSigmaDe()) > desetting_nsigmatof) {
      return;
    }
    m_hAnalysis.fill(HIST("de_selections"), DeSelections::kDePIDtof);
    m_hAnalysis.fill(HIST("nSigmaTPCDe"), track.p() * track.sign(), computeNSigmaDe(track));
    m_hAnalysis.fill(HIST("nSigmaITSDe"), track.p() * track.sign(), m_responseITS.nSigmaITS<o2::track::PID::Deuteron>(track.itsClusterSizes(), track.p(), track.eta()));
    m_hAnalysis.fill(HIST("nSigmaTOFDe"), track.p() * track.sign(), track.tofNSigmaDe());
    m_hAnalysis.fill(HIST("TOFmassDe"), track.p() * track.sign(), computeTOFmassDe<isMC>(track));
    m_hAnalysis.fill(HIST("pmatchingDe"), track.sign() * track.tpcInnerParam(), (track.tpcInnerParam() - track.p()) / track.tpcInnerParam());

    uint8_t partID = PartID::de;

    if constexpr (isMC) {
      if (!track.has_mcParticle() || track.sign() > 0) {
        return;
      }

      auto mcParticle = track.mcParticle();

      if (setting_smallTable) {
        m_ClusterStudiesTableMc(
          track.p() * track.sign(), // p_De,
          track.eta(),              // eta_De,
          track.phi(),              // phi_De,
          track.itsClusterSizes(),  // itsClSize_De,
          partID,                   // pdgCode_De,
          mcParticle.pdgCode());    // pdgCodeMc_De
      } else {
        m_ClusterStudiesTableMcExtra(
          track.p() * track.sign(),             // p_De,
          track.eta(),                          // eta_De,
          track.phi(),                          // phi_De,
          track.itsClusterSizes(),              // itsClSize_De,
          partID,                               // pdgCode_De,
          mcParticle.pdgCode(),                 // pdgCodeMc_De
          track.tpcInnerParam() * track.sign(), // pTPC_De,
          track.pidForTracking(),               // PIDinTrk_De,
          computeNSigmaDe(track),               // TpcNSigma_De,
          track.tofNSigmaDe(),                  // TofNSigma_De,
          computeTOFmassDe<isMC>(track),        // TofMass_De,
          -999.f,                               // cosPA,
          -999.f);                              // massMother
      }
    } else {
      if (setting_smallTable) {
        m_ClusterStudiesTable(
          track.p() * track.sign(), // p_De,
          track.eta(),              // eta_De,
          track.phi(),              // phi_De,
          track.itsClusterSizes(),  // itsClSize_De,
          partID);                  // pdgCode_De
      } else {
        m_ClusterStudiesTableExtra(
          track.p() * track.sign(),             // p_De,
          track.eta(),                          // eta_De,
          track.phi(),                          // phi_De,
          track.itsClusterSizes(),              // itsClSize_De,
          partID,                               // pdgCode_De,
          track.tpcInnerParam() * track.sign(), // pTPC_De,
          track.pidForTracking(),               // PIDinTrk_De,
          computeNSigmaDe(track),               // TpcNSigma_De,
          track.tofNSigmaDe(),                  // TofNSigma_De,
          computeTOFmassDe<isMC>(track),        // TofMass_De,
          -999.f,                               // cosPA,
          -999.f);                              // massMother
      }
    }

    m_hAnalysis.fill(HIST("isPositive"), track.sign() > 0);
  }

  template <bool isMC = false, typename Track>
  void fillHe3Table(const Track& track)
  {
    m_hAnalysis.fill(HIST("he3_selections"), He3Selections::kHe3NoCut);

    if (track.itsNCls() < he3setting_nClsIts) {
      return;
    }
    m_hAnalysis.fill(HIST("he3_selections"), He3Selections::kHe3NClsIts);
    if (!selectionPIDtpcHe3(track)) {
      return;
    }
    m_hAnalysis.fill(HIST("he3_selections"), He3Selections::kHe3PIDtpc);
    float tofMass = track.hasTOF() ? computeTOFmassHe3<isMC>(track) : -999.f;
    if (track.hasTOF() && (tofMass < he3setting_tofmasslow || tofMass > he3setting_tofmasshigh)) {
      return;
    }
    uint8_t partID = PartID::he;
    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    float correctedTPCinnerParam = (heliumPID && he3setting_compensatePIDinTracking) ? track.tpcInnerParam() / 2.f : track.tpcInnerParam();

    m_hAnalysis.fill(HIST("he3_selections"), He3Selections::kHe3PIDtof);
    m_hAnalysis.fill(HIST("nSigmaTPCHe"), track.p() * track.sign(), computeNSigmaHe3(track));
    m_hAnalysis.fill(HIST("nSigmaITSHe"), track.p() * track.sign(), m_responseITS.nSigmaITS<o2::track::PID::Helium3>(track.itsClusterSizes(), track.p(), track.eta()));
    m_hAnalysis.fill(HIST("TOFmassHe"), track.p() * track.sign(), tofMass);
    m_hAnalysis.fill(HIST("pmatchingHe"), track.sign() * correctedTPCinnerParam, (correctedTPCinnerParam - track.p()) / correctedTPCinnerParam);

    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return;
      }
      auto mcParticle = track.mcParticle();

      if (setting_smallTable) {
        m_ClusterStudiesTableMc(
          track.p() * track.sign(), // p_He3,
          track.eta(),              // eta_He3,
          track.phi(),              // phi_He3,
          track.itsClusterSizes(),  // itsClSize_He3,
          partID,                   // pdgCode_He3,
          mcParticle.pdgCode());    // pdgCodeMc_He3
      } else {
        m_ClusterStudiesTableMcExtra(
          track.p() * track.sign(),              // p_He3
          track.eta(),                           // eta_He3
          track.phi(),                           // phi_He3
          track.itsClusterSizes(),               // itsClSize_He3
          partID,                                // pdgCode_He3
          mcParticle.pdgCode(),                  // pdgCodeMc_He3
          correctedTPCinnerParam * track.sign(), // pTPC_He3
          track.pidForTracking(),                // PIDinTrk_He3
          computeNSigmaHe3(track),               // TpcNSigma_He3
          -999.f,                                // TofNSigma_He3
          tofMass,                               // TofMass_He3
          -999.f,                                // cosPA_He3
          -999.f);                               // massMother_He3
      }

    } else {
      if (setting_smallTable) {
        m_ClusterStudiesTable(
          track.p() * track.sign(), // p_He3,
          track.eta(),              // eta_He3,
          track.phi(),              // phi_He3,
          track.itsClusterSizes(),  // itsClSize_He3,
          partID);                  // pdgCode_He3
      } else {
        m_ClusterStudiesTableExtra(
          track.p() * track.sign(),              // p_He3,
          track.eta(),                           // eta_He3,
          track.phi(),                           // phi_He3,
          track.itsClusterSizes(),               // itsClSize_He3,
          partID,                                // pdgCode_He3,
          correctedTPCinnerParam * track.sign(), // pTPC_He3,
          track.pidForTracking(),                // PIDinTrk_He3,
          computeNSigmaHe3(track),               // TpcNSigma_He3,
          -999.f,                                // TofNSigma_He3,
          tofMass,                               // TofMass_He3,
          -999.f,                                // cosPA,
          -999.f);                               // massMother
      }
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

    if (setting_smallTable) {
      m_ClusterStudiesTable(
        track.p() * track.sign(),
        track.eta(),
        track.phi(),
        track.itsClusterSizes(),
        partID);
    } else {
      m_ClusterStudiesTableExtra(
        track.p() * track.sign(),             // p,
        track.eta(),                          // eta,
        track.phi(),                          // phi,
        track.itsClusterSizes(),              // itsClSize,
        partID,                               // pdgCode,
        track.tpcInnerParam() * track.sign(), // pTPC,
        track.pidForTracking(),               // PIDinTrk,
        tpcNSigma,                            // TpcNSigma,
        -999.f,                               // TofNSigma,
        -999.f,                               // TofMass,
        -999.f,                               // cosPA,
        -999.f);                              // massMother
    }
  }

  // =========================================================================================================

  void processDataV0Casc(CollisionsCustom const& collisions, TracksFullIU const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collisionSelection(collision)) {
        continue;
      }

      m_hAnalysis.fill(HIST("zVtx"), collision.posZ());
      std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

      const uint64_t collIdx = collision.globalIndex();
      auto v0Table_thisCollision = v0s.sliceBy(m_perCollisionV0, collIdx);
      auto cascTable_thisCollision = cascades.sliceBy(m_perCollisionCascade, collIdx);
      v0Table_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&v0s);

      if (setting_fillV0 || setting_fillK) {
        m_v0TrackParCovs.clear();
        for (auto& v0 : v0Table_thisCollision) {
          fillV0Cand</*isMC*/ false>(PV, v0, tracks);
        }
      }
      if (setting_fillK) { // the v0 loops are needed for the Ks
        for (auto& cascade : cascTable_thisCollision) {
          fillKCand</*isMC*/ false>(PV, cascade, tracks);
        }
      }
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataV0Casc, "process Data V0 and cascade", false);

  void processDataNuclei(CollisionsCustom const& collisions, TracksFullIU const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!collisionSelection(collision)) {
        continue;
      }

      m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(m_perCol, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track : TrackTable_thisCollision) {
        if (!nucleiTrackSelection(track)) {
          continue;
        }

        if (setting_fillDe)
          fillDeTable</*isMC*/ false>(track);
        if (setting_fillHe3)
          fillHe3Table</*isMC*/ false>(track);
      }
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataNuclei, "process Data Nuclei", false);

  /**
   * @brief Produce a dataset with high purity p, K, #pi
   */
  void processDataPKPi(CollisionsCustom const& collisions, TracksFullIU const& tracks)
  {
    for (const auto& collision : collisions) {
      if (!collisionSelection(collision)) {
        continue;
      }

      m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(m_perCol, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track : TrackTable_thisCollision) {
        if (!qualityTrackSelection(track)) {
          continue;
        }

        if (setting_fillPKPi)
          fillPKPiTable(track);
      }
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processDataPKPi, "process Data p, K, pi", false);

  void processMcV0Casc(CollisionsCustom const& collisions, TracksFullIUMc const& tracks, aod::V0s const& v0s, aod::Cascades const& cascades, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      if (!collisionSelection(collision)) {
        continue;
      }

      m_hAnalysis.fill(HIST("zVtx"), collision.posZ());
      std::array<float, 3> PV = {collision.posX(), collision.posY(), collision.posZ()};

      const uint64_t collIdx = collision.globalIndex();
      auto v0Table_thisCollision = v0s.sliceBy(m_perCollisionV0, collIdx);
      auto cascTable_thisCollision = cascades.sliceBy(m_perCollisionCascade, collIdx);
      v0Table_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&tracks);
      cascTable_thisCollision.bindExternalIndices(&v0s);

      m_v0TrackParCovs.clear();
      if (setting_fillV0 || setting_fillK) {
        for (auto& v0 : v0Table_thisCollision) {
          fillV0Cand</*isMC*/ true>(PV, v0, tracks);
        }
      }

      if (setting_fillK) { // the v0 loops are needed for the Ks
        for (auto& cascade : cascTable_thisCollision) {
          fillKCand</*isMC*/ true>(PV, cascade, tracks);
        }
      }
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcV0Casc, "process Mc V0 and cascade", false);

  void processMcNuclei(CollisionsCustom const& collisions, TracksFullIUMc const& tracks, aod::BCs const&, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      if (!collisionSelection(collision)) {
        continue;
      }

      m_hAnalysis.fill(HIST("zVtx"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto TrackTable_thisCollision = tracks.sliceBy(m_perColMC, collIdx);
      TrackTable_thisCollision.bindExternalIndices(&tracks);

      for (auto track : TrackTable_thisCollision) {
        if (!nucleiTrackSelection(track)) {
          continue;
        }

        if (setting_fillDe) {
          fillDeTable</*isMC*/ true>(track);
        }
        if (setting_fillHe3) {
          fillHe3Table</*isMC*/ true>(track);
        }
      }
    }
  }
  PROCESS_SWITCH(LfTreeCreatorClusterStudies, processMcNuclei, "process Mc Nuclei", false);

}; // LfTreeCreatorClusterStudies

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LfTreeCreatorClusterStudies>(cfgc)};
}
