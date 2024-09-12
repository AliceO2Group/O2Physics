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
// Nuclei spectra analysis task
// ========================
//
// Executable + dependencies:
//
// Data (run3):
// o2-analysis-lf-nuclei-spectra, o2-analysis-timestamp
// o2-analysis-pid-tof-base, o2-analysis-multiplicity-table, o2-analysis-event-selection
// (to add flow: o2-analysis-qvector-table, o2-analysis-centrality-table)

#include <cmath>

#include "Math/Vector4D.h"

#include "CCDB/BasicCCDBManager.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Tools/TrackTuner.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"

#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct NucleusCandidate {
  int globalIndex;
  int collTrackIndex;
  float pt;
  float eta;
  float phi;
  float tpcInnerParam;
  float beta;
  float zVertex;
  float DCAxy;
  float DCAz;
  float TPCsignal;
  float ITSchi2;
  float TPCchi2;
  std::array<float, 5> nSigmaTPC;
  std::array<float, 5> tofMasses;
  bool fillTree;
  bool fillDCAHist;
  bool correctPV;
  bool isSecondary;
  bool fromWeakDecay;
  uint16_t flags;
  uint8_t TPCfindableCls;
  uint8_t TPCcrossedRows;
  uint8_t ITSclsMap;
  uint8_t TPCnCls;
  uint8_t TPCnClsShared;
  uint8_t ITSnCls;
  uint32_t clusterSizesITS;
};

struct NucleusCandidateFlow {
  float centFV0A;
  float centFT0M;
  float centFT0A;
  float centFT0C;
  float psiFT0A;
  float multFT0A;
  float psiFT0C;
  float multFT0C;
  float psiTPC;
  float psiTPCl;
  float psiTPCr;
  int multTPC;
};

namespace nuclei
{
constexpr double bbMomScalingDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double betheBlochDefault[5][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-239.99, 1.155, 1.099, 1.137, 1.006, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09},
  {-586.66, 1.859, 4.435, 0.282, 3.201, 0.09}};
constexpr double nSigmaTPCdefault[5][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
// constexpr double nSigmaTOFdefault[5][2]{
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.}};
constexpr double DCAcutDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr int TreeConfigDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr int FlowHistDefault[5][1]{
  {0},
  {0},
  {0},
  {0},
  {0}};
constexpr int DCAHistDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr double DownscalingDefault[5][1]{
  {1.},
  {1.},
  {1.},
  {1.},
  {1.}};
// constexpr bool storeTreesDefault[5]{false, false, false, false, false};
constexpr int species{5};
constexpr int codes[5]{2212, 1000010020, 1000010030, 1000020030, 1000020040};
constexpr float charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};
constexpr float masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};
static const std::vector<std::string> matter{"M", "A"};
static const std::vector<std::string> pidName{"TPC", "TOF"};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> flowConfigNames{"Save flow hists"};
static const std::vector<std::string> DCAConfigNames{"Save DCA hist", "Matter/Antimatter"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};
static const std::vector<std::string> DownscalingConfigName{"Fraction of kept candidates"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> binnedVariableNames{"DCAxy", "DCAz", "TPCnsigma", "TOFnsigma", "TOFmass"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};

float pidCuts[2][5][2];
std::shared_ptr<TH3> hNsigma[2][5][2];
std::shared_ptr<TH3> hTOFmass[5][2];
std::shared_ptr<TH2> hGenNuclei[5][2];
std::shared_ptr<TH3> hMomRes[5][2];
std::shared_ptr<TH3> hNsigmaEta[2][5][2];
std::shared_ptr<TH3> hTOFmassEta[5][2];
std::shared_ptr<TH3> hDCAxy[2][5][2];
std::shared_ptr<TH3> hDCAz[2][5][2];
std::shared_ptr<TH2> hGloTOFtracks[2];
std::shared_ptr<TH2> hDeltaP[2][5];
std::shared_ptr<THnSparse> hFlowHists[2][5];
std::shared_ptr<THnSparse> hDCAHists[2][5];
o2::base::MatLayerCylSet* lut = nullptr;

std::vector<NucleusCandidate> candidates;
std::vector<NucleusCandidateFlow> candidates_flow;

enum centDetectors {
  kFV0A = 0,
  kFT0M = 1,
  kFT0A = 2,
  kFT0C = 3
};

static const std::vector<std::string> centDetectorNames{"FV0A", "FT0M", "FT0A", "FT0C"};
} // namespace nuclei

struct nucleiSpectra {
  enum {
    kProton = BIT(0),
    kDeuteron = BIT(1),
    kTriton = BIT(2),
    kHe3 = BIT(3),
    kHe4 = BIT(4),
    kHasTOF = BIT(5),
    kHasTRD = BIT(6),
    kIsAmbiguous = BIT(7), /// just a placeholder now
    kITSrof = BIT(8),
    kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
    kIsSecondaryFromMaterial = BIT(10),
    kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
  };

  Produces<o2::aod::NucleiTable> nucleiTable;
  Produces<o2::aod::NucleiTableMC> nucleiTableMC;
  Produces<o2::aod::NucleiTableFlow> nucleiTableFlow;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  TrackTuner trackTunerObj;

  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};

  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<float> cfgCMrapidity{"cfgCMrapidity", 0.f, "Rapidity of the center of mass (only for p-Pb)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutTpcMom{"cfgCutTpcMom", 0.2f, "Minimum TPC momentum for tracks"};
  Configurable<float> cfgCutRapidityMin{"cfgCutRapidityMin", -0.5, "Minimum rapidity for tracks"};
  Configurable<float> cfgCutRapidityMax{"cfgCutRapidityMax", 0.5, "Maximum rapidity for tracks"};
  Configurable<bool> cfgCutOnReconstructedRapidity{"cfgCutOnReconstructedRapidity", false, "Cut on reconstructed rapidity"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 5, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};

  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {nuclei::bbMomScalingDefault[0], 5, 2, nuclei::names, nuclei::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], 5, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTPC{"cfgNsigmaTPC", {nuclei::nSigmaTPCdefault[0], 5, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgDCAcut{"cfgDCAcut", {nuclei::DCAcutDefault[0], 5, 2, nuclei::names, nuclei::nDCAConfigName}, "Max DCAxy and DCAz for light nuclei"};
  Configurable<LabeledArray<double>> cfgDownscaling{"cfgDownscaling", {nuclei::DownscalingDefault[0], 5, 1, nuclei::names, nuclei::DownscalingConfigName}, "Fraction of kept candidates for light nuclei"};
  Configurable<LabeledArray<int>> cfgTreeConfig{"cfgTreeConfig", {nuclei::TreeConfigDefault[0], 5, 2, nuclei::names, nuclei::treeConfigNames}, "Filtered trees configuration"};
  Configurable<LabeledArray<int>> cfgDCAHists{"cfgDCAHists", {nuclei::DCAHistDefault[0], 5, 2, nuclei::names, nuclei::DCAConfigNames}, "DCA hist configuration"};
  Configurable<LabeledArray<int>> cfgFlowHist{"cfgFlowHist", {nuclei::FlowHistDefault[0], 5, 1, nuclei::names, nuclei::flowConfigNames}, "Flow hist configuration"};

  ConfigurableAxis cfgDCAxyBinsProtons{"cfgDCAxyBinsProtons", {1500, -1.5f, 1.5f}, "DCAxy binning for Protons"};
  ConfigurableAxis cfgDCAxyBinsDeuterons{"cfgDCAxyBinsDeuterons", {1500, -1.5f, 1.5f}, "DCAxy binning for Deuterons"};
  ConfigurableAxis cfgDCAxyBinsTritons{"cfgDCAxyBinsTritons", {1500, -1.5f, 1.5f}, "DCAxy binning for Tritons"};
  ConfigurableAxis cfgDCAxyBinsHe3{"cfgDCAxyBinsHe3", {1500, -1.5f, 1.5f}, "DCAxy binning for He3"};
  ConfigurableAxis cfgDCAxyBinsAlpha{"cfgDCAxyBinsAlpha", {1500, -1.5f, 1.5f}, "DCAxy binning for Alpha"};

  ConfigurableAxis cfgPtBinsProtons{"cfgPtBinsProtons", {100, 0., 10.}, "Pt binning for Protons"};
  ConfigurableAxis cfgPtBinsDeuterons{"cfgPtBinsDeuterons", {100, 0., 10.}, "Pt binning for Deuterons"};
  ConfigurableAxis cfgPtBinsTritons{"cfgPtBinsTritons", {100, 0., 10.}, "Pt binning for Tritons"};
  ConfigurableAxis cfgPtBinsHe3{"cfgPtBinsHe3", {100, 0., 10.}, "Pt binning for He3"};
  ConfigurableAxis cfgPtBinsAlpha{"cfgPtBinsAlpha", {100, 0., 10.}, "Pt binning for Alpha"};

  ConfigurableAxis cfgCentralityBins{"cfgCentralityBins", {100, 0., 100.}, "Centrality binning"};
  ConfigurableAxis cfgMomResBins{"cfgMomResBins", {200, -1., 1.}, "Momentum resolution binning"};
  ConfigurableAxis cfgNsigmaTPCbins{"cfgNsigmaTPCbins", {100, -5., 5.}, "nsigma_TPC binning"};
  ConfigurableAxis cfgNsigmaTOFbins{"cfgNsigmaTOFbins", {100, -5., 5.}, "nsigma_TOF binning"};
  ConfigurableAxis cfgTOFmassBins{"cfgTOFmassBins", {200, -5., 5.}, "TOF mass binning"};
  ConfigurableAxis cfgV2Bins{"cfgV2Bins", {100, -1.f, 1.f}, "Binning for v2"};
  ConfigurableAxis cfgNITSClusBins{"cfgNITSClusBins", {3, 4.5, 7.5}, "N ITS clusters binning"};
  ConfigurableAxis cfgNTPCClusBins{"cfgNTPCClusBins", {3, 89.5, 159.5}, "N TPC clusters binning"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};

  // configurables for track tuner
  Configurable<bool> cfgUseTrackTuner{"cfgUseTrackTuner", false, "Apply track tuner corrections to MC tracks"};
  Configurable<std::string> cfgTrackTunerParams{"cfgTrackTunerParams", "debugInfo=0|updateTrackDCAs=1|updateTrackCovMat=1|updateCurvature=0|updateCurvatureIU=0|updatePulls=1|isInputFileFromCCDB=1|pathInputFile=Users/m/mfaggin/test/inputsTrackTuner/pp2023/smoothHighPtMC|nameInputFile=trackTuner_DataLHC23fPass1_McLHC23k4b_run535085.root|pathFileQoverPt=Users/h/hsharma/qOverPtGraphs|nameFileQoverPt=D0sigma_Data_removal_itstps_MC_LHC22b1b.root|usePvRefitCorrections=0|qOverPtMC=-1.|qOverPtData=-1.", "TrackTuner parameter initialization (format: <name>=<value>|<name>=<value>)"};
  // running variables for track tuner
  o2::dataformats::DCA mDcaInfoCov;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // CCDB options
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  int mRunNumber = 0;
  float mBz = 0.f;

  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta && aod::track::tpcInnerParam > cfgCutTpcMom;

  using TrackCandidates = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>>;

  // Collisions with chentrality
  using CollWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>::iterator;

  // Flow analysis
  using CollWithEP = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::EPCalibrationTables>::iterator;

  using CollWithQvec = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs>::iterator;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  double getPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  double computeAbsoDecL(aod::McParticles::iterator particle)
  {
    if (!particle.has_daughters())
      return -1.f;

    float mothVtx[3]{particle.vx(), particle.vy(), particle.vz()};
    float dauVtx[3]{0.f, 0.f, 0.f};
    auto daughters = particle.daughters_as<aod::McParticles>();
    for (const auto& dau : daughters) {
      if (abs(dau.pdgCode()) != 22 && abs(dau.pdgCode()) != 11) {
        dauVtx[0] = dau.vx();
        dauVtx[1] = dau.vy();
        dauVtx[2] = dau.vz();
        break;
      }
    }
    return std::hypot(mothVtx[0] - dauVtx[0], mothVtx[1] - dauVtx[1], mothVtx[2] - dauVtx[2]);
  }

  float computeEventPlane(float y, float x)
  {
    return 0.5 * std::atan2(y, x);
  }

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.selection_bit(aod::evsel::kIsTriggerTVX) && collision.posZ() > -cfgCutVertex && collision.posZ() < cfgCutVertex && collision.selection_bit(aod::evsel::kNoTimeFrameBorder);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fHe");
      zorro.populateHistRegistry(spectra, bc.runNumber());
    }
    auto timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(nuclei::lut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  void init(o2::framework::InitContext&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)nuclei::centDetectorNames[cfgCentralityEstimator])};
    const AxisSpec nSigmaAxes[2]{{cfgNsigmaTPCbins, "n#sigma_{TPC}"}, {cfgNsigmaTOFbins, "n#sigma_{TOF}"}};
    const AxisSpec tofMassAxis{cfgTOFmassBins, "TOF mass - PDG mass"};
    const AxisSpec ptResAxis{cfgMomResBins, "#Delta#it{p}_{T}/#it{p}_{T}"};
    const AxisSpec v2Axis{cfgV2Bins, "cos(2(#phi - #Psi_{2}))"};
    const AxisSpec nITSClusAxis{cfgNITSClusBins, "N ITS clusters"};
    const AxisSpec nTPCClusAxis{cfgNTPCClusBins, "N TPC clusters"};
    const AxisSpec hasTRDAxis{2, -0.5, 1.5, "Has TRD"};
    const AxisSpec correctPVAxis{2, -0.5, 1.5, "Correct PV"};
    const AxisSpec isSecondaryAxis{2, -0.5, 1.5, "Is secondary"};
    const AxisSpec fromWeakDecayAxis{2, -0.5, 1.5, "From weak decay"};

    const AxisSpec ptAxes[5]{
      {cfgPtBinsProtons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsDeuterons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsTritons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsHe3, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsAlpha, "#it{p}_{T} (GeV/#it{c})"}};
    const AxisSpec dcaxyAxes[5]{
      {cfgDCAxyBinsProtons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsDeuterons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsTritons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsHe3, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsAlpha, "DCA_{xy} (cm)"}};
    const AxisSpec dcazAxes[5]{
      {cfgDCAxyBinsProtons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsDeuterons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsTritons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsHe3, "DCA_{z} (cm)"},
      {cfgDCAxyBinsAlpha, "DCA_{z} (cm)"}};
    const AxisSpec etaAxis{40, -1., 1., "#eta"};

    spectra.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    if (doprocessMC) {
      spectra.add("hGenVtxZ", " generated collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    }
    spectra.add("hRecVtxZDataITSrof", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    spectra.add("hTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTpcSignalDataSelected", "Specific energy loss for selected particles", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTofSignalData", "TOF beta", HistType::kTH2F, {{500, 0., 5., "#it{p} (GeV/#it{c})"}, {750, 0, 1.5, "TOF #beta"}});
    for (int iC{0}; iC < 2; ++iC) {
      nuclei::hGloTOFtracks[iC] = spectra.add<TH2>(fmt::format("hTPCTOFtracks{}", nuclei::matter[iC]).data(), fmt::format("Global vs TOF matched {} tracks in a collision", nuclei::chargeLabelNames[iC]).data(), HistType::kTH2D, {{300, -0.5, 300.5, "Number of global tracks"}, {300, -0.5, 300.5, "Number of TOF matched tracks"}});

      for (int iS{0}; iS < nuclei::species; ++iS) {
        nuclei::hNsigma[0][iS][iC] = spectra.add<TH3>(fmt::format("h{}nsigma{}_{}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], nSigmaAxes[0]});
        nuclei::hNsigmaEta[0][iS][iC] = spectra.add<TH3>(fmt::format("h{}nsigmaEta{}_{}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {} vs #eta", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {etaAxis, ptAxes[iS], nSigmaAxes[0]});

        for (int iPID{0}; iPID < 2; ++iPID) {
          nuclei::hDCAxy[iPID][iS][iC] = spectra.add<TH3>(fmt::format("hDCAxy{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAxy {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], dcaxyAxes[iS]});
          nuclei::hDCAz[iPID][iS][iC] = spectra.add<TH3>(fmt::format("hDCAz{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAz {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], dcazAxes[iS]});
        }
        nuclei::hTOFmass[iS][iC] = spectra.add<TH3>(fmt::format("h{}TOFmass{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], tofMassAxis});
        nuclei::hTOFmassEta[iS][iC] = spectra.add<TH3>(fmt::format("h{}TOFmassEta{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH3D, {etaAxis, ptAxes[iS], tofMassAxis});
        nuclei::hDeltaP[iC][iS] = spectra.add<TH2>(fmt::format("hDeltaP{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("#Delta#it{{p}}/#it{{p}} {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {{232, 0.2, 6., "#it{{p}} (GeV/#it{{c}})"}, {200, -1.0, 1.0, "(#it{p}_{IU} - #it{p}_{TPC}) / #it{p}"}});
        if (doprocessMC) {
          nuclei::hMomRes[iS][iC] = spectra.add<TH3>(fmt::format("h{}MomRes{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Momentum resolution {}", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], ptResAxis});
          nuclei::hGenNuclei[iS][iC] = spectra.add<TH2>(fmt::format("h{}Gen{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Generated {}", nuclei::names[iS]).data(), HistType::kTH2D, {centAxis, ptAxes[iS]});
          nuclei::hDCAHists[iC][iS] = spectra.add<THnSparse>(fmt::format("hDCAHists{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCA histograms {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTHnSparseF, {ptAxes[iS], dcaxyAxes[iS], dcazAxes[iS], nSigmaAxes[0], tofMassAxis, nITSClusAxis, nTPCClusAxis, correctPVAxis, isSecondaryAxis, fromWeakDecayAxis});
        } else {
          if (doprocessDataFlow) {
            if (cfgFlowHist->get(iS)) {
              nuclei::hFlowHists[iC][iS] = spectra.add<THnSparse>(fmt::format("hFlowHists{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Flow histograms {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTHnSparseF, {centAxis, ptAxes[iS], nSigmaAxes[0], tofMassAxis, v2Axis, nITSClusAxis, nTPCClusAxis});
            }
          }
          nuclei::hDCAHists[iC][iS] = spectra.add<THnSparse>(fmt::format("hDCAHists{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCA histograms {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTHnSparseF, {ptAxes[iS], dcaxyAxes[iS], dcazAxes[iS], nSigmaAxes[0], tofMassAxis, nITSClusAxis, nTPCClusAxis});
        }
      }
    }

    for (int iS{0}; iS < nuclei::species; ++iS) {
      for (int iMax{0}; iMax < 2; ++iMax) {
        nuclei::pidCuts[0][iS][iMax] = cfgNsigmaTPC->get(iS, iMax);
      }
    }

    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    // TrackTuner initialization
    if (cfgUseTrackTuner) {
      std::string outputStringParams = trackTunerObj.configParams(cfgTrackTunerParams);
      spectra.add("hTrackTunedTracks", outputStringParams.c_str(), HistType::kTH1F, {{1, 0.5, 1.5, ""}});
      trackTunerObj.getDcaGraphs();
    }
  }

  template <typename Tcoll>
  float getCentrality(Tcoll const& collision)
  {
    float centrality = 1.;
    if constexpr (std::is_same<Tcoll, CollWithCent>::value || std::is_same<Tcoll, CollWithEP>::value || std::is_same<Tcoll, CollWithQvec>::value) {
      if (cfgCentralityEstimator == nuclei::centDetectors::kFV0A) {
        centrality = collision.centFV0A();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0M) {
        centrality = collision.centFT0M();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0A) {
        centrality = collision.centFT0A();
      } else if (cfgCentralityEstimator == nuclei::centDetectors::kFT0C) {
        centrality = collision.centFT0C();
      } else {
        LOG(warning) << "Centrality estimator not valid. Possible values: (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3). Centrality set to 1.";
      }
    }
    return centrality;
  }

  template <typename Tcoll, typename Ttrks>
  void fillDataInfo(Tcoll const& collision, Ttrks const& tracks)
  {
    o2::pid::tof::Beta<typename Ttrks ::iterator> responseBeta;
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    gRandom->SetSeed(bc.timestamp());

    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      spectra.fill(HIST("hRecVtxZDataITSrof"), collision.posZ());
    }

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    float centrality = getCentrality(collision);

    const double bgScalings[5][2]{
      {nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / nuclei::masses[0], nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / nuclei::masses[0]},
      {nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / nuclei::masses[1], nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / nuclei::masses[1]},
      {nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / nuclei::masses[2], nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / nuclei::masses[2]},
      {nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[3], nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[3]},
      {nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[4], nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[4]}};

    int nGloTracks[2]{0, 0}, nTOFTracks[2]{0, 0};
    for (auto& track : tracks) { // start loop over tracks
      if (track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC ||
          track.tpcNClsCrossedRows() < 70 ||
          track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
          track.tpcChi2NCl() > 4.f ||
          track.itsChi2NCl() > 36.f) {
        continue;
      }

      // temporary fix: tpcInnerParam() returns the momentum in all the software tags before
      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();

      spectra.fill(HIST("hTpcSignalData"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      float nSigma[2][5]{
        {-10., -10., -10., -10., -10.},
        {0.f, 0.f, 0.f, 0.f, 0.f}}; /// then we will calibrate the TOF mass for the He3 and Alpha
      const int iC{track.sign() < 0};

      /// Checking if we have outliers in the TPC-TOF correlation
      nGloTracks[iC]++;
      if (track.hasTOF()) {
        nTOFTracks[iC]++;
      }

      bool selectedTPC[5]{false}, goodToAnalyse{false};
      std::array<float, 5> nSigmaTPC;
      for (int iS{0}; iS < nuclei::species; ++iS) {
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScalings[iS][iC]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};
        nSigma[0][iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        nSigmaTPC[iS] = nSigma[0][iS];
        selectedTPC[iS] = (nSigma[0][iS] > nuclei::pidCuts[0][iS][0] && nSigma[0][iS] < nuclei::pidCuts[0][iS][1]);
        goodToAnalyse = goodToAnalyse || selectedTPC[iS];
        if (selectedTPC[iS] && track.p() > 0.2) {
          nuclei::hDeltaP[iC][iS]->Fill(track.p(), 1 - correctedTpcInnerParam / track.p());
        }
      }
      if (!goodToAnalyse) {
        continue;
      }

      mDcaInfoCov.set(999, 999, 999, 999, 999);
      setTrackParCov(track, mTrackParCov);
      mTrackParCov.setPID(track.pidForTracking());
      if constexpr (
        requires {
          track.has_mcParticle();
        }) {
        if (cfgUseTrackTuner) {
          bool hasMcParticle = track.has_mcParticle();
          if (hasMcParticle) {
            spectra.get<TH1>(HIST("hTrackTunedTracks"))->Fill(1); // all tracks
            auto mcParticle = track.mcParticle();
            trackTunerObj.tuneTrackParams(mcParticle, mTrackParCov, matCorr, &mDcaInfoCov, spectra.get<TH1>(HIST("hTrackTunedTracks")));
          }
        }
      }

      gpu::gpustd::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, mTrackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

      float beta{responseBeta.GetBeta(track)};
      spectra.fill(HIST("hTpcSignalDataSelected"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      spectra.fill(HIST("hTofSignalData"), correctedTpcInnerParam, beta);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      uint16_t flag = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12);
      std::array<float, 5> tofMasses{-3.f, -3.f, -3.f, -3.f, -3.f};
      bool fillTree{false};
      bool fillDCAHist{false};
      bool correctPV{false};
      bool isSecondary{false};
      bool fromWeakDecay{false};

      if (track.hasTOF()) {
        flag |= kHasTOF;
      }
      if (track.hasTRD()) {
        flag |= kHasTRD;
      }
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        flag |= kITSrof;
      }
      for (int iS{0}; iS < nuclei::species; ++iS) {
        bool selectedTOF{false};
        if (std::abs(dcaInfo[1]) > cfgDCAcut->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{mTrackParCov.getPt() * nuclei::charges[iS], mTrackParCov.getEta(), mTrackParCov.getPhi(), nuclei::masses[iS]};
        float y{fvector.Rapidity() + cfgCMrapidity};
        for (int iPID{0}; iPID < 2; ++iPID) {
          if (selectedTPC[iS]) {
            if (iPID && !track.hasTOF()) {
              continue;
            } else if (iPID) {
              selectedTOF = true;
            }
            if (!cfgCutOnReconstructedRapidity || (y > cfgCutRapidityMin && y < cfgCutRapidityMax)) {
              nuclei::hDCAxy[iPID][iS][iC]->Fill(centrality, fvector.pt(), dcaInfo[0]);
              nuclei::hDCAz[iPID][iS][iC]->Fill(centrality, fvector.pt(), dcaInfo[1]);
              if (std::abs(dcaInfo[0]) < cfgDCAcut->get(iS, 0u)) {
                if (!iPID) { /// temporary exclusion of the TOF nsigma PID for the He3 and Alpha
                  nuclei::hNsigma[iPID][iS][iC]->Fill(centrality, fvector.pt(), nSigma[iPID][iS]);
                  nuclei::hNsigmaEta[iPID][iS][iC]->Fill(fvector.eta(), fvector.pt(), nSigma[iPID][iS]);
                }
                if (iPID) {
                  float charge{1};
                  if (iS == 3 || iS == 4)
                    charge = 2;
                  tofMasses[iS] = correctedTpcInnerParam * charge * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS];
                  nuclei::hTOFmass[iS][iC]->Fill(centrality, fvector.pt(), tofMasses[iS]);
                  nuclei::hTOFmassEta[iS][iC]->Fill(fvector.eta(), fvector.pt(), tofMasses[iS]);
                }

                if (cfgFlowHist->get(iS) && doprocessDataFlow) {
                  if constexpr (std::is_same<Tcoll, CollWithEP>::value) {
                    auto deltaPhiInRange = getPhiInRange(fvector.phi() - collision.psiFT0C());
                    auto v2 = std::cos(2.0 * deltaPhiInRange);
                    nuclei::hFlowHists[iC][iS]->Fill(collision.centFT0C(), fvector.pt(), nSigma[0][iS], tofMasses[iS], v2, track.itsNCls(), track.tpcNClsFound());
                  }
                } else if (cfgFlowHist->get(iS) && doprocessDataFlowAlternative) {
                  if constexpr (std::is_same<Tcoll, CollWithQvec>::value) {
                    auto deltaPhiInRange = getPhiInRange(fvector.phi() - computeEventPlane(collision.qvecFT0CIm(), collision.qvecFT0CRe()));
                    auto v2 = std::cos(2.0 * deltaPhiInRange);
                    nuclei::hFlowHists[iC][iS]->Fill(collision.centFT0C(), fvector.pt(), nSigma[0][iS], tofMasses[iS], v2, track.itsNCls(), track.tpcNClsFound());
                  }
                }
              }
            }
          }
        }
        if (selectedTPC[iS]) {
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          !fillTree && cfgTreeConfig->get(iS, 0u) ? fillTree = true : fillTree;
          !fillDCAHist && cfgDCAHists->get(iS, iC) ? fillDCAHist = true : fillDCAHist;
          bool setPartFlag = cfgTreeConfig->get(iS, 0u) || cfgDCAHists->get(iS, iC);
          if (setPartFlag) {
            if (cfgDownscaling->get(iS) < 1. && gRandom->Rndm() > cfgDownscaling->get(iS)) {
              continue;
            }
            flag |= BIT(iS);
          }
        }
      }
      if (flag & (kProton | kDeuteron | kTriton | kHe3 | kHe4) || doprocessMC) { /// ignore PID pre-selections for the MC
        if constexpr (std::is_same<Tcoll, CollWithEP>::value) {
          nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
            collision.centFV0A(),
            collision.centFT0M(),
            collision.centFT0A(),
            collision.centFT0C(),
            collision.psiFT0A(),
            collision.multFT0A(),
            collision.psiFT0C(),
            collision.multFT0C(),
            collision.psiTPC(),
            collision.psiTPCL(),
            collision.psiTPCR(),
            collision.multTPC()});
        } else if constexpr (std::is_same<Tcoll, CollWithQvec>::value) {
          nuclei::candidates_flow.emplace_back(NucleusCandidateFlow{
            collision.centFV0A(),
            collision.centFT0M(),
            collision.centFT0A(),
            collision.centFT0C(),
            computeEventPlane(collision.qvecFT0AIm(), collision.qvecFT0ARe()),
            collision.multFT0A(),
            computeEventPlane(collision.qvecFT0CIm(), collision.qvecFT0CRe()),
            collision.multFT0C(),
            -999.,
            computeEventPlane(collision.qvecBNegIm(), collision.qvecBNegRe()),
            computeEventPlane(collision.qvecBPosIm(), collision.qvecBPosRe()),
            collision.multTPC()});
        }
        nuclei::candidates.emplace_back(NucleusCandidate{
          static_cast<int>(track.globalIndex()), static_cast<int>(track.collisionId()), (1 - 2 * iC) * mTrackParCov.getPt(), mTrackParCov.getEta(), mTrackParCov.getPhi(),
          correctedTpcInnerParam, beta, collision.posZ(), dcaInfo[0], dcaInfo[1], track.tpcSignal(), track.itsChi2NCl(), track.tpcChi2NCl(),
          nSigmaTPC, tofMasses, fillTree, fillDCAHist, correctPV, isSecondary, fromWeakDecay, flag, track.tpcNClsFindable(), static_cast<uint8_t>(track.tpcNClsCrossedRows()), track.itsClusterMap(),
          static_cast<uint8_t>(track.tpcNClsFound()), static_cast<uint8_t>(track.tpcNClsShared()), static_cast<uint8_t>(track.itsNCls()), static_cast<uint32_t>(track.itsClusterSizes())});
      }
    } // end loop over tracks

    nuclei::hGloTOFtracks[0]->Fill(nGloTracks[0], nTOFTracks[0]);
    nuclei::hGloTOFtracks[1]->Fill(nGloTracks[1], nTOFTracks[1]);
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    if (!eventSelection(collision)) {
      return;
    }
    if (cfgSkimmedProcessing) {
      zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
    }

    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      if (c.fillTree) {
        nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS);
      }
      if (c.fillDCAHist) {
        for (int iS{0}; iS < nuclei::species; ++iS) {
          if (c.flags & BIT(iS)) {
            nuclei::hDCAHists[c.pt < 0][iS]->Fill(std::abs(c.pt), c.DCAxy, c.DCAz, c.nSigmaTPC[iS], c.tofMasses[iS], c.ITSnCls, c.TPCnCls);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processData, "Data analysis", true);

  void processDataFlow(CollWithEP const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    nuclei::candidates_flow.clear();
    if (!eventSelection(collision)) {
      return;
    }
    if (!collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      if (c.fillTree) {
        nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS);
      }
      if (c.fillDCAHist) {
        for (int iS{0}; iS < nuclei::species; ++iS) {
          if (c.flags & BIT(iS)) {
            nuclei::hDCAHists[c.pt < 0][iS]->Fill(std::abs(c.pt), c.DCAxy, c.DCAz, c.nSigmaTPC[iS], c.tofMasses[iS], c.ITSnCls, c.TPCnCls);
          }
        }
      }
    }
    for (auto& c : nuclei::candidates_flow) {
      nucleiTableFlow(c.centFV0A, c.centFT0M, c.centFT0A, c.centFT0C, c.psiFT0A, c.multFT0A, c.psiFT0C, c.multFT0C, c.psiTPC, c.psiTPCl, c.psiTPCr, c.multTPC);
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processDataFlow, "Data analysis with flow", false);

  void processDataFlowAlternative(CollWithQvec const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    nuclei::candidates_flow.clear();
    if (!eventSelection(collision)) {
      return;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      if (c.fillTree) {
        nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS);
      }
      if (c.fillDCAHist) {
        for (int iS{0}; iS < nuclei::species; ++iS) {
          if (c.flags & BIT(iS)) {
            nuclei::hDCAHists[c.pt < 0][iS]->Fill(std::abs(c.pt), c.DCAxy, c.DCAz, c.nSigmaTPC[iS], c.tofMasses[iS], c.ITSnCls, c.TPCnCls);
          }
        }
      }
    }
    for (auto& c : nuclei::candidates_flow) {
      nucleiTableFlow(c.centFV0A, c.centFT0M, c.centFT0A, c.centFT0C, c.psiFT0A, c.multFT0A, c.psiFT0C, c.multFT0C, c.psiTPC, c.psiTPCl, c.psiTPCr, c.multTPC);
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processDataFlowAlternative, "Data analysis with flow - alternative framework", false);

  Preslice<TrackCandidates> tracksPerCollisions = aod::track::collisionId;
  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mcCollisions, soa::Join<TrackCandidates, aod::McTrackLabels> const& tracks, aod::McParticles const& particlesMC, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    for (auto& c : mcCollisions) {
      spectra.fill(HIST("hGenVtxZ"), c.posZ());
    }
    std::vector<bool> goodCollisions(mcCollisions.size(), false);
    for (auto& collision : collisions) {
      if (!eventSelection(collision)) {
        continue;
      }
      goodCollisions[collision.mcCollisionId()] = true;
      const auto& slicedTracks = tracks.sliceBy(tracksPerCollisions, collision.globalIndex());
      fillDataInfo(collision, slicedTracks);
    }
    std::vector<bool> isReconstructed(particlesMC.size(), false);
    for (auto& c : nuclei::candidates) {
      auto label = tracks.iteratorAt(c.globalIndex);
      if (label.mcParticleId() < -1 || label.mcParticleId() >= particlesMC.size()) {
        continue;
      }
      auto particle = particlesMC.iteratorAt(label.mcParticleId());
      bool storeIt{false};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (std::abs(particle.pdgCode()) == nuclei::codes[iS]) {
          if (c.fillTree && !storeIt) {
            nuclei::hMomRes[iS][particle.pdgCode() < 0]->Fill(1., std::abs(c.pt * nuclei::charges[iS]), 1. - std::abs(c.pt * nuclei::charges[iS]) / particle.pt());
            storeIt = cfgTreeConfig->get(iS, 0u) || cfgTreeConfig->get(iS, 1u); /// store only the particles of interest
          }
          auto coll = collisions.iteratorAt(c.collTrackIndex);
          int collMCGlobId = coll.mcCollisionId();
          if (particle.mcCollisionId() == collMCGlobId) {
            c.correctPV = true;
          }
          if (!particle.isPhysicalPrimary()) {
            c.isSecondary = true;
            if (particle.getProcess() == 4) {
              c.fromWeakDecay = true;
            }
          }
          if (c.fillDCAHist && cfgDCAHists->get(iS, c.pt < 0)) {
            nuclei::hDCAHists[c.pt < 0][iS]->Fill(std::abs(c.pt), c.DCAxy, c.DCAz, c.nSigmaTPC[iS], c.tofMasses[iS], c.ITSnCls, c.TPCnCls, c.correctPV, c.isSecondary, c.fromWeakDecay);
          }
        }
      }
      if (!storeIt) {
        continue;
      }
      isReconstructed[particle.globalIndex()] = true;
      if (particle.isPhysicalPrimary()) {
        c.flags |= kIsPhysicalPrimary;
      } else if (particle.has_mothers()) {
        c.flags |= kIsSecondaryFromWeakDecay;
      } else {
        c.flags |= kIsSecondaryFromMaterial;
      }
      float absoDecL = computeAbsoDecL(particle);
      nucleiTableMC(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS, particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), goodCollisions[particle.mcCollisionId()], absoDecL);
    }

    int index{0};
    for (auto& particle : particlesMC) {
      int pdg{std::abs(particle.pdgCode())};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (pdg != nuclei::codes[iS]) {
          continue;
        }
        uint16_t flags{kIsPhysicalPrimary};
        if (particle.isPhysicalPrimary()) {
          if (particle.y() > cfgCutRapidityMin && particle.y() < cfgCutRapidityMax) {
            nuclei::hGenNuclei[iS][particle.pdgCode() < 0]->Fill(1., particle.pt());
          }
        } else {
          continue; /// for not-reconstructed particles we store only the primaries
        }

        if (!isReconstructed[index] && (cfgTreeConfig->get(iS, 0u) || cfgTreeConfig->get(iS, 1u))) {
          float absDecL = computeAbsoDecL(particle);
          nucleiTableMC(999., 999., 999., 0., 0., 999., 999., 999., -1, -1, -1, flags, 0, 0, 0, 0, 0, 0, particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), goodCollisions[particle.mcCollisionId()], absDecL);
        }
        break;
      }
      index++;
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processMC, "MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiSpectra>(cfgc, TaskName{"nuclei-spectra"})};
}
