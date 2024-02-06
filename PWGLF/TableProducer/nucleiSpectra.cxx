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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/Qvectors.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/Track.h"

#include "PWGLF/DataModel/LFSlimNucleiTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct NucleusCandidate {
  int globalIndex;
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
  uint16_t flags;
  uint8_t TPCfindableCls;
  uint8_t TPCcrossedRows;
  uint8_t ITSclsMap;
  uint8_t TPCnCls;
  uint32_t clusterSizesITS;
  int selCollIndex;
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
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
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
// constexpr bool storeTreesDefault[5]{false, false, false, false, false};
constexpr int species{5};
constexpr int codes[5]{2212, 1000010020, 1000010030, 1000020030, 1000020040};
constexpr float charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};
constexpr float masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};
static const std::vector<std::string> matter{"M", "A"};
static const std::vector<std::string> pidName{"TPC", "TOF"};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};
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
o2::base::MatLayerCylSet* lut = nullptr;

std::vector<NucleusCandidate> candidates;

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
    kIsReconstructed = BIT(6),
    kIsAmbiguous = BIT(7), /// just a placeholder now
    kPositive = BIT(8),
    kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
    kIsSecondaryFromMaterial = BIT(10),
    kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
  };

  Produces<o2::aod::NucleiTable> nucleiTable;
  Produces<o2::aod::NucleiTableMC> nucleiTableMC;
  Produces<o2::aod::NucleiFlowColls> nucleiFlowTable;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

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
  Configurable<LabeledArray<int>> cfgTreeConfig{"cfgTreeConfig", {nuclei::TreeConfigDefault[0], 5, 2, nuclei::names, nuclei::treeConfigNames}, "Filtered trees configuration"};

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
  ConfigurableAxis cfgSpBins{"cfgSpBins", {100, -1.f, 1.f}, "Binning for scalar product"};

  // CCDB options
  Configurable<double> cfgBz{"cfgBz", -999, "bz field, -999 is automatic"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgGRPpath{"cfgGRPpath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> cfgGRPmagPath{"cfgGRPmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> cfgLUTpath{"cfgLUTpath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> cfgGeoPath{"cfgGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  int mRunNumber = 0;
  float mBz = 0.f;

  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta && aod::track::tpcInnerParam > cfgCutTpcMom;

  using TrackCandidates = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>>;

  // Collisions with chentrality
  using CollWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>::iterator;

  // Flow analysis
  using CollWithQvec = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>::iterator;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  o2::pid::tof::Beta<TrackCandidates::iterator> responseBeta;

  template <class collision_t>
  bool eventSelection(collision_t& collision)
  {
    return collision.sel8() && collision.posZ() > -cfgCutVertex && collision.posZ() < cfgCutVertex;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    if (cfgBz > -990) {
      mBz = cfgBz;
    } else {
      o2::parameters::GRPObject* grpo{ccdb->getForTimeStamp<o2::parameters::GRPObject>(cfgGRPpath, run3grp_timestamp)};
      o2::parameters::GRPMagField* grpmag{nullptr};
      if (grpo) {
        mBz = grpo->getNominalL3Field();
      } else {
        grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(cfgGRPmagPath, run3grp_timestamp);
        if (!grpmag) {
          LOG(fatal) << "Got nullptr from CCDB for path " << cfgGRPmagPath << " of object GRPMagField and " << cfgGRPpath << " of object GRPObject for timestamp " << run3grp_timestamp;
        }
        mBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      }
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << mBz << " kZG";
    }
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)nuclei::centDetectorNames[cfgCentralityEstimator])};
    const AxisSpec nSigmaAxes[2]{{cfgNsigmaTPCbins, "n#sigma_{TPC}"}, {cfgNsigmaTOFbins, "n#sigma_{TOF}"}};
    const AxisSpec tofMassAxis{cfgTOFmassBins, "TOF mass - PDG mass"};
    const AxisSpec ptResAxis{cfgMomResBins, "#Delta#it{p}_{T}/#it{p}_{T}"};

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

    const AxisSpec ft0Aft0CspAxis{cfgSpBins, "#vec{Q}_{2}^{FT0A} #upoint #vec{Q}_{2}^{FT0C}"};
    const AxisSpec fv0Aft0CspAxis{cfgSpBins, "#vec{Q}_{2}^{FV0A} #upoint #vec{Q}_{2}^{FT0C}"};
    const AxisSpec fv0Aft0AspAxis{cfgSpBins, "#vec{Q}_{2}^{FV0A} #upoint #vec{Q}_{2}^{FT0A}"};

    spectra.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
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
        nuclei::hMomRes[iS][iC] = spectra.add<TH3>(fmt::format("h{}MomRes{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Momentum resolution {}", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], ptResAxis});
        nuclei::hGenNuclei[iS][iC] = spectra.add<TH2>(fmt::format("h{}Gen{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Generated {}", nuclei::names[iS]).data(), HistType::kTH2D, {centAxis, ptAxes[iS]});
        nuclei::hDeltaP[iC][iS] = spectra.add<TH2>(fmt::format("hDeltaP{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("#Delta#it{{p}}/#it{{p}} {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {{232, 0.2, 6., "#it{{p}} (GeV/#it{{c}})"}, {200, -1.0, 1.0, "(#it{p}_{IU} - #it{p}_{TPC}) / #it{p}"}});
      }
    }
    if (doprocessDataFlow) {
      spectra.add("hScalarProductFT0AvsFT0C", "", HistType::kTH2F, {ft0Aft0CspAxis, centAxis});
      spectra.add("hScalarProductFV0AvsFT0C", "", HistType::kTH2F, {fv0Aft0CspAxis, centAxis});
      spectra.add("hScalarProductFV0AvsFT0A", "", HistType::kTH2F, {fv0Aft0AspAxis, centAxis});
    }

    for (int iS{0}; iS < nuclei::species; ++iS) {
      for (int iMax{0}; iMax < 2; ++iMax) {
        nuclei::pidCuts[0][iS][iMax] = cfgNsigmaTPC->get(iS, iMax);
      }
    }

    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    o2::base::Propagator::Instance(true)->setMatLUT(nuclei::lut);
  }

  template <typename Tcoll>
  float getCentrality(Tcoll const& collision)
  {
    float centrality = 1.;
    if constexpr (std::is_same<Tcoll, CollWithCent>::value || std::is_same<Tcoll, CollWithQvec>::value) {
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
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    float centrality = getCentrality(collision);

    if constexpr (std::is_same<Tcoll, CollWithQvec>::value) {
      if (doprocessDataFlow) {
        spectra.fill(HIST("hScalarProductFT0AvsFT0C"), collision.qvecFT0ARe() * collision.qvecFT0CRe() + collision.qvecFT0AIm() * collision.qvecFT0CIm(), centrality);
        spectra.fill(HIST("hScalarProductFV0AvsFT0C"), collision.qvecFV0ARe() * collision.qvecFT0CRe() + collision.qvecFV0AIm() * collision.qvecFT0CIm(), centrality);
        spectra.fill(HIST("hScalarProductFV0AvsFT0A"), collision.qvecFT0ARe() * collision.qvecFV0ARe() + collision.qvecFT0AIm() * collision.qvecFV0AIm(), centrality);
      }
    }

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
      spectra.fill(HIST("hTpcSignalData"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
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
      for (int iS{0}; iS < nuclei::species; ++iS) {
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScalings[iS][iC]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};
        nSigma[0][iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        selectedTPC[iS] = (nSigma[0][iS] > nuclei::pidCuts[0][iS][0] && nSigma[0][iS] < nuclei::pidCuts[0][iS][1]);
        goodToAnalyse = goodToAnalyse || selectedTPC[iS];
        if (selectedTPC[iS] && track.p() > 0.2) {
          nuclei::hDeltaP[iC][iS]->Fill(track.p(), 1 - track.tpcInnerParam() / track.p());
        }
      }
      if (!goodToAnalyse) {
        continue;
      }

      auto trackParCov = getTrackParCov(track); // should we set the charge according to the nucleus?
      gpu::gpustd::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

      float beta{responseBeta.GetBeta(track)};
      spectra.fill(HIST("hTpcSignalDataSelected"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      spectra.fill(HIST("hTofSignalData"), track.tpcInnerParam(), beta);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta)); /// sometimes beta > 1 or < 0, to be checked
      uint16_t flag = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12);
      if (track.hasTOF()) {
        flag |= kHasTOF;
      }
      if (!iC) {
        flag |= kPositive;
      }
      for (int iS{0}; iS < nuclei::species; ++iS) {
        bool selectedTOF{false};
        if (std::abs(dcaInfo[1]) > cfgDCAcut->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{trackParCov.getPt() * nuclei::charges[iS], trackParCov.getEta(), trackParCov.getPhi(), nuclei::masses[iS]};
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
                  float mass{track.tpcInnerParam() * nuclei::charges[iS] * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS]};
                  nuclei::hTOFmass[iS][iC]->Fill(centrality, fvector.pt(), mass);
                  nuclei::hTOFmassEta[iS][iC]->Fill(fvector.eta(), fvector.pt(), mass);
                }
              }
            }
          }
        }

        if (cfgTreeConfig->get(iS, 0u) && selectedTPC[iS]) {
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          flag |= BIT(iS);
        }
      }
      if (flag & (kProton | kDeuteron | kTriton | kHe3 | kHe4)) {
        if constexpr (std::is_same<Tcoll, CollWithQvec>::value) {
          if (nuclei::candidates.empty()) {
            nucleiFlowTable(collision.centFV0A(),
                            collision.centFT0M(),
                            collision.centFT0A(),
                            collision.centFT0C(),
                            collision.qvecFV0ARe(),
                            collision.qvecFV0AIm(),
                            collision.sumAmplFV0A(),
                            collision.qvecFT0MRe(),
                            collision.qvecFT0MIm(),
                            collision.sumAmplFT0M(),
                            collision.qvecFT0ARe(),
                            collision.qvecFT0AIm(),
                            collision.sumAmplFT0A(),
                            collision.qvecFT0CRe(),
                            collision.qvecFT0CIm(),
                            collision.sumAmplFT0C(),
                            collision.qvecBPosRe(),
                            collision.qvecBPosIm(),
                            collision.nTrkBPos(),
                            collision.qvecBNegRe(),
                            collision.qvecBNegIm(),
                            collision.nTrkBNeg());
          }
        }
        nuclei::candidates.emplace_back(NucleusCandidate{static_cast<int>(track.globalIndex()), (1 - 2 * iC) * trackParCov.getPt(), trackParCov.getEta(), trackParCov.getPhi(), track.tpcInnerParam(), beta, collision.posZ(), dcaInfo[0], dcaInfo[1], track.tpcSignal(), track.itsChi2NCl(),
                                                         track.tpcChi2NCl(), flag, track.tpcNClsFindable(), static_cast<uint8_t>(track.tpcNClsCrossedRows()), track.itsClusterMap(), static_cast<uint8_t>(track.tpcNClsFound()), static_cast<uint32_t>(track.itsClusterSizes()), static_cast<int>(nucleiFlowTable.lastIndex())});
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
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.clusterSizesITS, c.selCollIndex);
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processData, "Data analysis", true);

  void processDataFlow(CollWithQvec const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    if (!eventSelection(collision)) {
      return;
    }
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      nucleiTable(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.clusterSizesITS, c.selCollIndex);
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processDataFlow, "Data analysis with flow", false);

  Preslice<TrackCandidates> tracksPerCollisions = aod::track::collisionId;
  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mcCollisions, TrackCandidates const& tracks, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
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
      auto label = trackLabelsMC.iteratorAt(c.globalIndex);
      if (label.mcParticleId() < -1 || label.mcParticleId() >= particlesMC.size()) {
        continue;
      }
      auto particle = particlesMC.iteratorAt(label.mcParticleId());
      isReconstructed[particle.globalIndex()] = true;
      if (particle.isPhysicalPrimary()) {
        c.flags |= kIsPhysicalPrimary;
      } else if (particle.has_mothers()) {
        c.flags |= kIsSecondaryFromWeakDecay;
      } else {
        c.flags |= kIsSecondaryFromMaterial;
      }

      nucleiTableMC(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.clusterSizesITS, particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), goodCollisions[particle.mcCollisionId()]);
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (std::abs(particle.pdgCode()) == nuclei::codes[iS]) {
          nuclei::hMomRes[iS][particle.pdgCode() < 0]->Fill(1., std::abs(c.pt * nuclei::charges[iS]), 1. - std::abs(c.pt * nuclei::charges[iS]) / particle.pt());
          break;
        }
      }
    }

    int index{0};
    for (auto& particle : particlesMC) {
      int pdg{std::abs(particle.pdgCode())};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (pdg != nuclei::codes[iS]) {
          continue;
        }
        uint16_t flags{0u};
        if (particle.isPhysicalPrimary()) {
          flags |= kIsPhysicalPrimary;
          if (particle.y() > cfgCutRapidityMin && particle.y() < cfgCutRapidityMax) {
            nuclei::hGenNuclei[iS][particle.pdgCode() < 0]->Fill(1., particle.pt());
          }
        } else if (particle.has_mothers()) {
          flags |= kIsSecondaryFromWeakDecay;
        } else {
          flags |= kIsSecondaryFromMaterial;
        }

        if (!isReconstructed[index] && (cfgTreeConfig->get(iS, 0u) || cfgTreeConfig->get(iS, 1u))) {
          nucleiTableMC(999., 999., 999., 0., 0., 999., 999., 999., -1, -1, -1, flags, 0, 0, 0, 0, 0, particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), goodCollisions[particle.mcCollisionId()]);
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
