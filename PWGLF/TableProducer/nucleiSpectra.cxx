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

uint8_t getBinnedValue(double val, double max)
{
  if (val >= max) {
    return 255u;
  } else if (val < -max) {
    return 0u;
  } else {
    return 1u + static_cast<uint8_t>(254 * (val - max) / (2 * max));
  }
}

float getBinCenter(uint8_t bin, double max)
{
  if (bin == 0u) {
    return -max;
  } else if (bin == 255u) {
    return max;
  } else {
    return -max + (bin + 0.5) / (2 * max);
  }
}

struct NucleusCandidate {
  NucleusCandidate(int idx, float aPt, float aEta, uint8_t z, uint8_t aITSclsMap, uint8_t aTPCnCls, int8_t aDCAxy, int8_t aDCAz, uint16_t aFlags, uint8_t aTPCnsigma, uint8_t aTOFmass) : globalIndex(idx), pt(aPt), eta(aEta), zVertex(z), ITSclsMap(aITSclsMap), TPCnCls(aTPCnCls), DCAxy(aDCAxy), DCAz(aDCAz), flags(aFlags), TPCnsigma(aTPCnsigma), TOFmass(aTOFmass)
  {
  }
  int globalIndex;
  float pt;
  float eta;
  uint8_t zVertex;
  uint8_t ITSclsMap;
  uint8_t TPCnCls;
  int8_t DCAxy;
  int8_t DCAz;
  uint16_t flags;
  uint8_t TPCnsigma;
  uint8_t TOFmass;
};

namespace nuclei
{
constexpr double bbMomScalingDefault[4][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double betheBlochDefault[4][6]{
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
constexpr double nSigmaTPCdefault[4][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
// constexpr double nSigmaTOFdefault[4][2]{
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.},
//   {-5., 5.}};
constexpr double DCAcutDefault[4][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr int TreeConfigDefault[4][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr double BinnedVariablesDefaultMax[5][1]{
  {1.27},
  {2.54},
  {5.08},
  {5.08},
  {5.08}};
// constexpr bool storeTreesDefault[4]{false, false, false, false};
constexpr int species{4};
constexpr int codes[4]{1000010020, 1000010030, 1000020030, 1000020040};
constexpr float charges[4]{1.f, 1.f, 2.f, 2.f};
constexpr float masses[4]{MassDeuteron, MassTriton, MassHelium3, MassAlpha};
static const std::vector<std::string> matter{"M", "A"};
static const std::vector<std::string> pidName{"TPC", "TOF"};
static const std::vector<std::string> names{"deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> binnedVariableNames{"DCAxy", "DCAz", "TPCnsigma", "TOFnsigma", "TOFmass"};
static const std::vector<std::string> binnedLabelNames{"Maximum value of binned variables"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};

float pidCuts[2][4][2];
std::shared_ptr<TH3> hNsigma[2][4][2];
std::shared_ptr<TH3> hTOFmass[4][2];
std::shared_ptr<TH2> hGenNuclei[4][2];
std::shared_ptr<TH3> hMomRes[4][2];
std::shared_ptr<TH3> hDCAxy[2][4][2];
o2::base::MatLayerCylSet* lut = nullptr;

std::vector<NucleusCandidate> candidates;
} // namespace nuclei

struct nucleiSpectra {
  enum {
    kDeuteron = BIT(0),
    kTriton = BIT(1),
    kHe3 = BIT(2),
    kHe4 = BIT(3),
    kHasTOF = BIT(4),
    kIsReconstructed = BIT(5),
    kIsAmbiguous = BIT(6),       /// just a placeholder now
    kIsPhysicalPrimary = BIT(8), /// MC flags starting from the second half of the short
    kIsSecondaryFromMaterial = BIT(9),
    kIsSecondaryFromWeakDecay = BIT(10)
  };

  Produces<o2::aod::NucleiTable> nucleiTable;
  Produces<o2::aod::NucleiTableMC> nucleiTableMC;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> cfgCentralityEstimator{"cfgCentralityEstimator", "V0A", "Centrality estimator name"};
  Configurable<float> cfgCMrapidity{"cfgCMrapidity", 0.f, "Rapidity of the center of mass (only for p-Pb)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutRapidityMin{"cfgCutRapidityMin", -0.5, "Minimum rapidity for tracks"};
  Configurable<float> cfgCutRapidityMax{"cfgCutRapidityMax", 0.5, "Maximum rapidity for tracks"};
  Configurable<bool> cfgCutOnReconstructedRapidity{"cfgCutOnReconstructedRapidity", false, "Cut on reconstructed rapidity"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 5, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};

  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {nuclei::bbMomScalingDefault[0], 4, 2, nuclei::names, nuclei::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], 4, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTPC{"cfgNsigmaTPC", {nuclei::nSigmaTPCdefault[0], 4, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgDCAcut{"cfgDCAcut", {nuclei::DCAcutDefault[0], 4, 2, nuclei::names, nuclei::nDCAConfigName}, "Max DCAxy and DCAz for light nuclei"};
  Configurable<LabeledArray<int>> cfgTreeConfig{"cfgTreeConfig", {nuclei::TreeConfigDefault[0], 4, 2, nuclei::names, nuclei::treeConfigNames}, "Filtered trees configuration"};
  Configurable<LabeledArray<double>> cfgBinnedVariables{"cfgBinnedVariables", {nuclei::BinnedVariablesDefaultMax[0], 5, 1, nuclei::binnedVariableNames, nuclei::binnedLabelNames}, "Maximum value for the binned variables"};

  ConfigurableAxis cfgDCAxyBinsDeuterons{"cfgDCAxyBinsDeuterons", {300, -3.f, 3.f}, "DCAxy binning for Deuterons"};
  ConfigurableAxis cfgDCAxyBinsTritons{"cfgDCAxyBinsTritons", {300, -3.f, 3.f}, "DCAxy binning for Tritons"};
  ConfigurableAxis cfgDCAxyBinsHe3{"cfgDCAxyBinsHe3", {300, -3.f, 3.f}, "DCAxy binning for He3"};
  ConfigurableAxis cfgDCAxyBinsAlpha{"cfgDCAxyBinsAlpha", {300, -3.f, 3.f}, "DCAxy binning for Alpha"};

  ConfigurableAxis cfgPtBinsDeuterons{"cfgPtBinsDeuterons", {100, 0., 10.}, "Pt binning for Deuterons"};
  ConfigurableAxis cfgPtBinsTritons{"cfgPtBinsTritons", {100, 0., 10.}, "Pt binning for Tritons"};
  ConfigurableAxis cfgPtBinsHe3{"cfgPtBinsHe3", {100, 0., 10.}, "Pt binning for He3"};
  ConfigurableAxis cfgPtBinsAlpha{"cfgPtBinsAlpha", {100, 0., 10.}, "Pt binning for Alpha"};

  ConfigurableAxis cfgCentralityBins{"cfgCentralityBins", {100, 0., 100.}, "Centrality binning"};
  ConfigurableAxis cfgMomResBins{"cfgMomResBins", {200, -1., 1.}, "Momentum resolution binning"};
  ConfigurableAxis cfgNsigmaTPCbins{"cfgNsigmaTPCbins", {100, -5., 5.}, "nsigma_TPC binning"};
  ConfigurableAxis cfgNsigmaTOFbins{"cfgNsigmaTOFbins", {100, -5., 5.}, "nsigma_TOF binning"};
  ConfigurableAxis cfgTOFmassBins{"cfgTOFmassBins", {200, -5., 5.}, "TOF mass binning"};

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

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta;

  using TrackCandidates = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime>>;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  o2::pid::tof::Beta<TrackCandidates::iterator> responseBeta;

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

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)cfgCentralityEstimator)};
    const AxisSpec nSigmaAxes[2]{{cfgNsigmaTPCbins, "n#sigma_{TPC}"}, {cfgNsigmaTOFbins, "n#sigma_{TOF}"}};
    const AxisSpec tofMassAxis{cfgTOFmassBins, "TOF mass - PDG mass"};
    const AxisSpec ptResAxis{cfgMomResBins, "#Delta#it{p}_{T}/#it{p}_{T}"};

    const AxisSpec ptAxes[4]{
      {cfgPtBinsDeuterons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsTritons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsHe3, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsAlpha, "#it{p}_{T} (GeV/#it{c})"}};
    const AxisSpec dcaAxes[4]{
      {cfgDCAxyBinsDeuterons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsTritons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsHe3, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsAlpha, "DCA_{xy} (cm)"}};

    spectra.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    spectra.add("hTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTpcSignalDataSelected", "Specific energy loss for selected particles", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTofSignalData", "TOF beta", HistType::kTH2F, {{500, 0., 5., "#it{p} (GeV/#it{c})"}, {750, 0, 1.5, "TOF #beta"}});
    for (int iC{0}; iC < 2; ++iC) {
      for (int iS{0}; iS < nuclei::species; ++iS) {
        nuclei::hNsigma[0][iS][iC] = spectra.add<TH3>(fmt::format("h{}nsigma{}_{}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], nSigmaAxes[0]});

        for (int iPID{0}; iPID < 2; ++iPID) {
          nuclei::hDCAxy[iPID][iS][iC] = spectra.add<TH3>(fmt::format("hDCAxy{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAxy {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], dcaAxes[iS]});
        }
        nuclei::hTOFmass[iS][iC] = spectra.add<TH3>(fmt::format("h{}TOFmass{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], tofMassAxis});
        nuclei::hMomRes[iS][iC] = spectra.add<TH3>(fmt::format("h{}MomRes{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Momentum resolution {}", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], ptResAxis});
        nuclei::hGenNuclei[iS][iC] = spectra.add<TH2>(fmt::format("h{}Gen{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Generated {}", nuclei::names[iS]).data(), HistType::kTH2D, {centAxis, ptAxes[iS]});
      }
    }

    for (int iS{0}; iS < 4; ++iS) {
      for (int iMax{0}; iMax < 2; ++iMax) {
        nuclei::pidCuts[0][iS][iMax] = cfgNsigmaTPC->get(iS, iMax);
      }
    }

    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
    o2::base::Propagator::Instance(true)->setMatLUT(nuclei::lut);
  }

  template <typename TC>
  void fillDataInfo(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TC const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    nuclei::candidates.clear();
    // collision process loop
    if (!collision.sel8()) {
      return;
    }
    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};
    uint8_t zVert{getBinnedValue(collision.posZ(), 10.16)};

    const double bgScalings[4][2]{
      {nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / nuclei::masses[0], nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / nuclei::masses[0]},
      {nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / nuclei::masses[1], nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / nuclei::masses[1]},
      {nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / nuclei::masses[2], nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / nuclei::masses[2]},
      {nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[3], nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[3]}};

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
      float nSigma[2][4]{
        {-10., -10., -10., -10.},
        {0.f, 0.f, 0.f, 0.f}}; /// then we will calibrate the TOF mass for the He3 and Alpha
      const int iC{track.sign() < 0};
      bool selectedTPC[4]{false}, goodToAnalyse{false};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScalings[iS][iC]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};
        nSigma[0][iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        selectedTPC[iS] = (nSigma[0][iS] > nuclei::pidCuts[0][iS][0] && nSigma[0][iS] < nuclei::pidCuts[0][iS][1]);
        goodToAnalyse = goodToAnalyse || selectedTPC[iS];
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
      for (int iS{0}; iS < nuclei::species; ++iS) {
        bool selectedTOF{false};
        if (std::abs(dcaInfo[1]) > cfgDCAcut->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{trackParCov.getPt() * nuclei::charges[iS], trackParCov.getEta(), trackParCov.getPhi(), nuclei::masses[iS]};
        float y{fvector.Rapidity() + cfgCMrapidity};
        if (cfgCutOnReconstructedRapidity && (y < cfgCutRapidityMin || y > cfgCutRapidityMax)) {
          continue;
        }

        for (int iPID{0}; iPID < 2; ++iPID) {
          if (selectedTPC[iS]) {
            if (iPID && !track.hasTOF()) {
              continue;
            } else if (iPID) {
              selectedTOF = true;
            }
            nuclei::hDCAxy[iPID][iS][iC]->Fill(1., fvector.pt(), dcaInfo[0]);
            if (std::abs(dcaInfo[0]) < cfgDCAcut->get(iS, 0u)) {
              if (!iPID) { /// temporary exclusion of the TOF nsigma PID for the He3 and Alpha
                nuclei::hNsigma[iPID][iS][iC]->Fill(1., fvector.pt(), nSigma[iPID][iS]);
              }
              if (iPID) {
                float mass{track.tpcInnerParam() * nuclei::charges[iS] * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS]};
                nuclei::hTOFmass[iS][iC]->Fill(1., fvector.pt(), mass);
              }
            }
          }
        }
        uint16_t flag{kIsReconstructed};
        if (cfgTreeConfig->get(iS, 0u) && selectedTPC[iS]) {
          int8_t massTOF{0u};
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          if (track.hasTOF()) {
            flag |= kHasTOF;
            massTOF = getBinnedValue(beta > 1.e-6f ? track.tpcInnerParam() * nuclei::charges[iS] * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS] : -999.f, cfgBinnedVariables->get(4u, 1u));
          }
          flag |= BIT(iS);
          uint8_t dcaxy = getBinnedValue(dcaInfo[0], cfgBinnedVariables->get(0u, 1u));
          uint8_t dcaz = getBinnedValue(dcaInfo[1], cfgBinnedVariables->get(1u, 1u));
          uint8_t nsigmaTPC = getBinnedValue(nSigma[0][iS], cfgBinnedVariables->get(2u, 1u));

          nuclei::candidates.emplace_back(track.globalIndex(), fvector.pt(), fvector.eta(), zVert, track.itsClusterMap(), track.tpcNClsFound(), dcaxy, dcaz, flag, nsigmaTPC, massTOF);
        }
      }
    } // end loop over tracks
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    fillDataInfo(collision, tracks);
    for (auto& c : nuclei::candidates) {
      nucleiTable(c.pt, c.eta, c.zVertex, c.ITSclsMap, c.TPCnCls, c.DCAxy, c.DCAz, c.flags, c.TPCnsigma, c.TOFmass);
    }
  }
  PROCESS_SWITCH(nucleiSpectra, processData, "Data analysis", true);

  Preslice<TrackCandidates> tracksPerCollisions = aod::track::collisionId;
  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>> const& collisions, TrackCandidates const& tracks, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, aod::BCsWithTimestamps const&)
  {
    for (auto& collision : collisions) {
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

      nucleiTableMC(c.pt, c.eta, c.zVertex, c.ITSclsMap, c.TPCnCls, c.DCAxy, c.DCAz, c.flags, c.TPCnsigma, c.TOFmass, particle.pt(), particle.eta(), particle.pdgCode());
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (std::abs(particle.pdgCode()) == nuclei::codes[iS]) {
          nuclei::hMomRes[iS][particle.pdgCode() < 0]->Fill(1., std::abs(c.pt), 1. - std::abs(c.pt) / particle.pt());
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
          nucleiTableMC(0, 0, 0, 0, 0, 0, 0, flags, 0, 0, particle.pt(), particle.eta(), particle.pdgCode());
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
