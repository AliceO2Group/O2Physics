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
// o2-analysis-lf-nuclei-spectra, o2-analysis-track-propagation, o2-analysis-timestamp
// o2-analysis-trackselection, o2-analysis-pid-tof-base, o2-analysis-pid-tof-full
// o2-analysis-pid-tpc-full, o2-analysis-multiplicity-table, o2-analysis-event-selection

#include <cmath>

#include "Math/Vector4D.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsTPC/BetheBlochAleph.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/Track.h"

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

float getBinCenter(uint8_t bin, double max) {
  if (bin == 0u) {
    return -max;
  } else if (bin == 255u) {
    return max;
  } else {
    return -max + (bin + 0.5) / (2 * max);
  }
}

namespace nuclei
{
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
constexpr double nSigmaTOFdefault[4][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
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
// constexpr int codes[4]{1000010020, 1000010030, 1000020030, 1000020040};
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

float pidCuts[2][4][2];
std::shared_ptr<TH3> hNsigma[2][4][2];
std::shared_ptr<TH3> hTOFmass[4][2];
std::shared_ptr<TH3> hDCAxy[2][4][2];
} // namespace nuclei

namespace o2::aod
{
namespace NucleiTableNS
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(ITSclsMap, itsClsMap, uint8_t);
DECLARE_SOA_COLUMN(TPCnCls, tpcNCls, uint8_t);
DECLARE_SOA_COLUMN(DCAxy, dcaxy, int8_t);
DECLARE_SOA_COLUMN(DCAz, dcaz, int8_t);
DECLARE_SOA_COLUMN(Flags, flags, uint16_t);
DECLARE_SOA_COLUMN(TPCnsigma, tpcnsigma, uint8_t);
DECLARE_SOA_COLUMN(TOFnsigma, tofnsigma, uint8_t);
DECLARE_SOA_COLUMN(TOFmass, tofmass, uint8_t);
} // namespace NucleiTableNS
DECLARE_SOA_TABLE(NucleiTable, "AOD", "NUCLEITABLE",
                  NucleiTableNS::Pt,
                  NucleiTableNS::Eta,
                  NucleiTableNS::ITSclsMap,
                  NucleiTableNS::TPCnCls,
                  NucleiTableNS::DCAxy,
                  NucleiTableNS::DCAz,
                  NucleiTableNS::Flags,
                  NucleiTableNS::TPCnsigma,
                  NucleiTableNS::TOFnsigma,
                  NucleiTableNS::TOFmass)
} //namespace o2::aod

struct NucleiSpectraTask {
  enum {
    kDeuteron = BIT(0),
    kTriton = BIT(1),
    kHe3 = BIT(2),
    kHe4 = BIT(3),
    kHasTOF = BIT(4)
  };

  Produces<o2::aod::NucleiTable> nucleiTable;

  Configurable<std::string> cfgCentralityEstimator{"cfgCentralityEstimator", "V0A", "Centrality estimator name"};
  Configurable<float> cfgCMrapidity{"cfgCMrapidity", 0.f, "Rapidity of the center of mass (only for p-Pb)"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgCutRapidityMin{"cfgCutRapidityMin", -0.5, "Minimum rapidity for tracks"};
  Configurable<float> cfgCutRapidityMax{"cfgCutRapidityMax", 0.5, "Maximum rapidity for tracks"};
  Configurable<float> cfgCutNclusITS{"cfgCutNclusITS", 4, "Minimum number of ITS clusters"};
  Configurable<float> cfgCutNclusTPC{"cfgCutNclusTPC", 70, "Minimum number of TPC clusters"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], 4, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTPC{"cfgNsigmaTPC", {nuclei::nSigmaTPCdefault[0], 4, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgNsigmaTOF{"cfgNsigmaTOF", {nuclei::nSigmaTOFdefault[0], 4, 2, nuclei::names, nuclei::nSigmaConfigName}, "TOF nsigma selection for light nuclei"};
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
  ConfigurableAxis cfgNsigmaTPCbins{"cfgNsigmaTPCbins", {100, -5., 5.}, "nsigma_TPC binning"};
  ConfigurableAxis cfgNsigmaTOFbins{"cfgNsigmaTOFbins", {100, -5., 5.}, "nsigma_TOF binning"};
  ConfigurableAxis cfgTOFmassBins{"cfgTOFmassBins", {200, -5., 5.}, "TOF mass binning"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (requireGlobalTrackInFilter());

  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>>;
  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(o2::framework::InitContext&)
  {

    const AxisSpec centAxis{cfgCentralityBins, fmt::format("{} percentile", (std::string)cfgCentralityEstimator)};
    const AxisSpec nSigmaAxes[2]{{cfgNsigmaTPCbins, "n#sigma_{TPC}"}, {cfgNsigmaTOFbins, "n#sigma_{TOF}"}};
    const AxisSpec tofMassAxis{cfgTOFmassBins, "TOF mass - PDG mass"};

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
    for (int iC{0}; iC < 2; ++iC) {
      for (int iS{0}; iS < nuclei::species; ++iS) {
        for (int iPID{0}; iPID < 2; ++iPID) {
          nuclei::hNsigma[iPID][iS][iC] = spectra.add<TH3>(fmt::format("h{}nsigma{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], nSigmaAxes[iPID]});
          nuclei::hDCAxy[iPID][iS][iC] = spectra.add<TH3>(fmt::format("hDCAxy{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAxy {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], dcaAxes[iS]});
        }
        nuclei::hTOFmass[iS][iC] = spectra.add<TH3>(fmt::format("h{}TOFmass{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH3D, {centAxis, ptAxes[iS], tofMassAxis});
      }
    }

    std::cout << "CONFIGURATION" << std::endl;
    for (int iS{0}; iS < 4; ++iS) {
      std::cout << cfgTreeConfig->get(iS, 0u) << std::endl;
      for (int iMax{0}; iMax < 2; ++iMax) {
        nuclei::pidCuts[0][iS][iMax] = cfgNsigmaTPC->get(iS, iMax);
        nuclei::pidCuts[1][iS][iMax] = cfgNsigmaTOF->get(iS, iMax);
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, TrackCandidates const& tracks)
  {
    // collision process loop
    if (!collision.sel8()) {
      return;
    }
    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    for (auto& track : tracks) { // start loop over tracks
      if (track.itsNCls() < cfgCutNclusITS ||
          track.tpcNClsFound() < cfgCutNclusTPC ||
          std::abs(track.eta()) > cfgCutEta) {
        continue;
      }
      const int iC{track.sign() < 0};
      spectra.fill(HIST("hTpcSignalData"), track.tpcInnerParam() * track.sign(), track.tpcSignal());
      float nSigma[2][4]{
        {track.tpcNSigmaDe(), track.tpcNSigmaTr(), track.tpcNSigmaHe(), track.tpcNSigmaAl()},
        {track.tofNSigmaDe(), track.tofNSigmaTr(), track.tofNSigmaHe(), track.tofNSigmaAl()}};
      float beta{track.beta()};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        bool selectedTPC{false}, selectedTOF{false};
        if (std::abs(track.dcaZ()) > cfgDCAcut->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{track.pt() * nuclei::charges[iS], track.eta(), track.phi(), nuclei::masses[iS]};
        float y{fvector.Rapidity() + cfgCMrapidity};
        if (y < cfgCutRapidityMin || y > cfgCutRapidityMax) {
          continue;
        }

        if (cfgBetheBlochParams->get(iS, 5u) > 0.f) {
          double expBethe{tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() / nuclei::masses[iS]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};
          double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};
          nSigma[0][iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        }
        for (int iPID{0}; iPID < 2; ++iPID) {
          if (nSigma[0][iS] > nuclei::pidCuts[0][iS][0] && nSigma[0][iS] < nuclei::pidCuts[0][iS][1]) {
            selectedTPC = true;
            if (iPID && (!track.hasTOF() || nSigma[1][iS] < nuclei::pidCuts[1][iS][0] || nSigma[1][iS] > nuclei::pidCuts[1][iS][1])) {
              continue;
            } else if (iPID) {
              selectedTOF = true;
            }
            nuclei::hDCAxy[iPID][iS][iC]->Fill(1., fvector.pt(), track.dcaXY());
            if (std::abs(track.dcaXY()) < cfgDCAcut->get(iS, 0u)) {
              nuclei::hNsigma[iPID][iS][iC]->Fill(1., fvector.pt(), nSigma[iPID][iS]);
              if (iPID) {
                float mass{track.tpcInnerParam() * nuclei::charges[iS] * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS]};
                nuclei::hTOFmass[iS][iC]->Fill(1., fvector.pt(), mass);
              }
            }
          }
        }
        uint16_t flag{0u};
        if (cfgTreeConfig->get(iS, 0u) && selectedTPC) {
          int8_t massTOF{0u};
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          if (track.hasTOF()) {
            flag |= kHasTOF;
            massTOF = getBinnedValue(beta > 1.e-6f ? track.tpcInnerParam() * nuclei::charges[iS] * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS] : -999.f, cfgBinnedVariables->get(4u, 1u));
          }
          flag |= BIT(iS);
          int8_t dcaxy = getBinnedValue(track.dcaXY(), cfgBinnedVariables->get(0u, 1u));
          int8_t dcaz = getBinnedValue(track.dcaZ(), cfgBinnedVariables->get(1u, 1u));
          int8_t nsigmaTPC = getBinnedValue(nSigma[0][iS], cfgBinnedVariables->get(2u, 1u));
          int8_t nsigmaTOF = getBinnedValue(nSigma[1][iS], cfgBinnedVariables->get(3u, 1u));

          nucleiTable(track.sign() * track.pt() * nuclei::charges[iS], track.eta(), track.itsClusterMap(), track.tpcNClsFound(), dcaxy, dcaz, flag, nsigmaTPC, nsigmaTOF, massTOF);
        }
      }
    } // end loop over tracks
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiSpectraTask>(cfgc, TaskName{"nuclei-spectra"})};
}
