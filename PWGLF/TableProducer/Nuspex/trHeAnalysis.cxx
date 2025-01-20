// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file trHeAnalysis.cxx
/// \brief Triton and Helion Analysis on pp Data
/// \author Matthias Herzer <matthias.herzer@cern.ch>, Goethe University
/// Frankfurt

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFNucleiTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <TF1.h>
#include <TLorentzVector.h>

namespace o2::aod
{
namespace h3_data
{
DECLARE_SOA_COLUMN(TPt, tPt, float);
DECLARE_SOA_COLUMN(TEta, tEta, float);
DECLARE_SOA_COLUMN(TPhi, tPhi, float);
DECLARE_SOA_COLUMN(TCharge, tCharge, int8_t);
DECLARE_SOA_COLUMN(TP, tP, float);
DECLARE_SOA_COLUMN(TH3DeDx, tH3DeDx, float);
DECLARE_SOA_COLUMN(TnSigmaTpc, tnSigmaTpc, float);
DECLARE_SOA_COLUMN(TTofSignalH3, tTofSignalH3, float);
DECLARE_SOA_COLUMN(TDcaXY, tDcaXY, float);
DECLARE_SOA_COLUMN(TDcaZ, tDcaZ, float);
DECLARE_SOA_COLUMN(TSigmaYX, tSigmaYX, float);
DECLARE_SOA_COLUMN(TSigmaXYZ, tSigmaXYZ, float);
DECLARE_SOA_COLUMN(TSigmaZ, tSigmaZ, float);
DECLARE_SOA_COLUMN(TnTpcCluster, tnTpcCluster, int);
DECLARE_SOA_COLUMN(TnItsCluster, tnItsCluster, int);
DECLARE_SOA_COLUMN(TTpcChi2NCl, tTpcChi2NCl, float);
DECLARE_SOA_COLUMN(TItsChi2NCl, tItsChi2NCl, float);
DECLARE_SOA_COLUMN(TPassedTpcRefit, tPassedTpcRefit, bool);
DECLARE_SOA_COLUMN(TPassedItsRefit, tPassedItsRefit, bool);
DECLARE_SOA_COLUMN(TRigidity, tRigidity, float);
DECLARE_SOA_COLUMN(TItsClusterSize, tItsClusterSize, float);
} // namespace h3_data
DECLARE_SOA_TABLE(H3Data, "AOD", "h3_data", h3_data::TPt, h3_data::TEta,
                  h3_data::TPhi, h3_data::TCharge, h3_data::TH3DeDx,
                  h3_data::TnSigmaTpc, h3_data::TTofSignalH3, h3_data::TDcaXY,
                  h3_data::TDcaZ, h3_data::TSigmaYX, h3_data::TSigmaXYZ,
                  h3_data::TSigmaZ, h3_data::TnTpcCluster,
                  h3_data::TnItsCluster, h3_data::TTpcChi2NCl,
                  h3_data::TItsChi2NCl, h3_data::TPassedTpcRefit,
                  h3_data::TPassedItsRefit, h3_data::TRigidity,
                  h3_data::TItsClusterSize);
namespace he_data
{
DECLARE_SOA_COLUMN(TPt, tPt, float);
DECLARE_SOA_COLUMN(TEta, tEta, float);
DECLARE_SOA_COLUMN(TPhi, tPhi, float);
DECLARE_SOA_COLUMN(TCharge, tCharge, int8_t);
DECLARE_SOA_COLUMN(TP, tP, float);
DECLARE_SOA_COLUMN(THeDeDx, tHeDeDx, float);
DECLARE_SOA_COLUMN(TnSigmaTpc, tnSigmaTpc, float);
DECLARE_SOA_COLUMN(TTofSignalHe, tTofSignalHe, float);
DECLARE_SOA_COLUMN(TDcaXY, tDcaXY, float);
DECLARE_SOA_COLUMN(TDcaZ, tDcaZ, float);
DECLARE_SOA_COLUMN(TSigmaYX, tSigmaYX, float);
DECLARE_SOA_COLUMN(TSigmaXYZ, tSigmaXYZ, float);
DECLARE_SOA_COLUMN(TSigmaZ, tSigmaZ, float);
DECLARE_SOA_COLUMN(TnTpcCluster, tnTpcCluster, int);
DECLARE_SOA_COLUMN(TnItsCluster, tnItsCluster, int);
DECLARE_SOA_COLUMN(TTpcChi2NCl, tTpcChi2NCl, float);
DECLARE_SOA_COLUMN(TItsChi2NCl, tItsChi2NCl, float);
DECLARE_SOA_COLUMN(TPassedTpcRefit, tPassedTpcRefit, bool);
DECLARE_SOA_COLUMN(TPassedItsRefit, tPassedItsRefit, bool);
DECLARE_SOA_COLUMN(TRigidity, tRigidity, float);
DECLARE_SOA_COLUMN(TItsClusterSize, tItsClusterSize, float);
} // namespace he_data
DECLARE_SOA_TABLE(HeData, "AOD", "he_data", he_data::TPt, he_data::TEta,
                  he_data::TPhi, he_data::TCharge, he_data::THeDeDx,
                  he_data::TnSigmaTpc, he_data::TTofSignalHe, he_data::TDcaXY,
                  he_data::TDcaZ, he_data::TSigmaYX, he_data::TSigmaXYZ,
                  he_data::TSigmaZ, he_data::TnTpcCluster,
                  he_data::TnItsCluster, he_data::TTpcChi2NCl,
                  he_data::TItsChi2NCl, he_data::TPassedTpcRefit,
                  he_data::TPassedItsRefit, he_data::TRigidity,
                  he_data::TItsClusterSize);
} // namespace o2::aod
namespace
{
const int nBetheParams = 6;
const int nParticles = 2;
static const std::vector<std::string> particleNames{"triton", "helion"};
static const std::vector<int> particlePdgCodes{
  o2::constants::physics::kTriton, o2::constants::physics::kHelium3};
static const std::vector<float> particleMasses{
  o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3};
static const std::vector<int> particleCharge{1, 2};
static const std::vector<std::string> betheBlochParNames{
  "p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr float betheBlochDefault[nParticles][nBetheParams]{
  {0.313129, 181.664226, 2779397163087.684082, 2.130773, 29.609643,
   0.09},                                                     // triton
  {70.584685, 3.196364, 0.133878, 2.731736, 1.675617, 0.09}}; // Helion

} // namespace
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using TracksFull =
  soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
            o2::aod::TracksDCA, aod::pidTOFmass, aod::pidTOFbeta,
            aod::pidTPCLfFullTr, aod::pidTPCLfFullHe,
            aod::TOFSignal, aod::TrackSelectionExtension>;

class Particle
{
 public:
  TString name;
  int pdgCode;
  float mass;
  int charge;
  float resolution;
  std::vector<float> betheParams;

  Particle(const std::string name_, int pdgCode_, float mass_, int charge_,
           LabeledArray<float> bethe)
  {
    name = TString(name_);
    pdgCode = pdgCode_;
    mass = mass_;
    charge = charge_;

    resolution =
      bethe.get(name, "resolution"); // Access the "resolution" parameter

    betheParams.clear();
    for (int i = 0; i < 5; ++i) {
      betheParams.push_back(bethe.get(name, i));
    }
  }
};

struct trHeAnalysis {
  Produces<o2::aod::H3Data> H3Data;
  Produces<o2::aod::HeData> HeData;
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  std::vector<Particle> particles;
  Configurable<bool> enableTr{"enableTr", true,
                              "Flag to enable triton analysis."};
  Configurable<bool> enableHe{"enableHe", true,
                              "Flag to enable helium-3 analysis."};

  ConfigurableAxis binsPt{"binsPt",
                          {VARIABLE_WIDTH,
                           0.0,
                           0.05,
                           0.1,
                           0.15,
                           0.2,
                           0.25,
                           0.3,
                           0.35,
                           0.4,
                           0.425,
                           0.45,
                           0.475,
                           0.5,
                           0.5125,
                           0.525,
                           0.5375,
                           0.55,
                           0.5625,
                           0.575,
                           0.5875,
                           0.6,
                           0.6125,
                           0.625,
                           0.6375,
                           0.65,
                           0.6625,
                           0.675,
                           0.6875,
                           0.7,
                           0.7125,
                           0.725,
                           0.7375,
                           0.75,
                           0.7625,
                           0.775,
                           0.7875,
                           0.8,
                           0.8125,
                           0.825,
                           0.8375,
                           0.85,
                           0.8625,
                           0.875,
                           0.8875,
                           0.9,
                           0.9125,
                           0.925,
                           0.9375,
                           0.95,
                           0.9625,
                           0.975,
                           0.9875,
                           1.0,
                           1.0125,
                           1.025,
                           1.0375,
                           1.05,
                           1.0625,
                           1.075,
                           1.0875,
                           1.1,
                           1.1125,
                           1.125,
                           1.1375,
                           1.15,
                           1.1625,
                           1.175,
                           1.1875,
                           1.2,
                           1.2125,
                           1.225,
                           1.2375,
                           1.25,
                           1.2625,
                           1.275,
                           1.2875,
                           1.3,
                           1.3125,
                           1.325,
                           1.3375,
                           1.35,
                           1.3625,
                           1.375,
                           1.3875,
                           1.4,
                           1.4125,
                           1.425,
                           1.4375,
                           1.45,
                           1.4625,
                           1.475,
                           1.4875,
                           1.5,
                           1.5125,
                           1.525,
                           1.5375,
                           1.55,
                           1.5625,
                           1.575,
                           1.5875,
                           1.6,
                           1.6125,
                           1.625,
                           1.6375,
                           1.65,
                           1.6625,
                           1.675,
                           1.6875,
                           1.7,
                           1.7125,
                           1.725,
                           1.7375,
                           1.75,
                           1.7625,
                           1.775,
                           1.7875,
                           1.8,
                           1.8125,
                           1.825,
                           1.8375,
                           1.85,
                           1.8625,
                           1.875,
                           1.8875,
                           1.9,
                           1.9125,
                           1.925,
                           1.9375,
                           1.95,
                           1.9625,
                           1.975,
                           1.9875,
                           2.0,
                           2.0625,
                           2.125,
                           2.1875,
                           2.25,
                           2.3125,
                           2.375,
                           2.4375,
                           2.5,
                           2.625,
                           2.75,
                           2.875,
                           3.0,
                           3.25,
                           3.5,
                           3.75,
                           4.0,
                           4.5,
                           5.0,
                           6.0,
                           7.0,
                           8.0},
                          ""};
  ConfigurableAxis binsPtHe{"binsPtHe",
                            {VARIABLE_WIDTH, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25,
                             2.50, 2.75, 3.0, 3.25, 3.50, 3.75, 4.0, 4.50, 5.0,
                             6.0, 7.0, 8.0},
                            ""};
  ConfigurableAxis binsPtZHe{"binsPtZHe",
                             {VARIABLE_WIDTH, 0.5, 0.625, 0.75, 0.875, 1.0,
                              1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0,
                              2.25, 2.5, 3.0, 3.5, 4.0},
                             ""};
  ConfigurableAxis binsPtTr{"binsPtTr",
                            {VARIABLE_WIDTH, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25,
                             2.50, 2.75, 3.0, 3.25, 3.50, 3.75, 4.0, 4.50, 5.0,
                             6.0, 7.0, 8.0},
                            ""};
  ConfigurableAxis binsPtZTr{"binsPtZTr",
                             {VARIABLE_WIDTH, 0.5, 0.625, 0.75, 0.875, 1.0,
                              1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2.0,
                              2.25, 2.5, 3.0, 3.5, 4.0},
                             ""};
  ConfigurableAxis binsDeDx{"binsDeDx", {600, 0.f, 3000.f}, ""};
  ConfigurableAxis binsBeta{"binsBeta", {120, 0.0, 1.2}, ""};
  ConfigurableAxis binsDca{"binsDca", {400, -1.f, 1.f}, ""};
  ConfigurableAxis binsSigmaTpc{"binsSigmaTpc", {1000, -100, 100}, ""};
  ConfigurableAxis binsSigmaTof{"binsSigmaTof", {1000, -100, 100}, ""};
  ConfigurableAxis binsMassTr{"binsMassTr", {250, -2.5, 2.5f}, ""};
  ConfigurableAxis binsMassHe{"binsMassHe", {300, -3., 3.f}, ""};
  // Set the event selection cuts
  struct : ConfigurableGroup {
    Configurable<bool> useSel8{"useSel8", true,
                               "Use Sel8 for run3 Event Selection"};
    Configurable<bool> tvxTrigger{
      "tvxTrigger", false, "Use TVX for Event Selection (default w/ Sel8)"};
    Configurable<bool> removeTFBorder{
      "removeTFBorder", false, "Remove TimeFrame border (default w/ Sel8)"};
    Configurable<bool> removeITSROFBorder{
      "removeITSROFBorder", false,
      "Remove ITS Read-Out Frame border (default w/ Sel8)"};
  } evselOptions;

  Configurable<bool> cfgTPCPidMethod{
    "cfgTPCPidMethod", false,
    "Using own or built in bethe parametrization"}; // false for built in
                                                    // method
  // Set the multiplity event limits
  Configurable<float> cfgLowMultCut{
    "cfgLowMultCut", 0.0f, "Accepted multiplicity percentage lower limit"};
  Configurable<float> cfgHighMultCut{
    "cfgHighMultCut", 100.0f,
    "Accepted multiplicity percentage higher limit"};

  // Set the z-vertex event cut limits
  Configurable<float> cfgHighCutVertex{"cfgHighCutVertex", 10.0f,
                                       "Accepted z-vertex upper limit"};
  Configurable<float> cfgLowCutVertex{"cfgLowCutVertex", -10.0f,
                                      "Accepted z-vertex lower limit"};

  // Set the quality cuts for tracks
  Configurable<bool> rejectFakeTracks{
    "rejectFakeTracks", false, "Flag to reject ITS-TPC fake tracks (for MC)"};
  Configurable<float> cfgCutITSClusters{"cfgCutITSClusters", -1.f,
                                        "Minimum number of ITS clusters"};
  Configurable<float> cfgCutTPCXRows{"cfgCutTPCXRows", -1.f,
                                     "Minimum number of crossed TPC rows"};
  Configurable<float> cfgCutTPCClusters{"cfgCutTPCClusters", 40.f,
                                        "Minimum number of found TPC clusters"};
  Configurable<int> nITSLayer{"nITSLayer", 0, "ITS Layer (0-6)"};
  Configurable<float> cfgCutTPCcrRowToFindableCl{
    "cfgCutTPCcrRowToFindableCl", 0.8f,
    "Minimum ratio of crossed rows to findable cluster in TPC"};
  Configurable<float> cfgCutmaxChi2TPC{"cfgCutmaxChi2TPC", 4.f,
                                       "Maximum chi2 per cluster for TPC"};
  Configurable<float> cfgCutmaxChi2ITS{"cfgCutmaxChi2ITS", 36.f,
                                       "Maximum chi2 per cluster for ITS"};

  // Set the kinematic and PID cuts for tracks
  struct : ConfigurableGroup {
    Configurable<float> pCut{
      "pCut", 0.3f, "Value of the p selection for spectra (default 0.3)"};
    Configurable<float> etaCut{
      "etaCut", 0.8f, "Value of the eta selection for spectra (default 0.8)"};
    Configurable<float> yLowCut{
      "yLowCut", -1.0f,
      "Value of the low rapidity selection for spectra (default -1.0)"};
    Configurable<float> yHighCut{
      "yHighCut", 1.0f,
      "Value of the high rapidity selection for spectra (default 1.0)"};
  } kinemOptions;

  struct : ConfigurableGroup {
    Configurable<float> nsigmaTPCTr{"nsigmaTPCTr", 5.f,
                                    "Value of the Nsigma TPC cut for tritons"};
    Configurable<float> nsigmaTPCHe{"nsigmaTPCHe", 5.f,
                                    "Value of the Nsigma TPC cut for helium-3"};
  } nsigmaTPCvar;
  Configurable<LabeledArray<float>> cfgBetheBlochParams{
    "cfgBetheBlochParams",
    {betheBlochDefault[0], nParticles, nBetheParams, particleNames,
     betheBlochParNames},
    "TPC Bethe-Bloch parameterisation for light nuclei"};

  void init(o2::framework::InitContext&)
  {
    const AxisSpec pAxis{binsPt, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptHeAxis{binsPtHe, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptZHeAxis{binsPtZHe, "#it{p}_{T}/z (GeV/#it{c})"};
    const AxisSpec ptTrAxis{binsPtTr, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptZTrAxis{binsPtZTr, "#it{p}_{T}/z (GeV/#it{c})"};
    const AxisSpec dedxAxis{binsDeDx, "d#it{E}/d#it{x} A.U."};
    const AxisSpec betaAxis{binsBeta, "TOF #beta"};
    const AxisSpec dcaxyAxis{binsDca, "DCAxy (cm)"};
    const AxisSpec dcazAxis{binsDca, "DCAz (cm)"};
    const AxisSpec massTrAxis{binsMassTr, ""};
    const AxisSpec massHeAxis{binsMassHe, ""};
    const AxisSpec sigmaTPCAxis{binsSigmaTpc, ""};
    const AxisSpec sigmaTOFAxis{binsSigmaTof, ""};

    histos.add<TH1>("histogram/pT",
                    "Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); counts",
                    HistType::kTH1F, {{500, 0., 10.}});
    histos.add<TH1>("histogram/p", "Track momentum; p (GeV/#it{c}); counts",
                    HistType::kTH1F, {{500, 0., 10.}});
    histos.add<TH2>("histogram/TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH2>("histogram/H3/H3-TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/H3/H3-TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH2>("histogram/He/He-TPCsignVsTPCmomentum",
                    "TPC <-dE/dX> vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TPC "
                    "<-dE/dx> (a.u.)",
                    HistType::kTH2F, {{400, -8.f, 8.f}, {dedxAxis}});
    histos.add<TH2>(
      "histogram/He/He-TOFbetaVsP",
      "TOF #beta vs #it{p}/Z; Signed #it{p} (GeV/#it{c}); TOF #beta",
      HistType::kTH2F, {{250, -5.f, 5.f}, {betaAxis}});
    histos.add<TH1>("event/eventSelection", "eventSelection", HistType::kTH1D,
                    {{7, -0.5, 6.5}});
    auto h = histos.get<TH1>(HIST("event/eventSelection"));
    h->GetXaxis()->SetBinLabel(1, "Total");
    h->GetXaxis()->SetBinLabel(2, "TVX trigger cut");
    h->GetXaxis()->SetBinLabel(3, "TF border cut");
    h->GetXaxis()->SetBinLabel(4, "ITS ROF cut");
    h->GetXaxis()->SetBinLabel(5, "TVX + TF + ITS ROF");
    h->GetXaxis()->SetBinLabel(6, "Sel8 cut");
    h->GetXaxis()->SetBinLabel(7, "Z-vert Cut");
    histos.add<TH1>("histogram/cuts", "cuts", HistType::kTH1D,
                    {{9, -0.5, 8.5}});
    auto h_cuts = histos.get<TH1>(HIST("histogram/cuts"));
    h_cuts->GetXaxis()->SetBinLabel(1, "total");
    h_cuts->GetXaxis()->SetBinLabel(2, "p cut");
    h_cuts->GetXaxis()->SetBinLabel(3, "eta cut");
    h_cuts->GetXaxis()->SetBinLabel(4, "TPC cluster");
    h_cuts->GetXaxis()->SetBinLabel(5, "ITS clsuter");
    h_cuts->GetXaxis()->SetBinLabel(6, "TPC crossed rows");
    h_cuts->GetXaxis()->SetBinLabel(7, "max chi2 ITS");
    h_cuts->GetXaxis()->SetBinLabel(8, "max chi2 TPC");
    h_cuts->GetXaxis()->SetBinLabel(9, "crossed rows over findable cluster");
    for (int i = 0; i < nParticles; i++) {
      particles.push_back(Particle(particleNames.at(i), particlePdgCodes.at(i),
                                   particleMasses.at(i), particleCharge.at(i),
                                   cfgBetheBlochParams));
    }
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& event,
               TracksFull const& tracks)
  {
    bool trRapCut = kFALSE;
    bool heRapCut = kFALSE;
    histos.fill(HIST("event/eventSelection"), 0);
    if ((event.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
        (event.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
        (event.selection_bit(aod::evsel::kIsTriggerTVX))) {
      histos.fill(HIST("event/eventSelection"), 4);
    }
    if (evselOptions.useSel8 && !event.sel8())
      return;
    histos.fill(HIST("event/eventSelection"), 5);
    if (event.posZ() < cfgLowCutVertex || event.posZ() > cfgHighCutVertex)
      return;
    histos.fill(HIST("event/eventSelection"), 6);
    if (cfgTPCPidMethod) {
      for (const auto& track : tracks) {
        trRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) <
            kinemOptions.yHighCut;
        heRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) <
            kinemOptions.yHighCut;
        histos.fill(HIST("histogram/cuts"), 0);
        if (std::abs(track.tpcInnerParam()) < kinemOptions.pCut) {
          histos.fill(HIST("histogram/cuts"), 1);
          continue;
        }
        if (std::abs(track.eta()) >= kinemOptions.etaCut) {
          histos.fill(HIST("histogram/cuts"), 2);
          continue;
        }
        if (track.tpcNClsFound() < cfgCutTPCClusters) {
          histos.fill(HIST("histogram/cuts"), 3);
          continue;
        }
        if (track.itsNCls() < cfgCutITSClusters) {
          histos.fill(HIST("histogram/cuts"), 4);
          continue;
        }
        if (track.tpcNClsCrossedRows() < cfgCutTPCXRows) {
          histos.fill(HIST("histogram/cuts"), 5);
          continue;
        }
        if (track.itsChi2NCl() > cfgCutmaxChi2ITS) {
          histos.fill(HIST("histogram/cuts"), 6);
          continue;
        }
        if (track.tpcChi2NCl() > cfgCutmaxChi2TPC) {
          histos.fill(HIST("histogram/cuts"), 7);
          continue;
        }
        if (track.tpcCrossedRowsOverFindableCls() <= cfgCutTPCcrRowToFindableCl) {
          histos.fill(HIST("histogram/cuts"), 8);
          continue;
        }
        histos.fill(HIST("histogram/pT"), track.pt());
        histos.fill(HIST("histogram/p"), track.p());
        histos.fill(HIST("histogram/TPCsignVsTPCmomentum"),
                    track.tpcInnerParam() / (1.f * track.sign()),
                    track.tpcSignal());
        histos.fill(HIST("histogram/TOFbetaVsP"),
                    track.p() / (1.f * track.sign()), track.beta());
        if (enableTr && trRapCut) {
          if (std::abs(getTPCnSigma(track, particles.at(0))) <
              nsigmaTPCvar.nsigmaTPCTr) {
            histos.fill(HIST("histogram/H3/H3-TPCsignVsTPCmomentum"),
                        track.tpcInnerParam() / (1.f * track.sign()),
                        track.tpcSignal());
            histos.fill(HIST("histogram/H3/H3-TOFbetaVsP"),
                        track.p() / (1.f * track.sign()), track.beta());
            float TPt = track.pt();
            float TEta = track.eta();
            float TPhi = track.phi();
            int8_t TCharge = track.sign();
            float TH3DeDx = track.tpcSignal();
            float TnSigmaTpc = track.tpcNSigmaTr();
            float TTofSignalH3 = track.mass();
            float TDcaXY = track.dcaXY();
            float TDcaZ = track.dcaZ();
            float TSigmaYX = track.sigmaY();
            float TSigmaXYZ = track.sigmaSnp();
            float TSigmaZ = track.sigmaZ();
            int TnTpcCluster = track.tpcNClsFound();
            int TnItsCluster = track.itsNCls();
            float TTpcChi2NCl = track.tpcChi2NCl();
            float TItsChi2NCl = track.itsChi2NCl();
            bool TPassedTpcRefit = track.passedTPCRefit();
            bool TPassedItsRefit = track.passedITSRefit();
            float TRigidity = track.tpcInnerParam();
            float TItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            H3Data(TPt, TEta, TPhi, TCharge, TH3DeDx, TnSigmaTpc, TTofSignalH3,
                   TDcaXY, TDcaZ, TSigmaYX, TSigmaXYZ, TSigmaZ, TnTpcCluster,
                   TnItsCluster, TTpcChi2NCl, TItsChi2NCl, TPassedTpcRefit,
                   TPassedItsRefit, TRigidity, TItsClusterSize);
          }
        }
        if (enableHe && heRapCut) {
          if (std::abs(getTPCnSigma(track, particles.at(1))) <
              nsigmaTPCvar.nsigmaTPCHe) {
            histos.fill(HIST("histogram/He/He-TPCsignVsTPCmomentum"),
                        track.tpcInnerParam() / (2.f * track.sign()),
                        track.tpcSignal());
            histos.fill(HIST("histogram/He/He-TOFbetaVsP"),
                        track.p() / (2.f * track.sign()), track.beta());
            float TPt = track.pt();
            float TEta = track.eta();
            float TPhi = track.phi();
            int8_t TCharge = track.sign();
            float THeDeDx = track.tpcSignal();
            float TnSigmaTpc = track.tpcNSigmaHe();
            float TTofSignalHe = track.mass();
            float TDcaXY = track.dcaXY();
            float TDcaZ = track.dcaZ();
            float TSigmaYX = track.sigmaY();
            float TSigmaXYZ = track.sigmaSnp();
            float TSigmaZ = track.sigmaZ();
            int TnTpcCluster = track.tpcNClsFound();
            int TnItsCluster = track.itsNCls();
            float TTpcChi2NCl = track.tpcChi2NCl();
            float TItsChi2NCl = track.itsChi2NCl();
            bool TPassedTpcRefit = track.passedTPCRefit();
            bool TPassedItsRefit = track.passedITSRefit();
            float TRigidity = track.tpcInnerParam();
            float TItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            HeData(TPt, TEta, TPhi, TCharge, THeDeDx, TnSigmaTpc, TTofSignalHe,
                   TDcaXY, TDcaZ, TSigmaYX, TSigmaXYZ, TSigmaZ, TnTpcCluster,
                   TnItsCluster, TTpcChi2NCl, TItsChi2NCl, TPassedTpcRefit,
                   TPassedItsRefit, TRigidity, TItsClusterSize);
          }
        }
      }
    }
    if (!cfgTPCPidMethod) {
      for (const auto& track : tracks) {
        trRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Triton)) <
            kinemOptions.yHighCut;
        heRapCut =
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) >
            kinemOptions.yLowCut &&
          track.rapidity(o2::track::PID::getMass2Z(o2::track::PID::Helium3)) <
            kinemOptions.yHighCut;
        histos.fill(HIST("histogram/cuts"), 0);
        if (std::abs(track.tpcInnerParam()) < kinemOptions.pCut) {
          histos.fill(HIST("histogram/cuts"), 1);
          continue;
        }
        if (std::abs(track.eta()) < kinemOptions.etaCut) {
          histos.fill(HIST("histogram/cuts"), 2);
          continue;
        }
        if (track.tpcNClsFound() < cfgCutTPCClusters) {
          histos.fill(HIST("histogram/cuts"), 3);
          continue;
        }
        if (track.itsNCls() < cfgCutITSClusters) {
          histos.fill(HIST("histogram/cuts"), 4);
          continue;
        }
        if (track.tpcNClsCrossedRows() < cfgCutTPCXRows) {
          histos.fill(HIST("histogram/cuts"), 5);
          continue;
        }
        if (track.itsChi2NCl() > cfgCutmaxChi2ITS) {
          histos.fill(HIST("histogram/cuts"), 6);
          continue;
        }
        if (track.tpcChi2NCl() > cfgCutmaxChi2TPC) {
          histos.fill(HIST("histogram/cuts"), 7);
          continue;
        }
        if (track.tpcCrossedRowsOverFindableCls() <= cfgCutTPCcrRowToFindableCl) {
          histos.fill(HIST("histogram/cuts"), 8);
          continue;
        }
        histos.fill(HIST("histogram/pT"), track.pt());
        histos.fill(HIST("histogram/p"), track.p());
        histos.fill(HIST("histogram/TPCsignVsTPCmomentum"),
                    track.tpcInnerParam() / (1.f * track.sign()),
                    track.tpcSignal());
        histos.fill(HIST("histogram/TOFbetaVsP"),
                    track.p() / (1.f * track.sign()), track.beta());
        if (enableTr && trRapCut) {
          if (std::abs(track.tpcNSigmaTr()) < nsigmaTPCvar.nsigmaTPCTr) {
            histos.fill(HIST("histogram/H3/H3-TPCsignVsTPCmomentum"),
                        track.tpcInnerParam() / (1.f * track.sign()),
                        track.tpcSignal());
            histos.fill(HIST("histogram/H3/H3-TOFbetaVsP"),
                        track.p() / (1.f * track.sign()), track.beta());
            float TPt = track.pt();
            float TEta = track.eta();
            float TPhi = track.phi();
            int8_t TCharge = track.sign();
            float TH3DeDx = track.tpcSignal();
            float TnSigmaTpc = track.tpcNSigmaTr();
            float TTofSignalH3 = track.mass();
            float TDcaXY = track.dcaXY();
            float TDcaZ = track.dcaZ();
            float TSigmaYX = track.sigmaY();
            float TSigmaXYZ = track.sigmaSnp();
            float TSigmaZ = track.sigmaZ();
            int TnTpcCluster = track.tpcNClsFound();
            int TnItsCluster = track.itsNCls();
            float TTpcChi2NCl = track.tpcChi2NCl();
            float TItsChi2NCl = track.itsChi2NCl();
            bool TPassedTpcRefit = track.passedTPCRefit();
            bool TPassedItsRefit = track.passedITSRefit();
            float TRigidity = track.tpcInnerParam();
            float TItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            H3Data(TPt, TEta, TPhi, TCharge, TH3DeDx, TnSigmaTpc, TTofSignalH3,
                   TDcaXY, TDcaZ, TSigmaYX, TSigmaXYZ, TSigmaZ, TnTpcCluster,
                   TnItsCluster, TTpcChi2NCl, TItsChi2NCl, TPassedTpcRefit,
                   TPassedItsRefit, TRigidity, TItsClusterSize);
          }
        }
        if (enableHe && heRapCut) {
          if (std::abs(track.tpcNSigmaHe()) < nsigmaTPCvar.nsigmaTPCHe) {
            histos.fill(HIST("histogram/He/He-TPCsignVsTPCmomentum"),
                        track.tpcInnerParam() / (2.f * track.sign()),
                        track.tpcSignal());
            histos.fill(HIST("histogram/He/He-TOFbetaVsP"),
                        track.p() / (2.f * track.sign()), track.beta());
            float TPt = track.pt();
            float TEta = track.eta();
            float TPhi = track.phi();
            int8_t TCharge = track.sign();
            float THeDeDx = track.tpcSignal();
            float TnSigmaTpc = track.tpcNSigmaHe();
            float TTofSignalHe = track.mass();
            float TDcaXY = track.dcaXY();
            float TDcaZ = track.dcaZ();
            float TSigmaYX = track.sigmaY();
            float TSigmaXYZ = track.sigmaSnp();
            float TSigmaZ = track.sigmaZ();
            int TnTpcCluster = track.tpcNClsFound();
            int TnItsCluster = track.itsNCls();
            float TTpcChi2NCl = track.tpcChi2NCl();
            float TItsChi2NCl = track.itsChi2NCl();
            bool TPassedTpcRefit = track.passedTPCRefit();
            bool TPassedItsRefit = track.passedITSRefit();
            float TRigidity = track.tpcInnerParam();
            float TItsClusterSize =
              getMeanItsClsSize(track) / std::cosh(track.eta());
            HeData(TPt, TEta, TPhi, TCharge, THeDeDx, TnSigmaTpc, TTofSignalHe,
                   TDcaXY, TDcaZ, TSigmaYX, TSigmaXYZ, TSigmaZ, TnTpcCluster,
                   TnItsCluster, TTpcChi2NCl, TItsChi2NCl, TPassedTpcRefit,
                   TPassedItsRefit, TRigidity, TItsClusterSize);
          }
        }
      }
    }
  }

  template <class T>
  float getTPCnSigma(T const& track, Particle const& particle)
  {
    const float rigidity = track.tpcInnerParam();
    if (!track.hasTPC())
      return -999;

    float expBethe{tpc::BetheBlochAleph(
      static_cast<float>(particle.charge * rigidity / particle.mass),
      particle.betheParams[0], particle.betheParams[1],
      particle.betheParams[2], particle.betheParams[3],
      particle.betheParams[4])};
    float expSigma{expBethe * particle.resolution};
    float sigmaTPC =
      static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }

  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    for (int i = 0; i < 8; i++) {
      sum += (track.itsClusterSizes() >> (4 * i) & 15);
      if (track.itsClusterSizes() >> (4 * i) & 15)
        n++;
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<trHeAnalysis>(cfgc),
  };
}
