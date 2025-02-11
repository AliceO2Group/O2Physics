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

/// \file uccZdc.cxx
///
/// \brief task for analysis of UCC with the ZDC
/// \author Omar Vazquez (omar.vazquez.rueda@cern.ch)
/// \since January 29, 2025

#include <cmath>
#include <cstdlib>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h" // required for Filter op.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "TPDGCode.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::evsel;
using namespace o2::constants::physics;
using namespace o2::constants::math;

namespace o2::aod
{
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::FT0MultZeqs,
                            o2::aod::CentFT0Cs, aod::TPCMults>;
using BCsRun3 =
  soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using TracksSel =
  soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA>;
} // namespace o2::aod

using SimCollisions = soa::Join<aod::Collisions, aod::EvSels,
                                aod::McCollisionLabels, o2::aod::CentFT0Cs>;
using SimTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                            aod::TracksDCA, aod::McTrackLabels>;

struct UccZdc {
  // Event selection
  Configurable<float> posZcut{"posZcut", +10.0, "z-vertex position cut"};
  Configurable<float> minT0CcentCut{"minT0CcentCut", 0.0, "Min T0C Cent. cut"};
  Configurable<float> maxT0CcentCut{"maxT0CcentCut", 90.0, "Max T0C Cent. cut"};

  // Track selection settings
  // Configurable<int> minItsNclusters{"minItsNclusters", 5,
  //                                   "minimum number of ITS clusters"};
  // Configurable<int> minTpcNclusters{"minTpcNclusters", 70,
  //                                   "minimum number of TPC clusters"};
  // Configurable<int> minTpcNcrossedRows{
  //     "minTpcNcrossedRows", 70, "minimum number of TPC crossed pad rows"};
  // Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0,
  //                                      "maximum TPC chi^2/Ncls"};
  // Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0,
  //                                      "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.1, "minimum pt of the tracks"};
  Configurable<double> maxPt{"maxPt", 50., "maximum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  // Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  // Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  // Configurable<bool> setDCAselectionPtDep{"setDCAselectionPtDep", true,
  //                                         "require pt dependent selection"};
  // Configurable<double> par0{"par0", 0.0105, "par 0"};
  // Configurable<double> par1{"par1", 0.035, "par 1"};
  // Configurables, binning
  Configurable<int> nBinsAmpFV0{"nBinsAmpFV0", 1000,
                                "Number of bins FV0 amplitude"};
  Configurable<float> maxAmpFV0{"maxAmpFV0", 3000, "Max FV0 amplitude"};
  Configurable<int> nBinsAmpFT0{"nBinsAmpFT0", 1000,
                                "Number of bins FT0 amplitude"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 3000, "Max FT0 amplitude"};
  Configurable<int> nBinsNch{"nBinsNch", 2500, "# of bins for midrapidity Nch"};
  Configurable<float> maxNch{"maxNch", 2500, "Max Nch at midrapidity"};
  Configurable<int> nBinsZDC{"nBinsZDC", 1025, "nBinsZDC"};
  Configurable<int> nBinsZEM{"nBinsZEM", 100, "nBinsZEM"};
  Configurable<float> maxZN{"maxZN", 4099.5, "Max ZN signal"};
  Configurable<float> maxZP{"maxZP", 3099.5, "Max ZP signal"};
  Configurable<float> maxZEM{"maxZEM", 3099.5, "Max ZEM signal"};
  Configurable<int> nBinsTDC{"nBinsTDC", 480, "nbinsTDC"};
  Configurable<float> minTdc{"minTdc", -15.0, "minimum TDC"};
  Configurable<float> maxTdc{"maxTdc", 15.0, "maximum TDC"};
  Configurable<float> minMeanpT{"minMeanpT", 0.5, "minimum [pT]"};
  Configurable<float> maxMeanpT{"maxMeanpT", 1.1, "maximum [pT]"};
  Configurable<int> nBinsMeanpT{"nBinsMeanpT", 160, "# bins [pT]"};
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0., 0.1, 0.25, 0.5, 1., 2., 4., 6., 8., 10., 20.}, "pT binning"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.}, "T0C binning"};

  // Configurable flags ZDC
  Configurable<bool> isAmpZDC{"isAmpZDC", false, "Use amplitude ZDC?"};
  Configurable<bool> isCommPMT{"isCommPMT", false, "Use common PMT ZDC?"};
  Configurable<bool> isSumTowers{"isSumTowers", false, "Use sum of Tow ZDC?"};
  Configurable<bool> isTDCcut{"isTDCcut", false, "Use TDC cut?"};
  Configurable<bool> isZEMcut{"isZEMcut", true, "Use ZEM cut?"};
  Configurable<float> zemCut{"zemCut", 1000.0, "ZEM cut"};
  Configurable<float> tdcCut{"tdcCut", 1.0, "TDC cut"};

  // Configurable<float> tdcZNmincut{"tdcZNmincut", -4.0, "Min ZN TDC cut"};
  // Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", -4.0, "Max ZN TDC cut"};
  // Configurable<float> tdcZPmincut{"tdcZPmincut", -4.0, "Min ZP TDC cut"};
  // Configurable<float> tdcZPmaxcut{"tdcZPmaxcut", -4.0, "Max ZP TDC cut"};

  // Filters
  Filter collFilter = (nabs(aod::collision::posZ) < posZcut);
  Filter trackFilter = (requireGlobalTrackInFilter());

  // Apply Filters
  using TheFilteredCollisions = soa::Filtered<o2::aod::ColEvSels>;
  using TheFilteredCollision = TheFilteredCollisions::iterator;
  using TheFilteredTracks = soa::Filtered<o2::aod::TracksSel>;
  // using TheFilteredTrack = TheFilteredTracks::iterator;

  // Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Histograms: Sim
  HistogramRegistry registrySim{
    "registrySim",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisEvent{10, 0., +10.0, ""};
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisCent{binsCent, "T0C centrality"};

    //  Histograms: paritcle-level info
    registryData.add("EtaVsPhi", ";#eta;#varphi", kTH2F,
                     {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});
    registryData.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    registryData.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});
    registryData.add("dcaXYvspT", ";DCA_{xy};p_{T};",
                     {HistType::kTH2F, {{{50, -1., 1.}, {axisPt}}}});

    registryData.add("hT0C_cent", ";T0C centrality;Entries", kTH1F, {axisCent});
    registryData.add("hEventCounter", "Event counter", kTH1F, {axisEvent});
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(1, "total");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(2, "sel8?");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(3, "zdc?");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(4, "t0?");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(5, "v0a?");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(6, "TDC cut");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(7, "zem cut");
    registryData.get<TH1>(HIST("hEventCounter"))
      ->GetXaxis()
      ->SetBinLabel(8, "min < t0c < max");
    registryData.add("ZNA", "ZNA; ZNA amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZDC, -0.5, maxZN}}});
    registryData.add("ZPA", "ZPA; ZPA amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZDC, -0.5, maxZP}}});
    registryData.add("ZNC", "ZNC; ZNC amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZDC, -0.5, maxZN}}});
    registryData.add("ZPC", "ZPC; ZPC amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZDC, -0.5, maxZP}}});
    registryData.add("ZEM1", "ZEM1; ZEM1 amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZEM, -0.5, maxZEM}}});
    registryData.add("ZEM2", "ZEM2; ZEM2 amplitude; Entries",
                     {HistType::kTH1F, {{nBinsZEM, -0.5, maxZEM}}});
    registryData.add("ZNvsZEM", "ZNvsZEM; ZEM; ZNA+ZNC",
                     {HistType::kTH2F,
                      {{{nBinsZDC, -0.5, maxZEM}, {nBinsZEM, -0.5, maxZN}}}});
    registryData.add("ZNAvsZNC", "ZNAvsZNC; ZNC; ZNA",
                     {HistType::kTH2F,
                      {{{nBinsZDC, -0.5, maxZN}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZPAvsZPC", "ZPAvsZPC; ZPA; ZPC",
                     {HistType::kTH2F,
                      {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZP}}}});
    registryData.add("ZNAvsZPA", "ZNAvsZPA; ZPA; ZNA",
                     {HistType::kTH2F,
                      {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNCvsZPC", "ZNCvsZPC; ZPC; ZNC",
                     {HistType::kTH2F,
                      {{{nBinsZDC, -0.5, maxZP}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNCvstdc", "ZNCvstdc; time ZNC; ZNC",
                     {HistType::kTH2F,
                      {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNAvstdc", "ZNAvstdc; time ZNA; ZNA",
                     {HistType::kTH2F,
                      {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZPCvstdc", "ZPCvstdc; time ZPC; ZPC",
                     {HistType::kTH2F,
                      {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZP}}}});
    registryData.add("ZPAvstdc", "ZPAvstdc; time ZPA; ZPA",
                     {HistType::kTH2F,
                      {{{nBinsTDC, minTdc, maxTdc}, {nBinsZDC, -0.5, maxZP}}}});
    registryData.add(
      "ZEM1vstdc", "ZEM1vstdc; time ZEM1; ZEM1",
      {HistType::kTH2F,
       {{{nBinsTDC, minTdc, maxTdc}, {nBinsZEM, -0.5, maxZEM}}}});
    registryData.add(
      "ZEM2vstdc", "ZEM2vstdc; time ZEM2; ZEM2",
      {HistType::kTH2F,
       {{{nBinsTDC, minTdc, maxTdc}, {nBinsZEM, -0.5, maxZEM}}}});
    registryData.add(
      "debunch", "ZN sum vs. ZN diff.; t_{ZDC}-t_{ZDA}; t_{ZDC}+t_{ZDA}",
      {HistType::kTH2F,
       {{{nBinsTDC, minTdc, maxTdc}, {nBinsTDC, minTdc, maxTdc}}}});
    registryData.add(
      "ZNvsFV0A", "ZNvsFV0A",
      {HistType::kTH2F,
       {{{nBinsAmpFV0, 0., maxAmpFV0}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add(
      "ZNvsFT0", "FT0",
      {HistType::kTH2F,
       {{{nBinsAmpFT0, 0., maxAmpFT0}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNCvsNch", ";Nch (|#eta|<0.8);ZNC",
                     {HistType::kTH2F,
                      {{{nBinsNch, -0.5, maxNch}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNAvsNch", ";Nch (|#eta|<0.8);ZNA",
                     {HistType::kTH2F,
                      {{{nBinsNch, -0.5, maxNch}, {nBinsZDC, -0.5, maxZN}}}});
    registryData.add("ZNCvsNchvspT", ";Nch (|#eta|<0.8);ZNC;[p_{T}]",
                     {HistType::kTH3F,
                      {{{nBinsNch, -0.5, maxNch},
                        {nBinsZDC, -0.5, maxZN},
                        {nBinsMeanpT, minMeanpT, maxMeanpT}}}});
    registryData.add("ZNAvsNchvspT", ";Nch (|#eta|<0.8);ZNA;[p_{T}]",
                     {HistType::kTH3F,
                      {{{nBinsNch, -0.5, maxNch},
                        {nBinsZDC, -0.5, maxZN},
                        {nBinsMeanpT, minMeanpT, maxMeanpT}}}});

    // MC Histograms
    registrySim.add("hEvent_MC_rec", "Event counter", kTH1F, {axisEvent});
    registrySim.add("hT0C_cent_rec", ";T0C centrality;Entries", kTH1F,
                    {axisCent});
    registrySim.add("Pt_MC_rec_ch", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_rec_pi", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_rec_ka", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_rec_pr", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_rec_sigpos", "#Sigma^{+};p_{T};Entries;", kTH1F,
                    {axisPt});
    registrySim.add("Pt_MC_rec_signeg", "#Sigma^{-};p_{T};Entries;", kTH1F,
                    {axisPt});
    registrySim.add("Pt_MC_rec_re", "Remaining ch particles;p_{T};Entries;",
                    kTH1F, {axisPt});
    registrySim.add("EtaVsPhi_MC_rec", ";#eta;#varphi", kTH2F,
                    {{{axisEta}, {100, -0.1 * PI, +2.1 * PI}}});

    registrySim.add("numberOfRecoCollisions", "",
                    {HistType::kTH1F, {{6, -0.5, 5.5}}});
    registrySim.add("hEvent_MC_tru", "Event counter", kTH1F, {axisEvent});
    registrySim.add("hZpos_MC_tru", "z_{vtx}",
                    {HistType::kTH1F, {{48, -12., 12}}});
    registrySim.add("hZpos_MC_rec", "z_{vtx}",
                    {HistType::kTH1F, {{48, -12., 12}}});
    registrySim.add(
      "aV0Avsb", ";V0A amplitude; Impact parameter",
      {HistType::kTH2F, {{{nBinsAmpFV0, 0., maxAmpFV0}, {19, 0., 18.}}}});
    registrySim.add(
      "aT0Avsb", ";T0A amplitude; Impact parameter",
      {HistType::kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {19, 0., 18.}}}});
    registrySim.add(
      "aT0Cvsb", ";T0C amplitude; Impact parameter",
      {HistType::kTH2F, {{{nBinsAmpFT0, 0., maxAmpFT0}, {19, -0.5, 18.5}}}});
    registrySim.add("Pt_MC_tru_ch", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_tru_pi", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_tru_ka", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_tru_pr", ";p_{T};Entries;", kTH1F, {axisPt});
    registrySim.add("Pt_MC_tru_sigpos", "#Sigma^{+};p_{T};Entries;", kTH1F,
                    {axisPt});
    registrySim.add("Pt_MC_tru_signeg", "#Sigma^{-};p_{T};Entries;", kTH1F,
                    {axisPt});
    registrySim.add("Pt_MC_tru_re", "Remaining ch particles;p_{T};Entries;",
                    kTH1F, {axisPt});
  }

  void processZdcCollAss(TheFilteredCollision const& collision,
                         o2::aod::BCsRun3 const& /*bcs*/,
                         aod::Zdcs const& /*zdcs*/, aod::FV0As const& /*fv0as*/,
                         aod::FT0s const& /*ft0s*/,
                         TheFilteredTracks const& tracks)
  {
    registryData.fill(HIST("hEventCounter"), 0.5);
    if (!collision.sel8()) {
      return;
    }
    registryData.fill(HIST("hEventCounter"), 1.5);

    const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
    // Does the event has ZDC?
    if (!foundBC.has_zdc()) {
      return;
    }
    registryData.fill(HIST("hEventCounter"), 2.5);

    float aT0A{0.0};
    float aT0C{0.0};
    float aV0A{0.0};
    float sumT0s{0.0};
    float znA{0.0};
    float znC{0.0};
    float zpA{0.0};
    float zpC{0.0};
    float aZEM1{0.0};
    float aZEM2{0.0};
    float tZEM1{0.0};
    float tZEM2{0.0};
    float tZNA{0.0};
    float tZNC{0.0};
    float tZPA{0.0};
    float tZPC{0.0};
    float sumZNs{0.0};
    float sumZEMs{0.0};
    float tZDCdif{0.0};
    float tZDCsum{0.0};

    aZEM1 = foundBC.zdc().amplitudeZEM1();
    aZEM2 = foundBC.zdc().amplitudeZEM2();
    tZEM1 = foundBC.zdc().timeZEM1();
    tZEM2 = foundBC.zdc().timeZEM2();
    tZNA = foundBC.zdc().timeZNA();
    tZNC = foundBC.zdc().timeZNC();
    tZPA = foundBC.zdc().timeZPA();
    tZPC = foundBC.zdc().timeZPC();
    tZDCdif = tZNC + tZPC - tZNA - tZPA;
    tZDCsum = tZNC + tZPC + tZNA + tZPA;

    if (isAmpZDC) {
      znA = foundBC.zdc().amplitudeZNA();
      znC = foundBC.zdc().amplitudeZNC();
      zpA = foundBC.zdc().amplitudeZPA();
      zpC = foundBC.zdc().amplitudeZPC();
    } else if (isCommPMT) {
      znA = foundBC.zdc().energyCommonZNA();
      znC = foundBC.zdc().energyCommonZNC();
      zpA = foundBC.zdc().energyCommonZPA();
      zpC = foundBC.zdc().energyCommonZPC();
    } else if (isSumTowers) {
      for (const auto& eZNA : foundBC.zdc().energySectorZNA())
        znA += eZNA;
      for (const auto& eZNC : foundBC.zdc().energySectorZNC())
        znC += eZNC;
      for (const auto& eZPA : foundBC.zdc().energySectorZPA())
        zpA += eZPA;
      for (const auto& eZPC : foundBC.zdc().energySectorZPC())
        zpC += eZPC;
    } else {
      znA = -999.;
      znC = -999.;
      zpA = -999.;
      zpC = -999.;
    }
    sumZNs = znA + znC;
    sumZEMs = aZEM1 + aZEM2;

    if (foundBC.has_ft0()) {
      for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
        aT0A += amplitude;
      }
      for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
        aT0C += amplitude;
      }
      sumT0s = aT0A + aT0C;
      registryData.fill(HIST("hEventCounter"), 3.5);
    } else {
      aT0A = aT0C = -999.;
      sumT0s = -999.;
    }

    if (foundBC.has_fv0a()) {
      for (const auto& amplitude : foundBC.fv0a().amplitude()) {
        aV0A += amplitude;
      }
      registryData.fill(HIST("hEventCounter"), 4.5);
    } else {
      aV0A = -999.;
    }

    // TDC cut
    if (isTDCcut) {
      if (std::sqrt(std::pow(tZDCdif, 2.) + std::pow(tZDCsum, 2.)) > tdcCut) {
        return;
      }
    }
    registryData.fill(HIST("hEventCounter"), 5.5);

    // ZEM cut
    if (isZEMcut) {
      if (sumZEMs < zemCut) {
        return;
      }
    }
    registryData.fill(HIST("hEventCounter"), 6.5);

    // T0C centrality cut
    if (collision.centFT0C() < minT0CcentCut ||
        collision.centFT0C() > maxT0CcentCut) {
      return;
    }

    registryData.fill(HIST("hEventCounter"), 7.5);
    registryData.fill(HIST("hT0C_cent"), collision.centFT0C());
    registryData.fill(HIST("ZNA"), znA);
    registryData.fill(HIST("ZNC"), znC);
    registryData.fill(HIST("ZPA"), zpA);
    registryData.fill(HIST("ZPC"), zpC);
    registryData.fill(HIST("ZNAvsZNC"), znC, znA);
    registryData.fill(HIST("ZNAvsZPA"), zpA, znA);
    registryData.fill(HIST("ZNCvsZPC"), zpC, znC);
    registryData.fill(HIST("ZPAvsZPC"), zpC, zpA);
    registryData.fill(HIST("ZNvsZEM"), sumZEMs, sumZNs);
    registryData.fill(HIST("ZNvsFV0A"), aV0A / 100., sumZNs);
    registryData.fill(HIST("ZNvsFT0"), sumT0s / 100., sumZNs);
    registryData.fill(HIST("ZEM1"), aZEM1);
    registryData.fill(HIST("ZEM2"), aZEM2);
    registryData.fill(HIST("ZNCvstdc"), tZNC, znC);
    registryData.fill(HIST("ZNAvstdc"), tZNA, znA);
    registryData.fill(HIST("ZPCvstdc"), tZPC, zpC);
    registryData.fill(HIST("ZPAvstdc"), tZPA, zpA);
    registryData.fill(HIST("ZEM1vstdc"), tZEM1, aZEM1);
    registryData.fill(HIST("ZEM2vstdc"), tZEM2, aZEM2);
    registryData.fill(HIST("debunch"), tZDCdif, tZDCsum);

    float meanpT{0.0};
    const int64_t nch{tracks.size()};
    for (const auto& track : tracks) {
      // Track Selection
      // if (!track.isGlobalTrack()) {
      //   continue;
      // }
      // if (track.pt() < minPt || track.pt() > maxPt) {
      //   continue;
      // }
      // if (!passedTrackSelection(track)) {
      //   continue;
      // }
      meanpT += track.pt();
      registryData.fill(HIST("EtaVsPhi"), track.eta(), track.phi());
      registryData.fill(HIST("etaHistogram"), track.eta());
      registryData.fill(HIST("ptHistogram"), track.pt());
      registryData.fill(HIST("dcaXYvspT"), track.dcaXY(), track.pt());
    }

    if (nch > 0) {
      meanpT /= nch;
    }
    registryData.fill(HIST("ZNAvsNch"), nch, znA);
    registryData.fill(HIST("ZNCvsNch"), nch, znC);
    registryData.fill(HIST("ZNCvsNchvspT"), nch, znC, meanpT);
    registryData.fill(HIST("ZNAvsNchvspT"), nch, znA, meanpT);
  }
  PROCESS_SWITCH(UccZdc, processZdcCollAss,
                 "Processing ZDC w. collision association", true);

  Preslice<aod::McParticles> perMCCollision = aod::mcparticle::mcCollisionId;
  Preslice<SimTracks> perCollision = aod::track::collisionId;
  void processMC(aod::McCollisions const& mcCollisions,
                 o2::aod::BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/,
                 aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/,
                 SimCollisions const& collisions,
                 aod::McParticles const& mcParticles,
                 SimTracks const& simTracks)
  {
    // Generated MC
    for (const auto& mccollision : mcCollisions) {
      registrySim.fill(HIST("hEvent_MC_tru"), 0.5);
      // Z-vtx position cut
      if (std::fabs(mccollision.posZ()) > posZcut) {
        continue;
      }
      registrySim.fill(HIST("hZpos_MC_tru"), mccollision.posZ());
      registrySim.fill(HIST("hEvent_MC_tru"), 1.5);

      auto mcParticlesPerColl =
        mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());

      for (const auto& particle : mcParticlesPerColl) {
        if (particle.eta() < minEta || particle.eta() > maxEta) {
          continue;
        }
        if (particle.pt() < minPt || particle.pt() > maxPt) {
          continue;
        }
        registrySim.fill(HIST("Pt_MC_tru_ch"), particle.pt());
        if (particle.pdgCode() == PDG_t::kPiPlus ||
            particle.pdgCode() == PDG_t::kPiMinus) { // pion
          registrySim.fill(HIST("Pt_MC_tru_pi"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kKPlus ||
                   particle.pdgCode() == PDG_t::kKMinus) { // kaon
          registrySim.fill(HIST("Pt_MC_tru_ka"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kProton ||
                   particle.pdgCode() == PDG_t::kProtonBar) { // proton
          registrySim.fill(HIST("Pt_MC_tru_pr"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaPlus ||
                   particle.pdgCode() ==
                     PDG_t::kSigmaBarMinus) { // positive sigma
          registrySim.fill(HIST("Pt_MC_tru_sigpos"), particle.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaMinus ||
                   particle.pdgCode() ==
                     PDG_t::kSigmaBarPlus) { // negative sigma
          registrySim.fill(HIST("Pt_MC_tru_signeg"), particle.pt());
        } else { // rest
          registrySim.fill(HIST("Pt_MC_tru_re"), particle.pt());
        }
      }
    }
    registrySim.fill(HIST("numberOfRecoCollisions"), collisions.size());
    // if (collisions.size() == 0 || collisions.size() > 1) {
    //   return;
    // }
    //----- MC reconstructed -----//
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mccollision = collision.mcCollision_as<aod::McCollisions>();
      registrySim.fill(HIST("hEvent_MC_rec"), 0.5);

      // Event Selection
      if (!collision.sel8()) {
        continue;
      }
      registrySim.fill(HIST("hEvent_MC_rec"), 1.5);

      // Z-vertex position cut
      if (std::fabs(collision.posZ()) > posZcut) {
        continue;
      }

      // T0C centrality cut
      if (collision.centFT0C() < minT0CcentCut ||
          collision.centFT0C() > maxT0CcentCut) {
        continue;
      }

      registrySim.fill(HIST("hEvent_MC_rec"), 2.5);
      registrySim.fill(HIST("hZpos_MC_rec"), collision.posZ());

      const auto& foundBC = collision.foundBC_as<o2::aod::BCsRun3>();
      if (foundBC.has_zdc()) {
        return;
      }
      float aT0A{0.0};
      float aT0C{0.0};
      float aV0A{0.0};
      float b{mccollision.impactParameter()};
      if (foundBC.has_ft0()) {
        for (const auto& amplitude : foundBC.ft0().amplitudeA()) {
          aT0A += amplitude;
        }
        for (const auto& amplitude : foundBC.ft0().amplitudeC()) {
          aT0C += amplitude;
        }
      } else {
        aT0A = aT0C = -999;
      }
      if (foundBC.has_fv0a()) {
        for (const auto& amplitude : foundBC.fv0a().amplitude()) {
          aV0A += amplitude;
        }
      } else {
        aV0A = -999;
      }

      registrySim.fill(HIST("hT0C_cent_rec"), collision.centFT0C());
      registrySim.fill(HIST("aT0Avsb"), aT0A / 100., b);
      registrySim.fill(HIST("aT0Cvsb"), aT0C / 100., b);
      registrySim.fill(HIST("aV0Avsb"), aV0A / 100., b);

      auto groupedTracks =
        simTracks.sliceBy(perCollision, collision.globalIndex());

      for (const auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        // Track Selection
        if (!track.isGlobalTrack()) {
          continue;
        }
        if (track.pt() < minPt || track.pt() > maxPt) {
          continue;
        }
        // if (!passedTrackSelection(track)) {
        //   continue;
        // }

        registrySim.fill(HIST("EtaVsPhi_MC_rec"), track.eta(), track.phi());

        const auto particle = track.mcParticle();
        registrySim.fill(HIST("Pt_MC_rec_ch"), track.pt());
        if (particle.pdgCode() == PDG_t::kPiPlus ||
            particle.pdgCode() == PDG_t::kPiMinus) {
          registrySim.fill(HIST("Pt_MC_rec_pi"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kKPlus ||
                   particle.pdgCode() == PDG_t::kKMinus) {
          registrySim.fill(HIST("Pt_MC_rec_ka"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kProton ||
                   particle.pdgCode() == PDG_t::kProtonBar) {
          registrySim.fill(HIST("Pt_MC_rec_pr"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaPlus ||
                   particle.pdgCode() == PDG_t::kSigmaBarMinus) {
          registrySim.fill(HIST("Pt_MC_rec_sigpos"), track.pt());
        } else if (particle.pdgCode() == PDG_t::kSigmaMinus ||
                   particle.pdgCode() == PDG_t::kSigmaBarPlus) {
          registrySim.fill(HIST("Pt_MC_rec_signeg"), track.pt());
        } else {
          registrySim.fill(HIST("Pt_MC_rec_re"), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(UccZdc, processMC, "process pure simulation", false);

  // Single-Track Selection
  // template <typename T2>
  // bool passedTrackSelection(const T2& track) {
  //   if (track.eta() < minEta || track.eta() > maxEta) return false;
  //   if (track.pt() < minPt) return false;
  //
  //   if (!track.hasITS()) return false;
  //   if (track.itsNCls() < minItsNclusters) return false;
  //   if (!track.hasTPC()) return false;
  //   if (track.tpcNClsFound() < minTpcNclusters) return false;
  //   if (track.tpcNClsCrossedRows() < minTpcNcrossedRows) return false;
  //   if (track.tpcChi2NCl() > maxChiSquareTpc) return false;
  //   if (track.itsChi2NCl() > maxChiSquareIts) return false;
  //   // pt-dependent selection
  //   if (setDCAselectionPtDep) {
  //     if (std::fabs(track.dcaXY()) > (par0 + par1 / track.pt())) return
  //     false; if (std::fabs(track.dcaZ()) > (par0 + par1 / track.pt())) return
  //     false;
  //   }
  //   // standard selection
  //   if (!setDCAselectionPtDep) {
  //     if (std::fabs(track.dcaXY()) > maxDcaxy) return false;
  //     if (std::fabs(track.dcaZ()) > maxDcaz) return false;
  //   }
  //   return true;
  // }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UccZdc>(cfgc)};
}
