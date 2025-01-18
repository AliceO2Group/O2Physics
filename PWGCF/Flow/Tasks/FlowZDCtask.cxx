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

/// \file   FlowZDCtask.cxx
/// \author Sabrina Hernandez
/// \since  10/01/2024
/// \brief  task to evaluate flow and neutron skin with information from ZDC
#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <complex>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TComplex.h>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::mult;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>;
using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using AodZDCs = soa::Join<aod::ZDCMults, aod::Zdcs>;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowZDCtask {
  SliceCache cache;

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.9f, "Eta range for tracks") // changed from .8 bc zdc is at eta 8.78
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};
  Configurable<int> eventSelection{"eventSelection", 1, "event selection"};
  Configurable<float> maxZp{"maxZp", 3099.5, "Max ZP signal"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> etaGap{"etaGap", 0.5, "Eta gap"};
  Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
  Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};
  Configurable<float> maxZem{"maxZem", 3099.5, "Max ZEM signal"};
  // for ZDC info and analysis
  Configurable<int> nBinsADC{"nBinsADC", 1000, "nbinsADC"};
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<int> nBinsFT0Amp{"nBinsFT0Amp", 250000, "nbinsAmp"};
  Configurable<float> maxZn{"maxZn", 4099.5, "Max ZN signal"};
  Configurable<float> acceptanceZna{"acceptanceZna", 0.92, "ZNA acceptance factor"};
  Configurable<float> acceptanceZnc{"acceptanceZnc", 0.90, "ZNC acceptance factor"};
  Configurable<float> acceptanceZpa{"acceptanceZpa", 0.52, "ZPA acceptance factor"};
  Configurable<float> acceptanceZpc{"acceptanceZpc", 0.50, "ZPC acceptance factor"};

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {3500, 0, 3500}, "centrality axis for histograms"};
  ConfigurableAxis axisEnergy{"axisEnergy", {100, 0, 700}, "energy axis for zdc histos"};
  ConfigurableAxis axisMultTpc{"axisMultTpc", {1000, -0.5f, 1999.5f}, "TPCmultiplicity"};
  ConfigurableAxis axisZN{"axisZN", {5000, 0, 500}, "axisZN"};
  ConfigurableAxis axisZP{"axisZP", {5000, 0, 500}, "axisZP"};
  ConfigurableAxis axisFT0CAmp{"axisFT0CAmp", {5000, 0, 5000}, "axisFT0CAmp"};
  ConfigurableAxis axisFT0AAmp{"axisFT0AAmp", {5000, 0, 5000}, "axisFT0AAmp"};
  ConfigurableAxis axisFT0MAmp{"axisFT0MAmp", {10000, 0, 10000}, "axisFT0MAmp"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  Partition<AodTracks> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  std::complex<double> qTPC;       // init q TPC
  std::complex<double> qZNA{0, 0}; // init qZNA
  std::complex<double> qZNC{0, 0}; // init qZNC

  // Begin Histogram Registry

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TProfile> q2RealMean{TProfile("q2_real_mean", "q2 real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> q2ImagMean{TProfile("q2_imag_mean", "q2 imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> q2After{TProfile("q2after", "q2 recentered vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> q2Before{TProfile("q2before", "q2 re vs imag", 10, 0, 100.)};
  OutputObj<TProfile> q2ZnaReal{TProfile("Q2_ZNA_real", "q2_ZNA real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> q2ZnaImag{TProfile("Q2_ZNA_imag", "q2_ZNA imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> q2ZncReal{TProfile("Q2_ZNC_real", "q2_ZNC real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> qZncImag{TProfile("Q2_ZNC_imag", "q2_ZNC imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> avgQ2TPCRe{TProfile("avgQ2TPCRe", "Average Q2 Real part vs Centrality", 10, 0, 100)};
  OutputObj<TProfile> avgQ2TPCIm{TProfile("avgQ2TPCIm", "Average Q2 Imaginary part vs Centrality", 10, 0, 100)};
  OutputObj<TProfile> zdcZemEnergy{TProfile("ZDC_ZEM_Energy", "ZDC vs ZEM Energy", 10, 0, 1000)};
  OutputObj<TProfile> pCosPsiDifferences{TProfile("pCosPsiDifferences", "Differences in cos(psi) vs Centrality;Centrality;Mean cos(psi) Difference", 200, 0, 100, -1, 1)};
  OutputObj<TProfile> pSinPsiDifferences{TProfile("pSinPsiDifferences", "Differences in sin(psi) vs Centrality;Centrality;Mean sin(psi) Difference", 200, 0, 100, -1, 1)};
  OutputObj<TProfile> pZNvsFT0Ccent{TProfile("pZNvsFT0Ccent", "ZN Energy vs FT0C Centrality", 100, 0, 100, 0, 500)};
  OutputObj<TProfile> pZPvsFT0Ccent{TProfile("pZPvsFT0Ccent", "ZP Energy vs FT0C Centrality", 100, 0, 100, 0, 500)};
  OutputObj<TProfile> pZNratiovscent{TProfile("pZNratiovscent", "Ratio ZNC/ZNA vs FT0C Centrality", 100, 0, 100, 0, 5)};
  OutputObj<TProfile> pZPratiovscent{TProfile("pZPratiovscent", "Ratio ZPC/ZPA vs FT0C Centrality", 100, 0, 100, 0, 5)};

  double sumCosPsiDiff = 0.0; // Sum of cos(psiZNC) - cos(psiZNA)
  int countEvents = 0;        // Count of processed events

  void init(InitContext const&)
  {
    // define axes
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispt{100, 0, 2, "#pt"};

    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisPhi{100, 0, o2::constants::math::TwoPI, "#phi"};
    const AxisSpec axisQ{100, -1, 1, "Q"};
    const AxisSpec axisZNA{100, 0, 200, "energy"};
    const AxisSpec axisQZNA{100, -1, 1, "Q"};
    const AxisSpec axisREQ{100, -1, 1, "real Q"};
    const AxisSpec axisIMQ{100, -1, 1, "imag Q"};

    AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
    AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
    AxisSpec axisCent{8, 0.f, 105.f, "centrality"};
    AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};
    AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});

    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("centHistogram", "centHistogram", kTH1F, {axisCent});
    histos.add("multHistogram", "multHistogram", kTH1F, {axisMultiplicity});
    histos.add("multvsCent", "centrality vs multiplicity", kTH2F, {axisCent, axisMultiplicity});
    histos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    histos.add("TPCmultiplicity", "TPCmultiplicity", kTH1F, {axisMultTpc});

    histos.add("REqHistogram", "REqHistogram", kTH1F, {axisQ});
    histos.add("IMqHistogram", "IMqHistogram", kTH1F, {axisQ});

    histos.add("REqHistogramZNA", "REqHistogramZNA", kTH1F, {axisQZNA});
    histos.add("IMqHistogramZNA", "IMqHistogramZNA", kTH1F, {axisQZNA});

    histos.add("REqHistogramZNC", "REqHistogramZNC", kTH1F, {axisQZNA});
    histos.add("IMqHistogramZNC", "IMqHistogramZNC", kTH1F, {axisQZNA});

    histos.add("EnergyZNA", "ZNA Sector Energy", kTH1F, {axisEnergy});
    histos.add("EnergyZNC", "ZNC Sector Energy", kTH1F, {axisEnergy});
    histos.add("hCentFT0C", "FT0C Centrality Distribution", kTH1F, {{100, 0, 105}});
    histos.add("hZNvsFT0Ccent",
               "ZN Energy vs FT0C Centrality",
               kTH2F,
               {AxisSpec{100, 0, 100, "Centrality [%]"}, AxisSpec{100, 0, 500, "ZN Energy"}});

    histos.add("hZPvsFT0Ccent",
               "ZP Energy vs FT0C Centrality;Centrality [%];ZP Energy",
               kTH2F,
               {AxisSpec{100, 0, 100, "Centrality [%]"}, AxisSpec{100, 0, 500, "ZP Energy"}});
    // for q vector recentering
    histos.add("revsimag", "revsimag", kTH2F, {axisREQ, axisIMQ});

    if (doprocessZdcCollAssoc) { // Check if the process function for ZDCCollAssoc is enabled
      histos.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZn}}});
      histos.add("ZEM1coll", "ZEM1coll; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZem}}});
      histos.add("ZEM2coll", "ZEM2coll; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, maxZem}}});
      histos.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZem}, {nBinsAmp, -0.5, 2. * maxZn}}}});
      histos.add("ZNAvsZNCcoll", "ZNAvsZNCcoll; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZn}, {nBinsAmp, -0.5, maxZn}}}});

      histos.add("RealQHistogramZNA", "RealQHistogramZNA", kTH1F, {axisQZNA});
      histos.add("ImagQHistogramZNA", "ImagQHistogramZNA", kTH1F, {axisQZNA});
      histos.add("RealQHistogramZNC", "RealQHistogramZNC", kTH1F, {axisQZNA});
      histos.add("ImagQHistogramZNC", "ImagQHistogramZNC", kTH1F, {axisQZNA});

      histos.add("Acorrelations", "Acorrelations", kTH2F, {{axisQZNA}, {axisQZNA}});
      histos.add("SPAngleZNA", "Spectator Plane Angle ZNA;Angle (radians);Entries", {HistType::kTH1F, {{100, -o2::constants::math::PI, o2::constants::math::PI}}});
      histos.add("SPAngleZNC", "Spectator Plane Angle ZNC;Angle (radians);Entries", {HistType::kTH1F, {{100, -o2::constants::math::PI, o2::constants::math::PI}}});

      histos.add("RunningAverageCosPsiDiff", "Running Average of cos(psi) Differences;Running Average;Entries", {HistType::kTH1F, {{100, -1, 1}}});

      histos.add("CosPsiDifferences", "Differences in cos(psi);cos(psiZNC) - cos(psiZNA);Entries", {HistType::kTH1F, {{100, -2, 2}}});
      histos.add("hSinDifferences", "Differences in sin(psi);sin(psiZNC) - sin(psiZNA);Entries", {HistType::kTH1F, {{100, -2, 2}}});
      histos.add("CosPsiDifferencesAvg", "Differences in cos(psi);cos(psiZNC) - cos(psiZNA);Entries", {HistType::kTH2F, {{axisCent}, {100, -2, 2}}});
      histos.add("ZDC_energy_vs_ZEM", "ZDCvsZEM; ZEM; ZNA+ZNC+ZPA+ZPC", {HistType::kTH2F, {{{nBinsAmp, -0.5, maxZem}, {nBinsAmp, -0.5, 2. * maxZn}}}});
      // common energies information for ZDC
      histos.add("ZNCenergy", "ZN energy side c", kTH1F, {axisEnergy});
      histos.add("ZNAenergy", "ZN energy side a", kTH1F, {axisEnergy});
      histos.add("ZPCenergy", "ZP energy side c", kTH1F, {axisEnergy});
      histos.add("ZPAenergy", "ZP energy side a", kTH1F, {axisEnergy});
      histos.add("ZNenergy", "common zn (a + c sides) energy", kTH1F, {axisEnergy});
      histos.add("ZPenergy", "common zp energy (a + c sides)", kTH1F, {axisEnergy});
      histos.add("hFT0CAmp", ";Amplitude;counts", kTH1F, {axisFT0CAmp});
      histos.add("hFT0AAmp", ";Amplitude;counts", kTH1F, {axisFT0AAmp});
      histos.add("hFT0MAmp", ";Amplitude;counts", kTH1F, {axisFT0MAmp});
      histos.add("hMultT0A", ";Amplitude;counts", kTH1F, {{nBinsFT0Amp, 0, 250000}});
      histos.add("hMultT0C", ";Amplitude;counts", kTH1F, {{nBinsFT0Amp, 0, 250000}});
      histos.add("hMultT0M", ";Amplitude;counts", kTH1F, {{nBinsFT0Amp, 0, 250000}});
      histos.add("hZNvsFT0CAmp", "ZN Energy vs FT0C Amplitude", kTH2F, {axisFT0CAmp, axisZN});
      histos.add("hZPvsFT0CAmp", "ZP Energy vs FT0C Amplitude", kTH2F, {axisFT0CAmp, axisZP});
      histos.add("hZNvsMult", "ZN Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZN});
      histos.add("hZPvsMult", "ZP Energy vs Multiplicity", kTH2F, {axisMultiplicity, axisZP});
    }
  }

  void processQVector(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::ZDCMults const& /*zdcMults*/)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    if (!collision.sel8())
      return;
    histos.fill(HIST("centHistogram"), collision.centFT0C());
    const auto& tracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    const int multTPC = tracksGrouped.size();
    const auto cent = collision.centFT0C();
    int nTot = tracks.size();

    // this is the q vector for the TPC data. it is a complex function
    double qTpcReal = 0.0; // Initialize qTPC_real
    double qTpcIm = 0.0;   // init qTPC_imaginary

    int multTrk = tracks.size(); // Tracks are already filtered with GlobalTrack || GlobalTrackSDD

    if (cent < 0.0 && cent > 70)
      return;
    std::complex<double> qTPC(0, 0); // Starting with a q-vector of zero

    for (const auto& track : tracks) {
      double phi = track.phi();

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("phiHistogram"), track.phi());
      histos.fill(HIST("ptHistogram"), track.pt());

      qTPC += std::complex<double>(std::cos(2.0 * phi), std::sin(2.0 * phi));

      histos.fill(HIST("multvsCent"), cent, nTot);

    } // end track loop

    qTpcReal = qTPC.real() / nTot; // normalize these vectors by the total number of particles
    qTpcIm = qTPC.imag() / nTot;

    histos.fill(HIST("REqHistogram"), qTpcReal);
    histos.fill(HIST("IMqHistogram"), qTpcIm);

    histos.fill(HIST("multHistogram"), multTrk);
    histos.fill(HIST("TPCmultiplicity"), multTPC);

    histos.fill(HIST("revsimag"), qTpcReal, qTpcIm);
  }
  void processZdcCollAssoc(
    AodCollisions::iterator const& collision,
    AodTracks const& tracks,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/,
    aod::FT0s const& /*ft0s*/)
  {
    int nTot = tracks.size();
    double sumCosPsiDiff = 0.0; // initialize Sum of cosPsiDiff for averaging
    double sumSinPsiDiff = 0.0; // initialize Sum of cosPsiDiff for averaging
    int countEvents = 0;        // initialize Counter for the number of events processed
    double ft0aAmp = 0;
    double ft0cAmp = 0;
    // collision-based event selection
    if (!collision.sel8())
      return;
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (const auto& amplitude : ft0.amplitudeA()) {
        histos.fill(HIST("hFT0AAmp"), amplitude);
        ft0aAmp += amplitude;
      }
      for (const auto& amplitude : ft0.amplitudeC()) {
        histos.fill(HIST("hFT0CAmp"), amplitude);
        ft0cAmp += amplitude;
      }
    }
    double ft0mAmp = ft0aAmp + ft0cAmp;
    histos.fill(HIST("hFT0MAmp"), ft0mAmp);
    if (foundBC.has_zdc()) {
      const auto& zdcread = foundBC.zdc();
      const auto cent = collision.centFT0C();

      // ZDC data and histogram filling
      histos.get<TH1>(HIST("ZNAcoll"))->Fill(zdcread.amplitudeZNA());
      histos.get<TH1>(HIST("ZNCcoll"))->Fill(zdcread.amplitudeZNC());
      histos.get<TH2>(HIST("ZNvsZEMcoll"))->Fill(zdcread.amplitudeZEM1() + zdcread.amplitudeZEM2(), zdcread.amplitudeZNA() + zdcread.amplitudeZNC());
      histos.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(zdcread.amplitudeZNC(), zdcread.amplitudeZNA());

      histos.get<TH1>(HIST("ZEM1coll"))->Fill(zdcread.amplitudeZEM1());
      histos.get<TH1>(HIST("ZEM2coll"))->Fill(zdcread.amplitudeZEM2());

      float sumZNC = (zdcread.energySectorZNC())[0] + (zdcread.energySectorZNC())[1] + (zdcread.energySectorZNC())[2] + (zdcread.energySectorZNC())[3];
      float sumZNA = (zdcread.energySectorZNA())[0] + (zdcread.energySectorZNA())[1] + (zdcread.energySectorZNA())[2] + (zdcread.energySectorZNA())[3];
      float sumZPC = (zdcread.energySectorZPC())[0] + (zdcread.energySectorZPC())[1] + (zdcread.energySectorZPC())[2] + (zdcread.energySectorZPC())[3];
      float sumZPA = (zdcread.energySectorZPA())[0] + (zdcread.energySectorZPA())[1] + (zdcread.energySectorZPA())[2] + (zdcread.energySectorZPA())[3];
      float sumZDC = sumZPA + sumZPC + sumZNA + sumZNC;
      float sumZEM = zdcread.amplitudeZEM1() + zdcread.amplitudeZEM2();

      // common energies
      float commonSumZnc = (zdcread.energyCommonZNC()) / acceptanceZnc;
      float commonSumZna = (zdcread.energyCommonZNA()) / acceptanceZna;
      float commonSumZpc = (zdcread.energyCommonZPC()) / acceptanceZpc;
      float commonSumZpa = (zdcread.energyCommonZPA()) / acceptanceZpa;
      float sumZN = (sumZNC) + (sumZNA);
      float sumZP = (sumZPC) + (sumZPA);

      histos.fill(HIST("ZNenergy"), sumZN);
      histos.fill(HIST("ZPenergy"), sumZP);
      histos.fill(HIST("ZNCenergy"), commonSumZnc);
      histos.fill(HIST("ZNAenergy"), commonSumZna);
      histos.fill(HIST("ZPAenergy"), commonSumZpa);
      histos.fill(HIST("ZPCenergy"), commonSumZpc);
      histos.fill(HIST("hZNvsFT0Ccent"), cent, sumZN);
      histos.fill(HIST("hZPvsFT0Ccent"), cent, sumZP);
      histos.fill(HIST("hZNvsFT0CAmp"), ft0cAmp, sumZN);
      histos.fill(HIST("hZPvsFT0CAmp"), ft0cAmp, sumZP);
      histos.fill(HIST("hZNvsMult"), nTot, sumZN);
      histos.fill(HIST("hZPvsMult"), nTot, sumZP);

      float ratioZN = sumZNC / sumZNA;
      float ratioZP = sumZPC / sumZPA;
      pZNratiovscent->Fill(cent, ratioZN);
      pZPratiovscent->Fill(cent, ratioZP);
      pZNvsFT0Ccent->Fill(cent, sumZN);
      pZPvsFT0Ccent->Fill(cent, sumZP);

      histos.get<TH2>(HIST("ZDC_energy_vs_ZEM"))->Fill(sumZEM, sumZDC);

      // Spectator plane angle calculations and histograms
      const auto nTotZna = zdcread.amplitudeZNA();
      const auto nTotZnc = zdcread.amplitudeZNC();
      double qZnaReal = 0.0;
      double qZnaIm = 0.0;
      double qZncReal = 0.0;
      double qZncIm = 0.0;
      const double phiRadians[4] = {45 * o2::constants::math::PI / 180, 135 * o2::constants::math::PI / 180, 225 * o2::constants::math::PI / 180, 315 * o2::constants::math::PI / 180};
      std::complex<double> qZNA = std::complex<double>(0.0, 0.0);
      std::complex<double> qZNC = std::complex<double>(0.0, 0.0);

      for (int sector = 0; sector < 4; ++sector) {
        float energyZNA = zdcread.energySectorZNA()[sector];
        float energyZNC = zdcread.energySectorZNC()[sector];

        qZNA += std::complex<double>(std::cos(2 * phiRadians[sector]) * energyZNA / sumZNA, std::sin(2 * phiRadians[sector]) * energyZNA / sumZNA);
        qZNC += std::complex<double>(std::cos(2 * phiRadians[sector]) * energyZNC / sumZNC, std::sin(2 * phiRadians[sector]) * energyZNC / sumZNC);
      }

      qZnaReal = qZNA.real() / nTotZna;
      qZnaIm = qZNA.imag() / nTotZna;
      qZncReal = qZNC.real() / nTotZnc;
      qZncIm = qZNC.imag() / nTotZnc;

      histos.fill(HIST("Acorrelations"), qZNA.real(), qZNA.imag());
      histos.fill(HIST("RealQHistogramZNA"), qZnaReal);
      histos.fill(HIST("ImagQHistogramZNA"), qZnaIm);
      histos.fill(HIST("RealQHistogramZNC"), qZncReal);
      histos.fill(HIST("ImagQHistogramZNC"), qZncIm);

      // Calculate the spectator plane angles for ZNA and ZNC
      double psiZNA = std::atan2(qZNA.imag(), qZNA.real()) / 2.0;
      double psiZNC = std::atan2(qZNC.imag(), qZNC.real()) / 2.0;

      // Fill the histograms with the calculated angles
      histos.fill(HIST("SPAngleZNA"), psiZNA);
      histos.fill(HIST("SPAngleZNC"), psiZNC);

      double cosPsiDiff = std::cos(psiZNA) - std::cos(psiZNC);
      double sinPsiDiff = std::sin(psiZNA) - std::sin(psiZNC);

      sumCosPsiDiff += cosPsiDiff;
      sumSinPsiDiff += sinPsiDiff;
      ++countEvents;

      if (countEvents > 0) {
        double runningAverageCosPsiDiff = sumCosPsiDiff / countEvents;
        double runningAverageSinPsiDiff = sumSinPsiDiff / countEvents;
        histos.fill(HIST("RunningAverageCosPsiDiff"), runningAverageCosPsiDiff);
        pCosPsiDifferences->Fill(cent, runningAverageCosPsiDiff);
        pSinPsiDifferences->Fill(cent, runningAverageSinPsiDiff);
      }
      histos.fill(HIST("CosPsiDifferences"), cosPsiDiff);
      histos.fill(HIST("hSinDifferences"), sinPsiDiff);
    }
  }

  PROCESS_SWITCH(FlowZDCtask, processZdcCollAssoc, "Processing ZDC w. collision association", true);
  PROCESS_SWITCH(FlowZDCtask, processQVector, "Process before recentering", true);

}; // end of struct function

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowZDCtask>(cfgc)};
}
