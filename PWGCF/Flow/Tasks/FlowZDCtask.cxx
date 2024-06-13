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

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
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
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using aodZDCs = soa::Join<aod::ZDCMults, aod::Zdcs>;

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
  Configurable<float> MaxZP{"MaxZP", 3099.5, "Max ZP signal"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Z vertex cut"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> etaGap{"etaGap", 0.5, "Eta gap"};
  Configurable<float> minPt{"minPt", 0.2, "Minimum pt"};
  Configurable<float> maxPt{"maxPt", 20.0, "Maximum pt"};
  Configurable<float> MaxZEM{"MaxZEM", 3099.5, "Max ZEM signal"};
  // for ZDC info and analysis
  Configurable<int> nBinsADC{"nBinsADC", 1000, "nbinsADC"};
  Configurable<int> nBinsAmp{"nBinsAmp", 1025, "nbinsAmp"};
  Configurable<float> MaxZN{"MaxZN", 4099.5, "Max ZN signal"};

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  Partition<aodTracks> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  TComplex qTPC;       // init q TPC
  TComplex qZNA{0, 0}; // init qZNA
  TComplex qZNC{0, 0}; // init qZNC

  // Begin Histogram Registry

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TProfile> Q2_real_mean{TProfile("Q2_real_mean", "q2 real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2_imag_mean{TProfile("Q2_imag_mean", "q2 imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2_after{TProfile("Q2_after", "q2 recentered vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2before{TProfile("Q2_imag_mean", "q2 re vs imag", 10, 0, 100.)};
  OutputObj<TProfile> Q2_ZNA_real{TProfile("Q2_ZNA_real", "q2_ZNA real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2_ZNA_imag{TProfile("Q2_ZNA_imag", "q2_ZNA imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2_ZNC_real{TProfile("Q2_ZNC_real", "q2_ZNC real vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> Q2_ZNC_imag{TProfile("Q2_ZNC_imag", "q2_ZNC imag vs centrality", 10, 0, 100.)};
  OutputObj<TProfile> avgQ2TPCRe{TProfile("avgQ2TPCRe", "Average Q2 Real part vs Centrality", 10, 0, 100)};
  OutputObj<TProfile> avgQ2TPCIm{TProfile("avgQ2TPCIm", "Average Q2 Imaginary part vs Centrality", 10, 0, 100)};
  OutputObj<TProfile> ZDC_ZEM_Energy{TProfile("ZDC_ZEM_Energy", "ZDC vs ZEM Energy", 10, 0, 1000)};
  OutputObj<TProfile> pCosPsiDifferences{TProfile("pCosPsiDifferences", "Differences in cos(psi) vs Centrality;Centrality;Mean cos(psi) Difference", 200, 0, 100, -1, 1)};
  OutputObj<TProfile> pSinPsiDifferences{TProfile("pSinPsiDifferences", "Differences in sin(psi) vs Centrality;Centrality;Mean sin(psi) Difference", 200, 0, 100, -1, 1)};

  double sumCosPsiDiff = 0.0; // Sum of cos(psiZNC) - cos(psiZNA)
  int countEvents = 0;        // Count of processed events

  void init(InitContext const&)
  {
    // define axes
    const AxisSpec axisEta{30, -1.5, +1.5, "#eta"};
    const AxisSpec axispt{100, 0, 2, "#pt"};

    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisPhi{100, 0, 2 * TMath::Pi(), "#phi"};
    const AxisSpec axisQ{100, -1, 1, "Q"};
    const AxisSpec axisZNA{100, 0, 200, "energy"};
    const AxisSpec axisQZNA{100, -1, 1, "Q"};
    const AxisSpec axisREQ{100, -1, 1, "real Q"};
    const AxisSpec axisIMQ{100, -1, 1, "imag Q"};
    const AxisSpec axisEnergy{100, 0, 50, "energy"};

    AxisSpec axisVtxcounts{2, -0.5f, 1.5f, "Vtx info (0=no, 1=yes)"};
    AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
    AxisSpec axisCent{8, 0.f, 100.f, "centrality"};
    AxisSpec axisMult{1000, -0.5f, 1500.5f, "multiplicity"};
    AxisSpec axisMultTPC{1000, -0.5f, 1999.5f, "TPCmultiplicity"};
    AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};
    AxisSpec axisPtBins{{0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10., 13., 16., 20.}, "p_{T} (GeV/c)"};

    // create histograms
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});

    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("centHistogram", "centHistogram", kTH1F, {axisCent});
    histos.add("multHistogram", "multHistogram", kTH1F, {axisMult});
    histos.add("multvsCent", "centrality vs multiplicity", kTH2F, {axisCent, axisMult});
    histos.add("phiHistogram", "phiHistogram", kTH1F, {axisPhi});
    histos.add("TPCmultiplicity", "TPCmultiplicity", kTH1F, {axisMultTPC});

    histos.add("REqHistogram", "REqHistogram", kTH1F, {axisQ});
    histos.add("IMqHistogram", "IMqHistogram", kTH1F, {axisQ});

    histos.add("REqHistogramZNA", "REqHistogramZNA", kTH1F, {axisQZNA});
    histos.add("IMqHistogramZNA", "IMqHistogramZNA", kTH1F, {axisQZNA});

    histos.add("REqHistogramZNC", "REqHistogramZNC", kTH1F, {axisQZNA});
    histos.add("IMqHistogramZNC", "IMqHistogramZNC", kTH1F, {axisQZNA});

    histos.add("EnergyZNA", "ZNA Sector Energy", kTH1F, {axisEnergy});
    histos.add("EnergyZNC", "ZNC Sector Energy", kTH1F, {axisEnergy});
    // for q vector recentering
    histos.add("revsimag", "revsimag", kTH2F, {axisREQ, axisIMQ});

    if (doprocessZdcCollAssoc) { // Check if the process function for ZDCCollAssoc is enabled
      histos.add("ZNAcoll", "ZNAcoll; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      histos.add("ZNCcoll", "ZNCcoll; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZN}}});
      histos.add("ZEM1coll", "ZEM1coll; ZEM1 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      histos.add("ZEM2coll", "ZEM2coll; ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmp, -0.5, MaxZEM}}});
      histos.add("ZNvsZEMcoll", "ZNvsZEMcoll; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZEM}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
      histos.add("ZNAvsZNCcoll", "ZNAvsZNCcoll; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZN}, {nBinsAmp, -0.5, MaxZN}}}});

      histos.add("RealQHistogramZNA", "RealQHistogramZNA", kTH1F, {axisQZNA});
      histos.add("ImagQHistogramZNA", "ImagQHistogramZNA", kTH1F, {axisQZNA});
      histos.add("RealQHistogramZNC", "RealQHistogramZNC", kTH1F, {axisQZNA});
      histos.add("ImagQHistogramZNC", "ImagQHistogramZNC", kTH1F, {axisQZNA});

      histos.add("Acorrelations", "Acorrelations", kTH2F, {{axisQZNA}, {axisQZNA}});

      histos.add("SPAngleZNA", "Spectator Plane Angle ZNA;Angle (radians);Entries", {HistType::kTH1F, {{100, -TMath::Pi(), TMath::Pi()}}});
      histos.add("SPAngleZNC", "Spectator Plane Angle ZNC;Angle (radians);Entries", {HistType::kTH1F, {{100, -TMath::Pi(), TMath::Pi()}}});

      histos.add("RunningAverageCosPsiDiff", "Running Average of cos(psi) Differences;Running Average;Entries", {HistType::kTH1F, {{100, -1, 1}}});

      histos.add("CosPsiDifferences", "Differences in cos(psi);cos(psiZNC) - cos(psiZNA);Entries", {HistType::kTH1F, {{100, -2, 2}}});
      histos.add("hSinDifferences", "Differences in sin(psi);sin(psiZNC) - sin(psiZNA);Entries", {HistType::kTH1F, {{100, -2, 2}}});
      histos.add("CosPsiDifferencesAvg", "Differences in cos(psi);cos(psiZNC) - cos(psiZNA);Entries", {HistType::kTH2F, {{axisCent}, {100, -2, 2}}});
      histos.add("ZDC_energy_vs_ZEM", "ZDCvsZEM; ZEM; ZNA+ZNC+ZPA+ZPC", {HistType::kTH2F, {{{nBinsAmp, -0.5, MaxZEM}, {nBinsAmp, -0.5, 2. * MaxZN}}}});
    }
  }

  void processQVector(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::ZDCMults const& /*zdcMults*/)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    histos.fill(HIST("centHistogram"), collision.centFT0C());
    const auto& tracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    const int multTPC = tracksGrouped.size();
    const auto cent = collision.centFT0C();
    int Ntot = tracks.size();

    // this is the q vector for the TPC data. it is a complex function
    double qTPC_real = 0.0; // Initialize qTPC_real
    double qTPC_im = 0.0;   // init qTPC_imaginary

    Int_t multTrk = tracks.size(); // Tracks are already filtered with GlobalTrack || GlobalTrackSDD

    if (cent < 0.0 && cent > 70)
      return;
    TComplex qTPC(0, 0); // Starting with a q-vector of zero

    for (auto& track : tracks) {
      double phi = track.phi();

      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("phiHistogram"), track.phi());
      histos.fill(HIST("ptHistogram"), track.pt());

      qTPC += TComplex(TMath::Cos(2.0 * phi), TMath::Sin(2.0 * phi));

      histos.fill(HIST("multvsCent"), cent, Ntot);

    } // end track loop

    qTPC_real = qTPC.Re() / Ntot; // normalize these vectors by the total number of particles
    qTPC_im = qTPC.Im() / Ntot;

    histos.fill(HIST("REqHistogram"), qTPC_real);
    histos.fill(HIST("IMqHistogram"), qTPC_im);

    histos.fill(HIST("multHistogram"), multTrk);
    histos.fill(HIST("TPCmultiplicity"), multTPC);

    histos.fill(HIST("revsimag"), qTPC_real, qTPC_im);
  }
  void processZdcCollAssoc(
    ColEvSels const& cols,
    BCsRun3 const& /*bcs*/,
    aod::Zdcs const& /*zdcs*/)
  {
    double sumCosPsiDiff = 0.0; // initialize Sum of cosPsiDiff for averaging
    double sumSinPsiDiff = 0.0; // initialize Sum of cosPsiDiff for averaging
    int countEvents = 0;        // initialize Counter for the number of events processed
    // collision-based event selection
    for (auto& collision : cols) {
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if (foundBC.has_zdc()) {
        const auto& zdcread = foundBC.zdc();
        const auto cent = collision.centFT0C();

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
        histos.get<TH2>(HIST("ZDC_energy_vs_ZEM"))->Fill(sumZEM, sumZDC);

        const auto Ntot_ZNA = zdcread.amplitudeZNA();
        const auto Ntot_ZNC = zdcread.amplitudeZNC();
        double qZNA_real = 0.0; // Initialize qZNA_real
        double qZNA_im = 0.0;   // Initialize
        double qZNC_real = 0.0; // Initialize qZNC_real
        double qZNC_im = 0.0;   // Initialize
        const double phiRadians[4] = {45 * TMath::Pi() / 180, 135 * TMath::Pi() / 180, 225 * TMath::Pi() / 180, 315 * TMath::Pi() / 180};
        TComplex qZNA(0, 0), qZNC(0, 0);

        for (int sector = 0; sector < 4; ++sector) {
          // energy for current sector for ZNA and ZNC
          float energyZNA = zdcread.energySectorZNA()[sector];
          float energyZNC = zdcread.energySectorZNC()[sector];

          // Calculate q-vector from current sector and add it to the total q-vector
          qZNA += TComplex(TMath::Cos(2 * phiRadians[sector]) * energyZNA / sumZNA, TMath::Sin(2 * phiRadians[sector]) * energyZNA / sumZNA);
          qZNC += TComplex(TMath::Cos(2 * phiRadians[sector]) * energyZNC / sumZNC, TMath::Sin(2 * phiRadians[sector]) * energyZNC / sumZNC);
        }

        qZNA_real = qZNA.Re() / Ntot_ZNA;
        qZNA_im = qZNA.Im() / Ntot_ZNA;
        qZNC_real = qZNC.Re() / Ntot_ZNC;
        qZNC_im = qZNC.Im() / Ntot_ZNC;

        histos.fill(HIST("Acorrelations"), qZNA.Re(), qZNA.Im());
        histos.fill(HIST("RealQHistogramZNA"), qZNA_real);
        histos.fill(HIST("ImagQHistogramZNA"), qZNA_im);
        histos.fill(HIST("RealQHistogramZNC"), qZNC_real);
        histos.fill(HIST("ImagQHistogramZNC"), qZNC_im);

        // Calculate the spectator plane angles for ZNA and ZNC
        double psiZNA = TMath::ATan2(qZNA.Im(), qZNA.Re()) / 2;
        double psiZNC = TMath::ATan2(qZNC.Im(), qZNC.Re()) / 2;

        // Fill the histograms with the calculated angles
        histos.fill(HIST("SPAngleZNA"), psiZNA);
        histos.fill(HIST("SPAngleZNC"), psiZNC);

        double cosPsiDiff = TMath::Cos(psiZNA) - TMath::Cos(psiZNC);
        double sinPsiDiff = TMath::Sin(psiZNA) - TMath::Sin(psiZNC);

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
  }
  PROCESS_SWITCH(FlowZDCtask, processZdcCollAssoc, "Processing ZDC w. collision association", true);
  PROCESS_SWITCH(FlowZDCtask, processQVector, "Process before recentering", true);

}; // end of struct function

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowZDCtask>(cfgc)};
}
