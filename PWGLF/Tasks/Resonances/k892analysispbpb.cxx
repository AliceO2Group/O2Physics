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
/// \file k892analysispbpb.cxx
/// \brief K*0 spectra in Pb-Pb
/// \author Marta Urioni <marta.urioni@cern.ch>

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using std::array;
struct K892analysispbpb {
  SliceCache cache;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // histos
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};

  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};
  Configurable<int> cPIDBins{"cPIDBins", 65, "PID binning"};
  Configurable<float> cPIDQALimit{"cPIDQALimit", 6.5, "PID QA limit"};
  Configurable<int> cDCABins{"cDCABins", 300, "DCA binning"};
  Configurable<int> cPDGbins{"cPDGbins", 5000, "number of PDG bins"};
  Configurable<float> cPDGMax{"cPDGMax", 9500000.0f, "PDG limit"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "NoSameBunchPileUp and IsGoodZvtxFT0vsPV"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "Additional evsel3"};

  // presel
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutMaxOccupancy{"cfgCutMaxOccupancy", 2000.0f, "Accepted maximum Occupancy"};
  Configurable<bool> cfgApplyOccupancyCut{"cfgApplyOccupancyCut", false, "Apply maximum Occupancy"};

  // Track selections
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.0f, "TPC Crossed Rows to Findable Clusters"};

  Configurable<bool> cfgIsPhysicalPrimary{"cfgIsPhysicalPrimary", true, "Primary track selection in MC"}; // for MC bkg study

  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  Configurable<bool> cfgUseITSTPCrefit{"cfgUseITSTPCrefit", true, "Use ITS and TPC refit"};
  Configurable<float> cfgITSChi2Ncl{"cfgITSChi2Ncl", 999.0, "ITS Chi2/NCl"};
  Configurable<float> cfgTPCChi2Ncl{"cfgTPCChi2Ncl", 999.0, "TPC Chi2/NCl"};

  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"}; // TPC
  Configurable<float> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"}; // TOF
  Configurable<float> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<float> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"}; // TOF
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};            // By pass TOF PID selection
  Configurable<bool> cTofBetaCut{"cTofBetaCut", false, "selection on TOF beta"};

  Configurable<bool> cTPClowpt{"cTPClowpt", true, "apply TPC at low pt"};
  Configurable<bool> cTOFonlyHighpt{"cTOFonlyHighpt", false, "apply TOF only at high pt"};
  Configurable<bool> cTOFandTPCHighpt{"cTOFandTPCHighpt", false, "apply TOF and TPC at high pt"};

  // rotational bkg
  Configurable<int> cfgNoRotations{"cfgNoRotations", 3, "Number of rotations per pair for rotbkg"};
  Configurable<int> rotationalCut{"rotationalCut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<bool> cfgRotPi{"cfgRotPi", true, "rotate Pion"};

  // event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "Mixing bins - z-vertex"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};

  // cuts on mother
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", false, "Enamble additional cuts on mother"};
  Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
  Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum Minv of mother cut"};

  // configurables for partitions
  Configurable<float> cMaxPtTPC{"cMaxPtTPC", 1.2, "maximum pt to apply TPC PID and TOF if available"};
  Configurable<float> cMinPtTOF{"cMinPtTOF", 0.8, "minimum pt to require TOF PID in addition to TPC"};

  // plots
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAeventPlots{"additionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> additionalMEPlots{"additionalMEPlots", false, "Additional Mixed event plots"};

  // MC
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};

  TRandom* rand = new TRandom();

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{z}} (cm)"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisMom = {binsPt, "Mom #it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisDau = {binsPtQA, "Dau #it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};
    AxisSpec pdgCodeAxis = {cPDGbins, 0, cPDGMax};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};

    if ((!doprocessMC && !doprocessMCRun2) || doprocessMixedEventMC || doprocessMixedEventMCRun2) {
      // event histograms
      histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
      histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});

      if (additionalQAeventPlots) {
        // Test on Mixed event
        histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
        histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
        histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
        histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
        histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
        histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});

        // event histograms
        histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
        histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
        histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
      }
    }

    // Mass QA (quick check)
    histos.add("k892invmassDS", "Invariant mass of K(892)0 different sign", kTH1F, {invMassAxis});
    histos.add("k892invmassDSAnti", "Invariant mass of Anti-K(892)0 different sign", kTH1F, {invMassAxis});
    histos.add("k892invmassLS", "Invariant mass of K(892)0 like sign", kTH1F, {invMassAxis});
    histos.add("k892invmassLSAnti", "Invariant mass of Anti-K(892)0 like sign", kTH1F, {invMassAxis});
    if (doprocessMixedEvent || doprocessMixedEventRun2 || doprocessMixedEventMC || doprocessMixedEventMCRun2) {
      histos.add("k892invmassME", "Invariant mass of K(892)0 mixed event", kTH1F, {invMassAxis});
      if (additionalMEPlots) {
        histos.add("k892invmassME_DS", "Invariant mass of K(892)0 mixed event DS", kTH1F, {invMassAxis});
        histos.add("k892invmassME_DSAnti", "Invariant mass of K(892)0 mixed event DSAnti", kTH1F, {invMassAxis});
      }
    }

    if (additionalQAplots) {
      // TPC ncluster distirbutions
      histos.add("Ncluster/TPCnclusterpi", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("Ncluster/TPCnclusterka", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("Ncluster/TPCnclusterPhipi", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
      histos.add("Ncluster/TPCnclusterPhika", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});

      histos.add("Ncluster/TPCChi2ncluster", "TPC Chi2ncluster distribution", kTH1F, {{100, 0, 10, "TPC Chi2nCluster"}});
      histos.add("Ncluster/ITSChi2ncluster", "ITS Chi2ncluster distribution", kTH1F, {{100, 0, 40, "ITS Chi2nCluster"}});
      histos.add("Ncluster/ITSncluster", "ITS  ncluster distribution", kTH1F, {{10, 0, 10, "ITS nCluster"}});

      histos.add("QA/h2k892ptMothervsptPiDS", "Pt of K(892)0 differnt sign vs pt pion daughter", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QA/h2k892ptMothervsptPiDSAnti", "Pt of Anti-K(892)0 differnt sign vs pt pion daughter", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QA/h2k892ptMothervsptKaDS", "Pt of K(892)0 differnt sign vs pt kaon daughter", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QA/h2k892ptMothervsptKaDSAnti", "Pt of Anti-K(892)0 differnt sign vs pt kaon daughter", kTH2F, {ptAxisMom, ptAxisDau});

      histos.add("QAME/h2k892ptMothervsptPiDS", "Pt of Mother vs pt pion daughter, Mixed Event", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QAME/h2k892ptMothervsptPiDSAnti", "Pt of Anti-Mother vs pt pion daughter, Mixed Event", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QAME/h2k892ptMothervsptKaDS", "Pt of Mother vs pt kaon daughter, Mixed Event", kTH2F, {ptAxisMom, ptAxisDau});
      histos.add("QAME/h2k892ptMothervsptKaDSAnti", "Pt of Anti-Mother vs pt pion daughter, Mixed Event", kTH2F, {ptAxisMom, ptAxisDau});
    }

    // DCA QA
    histos.add("QA/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    // pT QA
    histos.add("QA/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QA/trkpT_ka", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    // PID QA
    histos.add("QA/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QA/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    histos.add("QA/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    histos.add("QA/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QA/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    histos.add("QA/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

    // inv mass histograms
    histos.add("h3k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassDSAnti", "Invariant mass of Anti-K(892)0 differnt sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLS", "Invariant mass of K(892)0 same sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLSAnti", "Invariant mass of Anti-K(892)0 same sign", kTH3F, {centAxis, ptAxis, invMassAxis});

    if (doprocessRotationalBkg || doprocessRotationalBkgMC) {
      histos.add("k892invmassRotDS", "Invariant mass of K(892)0 RotBkg", kTH1F, {invMassAxis});
      histos.add("k892invmassRotDSAnti", "Invariant mass of Anti-K(892)0 RotBkg", kTH1F, {invMassAxis});

      histos.add("h3k892invmassRotDS", "Invariant mass of K(892)0 Rotational Bkg", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassRotDSAnti", "Invariant mass of Anti-K(892)0 Rotational Bkg", kTH3F, {centAxis, ptAxis, invMassAxis});
    }

    if (doprocessMixedEvent || doprocessMixedEventRun2 || doprocessMixedEventMC || doprocessMixedEventMCRun2) {
      histos.add("h3k892invmassME", "Invariant mass of K(892)0 mixed event", kTH3F, {centAxis, ptAxis, invMassAxis});

      if (additionalMEPlots) {
        histos.add("QAME/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAME/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAME/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAME/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAME/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAME/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

        histos.add("h3k892invmassME_DS", "Invariant mass of K(892)0 mixed event DS", kTH3F, {centAxis, ptAxis, invMassAxis});
        histos.add("h3k892invmassME_DSAnti", "Invariant mass of K(892)0 mixed event DSAnti", kTH3F, {centAxis, ptAxis, invMassAxis});
      }
    }

    if (doprocessMixedEventMC || doprocessMixedEventMCRun2) {
      histos.add("h3k892invmassWrongDaughtersME_DS", "Invariant mass ME with wrong daughters DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassWrongDaughtersME_DSAnti", "Invariant mass ME with wrong daughters DS anti", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassRightDaughtersME_DS", "Invariant mass ME with right daughters DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassRightDaughtersME_DSAnti", "Invariant mass ME with right daughters DS anti", kTH3F, {centAxis, ptAxis, invMassAxis});
    }

    if (doprocessMC || doprocessMCRun2) {
      histos.add("QAevent/hMCrecCollSels", "MC Event statistics", HistType::kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("QAevent/hMultiplicityPercentMC", "Multiplicity percentile of MCrec collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

      histos.add("h1k892Recsplit", "k892 Rec split", HistType::kTH1F, {{200, 0.0f, 20.0f}});
      // MC QA
      histos.add("QAMCTrue/hGlobalIndexMotherRec", "index of rec mothers", HistType::kTH1F, {{static_cast<int>(1e5), 0.0f, 1e5f}});
      histos.add("QAMCTrue/hGlobalIndexMotherGen", "index of gen mothers", HistType::kTH1F, {{static_cast<int>(1e5), 0.0f, 1e5f}});
      histos.add("QAMCTrue/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("k892Recinvmass", "Inv mass distribution of Reconstructed MC K(892)", kTH1F, {invMassAxis});
      histos.add("k892RecinvmassAnti", "Inv mass distribution of Reconstructed MC AntiK(892)", kTH1F, {invMassAxis});
      histos.add("k892GenInvmass", "Invariant mass of generated K(892)0", kTH1F, {invMassAxis});
      histos.add("k892GenInvmassAnti", "Invariant mass of generated Anti-K(892)0", kTH1F, {invMassAxis});

      histos.add("h3Reck892invmass", "Invariant mass of Reconstructed MC K(892)0", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3Reck892invmassAnti", "Invariant mass of Reconstructed MC Anti-K(892)0", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("k892Gen", "pT distribution of True MC K(892)0", kTH3F, {mcLabelAxis, ptAxis, centAxis});
      histos.add("k892GenAnti", "pT distribution of True MC Anti-K(892)0", kTH3F, {mcLabelAxis, ptAxis, centAxis});
      histos.add("k892Rec", "pT distribution of Reconstructed MC K(892)0", kTH2F, {ptAxis, centAxis});
      histos.add("k892RecAnti", "pT distribution of Reconstructed MC Anti-K(892)0", kTH2F, {ptAxis, centAxis});
      histos.add("h3k892GenInvmass", "Invariant mass of generated K(892)0", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892GenInvmassAnti", "Invariant mass of generated Anti-K(892)0", kTH3F, {centAxis, ptAxis, invMassAxis});

      histos.add("h3Reck892invmassPtGen", "Invariant mass of Reconstructed MC K(892)0 with Pt Gen", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3Reck892invmassAntiPtGen", "Invariant mass of Reconstructed MC Anti-K(892)0 with Pt Gen", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3PtRecvsPtGenAnti", "reconstructed K* Pt vs generated K* pt", kTH3F, {centAxis, ptAxis, ptAxis});
      histos.add("h3PtRecvsPtGen", "reconstructed Anti-K* Pt vs generated Anti-K* pt", kTH3F, {centAxis, ptAxis, ptAxis});

      histos.add("h3k892invmassWrongDaughters_DS", "Invariant mass of K*0 with wrong daughters DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassWrongDaughters_DSAnti", "Invariant mass of K*0 with wrong daughters DS anti", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassRightDaughters_DS", "Invariant mass of K*0 with right daughters DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassRightDaughters_DSAnti", "Invariant mass of K*0 with right daughters DS anti", kTH3F, {centAxis, ptAxis, invMassAxis});

      histos.add("h3k892invmassSameMother_DS", "Invariant mass same mother DS", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3k892invmassSameMother_DSAnti", "Invariant mass same mother DS anti", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("h3PdgCodeSameMother_DS", "PDG code same mother DS", kTH3F, {centAxis, ptAxis, pdgCodeAxis});
      histos.add("h3PdgCodeSameMother_DSAnti", "PDG code same mother DS anti", kTH3F, {centAxis, ptAxis, pdgCodeAxis});
    }

    if (doprocessEvtLossSigLossMC) {
      histos.add("QAevent/hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("QAevent/hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      histos.add("QAevent/hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {{120, 0.0f, 120.0f}, impactParAxis});
      histos.add("QAevent/k892genBeforeEvtSel", "K* before event selections", kTH2F, {ptAxis, impactParAxis});
      histos.add("QAevent/k892genBeforeEvtSelAnti", "K* before event selections", kTH2F, {ptAxis, impactParAxis});
      histos.add("QAevent/k892genAfterEvtSel", "K* after event selections", kTH2F, {ptAxis, impactParAxis});
      histos.add("QAevent/k892genAfterEvtSelAnti", "K* after event selections", kTH2F, {ptAxis, impactParAxis});
    }

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;

  template <typename CollType>
  bool myEventSelections(const CollType& coll)
  {
    if (!coll.sel8())
      return false;
    if (std::abs(coll.posZ()) > cfgCutVertex)
      return false;
    if (timFrameEvsel && (!coll.selection_bit(aod::evsel::kNoTimeFrameBorder) || !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)))
      return false;
    if (additionalEvSel2 && (!coll.selection_bit(aod::evsel::kNoSameBunchPileup) || !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)))
      return false;
    if (additionalEvSel3 && (!coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)))
      return false;
    auto centrality = coll.centFT0C();
    if (centrality > cfgCutCentrality)
      return false;
    auto occupancy = coll.trackOccupancyInTimeRange();
    if (cfgApplyOccupancyCut && (occupancy > cfgCutMaxOccupancy))
      return false;

    return true;
  }

  template <typename CollType, typename bcType>
  bool myEventSelectionsRun2(const CollType& coll, const bcType&)
  {
    auto bc = coll.template bc_as<bcType>();
    if (!(bc.eventCuts() & BIT(aod::Run2EventCuts::kAliEventCutsAccepted)))
      return false;
    if (std::abs(coll.posZ()) > cfgCutVertex)
      return false;
    auto centrality = coll.centRun2V0M();
    if (centrality > cfgCutCentrality)
      return false;

    return true;
  }

  template <typename TrackType>
  bool trackCut(const TrackType& track)
  {
    // basic track cuts
    if (track.itsNCls() < cfgITScluster)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (track.itsChi2NCl() > cfgITSChi2Ncl)
      return false;
    if (track.tpcChi2NCl() > cfgTPCChi2Ncl)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgRatioTPCRowsOverFindableCls)
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgUseITSTPCrefit && (!(o2::aod::track::ITSrefit) || !(o2::aod::track::TPCrefit)))
      return false;

    return true;
  }

  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    if (cTOFonlyHighpt) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) <= cMaxTOFnSigmaKaon) { // tof cut only
        return true;
      }

    } else if (cTOFandTPCHighpt) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) <= cMaxTOFnSigmaKaon && candidate.hasTPC() && std::abs(candidate.tpcNSigmaKa()) <= cMaxTPCnSigmaKaon) { // tof and tpc cut
        return true;
      }

    } else {

      if (candidate.hasTPC() && std::abs(candidate.tpcNSigmaKa()) <= cMaxTPCnSigmaKaon) { // tpc cut, tof when available

        if (cTofBetaCut && candidate.hasTOF() && (candidate.beta() + 3 * candidate.betaerror() > 1))
          return false;

        if (cByPassTOF) // skip tof selection
          return true;

        if (candidate.hasTOF()) {
          if (std::abs(candidate.tofNSigmaKa()) <= cMaxTOFnSigmaKaon) {
            return true;
          }
        } else {
          return true;
        }
      }
    }

    return false;
  }

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    if (cTOFonlyHighpt) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) <= cMaxTOFnSigmaPion) { // tof cut only
        return true;
      }

    } else if (cTOFandTPCHighpt) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) <= cMaxTOFnSigmaPion && candidate.hasTPC() && std::abs(candidate.tpcNSigmaPi()) <= cMaxTPCnSigmaPion) { // tof and tpc cut
        return true;
      }

    } else {

      if (candidate.hasTPC() && std::abs(candidate.tpcNSigmaPi()) <= cMaxTPCnSigmaPion) { // tpc cut, tof when available

        if (cTofBetaCut && candidate.hasTOF() && (candidate.beta() + 3 * candidate.betaerror() > 1))
          return false;

        if (cByPassTOF) // skip tof selection
          return true;

        if (candidate.hasTOF()) {
          if (std::abs(candidate.tofNSigmaPi()) <= cMaxTOFnSigmaPion) {
            return true;
          }
        } else {
          return true;
        }
      }
    }

    return false;
  }

  template <bool IsMC, bool IsMix, bool IsRot, bool IsRun2, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {

    auto multiplicity = -999;

    if constexpr (!IsRun2)
      multiplicity = collision.centFT0C();
    else
      multiplicity = collision.centRun2V0M();

    auto oldindex = -999;
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance, ldaughterRot, lResonanceRot;
    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (additionalQAeventPlots) {
        if constexpr (!IsMC && !IsRot) {
          if constexpr (!IsMix) {
            histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
          } else {
            histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
          }
        }
      }

      //// Initialize variables
      // Trk1: Pion, Trk2: Kaon
      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();
      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      if (!selectionPIDPion(trk1) || !selectionPIDKaon(trk2))
        continue;

      if constexpr (IsMC) {
        if (cTPClowpt) {
          if (trk1ptPi >= cMaxPtTPC || trk2ptKa >= cMaxPtTPC)
            continue;
        } else if (cTOFonlyHighpt || cTOFandTPCHighpt) {
          if (trk1ptPi <= cMinPtTOF || trk2ptKa <= cMinPtTOF)
            continue;
        }
      }

      if (additionalQAplots && !IsMix && !IsRot) {
        // TPCncluster distributions
        histos.fill(HIST("Ncluster/TPCnclusterpi"), trk1.tpcNClsFound());
        histos.fill(HIST("Ncluster/TPCnclusterka"), trk2.tpcNClsFound());
        histos.fill(HIST("Ncluster/TPCnclusterPhipi"), trk1.tpcNClsFound(), trk1.phi());
        histos.fill(HIST("Ncluster/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());
        histos.fill(HIST("Ncluster/TPCChi2ncluster"), trk1.tpcChi2NCl());
        histos.fill(HIST("Ncluster/ITSChi2ncluster"), trk1.itsChi2NCl());
        histos.fill(HIST("Ncluster/ITSncluster"), trk1.itsNCls());
      }

      if constexpr (!IsMix && !IsRot) {
        //// QA plots after the selection
        //  --- PID QA Pion
        histos.fill(HIST("QA/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
          histos.fill(HIST("QA/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
        }
        //  --- PID QA Kaon
        histos.fill(HIST("QA/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QA/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QA/trkpT_ka"), trk2ptKa);
        histos.fill(HIST("QA/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QA/trkDCAxy_ka"), trk2.dcaXY());
        histos.fill(HIST("QA/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QA/trkDCAz_ka"), trk2.dcaZ());
      } else if (IsMix && additionalMEPlots) {
        //  --- PID QA Pion
        histos.fill(HIST("QAME/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QAME/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
          histos.fill(HIST("QAME/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
        }
        //  --- PID QA Kaon
        histos.fill(HIST("QAME/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QAME/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QAME/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
      }

      int track1Sign = trk1.sign();
      int track2Sign = trk2.sign();

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (std::abs(lResonance.Rapidity()) >= 0.5)
        continue;
      if (cfgCutsOnMother && !IsRot) {
        if (lResonance.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonance.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      //// Un-like sign pair only
      if (track1Sign * track2Sign < 0) {
        if constexpr (IsRot) { // rotational background
          for (int i = 0; i < cfgNoRotations; i++) {
            float theta = rand->Uniform(o2::constants::math::PI - o2::constants::math::PI / rotationalCut, o2::constants::math::PI + o2::constants::math::PI / rotationalCut);
            if (cfgRotPi) {
              ldaughterRot.SetPtEtaPhiM(trk1.pt(), trk1.eta(), trk1.phi() + theta, massPi);
              lResonanceRot = lDecayDaughter2 + ldaughterRot;
            } else {
              ldaughterRot.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi() + theta, massKa);
              lResonanceRot = lDecayDaughter1 + ldaughterRot;
            }
            if (std::abs(lResonanceRot.Rapidity()) >= 0.5)
              continue;
            if (cfgCutsOnMother) {
              if (lResonanceRot.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
                continue;
              if (lResonanceRot.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
                continue;
            }

            if (track1Sign < 0) {
              histos.fill(HIST("k892invmassRotDS"), lResonanceRot.M());
              histos.fill(HIST("h3k892invmassRotDS"), multiplicity, lResonanceRot.Pt(), lResonanceRot.M());
            } else if (track1Sign > 0) {
              histos.fill(HIST("k892invmassRotDSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassRotDSAnti"), multiplicity, lResonanceRot.Pt(), lResonanceRot.M());
            }
          }

        } else if constexpr (!IsMix) { // same event
          if (track1Sign < 0) {
            histos.fill(HIST("k892invmassDS"), lResonance.M());
            histos.fill(HIST("h3k892invmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
            if (additionalQAplots) {
              histos.fill(HIST("QA/h2k892ptMothervsptPiDS"), lResonance.Pt(), lDecayDaughter1.Pt());
              histos.fill(HIST("QA/h2k892ptMothervsptKaDS"), lResonance.Pt(), lDecayDaughter2.Pt());
            }
          } else if (track1Sign > 0) {
            histos.fill(HIST("k892invmassDSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassDSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            if (additionalQAplots) {
              histos.fill(HIST("QA/h2k892ptMothervsptPiDSAnti"), lResonance.Pt(), lDecayDaughter1.Pt());
              histos.fill(HIST("QA/h2k892ptMothervsptKaDSAnti"), lResonance.Pt(), lDecayDaughter2.Pt());
            }
          }

        } else { // mixed event
          histos.fill(HIST("k892invmassME"), lResonance.M());
          histos.fill(HIST("h3k892invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (additionalMEPlots) {
            if (track1Sign < 0) {
              histos.fill(HIST("k892invmassME_DS"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
              if (additionalQAplots) {
                histos.fill(HIST("QAME/h2k892ptMothervsptPiDS"), lResonance.Pt(), lDecayDaughter1.Pt());
                histos.fill(HIST("QAME/h2k892ptMothervsptKaDS"), lResonance.Pt(), lDecayDaughter2.Pt());
              }
            } else if (track1Sign > 0) {
              histos.fill(HIST("k892invmassME_DSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
              if (additionalQAplots) {
                histos.fill(HIST("QAME/h2k892ptMothervsptPiDSAnti"), lResonance.Pt(), lDecayDaughter1.Pt());
                histos.fill(HIST("QAME/h2k892ptMothervsptKaDSAnti"), lResonance.Pt(), lDecayDaughter2.Pt());
              }
            }
          }
        }

        // MC
        if constexpr (IsMC && !IsRot) {

          if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
            continue;

          const auto mctrack1 = trk1.mcParticle();
          const auto mctrack2 = trk2.mcParticle();
          int track1PDG = std::abs(mctrack1.pdgCode());
          int track2PDG = std::abs(mctrack2.pdgCode());

          if (cfgIsPhysicalPrimary && (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary()))
            continue;

          if (track1PDG != 211 || track2PDG != 321) {

            if (track1Sign < 0) {
              if constexpr (IsMix)
                histos.fill(HIST("h3k892invmassWrongDaughtersME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
              else
                histos.fill(HIST("h3k892invmassWrongDaughters_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (track1Sign > 0) {
              if constexpr (IsMix)
                histos.fill(HIST("h3k892invmassWrongDaughtersME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
              else
                histos.fill(HIST("h3k892invmassWrongDaughters_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }

            continue;
          }

          if (track1Sign < 0) {
            if constexpr (IsMix)
              histos.fill(HIST("h3k892invmassRightDaughtersME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            else
              histos.fill(HIST("h3k892invmassRightDaughters_DS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (track1Sign > 0) {
            if constexpr (IsMix)
              histos.fill(HIST("h3k892invmassRightDaughtersME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            else
              histos.fill(HIST("h3k892invmassRightDaughters_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }

          if constexpr (!IsMix) {

            bool isSameMother = false;
            bool isMotherOk = false;
            int pdgCodeMother = -999;
            float ptMother = -9999.;
            for (const auto& mothertrack1 : mctrack1.template mothers_as<aod::McParticles>()) {
              for (const auto& mothertrack2 : mctrack2.template mothers_as<aod::McParticles>()) {
                if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
                  continue;
                if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
                  continue;

                if (std::abs(mothertrack1.pdgCode()) == 1000822080) // Pb PDG code
                  continue;

                pdgCodeMother = mothertrack1.pdgCode();
                ptMother = mothertrack1.pt();
                isSameMother = true;

                if (std::abs(mothertrack1.pdgCode()) != 313)
                  continue;

                if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                  histos.fill(HIST("h1k892Recsplit"), mothertrack1.pt());
                  continue;
                }
                oldindex = mothertrack1.globalIndex();
                isMotherOk = true;
              }
            }

            if (isSameMother) {
              if (track1Sign < 0) {
                histos.fill(HIST("h3k892invmassSameMother_DS"), multiplicity, lResonance.Pt(), lResonance.M());
                histos.fill(HIST("h3PdgCodeSameMother_DS"), multiplicity, lResonance.Pt(), pdgCodeMother);
              } else if (track1Sign > 0) {
                histos.fill(HIST("h3k892invmassSameMother_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
                histos.fill(HIST("h3PdgCodeSameMother_DSAnti"), multiplicity, lResonance.Pt(), pdgCodeMother);
              }
            }

            if (!isMotherOk)
              continue;

            histos.fill(HIST("QAMCTrue/hGlobalIndexMotherRec"), oldindex);
            // Track selection check.
            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
            if (isTrk1hasTOF) {
              histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            }
            histos.fill(HIST("QAMCTrue/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
            if (isTrk2hasTOF) {
              histos.fill(HIST("QAMCTrue/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
            }

            // MC histograms
            if (pdgCodeMother > 0) {
              histos.fill(HIST("k892Rec"), lResonance.Pt(), multiplicity);
              histos.fill(HIST("k892Recinvmass"), lResonance.M());
              histos.fill(HIST("h3Reck892invmass"), multiplicity, lResonance.Pt(), lResonance.M());
              histos.fill(HIST("h3Reck892invmassPtGen"), multiplicity, ptMother, lResonance.M());
              histos.fill(HIST("h3PtRecvsPtGen"), multiplicity, lResonance.Pt(), ptMother);
            } else {
              histos.fill(HIST("k892RecAnti"), lResonance.Pt(), multiplicity);
              histos.fill(HIST("k892RecinvmassAnti"), lResonance.M());
              histos.fill(HIST("h3Reck892invmassAnti"), multiplicity, lResonance.Pt(), lResonance.M());
              histos.fill(HIST("h3Reck892invmassAntiPtGen"), multiplicity, ptMother, lResonance.M());
              histos.fill(HIST("h3PtRecvsPtGenAnti"), multiplicity, lResonance.Pt(), ptMother);
            }
          }
        } // end of IsMC

      } else if (track1Sign * track2Sign > 0) {
        if constexpr (!IsMix) {
          if (track1Sign < 0) {
            histos.fill(HIST("k892invmassLS"), lResonance.M());
            histos.fill(HIST("h3k892invmassLS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (track1Sign > 0) {
            histos.fill(HIST("k892invmassLSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassLSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      } // end on DS or LS if
    } // end of loop on track combinations
  } // end of fill histograms

  Filter collisionFilter = nabs(aod::collision::posZ) <= cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) <= cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) >= cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) <= cfgCutDCAxy) && (nabs(aod::track::dcaZ) <= cfgCutDCAz);

  // Data
  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFbeta>>;
  // MC
  using EventCandidatesMCrec = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using TrackCandidatesMCrec = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFbeta, aod::McTrackLabels>>;
  // ME run 3
  using BinningTypeVtxCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;

  // Data Run 2
  using Run2Events = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s>; //, aod::TrackletMults>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  // MC Run2
  using EventCandidatesMCrecRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentRun2V0Ms, aod::CentRun2CL0s>; // aod::TrackletMults>;
  // ME run 2
  using BinningTypeVtxCentRun2 = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentRun2V0M>;

  // partitions tpc low pt
  Partition<TrackCandidates> negPitpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) < cMaxPtTPC);
  Partition<TrackCandidates> posKatpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) < cMaxPtTPC);

  Partition<TrackCandidates> posPitpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) < cMaxPtTPC);
  Partition<TrackCandidates> negKatpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) < cMaxPtTPC);

  // tpc & tof, high pt
  Partition<TrackCandidates> negPitoftpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> posKatoftpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  Partition<TrackCandidates> posPitoftpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> negKatoftpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  // tof only, high pt
  Partition<TrackCandidates> negPitof = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> posKatof = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  Partition<TrackCandidates> posPitof = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> negKatof = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  template <bool IsMC, bool IsMix, bool IsRot, bool IsRun2, typename CollisionType, typename TracksType>
  void callFillHistoswithPartitions(const CollisionType& collision1, const TracksType&, const CollisionType& collision2, const TracksType&)
  {
    if (cTPClowpt) {
      //+-
      auto candPosPitpc = posPitpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candNegKatpc = negKatpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candPosPitpc, candNegKatpc);

      //-+
      auto candNegPitpc = negPitpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candPosKatpc = posKatpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candNegPitpc, candPosKatpc);

    } else if (cTOFandTPCHighpt) {
      //+-
      auto candPosPitoftpc = posPitoftpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candNegKatoftpc = negKatoftpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candPosPitoftpc, candNegKatoftpc);

      //-+
      auto candNegPitoftpc = negPitoftpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candPosKatoftpc = posKatoftpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candNegPitoftpc, candPosKatoftpc);

    } else if (cTOFonlyHighpt) {
      //+-
      auto candPosPitof = posPitof->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candNegKatof = negKatof->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candPosPitof, candNegKatof);

      //-+
      auto candNegPitof = negPitof->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto candPosKatof = posKatof->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      fillHistograms<IsMC, IsMix, IsRot, IsRun2>(collision1, candNegPitof, candPosKatof);
    }
  }

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {

    if (!myEventSelections(collision))
      return;

    auto centrality = collision.centFT0C();

    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1);
    histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
    histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), centrality);

    if (additionalQAeventPlots) {
      histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
      histos.fill(HIST("TestME/hnTrksSameE"), tracks.size());
    }
    //                            <IsMC, IsMix, IsRot, IsRun2>
    callFillHistoswithPartitions<false, false, false, false>(collision, tracks, collision, tracks);
  }
  PROCESS_SWITCH(K892analysispbpb, processSameEvent, "Process Same event", true);

  void processSameEventRun2(Run2Events::iterator const& collision, TrackCandidates const& tracks, BCsWithRun2Info const& bcs)
  {

    if (!myEventSelectionsRun2(collision, bcs))
      return;

    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1);
    histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
    histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), collision.centRun2V0M());

    if (additionalQAeventPlots) {
      histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
      histos.fill(HIST("TestME/hnTrksSameE"), tracks.size());
    }

    //                            <IsMC, IsMix, IsRot, IsRun2>
    callFillHistoswithPartitions<false, false, false, true>(collision, tracks, collision, tracks);
  }
  PROCESS_SWITCH(K892analysispbpb, processSameEventRun2, "Process Same event  Run2", false);

  void processRotationalBkg(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {

    if (!myEventSelections(collision))
      return;

    //                            <IsMC, IsMix, IsRot, IsRun2>
    callFillHistoswithPartitions<false, false, true, false>(collision, tracks, collision, tracks);
  }
  PROCESS_SWITCH(K892analysispbpb, processRotationalBkg, "Process Rotational Background", false);

  void processRotationalBkgMC(EventCandidatesMCrec::iterator const& recCollision, TrackCandidatesMCrec const& RecTracks)
  {

    if (!myEventSelections(recCollision))
      return;

    //             <IsMC, IsMix, IsRot, IsRun2>
    fillHistograms<true, false, true, false>(recCollision, RecTracks, RecTracks);
  }
  PROCESS_SWITCH(K892analysispbpb, processRotationalBkgMC, "Process Rotational Background MC", false);

  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVtxCent colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVtxCent> pairs{colBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      if (!myEventSelections(collision1) || !myEventSelections(collision2))
        continue;

      auto centrality = collision1.centFT0C();

      if (additionalQAeventPlots) {
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), centrality);
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision1.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), tracks1.size());
      }

      //                          <IsMC, IsMix, IsRot, IsRun2>
      callFillHistoswithPartitions<false, true, false, false>(collision1, tracks1, collision2, tracks2);
    }
  }
  PROCESS_SWITCH(K892analysispbpb, processMixedEvent, "Process Mixed event", true);

  void processMixedEventRun2(Run2Events const& collisions, TrackCandidates const& tracks, BCsWithRun2Info const& bcs)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVtxCentRun2 colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<Run2Events, TrackCandidates, BinningTypeVtxCentRun2> pairs{colBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      if (!myEventSelectionsRun2(collision1, bcs) || !myEventSelectionsRun2(collision2, bcs))
        continue;

      if (additionalQAeventPlots) {
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision1.centRun2V0M());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision1.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), tracks1.size());
      }

      //                          <IsMC, IsMix, IsRot, IsRun2>
      callFillHistoswithPartitions<false, true, false, true>(collision1, tracks1, collision2, tracks2);
    }
  }
  PROCESS_SWITCH(K892analysispbpb, processMixedEventRun2, "Process Mixed event Run2", false);

  void processMixedEventMC(EventCandidatesMCrec const& recCollisions, TrackCandidatesMCrec const& RecTracks, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(RecTracks);
    BinningTypeVtxCent colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<EventCandidatesMCrec, TrackCandidatesMCrec, BinningTypeVtxCent> pairs{colBinning, cfgNoMixedEvents, -1, recCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      if (!myEventSelections(collision1) || !myEventSelections(collision2))
        continue;

      if (additionalQAeventPlots) {
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision1.centFT0C());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision1.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), tracks1.size());
      }

      //            <IsMC, IsMix, IsRot, IsRun2>
      fillHistograms<true, true, false, false>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(K892analysispbpb, processMixedEventMC, "Process Mixed event MC", false);

  void processMixedEventMCRun2(EventCandidatesMCrecRun2 const& recCollisions, TrackCandidatesMCrec const& RecTracks, BCsWithRun2Info const& bcs, aod::McParticles const&)
  {
    auto tracksTuple = std::make_tuple(RecTracks);
    BinningTypeVtxCentRun2 colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<EventCandidatesMCrecRun2, TrackCandidatesMCrec, BinningTypeVtxCentRun2> pairs{colBinning, cfgNoMixedEvents, -1, recCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      if (!myEventSelectionsRun2(collision1, bcs) || !myEventSelectionsRun2(collision2, bcs))
        continue;

      if (additionalQAeventPlots) {
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision1.centRun2V0M());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision1.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), tracks1.size());
      }

      //            <IsMC, IsMix, IsRot, IsRun2>
      fillHistograms<true, true, false, true>(collision1, tracks1, tracks2);
    }
  }
  PROCESS_SWITCH(K892analysispbpb, processMixedEventMCRun2, "Process Mixed event MC Run2", false);

  void processEvtLossSigLossMC(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMCrec>& recCollisions)
  {

    // Event loss estimation
    auto impactPar = mcCollision.impactParameter();
    histos.fill(HIST("QAevent/hImpactParameterGen"), impactPar);

    bool isSel = false;
    auto centrality = -999.;
    if (recCollisions.size() > 0) {
      auto numcontributors = -999;
      for (const auto& RecCollision : recCollisions) {
        if (!myEventSelections(RecCollision))
          continue;

        if (RecCollision.numContrib() <= numcontributors)
          continue;
        else
          numcontributors = RecCollision.numContrib();

        centrality = RecCollision.centFT0C();
        isSel = true;
      }
    }

    if (isSel) {
      histos.fill(HIST("QAevent/hImpactParameterRec"), impactPar);
      histos.fill(HIST("QAevent/hImpactParvsCentrRec"), centrality, impactPar);
    }

    // Generated MC
    for (const auto& mcPart : mcParticles) {
      if (std::abs(mcPart.y()) >= 0.5 || std::abs(mcPart.pdgCode()) != 313)
        continue;

      // signal loss estimation
      if (mcPart.pdgCode() > 0) // no cuts, purely generated
        histos.fill(HIST("QAevent/k892genBeforeEvtSel"), mcPart.pt(), impactPar);
      else
        histos.fill(HIST("QAevent/k892genBeforeEvtSelAnti"), mcPart.pt(), impactPar);

      if (isSel) {
        // signal loss estimation
        if (mcPart.pdgCode() > 0) // no cuts, purely generated
          histos.fill(HIST("QAevent/k892genAfterEvtSel"), mcPart.pt(), impactPar);
        else
          histos.fill(HIST("QAevent/k892genAfterEvtSelAnti"), mcPart.pt(), impactPar);
      }
    } // end loop on gen particles
  }
  PROCESS_SWITCH(K892analysispbpb, processEvtLossSigLossMC, "Process Signal Loss, Event Loss", false);

  void processMC(aod::McCollisions::iterator const& /*mcCollision*/, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMCrec>& recCollisions, TrackCandidatesMCrec const& RecTracks)
  {

    histos.fill(HIST("QAevent/hMCrecCollSels"), 0);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 1);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 2);
      return;
    }
    for (const auto& RecCollision : recCollisions) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 3);

      if (!myEventSelections(RecCollision))
        continue;

      histos.fill(HIST("QAevent/hMCrecCollSels"), 8);
      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("QAevent/hMultiplicityPercentMC"), centrality);

      auto tracks = RecTracks.sliceByCached(aod::track::collisionId, RecCollision.globalIndex(), cache);

      //            <IsMC, IsMix, IsRot, IsRun2>
      fillHistograms<true, false, false, false>(RecCollision, tracks, tracks);

      // Generated MC
      for (const auto& mcPart : mcParticles) {
        if (std::abs(mcPart.y()) >= 0.5 || std::abs(mcPart.pdgCode()) != 313)
          continue;

        auto kDaughters = mcPart.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }

        TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

        auto daughtp = false;
        auto daughtk = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary())
            break;

          if (std::abs(kCurrentDaughter.pdgCode()) == 211) {
            daughtp = true;
            lDecayDaughter1.SetXYZM(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          } else if (std::abs(kCurrentDaughter.pdgCode()) == 321) {
            daughtk = true;
            lDecayDaughter2.SetXYZM(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }

        if (!daughtp || !daughtk)
          continue;

        lResonance = lDecayDaughter1 + lDecayDaughter2;

        histos.fill(HIST("QAMCTrue/hGlobalIndexMotherGen"), mcPart.globalIndex());

        if (mcPart.pdgCode() > 0) { // no cuts, purely generated
          histos.fill(HIST("k892GenInvmass"), lResonance.M());
          histos.fill(HIST("h3k892GenInvmass"), centrality, lResonance.Pt(), lResonance.M());
          histos.fill(HIST("k892Gen"), 3, mcPart.pt(), centrality);
        } else {
          histos.fill(HIST("k892GenInvmassAnti"), lResonance.M());
          histos.fill(HIST("h3k892GenInvmassAnti"), centrality, lResonance.Pt(), lResonance.M());
          histos.fill(HIST("k892GenAnti"), 3, mcPart.pt(), centrality);
        }

      } // end loop on gen particles

    } // end loop on rec collisions
  }
  PROCESS_SWITCH(K892analysispbpb, processMC, "Process Monte Carlo", false);

  void processMCRun2(aod::McCollisions::iterator const& /*mcCollision*/, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMCrecRun2>& recCollisions, TrackCandidatesMCrec const& RecTracks, BCsWithRun2Info const& bcs)
  {
    histos.fill(HIST("QAevent/hMCrecCollSels"), 0);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 1);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 2);
      return;
    }
    for (const auto& RecCollision : recCollisions) {
      histos.fill(HIST("QAevent/hMCrecCollSels"), 3);

      if (!myEventSelectionsRun2(RecCollision, bcs))
        continue;
      histos.fill(HIST("QAevent/hMCrecCollSels"), 8);

      auto centrality = RecCollision.centRun2V0M();
      histos.fill(HIST("QAevent/hMultiplicityPercentMC"), centrality);
      auto tracks = RecTracks.sliceByCached(aod::track::collisionId, RecCollision.globalIndex(), cache);

      //            <IsMC, IsMix, IsRot, IsRun2>
      fillHistograms<true, false, false, true>(RecCollision, tracks, tracks);

      // Generated MC
      for (const auto& mcPart : mcParticles) {
        if (std::abs(mcPart.y()) >= 0.5 || std::abs(mcPart.pdgCode()) != 313)
          continue;

        auto kDaughters = mcPart.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }

        TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

        auto daughtp = false;
        auto daughtk = false;
        for (const auto& kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary())
            break;

          if (std::abs(kCurrentDaughter.pdgCode()) == 211) {
            daughtp = true;
            lDecayDaughter1.SetXYZM(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          } else if (std::abs(kCurrentDaughter.pdgCode()) == 321) {
            daughtk = true;
            lDecayDaughter2.SetXYZM(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }

        if (!daughtp || !daughtk)
          continue;

        lResonance = lDecayDaughter1 + lDecayDaughter2;

        histos.fill(HIST("QAMCTrue/hGlobalIndexMotherGen"), mcPart.globalIndex());

        if (mcPart.pdgCode() > 0) { // no cuts, purely generated
          histos.fill(HIST("k892GenInvmass"), lResonance.M());
          histos.fill(HIST("h3k892GenInvmass"), centrality, lResonance.Pt(), lResonance.M());
          histos.fill(HIST("k892Gen"), 3, mcPart.pt(), centrality);
        } else {
          histos.fill(HIST("k892GenInvmassAnti"), lResonance.M());
          histos.fill(HIST("h3k892GenInvmassAnti"), centrality, lResonance.Pt(), lResonance.M());
          histos.fill(HIST("k892GenAnti"), 3, mcPart.pt(), centrality);
        }

      } // end loop on gen particles

    } // end loop on rec collisions
  }
  PROCESS_SWITCH(K892analysispbpb, processMCRun2, "Process Monte Carlo Run2", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K892analysispbpb>(cfgc)};
}
