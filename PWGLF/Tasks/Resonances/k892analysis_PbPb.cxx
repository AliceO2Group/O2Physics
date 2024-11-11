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
/// \author Marta Urioni <marta.urioni@cern.ch>

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
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
struct k892analysis_PbPb {
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

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "Additional evsel3"};

  // presel
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};

  // Track selections
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.0f, "TPC Crossed Rows to Findable Clusters"};

  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor

  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"}; // TPC
  Configurable<float> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"}; // TOF
  Configurable<float> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<float> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"}; // TOF
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};            // By pass TOF PID selection

  Configurable<bool> TofandTpcPID{"TOFandTPCPID", false, "apply both TOF and TPC PID"};

  Configurable<bool> tpclowpt{"tpclowpt", true, "apply TPC at low pt"};
  Configurable<bool> tofhighpt{"tofhighpt", false, "apply TOF at high pt"};

  // event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f}, "Mixing bins - z-vertex"};

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
  Configurable<bool> genacceptancecut{"genacceptancecut", false, "Acceptance cut on generated MC particles"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{z}} (cm)"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};

    if (doprocessSameEvent || doprocessMixedEvent) {
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
    if (doprocessMixedEvent) {
      histos.add("k892invmassME", "Invariant mass of K(892)0 mixed event", kTH1F, {invMassAxis});
      if (additionalMEPlots) {
        histos.add("k892invmassME_DS", "Invariant mass of K(892)0 mixed event DS", kTH1F, {invMassAxis});
        histos.add("k892invmassME_DSAnti", "Invariant mass of K(892)0 mixed event DSAnti", kTH1F, {invMassAxis});
      }
    }

    if (additionalQAplots) {
      // TPC ncluster distirbutions
      histos.add("TPCncluster/TPCnclusterpi", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("TPCncluster/TPCnclusterka", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      histos.add("TPCncluster/TPCnclusterPhipi", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
      histos.add("TPCncluster/TPCnclusterPhika", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
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

    // 3d histogram
    histos.add("h3k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassDSAnti", "Invariant mass of Anti-K(892)0 differnt sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLS", "Invariant mass of K(892)0 same sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3k892invmassLSAnti", "Invariant mass of Anti-K(892)0 same sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    if (doprocessMixedEvent) {
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

    if (doprocessMC) {
      histos.add("hMCrecCollSels", "MC Event statistics", HistType::kTH1F, {{10, 0.0f, 10.0f}});
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
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (track.itsNCls() < cfgITScluster)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
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

    return true;
  }

  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {

    if (TofandTpcPID) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon && candidate.hasTPC() && std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) { // tof and tpc cut
        return true;
      }

    } else {

      if (candidate.hasTPC() && std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) { // tpc cut, tof when available

        if (cByPassTOF) // skip tof selection
          return true;

        if (candidate.hasTOF()) {
          if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
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

    if (TofandTpcPID) {

      if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion && candidate.hasTPC() && std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) { // tof and tpc cut
        return true;
      }

    } else {

      if (candidate.hasTPC() && std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) { // tpc cut, tof when available

        if (cByPassTOF) // skip tof selection
          return true;

        if (candidate.hasTOF()) {
          if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
            return true;
          }
        } else {
          return true;
        }
      }
    }

    return false;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {

    auto multiplicity = collision.centFT0C();

    auto oldindex = -999;
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {

      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (additionalQAeventPlots) {
        if constexpr (!IsMC) {
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
        if (tpclowpt) {
          if (trk1ptPi > cMaxPtTPC || trk2ptKa > cMaxPtTPC)
            continue;
        } else if (tofhighpt) {
          if (trk1ptPi < cMinPtTOF || trk2ptKa < cMinPtTOF)
            continue;
        }
      }

      if (additionalQAplots) {
        // TPCncluster distributions
        histos.fill(HIST("TPCncluster/TPCnclusterpi"), trk1.tpcNClsFound());
        histos.fill(HIST("TPCncluster/TPCnclusterka"), trk2.tpcNClsFound());
        histos.fill(HIST("TPCncluster/TPCnclusterPhipi"), trk1.tpcNClsFound(), trk1.phi());
        histos.fill(HIST("TPCncluster/TPCnclusterPhika"), trk2.tpcNClsFound(), trk2.phi());
      }

      if constexpr (!IsMix) {
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
      } else if (additionalMEPlots) {
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

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (abs(lResonance.Rapidity()) >= 0.5)
        continue;
      if (cfgCutsOnMother) {
        if (lResonance.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonance.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      int track1Sign = trk1.sign();
      int track2Sign = trk2.sign();
      //// Un-like sign pair only

      if (track1Sign * track2Sign < 0) {
        if constexpr (!IsMix) {
          if (track1Sign < 0) {
            histos.fill(HIST("k892invmassDS"), lResonance.M());
            histos.fill(HIST("h3k892invmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (track1Sign > 0) {
            histos.fill(HIST("k892invmassDSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassDSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        } else {
          histos.fill(HIST("k892invmassME"), lResonance.M());
          histos.fill(HIST("h3k892invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (additionalMEPlots) {
            if (track1Sign < 0) {
              histos.fill(HIST("k892invmassME_DS"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (track1Sign > 0) {
              histos.fill(HIST("k892invmassME_DSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          }
        }

        // MC
        if constexpr (IsMC) {

          if (!trk1.has_mcParticle() || !trk2.has_mcParticle())
            continue;

          const auto mctrack1 = trk1.mcParticle();
          const auto mctrack2 = trk2.mcParticle();
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());

          if (!mctrack1.isPhysicalPrimary() || !mctrack2.isPhysicalPrimary())
            continue;

          if (track1PDG != 211 || track2PDG != 321)
            continue;

          bool ismotherok = false;
          int pdgcodeMother = -999;
          for (auto& mothertrack1 : mctrack1.template mothers_as<aod::McParticles>()) {
            for (auto& mothertrack2 : mctrack2.template mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode())
                continue;
              if (mothertrack1.globalIndex() != mothertrack2.globalIndex())
                continue;
              if (TMath::Abs(mothertrack1.pdgCode()) != 313)
                continue;
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("h1k892Recsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              pdgcodeMother = mothertrack1.pdgCode();
              ismotherok = true;
            }
          }

          if (!ismotherok)
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
          if (pdgcodeMother > 0) {
            histos.fill(HIST("k892Rec"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("k892Recinvmass"), lResonance.M());
            histos.fill(HIST("h3Reck892invmass"), multiplicity, lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("k892RecAnti"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("k892RecinvmassAnti"), lResonance.M());
            histos.fill(HIST("h3Reck892invmassAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
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
      } // end on DS or LS if tenses
    } // end of loop on tracks combinations
  } // ennd on fill histograms

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  // partitions tpc low pt
  Partition<TrackCandidates> negPitpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) < cMaxPtTPC);
  Partition<TrackCandidates> posKatpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) < cMaxPtTPC);

  Partition<TrackCandidates> posPitpc = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) < cMaxPtTPC);
  Partition<TrackCandidates> negKatpc = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) < cMaxPtTPC);

  // tpc & tof, high pt
  Partition<TrackCandidates> negPitof = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> posKatof = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  Partition<TrackCandidates> posPitof = (aod::track::signed1Pt > static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::track::pt) > cMinPtTOF);
  Partition<TrackCandidates> negKatof = (aod::track::signed1Pt < static_cast<float>(0)) && (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::track::pt) > cMinPtTOF);

  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    if (additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (additionalEvSel3 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return;
    }
    //    int occupancy = collision.trackOccupancyInTimeRange();

    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1);
    histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
    histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), centrality);

    if (additionalQAeventPlots) {
      histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
      histos.fill(HIST("TestME/hnTrksSameE"), tracks.size());
    }

    if (tpclowpt) {
      //+-
      auto candPosPitpc = posPitpc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto candNegKatpc = negKatpc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      fillHistograms<false, false>(collision, candPosPitpc, candNegKatpc);

      //-+
      auto candNegPitpc = negPitpc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto candPosKatpc = posKatpc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      fillHistograms<false, false>(collision, candNegPitpc, candPosKatpc);

    } else if (tofhighpt) {
      //+-
      auto candPosPitof = posPitof->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto candNegKatof = negKatof->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      fillHistograms<false, false>(collision, candPosPitof, candNegKatof);

      //-+
      auto candNegPitof = negPitof->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto candPosKatof = posKatof->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      fillHistograms<false, false>(collision, candNegPitof, candPosKatof);
    }
  }
  PROCESS_SWITCH(k892analysis_PbPb, processSameEvent, "Process Same event TOF High Pt", false);

  using BinningTypeVtxCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVtxCent colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVtxCent> pairs{colBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!collision1.sel8() || !collision2.sel8()) {
        return;
      }
      auto centrality = collision1.centFT0C();

      if (timFrameEvsel && (!collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        return;
      }
      if (additionalEvSel2 && (!collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        return;
      }
      if (additionalEvSel3 && (!collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        return;
      }

      if (additionalQAeventPlots) {
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision1.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), centrality);
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision1.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), tracks1.size());
      }

      if (tpclowpt) {

        //+-
        auto candPosPitpc = posPitpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
        auto candNegKatpc = negKatpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

        fillHistograms<false, true>(collision1, candPosPitpc, candNegKatpc);

        //-+
        auto candNegPitpc = negPitpc->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
        auto candPosKatpc = posKatpc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

        fillHistograms<false, true>(collision1, candNegPitpc, candPosKatpc);

      } else if (tofhighpt) {

        //+-
        auto candPosPitof = posPitof->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
        auto candNegKatof = negKatof->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

        fillHistograms<false, true>(collision1, candPosPitof, candNegKatof);

        //-+
        auto candNegPitof = negPitof->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
        auto candPosKatof = posKatof->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

        fillHistograms<false, true>(collision1, candNegPitof, candPosKatof);
      }
    }
  }
  PROCESS_SWITCH(k892analysis_PbPb, processMixedEvent, "Process Mixed event TPC low pt", true);

  // MC

  using EventCandidatesMCrec = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using TrackCandidatesMCrec = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels>>;

  void processMC(aod::McCollisions::iterator const& /*mcCollision*/, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMCrec>& recCollisions, TrackCandidatesMCrec const& RecTracks)
  {
    histos.fill(HIST("hMCrecCollSels"), 0);
    if (recCollisions.size() == 0) {
      histos.fill(HIST("hMCrecCollSels"), 1);
      return;
    }
    if (recCollisions.size() > 1) {
      histos.fill(HIST("hMCrecCollSels"), 2);
      return;
    }
    for (auto& RecCollision : recCollisions) {
      histos.fill(HIST("hMCrecCollSels"), 3);
      if (!RecCollision.sel8()) {
        continue;
      }
      histos.fill(HIST("hMCrecCollSels"), 4);
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        continue;
      }
      histos.fill(HIST("hMCrecCollSels"), 5);
      if (timFrameEvsel && (!RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      histos.fill(HIST("hMCrecCollSels"), 6);
      if (additionalEvSel2 && (!RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup) || !RecCollision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      histos.fill(HIST("hMCrecCollSels"), 7);
      if (additionalEvSel3 && (!RecCollision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
        continue;
      }
      histos.fill(HIST("hMCrecCollSels"), 8);

      auto centrality = RecCollision.centFT0C();
      histos.fill(HIST("QAevent/hMultiplicityPercentMC"), centrality);
      auto tracks = RecTracks.sliceByCached(aod::track::collisionId, RecCollision.globalIndex(), cache);
      fillHistograms<true, false>(RecCollision, tracks, tracks);

      // Generated MC
      for (auto& mcPart : mcParticles) {
        if (abs(mcPart.y()) >= 0.5 || abs(mcPart.pdgCode()) != 313)
          continue;

        auto kDaughters = mcPart.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }

        TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

        auto daughtp = false;
        auto daughtk = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary())
            break;
          if (genacceptancecut && (kCurrentDaughter.pt() < cfgCutPT || TMath::Abs(kCurrentDaughter.eta()) > cfgCutEta))
            break;

          if (abs(kCurrentDaughter.pdgCode()) == 211) {
            daughtp = true;
            lDecayDaughter1.SetXYZM(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massPi);
          } else if (abs(kCurrentDaughter.pdgCode()) == 321) {
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
  PROCESS_SWITCH(k892analysis_PbPb, processMC, "Process Monte Carlo", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<k892analysis_PbPb>(cfgc, TaskName{"k892analysis_PbPb"})};
}
