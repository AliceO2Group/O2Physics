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

/// \file kstar892analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// adaped from k892analysis.cxx by Bong-Hwi Lim <bong-hwi.lim@cern.ch>, Sawan Sawan <sawan.sawan@cern.ch>
/// \author Marta Urioni <marta.urioni@cern.ch>

#include <TLorentzVector.h>
#include "TF1.h"
#include <vector>

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct kstar892analysis {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ///// Configurables
  /// Histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};

  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};
  Configurable<int> cPIDBins{"cPIDBins", 65, "PID binning"};
  Configurable<float> cPIDQALimit{"cPIDQALimit", 6.5, "PID QA limit"};
  Configurable<int> cDCABins{"cDCABins", 150, "DCA binning"};
  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  /// Pre-selection cuts
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  /// DCA Selections
  // DCAr to PV
  Configurable<float> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<float> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<float> cfgHighPt{"cfgHighPt", 999.0, "Minimum Pt to apply tighter PID selections on Mixed Event"};
  // Pion
  Configurable<float> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};                                            // TPC
  Configurable<float> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                                            // TOF
  Configurable<float> cMaxTPCnSigmaPionHighPtME{"cMaxTPCnSigmaPionHighPtME", 2.0, "TPC nSigma cut for Pion at high pt for Mixed Event"}; // TPC
  Configurable<float> cMaxTOFnSigmaPionHighPtME{"cMaxTOFnSigmaPionHighPtME", 2.0, "TOF nSigma cut for Pion at high pt for Mixed Event"}; // TOF
  Configurable<float> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};                              // Combined
  // Kaon
  Configurable<float> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};                                            // TPC
  Configurable<float> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};                                            // TOF
  Configurable<float> cMaxTPCnSigmaKaonHighPtME{"cMaxTPCnSigmaKaonHighPtME", 2.0, "TPC nSigma cut for Kaon at pt>3GeV for Mixed Event"}; // TPC
  Configurable<float> cMaxTOFnSigmaKaonHighPtME{"cMaxTOFnSigmaKaonHighPtME", 2.0, "TOF nSigma cut for Kaon at pt>3GeV for Mixed Event"}; // TOF
  Configurable<float> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"};                              // Combined

  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};      // By pass TOF PID selection
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"}; // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAeventPlots{"additionalQAeventPlots", false, "Additional QA event plots"};
  Configurable<bool> additionalMEPlots{"additionalMEPlots", false, "Additional Mixed event plots"};

  Configurable<bool> tof_at_high_pt{"tof_at_high_pt", false, "Use TOF at high pT"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};

  // MC Event selection
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  // cuts on mother
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", false, "Enamble additional cuts on mother"};
  Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
  Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum Minv of mother cut"};

  // configurables for partitions
  Configurable<float> cMaxPtTPC{"cMaxPtTPC", 1.2, "maximum pt to apply TPC PID and TOF if available"};
  Configurable<float> cMinPtTOF{"cMinPtTOF", 0.8, "minimum pt to require TOF PID in addition to TPC"};

  void init(o2::framework::InitContext&)
  {

    // LOG(info)<< "\n cfgTPCcluster ============>"<< static_cast<uint8_t>(cfgTPCcluster);

    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec ImpactParAxis = {binsImpactPar, "Impact Parameter"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{z}} (cm)"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};

    if (additionalQAeventPlots) {
      // Test on Mixed event
      histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
      histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});

      // event histograms
      histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

      histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
    }

    // Mass QA (quick check)
    histos.add("k892invmassDS", "Invariant mass of K(892)0 different sign", kTH1F, {invMassAxis});
    histos.add("k892invmassDSAnti", "Invariant mass of Anti-K(892)0 different sign", kTH1F, {invMassAxis});
    histos.add("k892invmassLS", "Invariant mass of K(892)0 like sign", kTH1F, {invMassAxis});
    histos.add("k892invmassLSAnti", "Invariant mass of Anti-K(892)0 like sign", kTH1F, {invMassAxis});
    if (doprocessMELight || doprocessMELightWithTof || doprocessMELightTPCLowPt || doprocessMELightTOFHighPt) {
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
    if (doprocessMELight || doprocessMELightWithTof || doprocessMELightTPCLowPt || doprocessMELightTOFHighPt) {
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

    if (doprocessMCLight || doprocessMCLightWithTof || doprocessMCLightTPCLowPt || doprocessMCLightTOFHighPt) {
      // MC QA
      histos.add("QAMCTrue/ImpactParameter", "Impact parameter of event", HistType::kTH1F, {ImpactParAxis});
      histos.add("QAMCTrue/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMCTrue/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMCTrue/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAMCTrue/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("k892Recinvmass", "Inv mass distribution of Reconstructed MC K(892)", kTH1F, {invMassAxis});
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

      histos.add("ImpactParPlots/h3Reck892invmass", "Invariant mass of Reconstructed MC K(892)0", kTH3F, {ImpactParAxis, ptAxis, invMassAxis});
      histos.add("ImpactParPlots/h3Reck892invmassAnti", "Invariant mass of Reconstructed MC Anti-K(892)0", kTH3F, {ImpactParAxis, ptAxis, invMassAxis});
      histos.add("ImpactParPlots/k892Gen", "pT distribution of True MC K(892)0", kTH3F, {mcLabelAxis, ptAxis, ImpactParAxis});
      histos.add("ImpactParPlots/k892GenAnti", "pT distribution of True MC Anti-K(892)0", kTH3F, {mcLabelAxis, ptAxis, ImpactParAxis});
      histos.add("ImpactParPlots/k892Rec", "pT distribution of Reconstructed MC K(892)0", kTH2F, {ptAxis, ImpactParAxis});
      histos.add("ImpactParPlots/k892RecAnti", "pT distribution of Reconstructed MC Anti-K(892)0", kTH2F, {ptAxis, ImpactParAxis});
      histos.add("ImpactParPlots/h3k892GenInvmass", "Invariant mass of generated K(892)0", kTH3F, {ImpactParAxis, ptAxis, invMassAxis});
      histos.add("ImpactParPlots/h3k892GenInvmassAnti", "Invariant mass of generated Anti-K(892)0", kTH3F, {ImpactParAxis, ptAxis, invMassAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  ///////////////////////////////////////////////////////

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Filters
  Filter acceptanceFilter = nabs(aod::resodaughter::pt) >= cMinPtcut;
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) <= cMaxDCArToPVcut) && (nabs(aod::track::dcaZ) <= cMaxDCAzToPVcut);
  Filter primarytrackFilter = aod::resodaughter::isPVContributor && aod::resodaughter::isPrimaryTrack && aod::resodaughter::isGlobalTrackWoDCA;

  // partitions for data
  Partition<soa::Filtered<aod::ResoTracks>> resoKaWithTof = (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (aod::resodaughter::hasTOF == true);
  Partition<soa::Filtered<aod::ResoTracks>> resoPiWithTof = (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (aod::resodaughter::hasTOF == true);
  Partition<soa::Filtered<aod::ResoTracks>> resoKa = (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon);
  Partition<soa::Filtered<aod::ResoTracks>> resoPi = (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion);
  Partition<soa::Filtered<aod::ResoTracks>> resoKaTPClowPt = (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::resodaughter::pt) < cMaxPtTPC);
  Partition<soa::Filtered<aod::ResoTracks>> resoPiTPClowPt = (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::resodaughter::pt) < cMaxPtTPC);
  Partition<soa::Filtered<aod::ResoTracks>> resoKaTOFhighPt = (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (aod::resodaughter::hasTOF == true) && (nabs(aod::resodaughter::pt) >= cMinPtTOF);
  Partition<soa::Filtered<aod::ResoTracks>> resoPiTOFhighPt = (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (aod::resodaughter::hasTOF == true) && (nabs(aod::resodaughter::pt) >= cMinPtTOF);

  // Partitions for mc
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecKaWithTof = (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (aod::resodaughter::hasTOF == true);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecPiWithTof = (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (aod::resodaughter::hasTOF == true);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecKa = (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecPi = (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecKaTPClowPt = (nabs(aod::pidtpc::tpcNSigmaKa) <= cMaxTPCnSigmaKaon) && (nabs(aod::resodaughter::pt) < cMaxPtTPC);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecPiTPClowPt = (nabs(aod::pidtpc::tpcNSigmaPi) <= cMaxTPCnSigmaPion) && (nabs(aod::resodaughter::pt) < cMaxPtTPC);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecKaTOFhighPt = (nabs(aod::pidtof::tofNSigmaKa) <= cMaxTOFnSigmaKaon) && (aod::resodaughter::hasTOF == true) && (nabs(aod::resodaughter::pt) >= cMinPtTOF);
  Partition<soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>>> resoMCrecPiTOFhighPt = (nabs(aod::pidtof::tofNSigmaPi) <= cMaxTOFnSigmaPion) && (aod::resodaughter::hasTOF == true) && (nabs(aod::resodaughter::pt) >= cMinPtTOF);

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  //******************************************

  float massKa = MassKaonCharged;
  float massPi = MassPionCharged;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;

    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;

    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;

    return true;
  }
  // PID selection tools

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    if (tof_at_high_pt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
        tpcPIDPassed = true;
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
          tofPIDPassed = true;
        }
      } else {
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    if (tof_at_high_pt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
        tpcPIDPassed = true;
      }
      if (cByPassTOF && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedKaon > 0) && (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa() < nsigmaCutCombinedKaon * nsigmaCutCombinedKaon)) {
          tofPIDPassed = true;
        }
      } else {
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPionME(const T& candidate)
  {
    float NsigmaTOF;
    float NsigmaTPC;

    if (std::abs(candidate.pt()) < cfgHighPt) {
      NsigmaTOF = cMaxTOFnSigmaPion;
      NsigmaTPC = cMaxTPCnSigmaPion;
    } else {
      NsigmaTOF = cMaxTOFnSigmaPionHighPtME;
      NsigmaTPC = cMaxTPCnSigmaPionHighPtME;
    }

    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < NsigmaTPC) {
      tpcPIDPassed = true;
    }
    if (cByPassTOF && tpcPIDPassed) {
      return true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < NsigmaTOF) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPIDKaonME(const T& candidate)
  {
    float NsigmaTOF;
    float NsigmaTPC;

    if (std::abs(candidate.pt()) < cfgHighPt) {
      NsigmaTOF = cMaxTOFnSigmaKaon;
      NsigmaTPC = cMaxTPCnSigmaKaon;
    } else {
      NsigmaTOF = cMaxTOFnSigmaKaonHighPtME;
      NsigmaTPC = cMaxTPCnSigmaKaonHighPtME;
    }

    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < NsigmaTPC) {
      tpcPIDPassed = true;
    }
    if (cByPassTOF && tpcPIDPassed) {
      return true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < NsigmaTOF) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    // auto multNTracksPV = collision.multNTracksPV();
    float impactpar = -999.0;
    if constexpr (IsMC) {
      impactpar = collision.impactParameter();
    }
    auto multiplicity = collision.cent();

    if (additionalQAeventPlots) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {

      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      if (additionalQAeventPlots) {
        if constexpr (!IsMix) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
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

      if constexpr (!IsMix) {
        if (!selectionPIDPion(trk1) || !selectionPIDKaon(trk2))
          continue;
      } else {
        if (!selectionPIDPionME(trk1) || !selectionPIDKaonME(trk2))
          continue;
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

      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          if (trk1.sign() < 0) {
            histos.fill(HIST("k892invmassDS"), lResonance.M());
            histos.fill(HIST("h3k892invmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (trk1.sign() > 0) {
            histos.fill(HIST("k892invmassDSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassDSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        } else {
          histos.fill(HIST("k892invmassME"), lResonance.M());
          histos.fill(HIST("h3k892invmassME"), multiplicity, lResonance.Pt(), lResonance.M());
          if (additionalMEPlots) {
            if (trk1.sign() < 0) {
              histos.fill(HIST("k892invmassME_DS"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M());
            } else if (trk1.sign() > 0) {
              histos.fill(HIST("k892invmassME_DSAnti"), lResonance.M());
              histos.fill(HIST("h3k892invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            }
          }
        }

        // MC
        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != 211 || abs(trk2.pdgCode()) != 321)
            continue;
          if (trk1.motherId() != trk2.motherId()) // Same mother
            continue;
          if (abs(trk1.motherPDG()) != 313)
            continue;

          // Track selection check.
          histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), trk1.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAxy_ka"), trk2.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAz_pi"), trk1.dcaZ());
          histos.fill(HIST("QAMCTrue/trkDCAz_ka"), trk2.dcaZ());
          histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
          }
          histos.fill(HIST("QAMCTrue/TPC_Nsigmaka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTPC);
          if (isTrk2hasTOF) {
            histos.fill(HIST("QAMCTrue/TOF_Nsigma_ka_all"), multiplicity, trk2ptKa, trk2NSigmaKaTOF);
          }

          // MC histograms
          if (trk1.motherPDG() > 0) {
            histos.fill(HIST("k892Rec"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("ImpactParPlots/k892Rec"), lResonance.Pt(), impactpar);
            histos.fill(HIST("k892Recinvmass"), lResonance.M());
            histos.fill(HIST("h3Reck892invmass"), multiplicity, lResonance.Pt(), lResonance.M());
            histos.fill(HIST("ImpactParPlots/h3Reck892invmass"), impactpar, lResonance.Pt(), lResonance.M());
          } else {
            histos.fill(HIST("k892RecAnti"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("ImpactParPlots/k892RecAnti"), lResonance.Pt(), impactpar);
            histos.fill(HIST("k892Recinvmass"), lResonance.M());
            histos.fill(HIST("h3Reck892invmassAnti"), multiplicity, lResonance.Pt(), lResonance.M());
            histos.fill(HIST("ImpactParPlots/h3Reck892invmassAnti"), impactpar, lResonance.Pt(), lResonance.M());
          }
        }
      } else if (trk1.sign() * trk2.sign() > 0) {
        if constexpr (!IsMix) {
          if (trk1.sign() < 0) {
            histos.fill(HIST("k892invmassLS"), lResonance.M());
            histos.fill(HIST("h3k892invmassLS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else if (trk1.sign() > 0) {
            histos.fill(HIST("k892invmassLSAnti"), lResonance.M());
            histos.fill(HIST("h3k892invmassLSAnti"), multiplicity, lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  } // end of fillHistograms

  ///////////////////////////////////////////////////////////////////////////////
  /////                                                                     /////
  /////                          PROCESS FUNCTIONS                          /////
  /////                                                                     /////
  ///////////////////////////////////////////////////////////////////////////////

  void processDataLight(aod::ResoCollision& collision,
                        soa::Filtered<aod::ResoTracks> const&)
  {
    // LOG(info) << "new collision, zvtx: " << collision.posZ();
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);

    auto candPi = resoPi->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoKa->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<false, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processDataLight, "Process Event for data", false);

  void processDataLightWithTof(aod::ResoCollision& collision,
                               soa::Filtered<aod::ResoTracks> const&)
  {
    // LOG(info) << "new collision, zvtx: " << collision.posZ();
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);

    auto candPiwithtof = resoPiWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKawithtof = resoKaWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<false, false>(collision, candPiwithtof, candKawithtof);
  }
  PROCESS_SWITCH(kstar892analysis, processDataLightWithTof, "Process Event for data, tracks with TOF", false);

  void processDataTPCLowPt(aod::ResoCollision& collision,
                           soa::Filtered<aod::ResoTracks> const&)
  {
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);

    auto candPi = resoPiTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoKaTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<false, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processDataTPCLowPt, "Process Event for data, TPC at low pt", false);

  void processDataTOFHighPt(aod::ResoCollision& collision,
                            soa::Filtered<aod::ResoTracks> const&)
  {
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);

    auto candPi = resoPiTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoKaTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<false, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processDataTOFHighPt, "Process Event for data, TOF at high pt", false);

  void processMCLight(ResoMCCols::iterator const& collision,
                      soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> const&)
  {
    if (!collision.isInAfterAllCuts() || (abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;

    auto candPi = resoMCrecPi->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoMCrecKa->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<true, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processMCLight, "Process Event for MC (Reconstructed)", false);

  void processMCLightWithTof(ResoMCCols::iterator const& collision,
                             soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> const&)
  {
    if (!collision.isInAfterAllCuts() || (abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;

    auto candPiwithtof = resoMCrecPiWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKawithtof = resoMCrecKaWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<true, false>(collision, candPiwithtof, candKawithtof);
  }
  PROCESS_SWITCH(kstar892analysis, processMCLightWithTof, "Process Event for MC (Reconstructed), tracks with TOF", false);

  void processMCLightTPCLowPt(ResoMCCols::iterator const& collision,
                              soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> const&)
  {
    if (!collision.isInAfterAllCuts() || (abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;

    auto candPi = resoMCrecPiTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoMCrecKaTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<true, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processMCLightTPCLowPt, "Process Event for MC (reconstructed), TPC at low pt", false);

  void processMCLightTOFHighPt(ResoMCCols::iterator const& collision,
                               soa::Filtered<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> const&)
  {
    if (!collision.isInAfterAllCuts() || (abs(collision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;

    auto candPi = resoMCrecPiTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
    auto candKa = resoMCrecKaTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);

    fillHistograms<true, false>(collision, candPi, candKa);
  }
  PROCESS_SWITCH(kstar892analysis, processMCLightTOFHighPt, "Process Event for MC (reconstructed), TOF at high pt", false);

  void processMCTrue(ResoMCCols::iterator const& collision, aod::ResoMCParents& resoParents, aod::McParticles& mcparticles)
  {
    auto multiplicity = collision.cent();
    float impactpar = collision.impactParameter();
    histos.fill(HIST("QAMCTrue/ImpactParameter"), impactpar);

    for (auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 313 || abs(part.y()) >= 0.5)
        continue;
      bool pass1 = abs(part.daughterPDG1()) == 211 || abs(part.daughterPDG2()) == 211;
      bool pass2 = abs(part.daughterPDG1()) == 321 || abs(part.daughterPDG2()) == 321;

      if (!pass1 || !pass2)
        continue;

      TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

      std::vector<int> listDaughters{};
      std::array<int, 2> dauPdgs = {211, 321};

      auto mcparticle = part.mcParticle_as<aod::McParticles>();

      RecoDecay::getDaughters(mcparticle, &listDaughters, dauPdgs, 1);

      for (const auto& dauIdx : listDaughters) {
        auto dauPart = mcparticles.rawIteratorAt(dauIdx - mcparticles.offset());

        if (std::abs(dauPart.pdgCode()) == 211) {
          lDecayDaughter1.SetXYZM(dauPart.px(), dauPart.py(), dauPart.pz(), massPi);
        } else if (std::abs(dauPart.pdgCode()) == 321) {
          lDecayDaughter2.SetXYZM(dauPart.px(), dauPart.py(), dauPart.pz(), massKa);
        }
      }
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (part.pdgCode() > 0) { // no cuts, purely generated
        histos.fill(HIST("k892GenInvmass"), lResonance.M());
        histos.fill(HIST("h3k892GenInvmass"), multiplicity, lResonance.Pt(), lResonance.M());
        histos.fill(HIST("ImpactParPlots/h3k892GenInvmass"), impactpar, lResonance.Pt(), lResonance.M());
      } else {
        histos.fill(HIST("k892GenInvmassAnti"), lResonance.M());
        histos.fill(HIST("h3k892GenInvmassAnti"), multiplicity, lResonance.Pt(), lResonance.M());
        histos.fill(HIST("ImpactParPlots/h3k892GenInvmassAnti"), impactpar, lResonance.Pt(), lResonance.M());
      }

      if (collision.isVtxIn10()) // INEL10
      {
        if (part.pdgCode() > 0) {
          histos.fill(HIST("k892Gen"), 0, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892Gen"), 0, part.pt(), impactpar);
        } else {
          histos.fill(HIST("k892GenAnti"), 0, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892GenAnti"), 0, part.pt(), impactpar);
        }
      }
      if (collision.isVtxIn10() && collision.isInSel8()) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0) {
          histos.fill(HIST("k892Gen"), 1, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892Gen"), 1, part.pt(), impactpar);
        } else {
          histos.fill(HIST("k892GenAnti"), 1, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892GenAnti"), 1, part.pt(), impactpar);
        }
      }
      if (collision.isVtxIn10() && collision.isTriggerTVX()) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0) {
          histos.fill(HIST("k892Gen"), 2, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892Gen"), 2, part.pt(), impactpar);
        } else {
          histos.fill(HIST("k892GenAnti"), 2, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892GenAnti"), 2, part.pt(), impactpar);
        }
      }
      if (collision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0) {
          histos.fill(HIST("k892Gen"), 3, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892Gen"), 3, part.pt(), impactpar);
        } else {
          histos.fill(HIST("k892GenAnti"), 3, part.pt(), multiplicity);
          histos.fill(HIST("ImpactParPlots/k892GenAnti"), 3, part.pt(), impactpar);
        }
      }
    }
  }
  PROCESS_SWITCH(kstar892analysis, processMCTrue, "Process Event for MC (Generated)", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMELight(o2::aod::ResoCollisions& collisions, soa::Filtered<aod::ResoTracks> const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, soa::Filtered<aod::ResoTracks>, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);

      auto candPi = resoPi->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex(), cache);
      auto candKa = resoKa->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex(), cache);

      fillHistograms<false, true>(collision1, candPi, candKa);
    }
  }
  PROCESS_SWITCH(kstar892analysis, processMELight, "Process EventMixing light without partition", false);

  void processMELightWithTof(o2::aod::ResoCollisions& collisions, soa::Filtered<aod::ResoTracks> const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, soa::Filtered<aod::ResoTracks>, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);

      auto candPiwithtof = resoPiWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex(), cache);
      auto candKawithtof = resoKaWithTof->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex(), cache);

      fillHistograms<false, true>(collision1, candPiwithtof, candKawithtof);
    }
  }
  PROCESS_SWITCH(kstar892analysis, processMELightWithTof, "Process EventMixing light with tof", false);

  void processMELightTPCLowPt(o2::aod::ResoCollisions& collisions, soa::Filtered<aod::ResoTracks> const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, soa::Filtered<aod::ResoTracks>, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);

      auto candPi = resoPiTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex(), cache);
      auto candKa = resoKaTPClowPt->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex(), cache);

      fillHistograms<false, true>(collision1, candPi, candKa);
    }
  }
  PROCESS_SWITCH(kstar892analysis, processMELightTPCLowPt, "Process EventMixing light with TPC at low pt", false);

  void processMELightTOFHighPt(o2::aod::ResoCollisions& collisions, soa::Filtered<aod::ResoTracks> const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, soa::Filtered<aod::ResoTracks>, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);

      auto candPi = resoPiTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex(), cache);
      auto candKa = resoKaTOFhighPt->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex(), cache);

      fillHistograms<false, true>(collision1, candPi, candKa);
    }
  }
  PROCESS_SWITCH(kstar892analysis, processMELightTOFHighPt, "Process EventMixing light with TOF at high pt", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstar892analysis>(cfgc, TaskName{"lf-kstar892analysis"})};
}
