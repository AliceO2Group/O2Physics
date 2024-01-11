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

/// \file EbyEmeanPtFluc.cxx
/// \brief Calculate EbyE <pt> fluctuations with cummulant method.
///        For charged particles and identified particles.
///
/// \author Tanu Gahlaut <tanu.gahlaut@cern.ch>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

struct meanPtFluc {
  Configurable<float> ptMax{"ptMax", 2.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.15, "minimum pT"};
  Configurable<float> etaCut{"etaCut", 0.8, "Eta cut"};
  Configurable<float> rapCut{"rapCut", 0.5, "Rapidity Cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.12, "DCAxy cut"};
  Configurable<float> dcaZCut{"dcaZCut", 1.0, "DCAz cut"};
  Configurable<float> posZCut{"posZCut", 7.0, "cut for vertex Z"};
  Configurable<float> nSigCut1{"nSigCut1", 1.0, "nSigma cut (1)"};
  Configurable<float> nSigCut2{"nSigCut2", 2.0, "nSigma cut (2)"};
  Configurable<float> nSigCut3{"nSigCut3", 3.0, "nSigma cut (3)"};
  Configurable<float> nSigCut4{"nSigCut4", 4.0, "nSigma cut (4)"};
  Configurable<float> nSigCut5{"nSigCut5", 5.0, "nSigma cut (5)"};
  Configurable<float> nSigCut15{"nSigCut15", 1.5, "nSigma cut (1.5)"};
  Configurable<float> nSigCut25{"nSigCut25", 2.5, "nSigma cut (2.5)"};
  Configurable<float> piP1{"piP1", 0.65, "pion p (1)"};
  Configurable<float> piP2{"piP2", 0.70, "pion p (2)"};
  Configurable<float> piP3{"piP3", 1.40, "pion p (3)"};
  Configurable<float> piP4{"piP4", 1.70, "pion p (4)"};
  Configurable<float> kaP1{"kaP1", 0.20, "min kaon p (1)"};
  Configurable<float> kaP2{"kaP2", 0.5, "kaon p (2)"};
  Configurable<float> kaP3{"kaP3", 0.55, "kaon p (3)"};
  Configurable<float> kaP4{"kaP4", 0.60, "kaon p (4)"};
  Configurable<float> kaP5{"kaP5", 0.65, "kaon p (5)"};
  Configurable<float> kaP6{"kaP6", 1.10, "kaon p (6)"};
  Configurable<float> kaP7{"kaP7", 1.28, "kaon p (7)"};
  Configurable<float> kaP8{"kaP8", 1.50, "kaon p (8)"};
  Configurable<float> prP1{"prP1", 0.40, "min proton p (1)"};
  Configurable<float> prP2{"prP2", 0.95, "proton p (2)"};
  Configurable<float> prP3{"prP3", 1.00, "proton p (3)"};
  Configurable<float> prP4{"prP4", 1.05, "proton p (4)"};
  Configurable<float> prP5{"prP5", 1.13, "proton p (5)"};
  Configurable<float> prP6{"prP6", 1.18, "proton p (6)"};

  using MyAllTracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTOFFullPi, aod::pidTPCFullPi, aod::pidTOFFullPr, aod::pidTPCFullPr,
                                aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullEl, aod::pidTPCFullEl,
                                aod::pidTOFbeta>;
  using MyAllCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;

  Filter collisionZFilter = nabs(aod::collision::posZ) < posZCut;
  Filter collisionTrigger = o2::aod::evsel::sel8 == true;
  Filter trackDCA = nabs(aod::track::dcaXY) < dcaXYCut && nabs(aod::track::dcaZ) < dcaZCut;
  Filter trackEta = nabs(aod::track::eta) < etaCut;
  Filter trackPt = aod::track::pt > ptMin&& aod::track::pt < ptMax;
  Filter trackGobal = requireGlobalTrackInFilter();
  // Filter trackGobal = requirePrimaryTracksInFilter();
  using MyFilteredTracks = soa::Filtered<MyAllTracks>;
  using MyFilteredCollisions = soa::Filtered<MyAllCollisions>;

  HistogramRegistry hist{"hist", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext const&)
  {
    const AxisSpec axisEvents{5, 0, 5, "Counts"};
    const AxisSpec axisEta{100, -1., +1., "#eta"};
    const AxisSpec axisY{100, -1., +1., "Rapidity"};
    const AxisSpec axisPt{300, 0., 3., "p_{T} (GeV/c)"};
    const AxisSpec axisP{300, 0., 3., "p (GeV/c)"};
    const AxisSpec axisPart{500, 0., 5., " "};
    const AxisSpec axisMeanPt{100, 0., 3., "M(p_{T}) (GeV/c)"};
    const AxisSpec axisMult{100, 0, 100, "N_{ch}"};
    const AxisSpec axisMultTPC{200, 0, 800, "N_{TPC} "};
    const AxisSpec axisMultFT0M{150, 0, 15000, "N_{FT0M}"};
    const AxisSpec axisCentFT0M{50, 0, 101, "FT0M (%)"};
    const AxisSpec axisVtxZ{80, -20., 20., "V_{Z} (cm)"};
    const AxisSpec axisDCAz{100, -1.2, 1.2, "DCA_{Z} (cm)"};
    const AxisSpec axisDCAxy{100, -0.15, 0.15, "DCA_{XY} (cm)"};
    const AxisSpec axisTPCNsigma{500, -5., 5., "n #sigma_{TPC}"};
    const AxisSpec axisTOFNsigma{500, -5., 5., "n #sigma_{TOF}"};
    const AxisSpec axisTPCSignal{180, 20., 200., "#frac{dE}{dx}"};
    const AxisSpec axisTOFSignal{100, 0.2, 1.2, "TOF #beta"};
    const AxisSpec axisChi2{50, 0., 50., "Chi2"};
    const AxisSpec axisCrossedTPC{500, 0, 500, "Crossed TPC"};

    // QA checks:
    hist.add("QA/before/h_Counts", "Counts before cuts", kTH1D, {axisEvents});
    hist.add("QA/before/h_VtxZ", "V_{Z}", kTH1D, {axisVtxZ});

    hist.add("QA/before/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster", kTH1D, {axisChi2});
    hist.add("QA/before/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});

    hist.add("QA/before/Charged/h_Eta_ch", "#eta Charged Particles", kTH1D, {axisEta});
    hist.add("QA/before/Charged/h_Pt_ch", "p_{T} Charged Particles", kTH1D, {axisPt});
    hist.add("QA/before/Charged/h2_DcaZ_ch", "DCA_{Z}", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/before/Charged/h2_DcaXY_ch", "DCA_{XY}", kTH2D, {{axisPt}, {axisDCAxy}});

    hist.add("QA/before/h2_TPCSignal_b", "TPC Signal (before)", kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/before/h2_TOFSignal_b", "TOF Signal (before)", kTH2D, {{axisP}, {axisTOFSignal}});

    hist.add("QA/before/Pion/h2_TPCNsigma_pi", "n #sigma_{TPC} (Pions)",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/before/Pion/h2_TOFNsigma_pi", "n #sigma_{TOF} (Pions)",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/before/Pion/h2_TpcTofNsigma_pi", "n #sigma_{TPC} vs n #sigma_{TOF} (Pions)",
             kTH2D, {{{axisTPCNsigma}, {axisTOFNsigma}}});
    hist.add("QA/before/Kaon/h2_TPCNsigma_ka", "n #sigma_{TPC} Kaons",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/before/Kaon/h2_TOFNsigma_ka", "n #sigma_{TOF} Kaons",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/before/Kaon/h2_TpcTofNsigma_ka", "N_{TPC} igma vs n #sigma_{TOF} Kaons",
             kTH2D, {{{axisTPCNsigma}, {axisTOFNsigma}}});
    hist.add("QA/before/Proton/h2_TPCNsigma_pr", "n #sigma_{TPC} Protons",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/before/Proton/h2_TOFNsigma_pr", "n #sigma_{TOF} Protons",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/before/Proton/h2_TpcTofNsigma_pr", "n #sigma_{TPC} vs n #sigma_{TOF} Protons",
             kTH2D, {{{axisTPCNsigma}, {axisTOFNsigma}}});

    // after
    hist.add("QA/after/h_Counts", "Counts after cuts", kTH1D, {axisEvents});
    hist.add("QA/after/h_VtxZ", "V_{Z} (after)", kTH1D, {axisVtxZ});

    hist.add("QA/after/h_NTPC", "N_{TPC}", kTH1D, {axisMultTPC});
    hist.add("QA/after/h_NFT0M", "FT0M Multiplicity", kTH1D, {axisMultFT0M});
    hist.add("QA/after/h_Cent", "FT0M (%)", kTH1D, {axisCentFT0M});
    hist.add("QA/after/h2_NTPC_NFT0M", "N_{TPC} vs N_{FT0M}", kTH2D, {{axisMultFT0M}, {axisMultTPC}});
    hist.add("QA/after/p_NTPC_NFT0M", "N_{TPC} vs N_{FT0M} (Profile)", kTProfile, {{axisMultFT0M}});
    hist.add("QA/after/p_NFT0M_NTPC", "N_{FT0M} vs N_{TPC} (Profile)", kTProfile, {{axisMultTPC}});
    hist.add("QA/after/h2_NTPC_Cent", "N_{TPC} vs FT0M(%)", kTH2D, {{axisCentFT0M}, {axisMultTPC}});
    hist.add("QA/after/p_NTPC_Cent", "N_{TPC} vs FT0M(%) (Profile)", kTProfile, {{axisCentFT0M}});
    hist.add("QA/after/h2_NTPC_Nch", "N_{ch} vs N_{TPC}", kTH2D, {{axisMultTPC}, {axisMult}});

    hist.add("QA/after/h_TPCChi2perCluster", "TPC #Chi^{2}/Cluster (after)", kTH1D, {axisChi2});
    hist.add("QA/after/h_ITSChi2perCluster", "ITS #Chi^{2}/Cluster (after)", kTH1D, {axisChi2});
    hist.add("QA/after/h_crossedTPC", "Crossed TPC", kTH1D, {axisCrossedTPC});

    hist.add("QA/after/Charged/h_Mult_ch", "Multiplicity Charged Prticles", kTH1D, {axisMult});
    hist.add("QA/after/Charged/h_Eta_ch", "#eta Charged Particles (after)", kTH1D, {axisEta});
    hist.add("QA/after/Charged/h_Pt_ch", "p_{T} Charged Particles (after)", kTH1D, {axisPt});
    hist.add("QA/after/Charged/h2_DcaZ_ch", "DCA_{Z} Charged Particles (after)",
             kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/after/Charged/h2_DcaXY_ch", "DCA_{XY} Charged Particles (after)",
             kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/after/Charged/h2_Pt_Eta_ch", "p_{T} vs #eta (Charged Particles)",
             kTH2D, {{axisEta}, {axisPt}});

    hist.add("QA/after/h2_TPCSignal_a", "TPC Signal (after)", kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/h2_TOFSignal_a", "TOF Signal (after)", kTH2D, {{axisP}, {axisTOFSignal}});

    hist.add("QA/after/TPC/Pion/h_Pt_pi_TPC", "p_{T} (Pions) TPC", kTH1D, {axisPt});
    hist.add("QA/after/TPC/Kaon/h_Pt_ka_TPC", "p_{T} (Kaons) TPC", kTH1D, {axisPt});
    hist.add("QA/after/TPC/Proton/h_Pt_pr_TPC", "p_{T} (Protons) TPC ", kTH1D, {axisPt});
    hist.add("QA/after/TPC/Pion/h_rap_pi_TPC", "y (Pions) TPC ", kTH1D, {axisY});
    hist.add("QA/after/TPC/Kaon/h_rap_ka_TPC", "y (Kaons) TPC", kTH1D, {axisY});
    hist.add("QA/after/TPC/Proton/h_rap_pr_TPC", "y (Protons) TPC", kTH1D, {axisY});
    hist.add("QA/after/TPC/Pion/h2_TPCSignal_pi_b", "TPC Signal Pions",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TPC/Pion/h2_ExpTPCSignal_pi_b", "Expected TPC Signal Pions",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TPC/Kaon/h2_TPCSignal_ka_b", "TPC Signal Kaons",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TPC/Kaon/h2_ExpTPCSignal_ka_b", "Expected TPC Signal Kaons",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TPC/Proton/h2_TPCSignal_pr_b", "TPC Signal Protons",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TPC/Proton/h2_ExpTPCSignal_pr_b", "Expected TPC Signal Protons",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/TOF/Pion/h_Pt_pi_TOF", "p_{T} (Pions) TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/TOF/Kaon/h_Pt_ka_TOF", "p_{T} (Kaons) TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/TOF/Proton/h_Pt_pr_TOF", "p_{T} (Protons) TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/TOF/Pion/h_rap_pi_TOF", "y (Pions) TPC+TOF ", kTH1D, {axisY});
    hist.add("QA/after/TOF/Kaon/h_rap_ka_TOF", "y (Kaons) TPC+TOF", kTH1D, {axisY});
    hist.add("QA/after/TOF/Proton/h_rap_pr_TOF", "y (Protons) TPC+TOF", kTH1D, {axisY});
    hist.add("QA/after/TOF/Pion/h2_TOFSignal_pi_b", "TOF Signal Pions",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/TOF/Pion/h2_ExpTOFSignal_pi_b", "Expected TOF Signal Pions",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/TOF/Kaon/h2_TOFSignal_ka_b", "TOF Signal Kaons",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/TOF/Kaon/h2_ExpTOFSignal_ka_b", "Expected TOF Signal Kaons",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/TOF/Proton/h2_TOFSignal_pr_b", "TOF Signal Protons",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/TOF/Proton/h2_ExpTOFSignal_pr_b", "Expected TOF Signal Protons",
             kTH2D, {{axisP}, {axisTOFSignal}});

    hist.add("QA/after/Pion/h_Mult_pi", "Multiplicity Pion", kTH1D, {axisMult});
    hist.add("QA/after/Pion/h_Pt_pi", "p_{T} (Pions) TPC and TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/Pion/h_rap_pi", "y (Pions) TPC and TPC+TOF", kTH1D, {axisY});
    hist.add("QA/after/Pion/h2_Pt_rap_pi", "p_{T} vs y (Pions)", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/after/Pion/h2_DcaZ_pi", "DCA_{z} (Pions)", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/after/Pion/h2_DcaXY_pi", "DCA_{xy} (Pions)", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/after/Pion/h2_TPCNsigma_pi", "n #sigma_{TPC} (Pions)",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/after/Pion/h2_TOFNsigma_pi", "n #sigma_{TOF} (Pions)",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/after/Pion/h2_TpcTofNsigma_pi", "n #sigma_{TPC} vs n #sigma_{TOF} (Pions)",
             kTH2D, {{{axisTPCNsigma}, {axisTOFNsigma}}});
    hist.add("QA/after/Pion/h2_TPCSignal_pi_a", "TPC Signal Pions (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Pion/h2_TOFSignal_pi_a", "TOF Signal Pions (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/Pion/h2_ExpTPCSignal_pi_a", "Expected TPC Signal Pions (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Pion/h2_ExpTOFSignal_pi_a", "Expected TOF Signal Pions (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});

    hist.add("QA/after/Kaon/h_Mult_ka", "Multiplicity Kaon", kTH1D, {axisMult});
    hist.add("QA/after/Kaon/h_Pt_ka", "p_{T} (Kaons) TPC and TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/Kaon/h_rap_ka", "y (Kaons) TPC and TPC+TOF", kTH1D, {axisY});
    hist.add("QA/after/Kaon/h2_Pt_rap_ka", "p_{T} vs y (Kaons)", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/after/Kaon/h2_DcaZ_ka", "DCA_{z} Kaons", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/after/Kaon/h2_DcaXY_ka", "DCA_{xy} Kaons", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/after/Kaon/h2_TPCNsigma_ka", "n #sigma_{TPC} Kaons",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/after/Kaon/h2_TOFNsigma_ka", "n #sigma_{TOF} Kaons",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/after/Kaon/h2_TpcTofNsigma_ka", "n #sigma_{TPC} vs n #sigma_{TOF} Kaons",
             kTH2D, {{axisTPCNsigma}, {axisTOFNsigma}});
    hist.add("QA/after/Kaon/h2_TPCSignal_ka_a", "TPC Signal Kaons (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Kaon/h2_TOFSignal_ka_a", "TOF Signal Kaons (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/Kaon/h2_ExpTPCSignal_ka_a", "Expected TPC Signal Kaons (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Kaon/h2_ExpTOFSignal_ka_a", "Expected TOF Signal Kaons (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});

    hist.add("QA/after/Proton/h_Mult_pr", "Multiplicity Proton", kTH1D, {axisMult});
    hist.add("QA/after/Proton/h_Pt_pr", "p_{T} (Protons) TPC and TPC+TOF", kTH1D, {axisPt});
    hist.add("QA/after/Proton/h_rap_pr", "y(Protons) TPC and TPC+TOF", kTH1D, {axisY});
    hist.add("QA/after/Proton/h2_Pt_rap_pr", "p_{T} vs y (Protons)", kTH2D, {{axisY}, {axisPt}});
    hist.add("QA/after/Proton/h2_DcaZ_pr", "DCA_{z} (Protons)", kTH2D, {{axisPt}, {axisDCAz}});
    hist.add("QA/after/Proton/h2_DcaXY_pr", "DCA_{xy} (Protons)", kTH2D, {{axisPt}, {axisDCAxy}});
    hist.add("QA/after/Proton/h2_TPCNsigma_pr", "n #sigma_{TPC} (Protons)",
             kTH2D, {{axisP}, {axisTPCNsigma}});
    hist.add("QA/after/Proton/h2_TOFNsigma_pr", "n #sigma_{TOF} (Protons)",
             kTH2D, {{axisP}, {axisTOFNsigma}});
    hist.add("QA/after/Proton/h2_TpcTofNsigma_pr", "n #sigma_{TPC} vs n #sigma_{TOF} (Protons)",
             kTH2D, {{{axisTPCNsigma}, {axisTOFNsigma}}});
    hist.add("QA/after/Proton/h2_TPCSignal_pr_a", "TPC Signal Protons (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Proton/h2_TOFSignal_pr_a", "TOF Signal Protons (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});
    hist.add("QA/after/Proton/h2_ExpTPCSignal_pr_a", "Expected TPC Signal Protons (after)",
             kTH2D, {{axisP}, {axisTPCSignal}});
    hist.add("QA/after/Proton/h2_ExpTOFSignal_pr_a", "Expected TOF Signal Protons (after)",
             kTH2D, {{axisP}, {axisTOFSignal}});

    // Analysis:
    // Charged Particles
    hist.add("Analysis/Charged/h_Mult_ch", "Multiplicity of Charged Prticles", kTH1D, {axisMult});
    hist.add("Analysis/Charged/h_mean_Q1_ch", "mean p_{T} (Charged particles)", kTH1D, {axisMeanPt});
    hist.add("Analysis/Charged/p_mean_Q1_ch", "mean p_{T} (Charged particles)", kTProfile, {axisMult});
    hist.add("Analysis/Charged/h_mean_Q1_Mult_ch", "Mean p_{T} vs N_{ch} (Charged Particles)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Charged/h_twopart_Mult_ch", "Twopart vs N_{ch} (Charged Particles)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Charged/h_threepart_Mult_ch", "Threepart vs N_{ch} (Charged Particles)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Charged/h_fourpart_Mult_ch", "Fourpart vs N_{ch} (Charged Particles)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});

    // Pions
    hist.add("Analysis/Pion/h_Mult_pi", "Multiplicity Pion", kTH1D, {axisMult});
    hist.add("Analysis/Pion/h_mean_Q1_pi", "mean p_{T} Pion", kTH1D, {axisMeanPt});
    hist.add("Analysis/Pion/p_mean_Q1_pi", "mean p_{T} Pion", kTProfile, {axisMult});
    hist.add("Analysis/Pion/h_mean_Q1_Mult_pi", "mean_Q1_Mult Pion",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_twopart_Mult_pi", "twopart_Mult Pion",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_threepart_Mult_pi", "threepart_Mult Pion",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_fourpart_Mult_pi", "fourpart_Mult Pion",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});

    // Kaons
    hist.add("Analysis/Kaon/h_Mult_ka", "Multiplicity Kaon", kTH1D, {axisMult});
    hist.add("Analysis/Kaon/h_mean_Q1_ka", "mean p_{T} Kaon", kTH1D, {axisMeanPt});
    hist.add("Analysis/Kaon/p_mean_Q1_ka", "mean p_{T} Kaon", kTProfile, {axisMult});
    hist.add("Analysis/Kaon/h_mean_Q1_Mult_ka", "mean_Q1_Mult Kaon",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_twopart_Mult_ka", "twopart_Mult Kaon",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_threepart_Mult_ka", "threepart_Mult Kaon",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_fourpart_Mult_ka", "fourpart_Mult Kaon",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});

    // Protons
    hist.add("Analysis/Proton/h_Mult_pr", "Multiplicity Proton", kTH1D, {axisMult});
    hist.add("Analysis/Proton/h_mean_Q1_pr", "mean p_{T} Proton", kTH1D, {axisMeanPt});
    hist.add("Analysis/Proton/p_mean_Q1_pr", "mean p_{T} Proton", kTProfile, {axisMult});
    hist.add("Analysis/Proton/h_mean_Q1_Mult_pr", "mean_Q1_Mult Proton",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_twopart_Mult_pr", "twopart_Mult Proton",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_threepart_Mult_pr", "threepart_Mult Proton",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_fourpart_Mult_pr", "fourpart_Mult Proton",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});

    // Additional QA
    hist.add("Analysis/Charged/h_mean_Q1_Mult_ch_tof", "mean_Q1_Mult (Charged Particles) (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Charged/h_twopart_Mult_ch_tof", "twopart_Mult (Charged Particles) (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_mean_Q1_Mult_pi_tof", "mean_Q1_Mult Pion (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_mean_Q1_Mult_pi_tpc", "mean_Q1_Mult Pion (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_twopart_Mult_pi_tof", "twopart_Mult Pion (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Pion/h_twopart_Mult_pi_tpc", "twopart_Mult Pion (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_mean_Q1_Mult_pr_tof", "mean_Q1_Mult Proton (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_mean_Q1_Mult_pr_tpc", "mean_Q1_Mult Proton (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_twopart_Mult_pr_tof", "twopart_Mult Proton (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Proton/h_twopart_Mult_pr_tpc", "twopart_Mult Proton (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_mean_Q1_Mult_ka_tof", "mean_Q1_Mult Kaon (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_mean_Q1_Mult_ka_tpc", "mean_Q1_Mult Kaon (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_twopart_Mult_ka_tof", "twopart_Mult Kaon (TOF+TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
    hist.add("Analysis/Kaon/h_twopart_Mult_ka_tpc", "twopart_Mult Kaon (TPC)",
             kTHnSparseD, {{axisMultTPC}, {axisPart}, {axisMultFT0M}});
  }

  void process_QA(MyAllCollisions::iterator const& myCol, MyAllTracks const& myTracks)
  {
    for (auto& myTrack : myTracks) {
      hist.fill(HIST("QA/before/Charged/h_Eta_ch"), myTrack.eta());
      hist.fill(HIST("QA/before/Charged/h_Pt_ch"), myTrack.pt());
      hist.fill(HIST("QA/before/h_TPCChi2perCluster"), myTrack.tpcChi2NCl());
      hist.fill(HIST("QA/before/h_ITSChi2perCluster"), myTrack.itsChi2NCl());
      hist.fill(HIST("QA/before/h_crossedTPC"), myTrack.tpcNClsCrossedRows());
      hist.fill(HIST("QA/before/Charged/h2_DcaXY_ch"), myTrack.pt(), myTrack.dcaXY());
      hist.fill(HIST("QA/before/Charged/h2_DcaZ_ch"), myTrack.pt(), myTrack.dcaZ());
    }
    hist.fill(HIST("QA/before/h_VtxZ"), myCol.posZ());
    hist.fill(HIST("QA/before/h_Counts"), 2);
  }
  PROCESS_SWITCH(meanPtFluc, process_QA, "process QA", true);

  void process(MyFilteredCollisions::iterator const& col, MyFilteredTracks const& tracks)
  {
    double Cent_FT0M = 0;
    int N_Pi = 0, N_Ka = 0, N_Pr = 0;
    int Nch = 0, NTPC = 0, N_FT0M = 0;
    int N_Ka_tpc = 0, N_Pr_tpc = 0, N_Pi_tpc = 0;
    int Nch_tof = 0, N_Ka_tof = 0, N_Pr_tof = 0, N_Pi_tof = 0;
    double pt_ch = 0, Q1_ch = 0, Q2_ch = 0, Q3_ch = 0, Q4_ch = 0;
    double pt_Pi = 0, Q1_Pi = 0, Q2_Pi = 0, Q3_Pi = 0, Q4_Pi = 0;
    double pt_Pr = 0, Q1_Pr = 0, Q2_Pr = 0, Q3_Pr = 0, Q4_Pr = 0;
    double pt_Ka = 0, Q1_Ka = 0, Q2_Ka = 0, Q3_Ka = 0, Q4_Ka = 0;
    double Q1_tof = 0, Q1_Pi_tof = 0, Q1_Pr_tof = 0, Q1_Ka_tof = 0;
    double Q1_Pi_tpc = 0, Q1_Pr_tpc = 0, Q1_Ka_tpc = 0;
    double Q2_tof = 0, Q2_Pi_tof = 0, Q2_Pr_tof = 0, Q2_Ka_tof = 0;
    double Q2_Pi_tpc = 0, Q2_Pr_tpc = 0, Q2_Ka_tpc = 0;

    for (auto& track : tracks) {
      Nch++;
      pt_ch = track.pt();
      Q1_ch += pt_ch;
      Q2_ch += pt_ch * pt_ch;
      Q3_ch += pt_ch * pt_ch * pt_ch;
      Q4_ch += pt_ch * pt_ch * pt_ch * pt_ch;

      hist.fill(HIST("QA/after/Charged/h_Eta_ch"), track.eta());
      hist.fill(HIST("QA/after/Charged/h_Pt_ch"), track.pt());
      hist.fill(HIST("QA/after/Charged/h2_Pt_Eta_ch"), track.eta(), track.pt());
      hist.fill(HIST("QA/after/Charged/h2_DcaXY_ch"), track.pt(), track.dcaXY());
      hist.fill(HIST("QA/after/Charged/h2_DcaZ_ch"), track.pt(), track.dcaZ());

      hist.fill(HIST("QA/after/h_TPCChi2perCluster"), track.tpcChi2NCl());
      hist.fill(HIST("QA/after/h_ITSChi2perCluster"), track.itsChi2NCl());
      hist.fill(HIST("QA/after/h_crossedTPC"), track.tpcNClsCrossedRows());

      hist.fill(HIST("QA/before/h2_TOFSignal_b"), track.p(), track.beta());
      hist.fill(HIST("QA/before/h2_TPCSignal_b"), track.p(), track.tpcSignal());

      hist.fill(HIST("QA/before/Pion/h2_TPCNsigma_pi"), track.p(), track.tpcNSigmaPi());
      hist.fill(HIST("QA/before/Pion/h2_TOFNsigma_pi"), track.p(), track.tofNSigmaPi());
      hist.fill(HIST("QA/before/Pion/h2_TpcTofNsigma_pi"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      hist.fill(HIST("QA/before/Proton/h2_TPCNsigma_pr"), track.p(), track.tpcNSigmaPr());
      hist.fill(HIST("QA/before/Proton/h2_TOFNsigma_pr"), track.p(), track.tofNSigmaPr());
      hist.fill(HIST("QA/before/Proton/h2_TpcTofNsigma_pr"), track.tpcNSigmaPr(), track.tofNSigmaPr());
      hist.fill(HIST("QA/before/Kaon/h2_TPCNsigma_ka"), track.p(), track.tpcNSigmaKa());
      hist.fill(HIST("QA/before/Kaon/h2_TOFNsigma_ka"), track.p(), track.tofNSigmaKa());
      hist.fill(HIST("QA/before/Kaon/h2_TpcTofNsigma_ka"), track.tpcNSigmaKa(), track.tofNSigmaKa());

      // ###################################################//
      //              TPC (Without p cuts)                  //
      // ###################################################//
      if (abs(track.tpcNSigmaPi()) < nSigCut3) {
        if (abs(track.rapidity(massPi)) >= 0.5)
          continue;
        N_Pi_tpc++;
        Q1_Pi_tpc += track.pt();
        Q2_Pi_tpc += track.pt() * track.pt();
        hist.fill(HIST("QA/after/TPC/Pion/h_Pt_pi_TPC"), track.pt());
        hist.fill(HIST("QA/after/TPC/Pion/h_rap_pi_TPC"), track.rapidity(massPi));
        hist.fill(HIST("QA/after/TPC/Pion/h2_TPCSignal_pi_b"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/TPC/Pion/h2_ExpTPCSignal_pi_b"), track.p(), track.tpcExpSignalPi(track.tpcSignal()));
      }
      if (abs(track.tpcNSigmaKa()) < nSigCut3) {
        if (abs(track.rapidity(massKa)) >= 0.5)
          continue;
        N_Ka_tpc++;
        Q1_Ka_tpc += track.pt();
        Q2_Ka_tpc += track.pt() * track.pt();
        hist.fill(HIST("QA/after/TPC/Kaon/h_Pt_ka_TPC"), track.pt());
        hist.fill(HIST("QA/after/TPC/Kaon/h_rap_ka_TPC"), track.rapidity(massKa));
        hist.fill(HIST("QA/after/TPC/Kaon/h2_TPCSignal_ka_b"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/TPC/Kaon/h2_ExpTPCSignal_ka_b"), track.p(), track.tpcExpSignalKa(track.tpcSignal()));
      }
      if (abs(track.tpcNSigmaPr()) < nSigCut3) {
        if (abs(track.rapidity(massPr)) >= 0.5)
          continue;
        N_Pr_tpc++;
        Q1_Pr_tpc += track.pt();
        Q2_Pr_tpc += track.pt() * track.pt();
        hist.fill(HIST("QA/after/TPC/Proton/h_Pt_pr_TPC"), track.pt());
        hist.fill(HIST("QA/after/TPC/Proton/h_rap_pr_TPC"), track.rapidity(massPr));
        hist.fill(HIST("QA/after/TPC/Proton/h2_TPCSignal_pr_b"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/TPC/Proton/h2_ExpTPCSignal_pr_b"), track.p(), track.tpcExpSignalPr(track.tpcSignal()));
      }

      // ###################################################//
      //               TPC + TOF (Without p cuts)           //
      // ###################################################//
      if (track.hasTOF()) {
        Nch_tof++;
        Q1_tof += track.pt();
        Q2_tof += track.pt() * track.pt();

        if ((std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2)) < 6.0) {
          if (abs(track.rapidity(massPi)) >= 0.5)
            continue;
          N_Pi_tof++;
          Q1_Pi_tof += track.pt();
          Q2_Pi_tof += track.pt() * track.pt();
          hist.fill(HIST("QA/after/TOF/Pion/h_Pt_pi_TOF"), track.pt());
          hist.fill(HIST("QA/after/TOF/Pion/h_rap_pi_TOF"), track.rapidity(massPi));
          hist.fill(HIST("QA/after/TOF/Pion/h2_TOFSignal_pi_b"), track.p(), track.beta());
          hist.fill(HIST("QA/after/TOF/Pion/h2_ExpTOFSignal_pi_b"), track.p(), track.tofExpSignalPi(track.beta()));
        }
        if ((std::pow(track.tpcNSigmaKa(), 2) + std::pow(track.tofNSigmaKa(), 2)) < 6.0) {
          if (abs(track.rapidity(massKa)) >= 0.5)
            continue;
          N_Ka_tof++;
          Q1_Ka_tof += track.pt();
          Q2_Ka_tof += track.pt() * track.pt();
          hist.fill(HIST("QA/after/TOF/Kaon/h_Pt_ka_TOF"), track.pt());
          hist.fill(HIST("QA/after/TOF/Kaon/h_rap_ka_TOF"), track.rapidity(massKa));
          hist.fill(HIST("QA/after/TOF/Kaon/h2_TOFSignal_ka_b"), track.p(), track.beta());
          hist.fill(HIST("QA/after/TOF/Kaon/h2_ExpTOFSignal_ka_b"), track.p(), track.tofExpSignalKa(track.beta()));
        }
        if ((std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2)) < 6.0) {
          if (abs(track.rapidity(massPr)) >= 0.5)
            continue;
          N_Pr_tof++;
          Q1_Pr_tof += track.pt();
          Q2_Pr_tof += track.pt() * track.pt();
          hist.fill(HIST("QA/after/TOF/Proton/h_Pt_pr_TOF"), track.pt());
          hist.fill(HIST("QA/after/TOF/Proton/h_rap_pr_TOF"), track.rapidity(massPr));
          hist.fill(HIST("QA/after/TOF/Proton/h2_TOFSignal_pr_b"), track.p(), track.beta());
          hist.fill(HIST("QA/after/TOF/Proton/h2_ExpTOFSignal_pr_b"), track.p(), track.tofExpSignalPr(track.beta()));
        }
      }

      // ####################################################//
      //      TPC and TPC+TOF nSigma Cuts (with p cuts)      //
      // ####################################################//
      // For Pions:
      if ((track.hasTOF() == false &&
           ((abs(track.tpcNSigmaPi()) < nSigCut3 && track.p() <= piP1) || (abs(track.tpcNSigmaPi()) < nSigCut2 && track.p() > piP1 && track.p() <= piP2))) ||
          (track.hasTOF() && abs(track.tpcNSigmaPi()) < nSigCut4 && abs(track.tofNSigmaEl()) > nSigCut1 &&
           ((abs(track.tofNSigmaPi()) < nSigCut3 && track.p() <= piP3) || (abs(track.tofNSigmaPi()) < nSigCut25 && track.p() > piP3 && track.p() <= piP4) || (abs(track.tofNSigmaPi()) < nSigCut2 && track.p() > piP4)))) {
        if (abs(track.rapidity(massPi)) >= 0.5)
          continue;
        N_Pi++;
        pt_Pi = track.pt();
        Q1_Pi += pt_Pi;
        Q2_Pi += pt_Pi * pt_Pi;
        Q3_Pi += pt_Pi * pt_Pi * pt_Pi;
        Q4_Pi += pt_Pi * pt_Pi * pt_Pi * pt_Pi;
        hist.fill(HIST("QA/after/Pion/h_Pt_pi"), track.pt());
        hist.fill(HIST("QA/after/Pion/h_rap_pi"), track.rapidity(massPi));
        hist.fill(HIST("QA/after/Pion/h2_Pt_rap_pi"), track.rapidity(massPi), track.pt());
        hist.fill(HIST("QA/after/Pion/h2_DcaXY_pi"), track.pt(), track.dcaXY());
        hist.fill(HIST("QA/after/Pion/h2_DcaZ_pi"), track.pt(), track.dcaZ());
        hist.fill(HIST("QA/after/Pion/h2_TPCNsigma_pi"), track.p(), track.tpcNSigmaPi());
        hist.fill(HIST("QA/after/Pion/h2_TOFNsigma_pi"), track.p(), track.tofNSigmaPi());
        hist.fill(HIST("QA/after/Pion/h2_TpcTofNsigma_pi"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        hist.fill(HIST("QA/after/h2_TOFSignal_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/Pion/h2_TOFSignal_pi_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/h2_TPCSignal_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/Pion/h2_TPCSignal_pi_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/Pion/h2_ExpTOFSignal_pi_a"), track.p(), track.tofExpSignalPi(track.beta()));
        hist.fill(HIST("QA/after/Pion/h2_ExpTPCSignal_pi_a"), track.p(), track.tpcExpSignalPi(track.tpcSignal()));
      }

      // For Kaons:
      if ((track.hasTOF() == false &&
           ((abs(track.tpcNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP2) || (abs(track.tpcNSigmaKa()) < nSigCut25 && track.p() > kaP2 && track.p() <= kaP3) || (abs(track.tpcNSigmaKa()) < nSigCut2 && track.p() > kaP3 && track.p() <= kaP4) || (abs(track.tpcNSigmaKa()) < nSigCut15 && track.p() > kaP4 && track.p() <= kaP5))) ||
          (track.hasTOF() && abs(track.tpcNSigmaKa()) < nSigCut4 && abs(track.tofNSigmaEl()) > nSigCut1 &&
           ((abs(track.tofNSigmaKa()) < nSigCut3 && track.pt() > kaP1 && track.p() <= kaP6) || (abs(track.tofNSigmaKa()) < nSigCut2 && track.p() > kaP6 && track.p() <= kaP7) || (abs(track.tofNSigmaKa()) < nSigCut15 && track.p() > kaP7 && track.p() <= kaP8) || (abs(track.tofNSigmaKa()) < nSigCut1 && track.p() > kaP8)))) {
        if (abs(track.rapidity(massKa)) >= 0.5)
          continue;
        pt_Ka = track.pt();
        Q1_Ka += pt_Ka;
        Q2_Ka += pt_Ka * pt_Ka;
        Q3_Ka += pt_Ka * pt_Ka * pt_Ka;
        Q4_Ka += pt_Ka * pt_Ka * pt_Ka * pt_Ka;
        N_Ka++;
        hist.fill(HIST("QA/after/Kaon/h_Pt_ka"), track.pt());
        hist.fill(HIST("QA/after/Kaon/h_rap_ka"), track.rapidity(massKa));
        hist.fill(HIST("QA/after/Kaon/h2_Pt_rap_ka"), track.rapidity(massKa), track.pt());
        hist.fill(HIST("QA/after/Kaon/h2_DcaXY_ka"), track.pt(), track.dcaXY());
        hist.fill(HIST("QA/after/Kaon/h2_DcaZ_ka"), track.pt(), track.dcaZ());
        hist.fill(HIST("QA/after/Kaon/h2_TPCNsigma_ka"), track.p(), track.tpcNSigmaKa());
        hist.fill(HIST("QA/after/Kaon/h2_TOFNsigma_ka"), track.p(), track.tofNSigmaKa());
        hist.fill(HIST("QA/after/Kaon/h2_TpcTofNsigma_ka"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        hist.fill(HIST("QA/after/h2_TOFSignal_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/Kaon/h2_TOFSignal_ka_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/h2_TPCSignal_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/Kaon/h2_TPCSignal_ka_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/Kaon/h2_ExpTOFSignal_ka_a"), track.p(), track.tofExpSignalKa(track.beta()));
        hist.fill(HIST("QA/after/Kaon/h2_ExpTPCSignal_ka_a"), track.p(), track.tpcExpSignalKa(track.tpcSignal()));
      }

      // For Protons:
      if ((track.hasTOF() == false &&
           ((abs(track.tpcNSigmaPr()) < nSigCut3 && track.pt() > prP1 && track.p() <= prP2) || (abs(track.tpcNSigmaPr()) < nSigCut25 && track.p() > prP2 && track.p() <= prP3) || (abs(track.tpcNSigmaPr()) < nSigCut2 && track.p() > prP3 && track.p() <= prP4) || (abs(track.tpcNSigmaPr()) < nSigCut15 && track.p() > prP4 && track.p() <= prP5) || (abs(track.tpcNSigmaPr()) < nSigCut1 && track.p() > prP5 && track.p() <= prP6))) ||
          (track.hasTOF() && abs(track.tpcNSigmaPr()) < nSigCut4 && abs(track.tofNSigmaEl()) > nSigCut1 && abs(track.tofNSigmaPr()) < nSigCut3 && track.pt() > prP1)) {
        if (abs(track.rapidity(massPr)) >= 0.5)
          continue;
        pt_Pr = track.pt();
        Q1_Pr += pt_Pr;
        Q2_Pr += pt_Pr * pt_Pr;
        Q3_Pr += pt_Pr * pt_Pr * pt_Pr;
        Q4_Pr += pt_Pr * pt_Pr * pt_Pr * pt_Pr;
        N_Pr++;
        hist.fill(HIST("QA/after/Proton/h_Pt_pr"), track.pt());
        hist.fill(HIST("QA/after/Proton/h_rap_pr"), track.rapidity(massPr));
        hist.fill(HIST("QA/after/Proton/h2_DcaZ_pr"), track.pt(), track.dcaZ());
        hist.fill(HIST("QA/after/Proton/h2_DcaXY_pr"), track.pt(), track.dcaXY());
        hist.fill(HIST("QA/after/Proton/h2_Pt_rap_pr"), track.rapidity(massPr), track.pt());
        hist.fill(HIST("QA/after/Proton/h2_TPCNsigma_pr"), track.p(), track.tpcNSigmaPr());
        hist.fill(HIST("QA/after/Proton/h2_TOFNsigma_pr"), track.p(), track.tofNSigmaPr());
        hist.fill(HIST("QA/after/Proton/h2_TpcTofNsigma_pr"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        hist.fill(HIST("QA/after/Proton/h2_TPCSignal_pr_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/h2_TPCSignal_a"), track.p(), track.tpcSignal());
        hist.fill(HIST("QA/after/h2_TOFSignal_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/Proton/h2_TOFSignal_pr_a"), track.p(), track.beta());
        hist.fill(HIST("QA/after/Proton/h2_ExpTOFSignal_pr_a"), track.p(), track.tofExpSignalPr(track.beta()));
        hist.fill(HIST("QA/after/Proton/h2_ExpTPCSignal_pr_a"), track.p(), track.tpcExpSignalPr(track.tpcSignal()));
      }
    }

    NTPC = col.multTPC();
    N_FT0M = col.multFT0M();
    Cent_FT0M = col.centFT0M();
    hist.fill(HIST("QA/after/h_VtxZ"), col.posZ());
    hist.fill(HIST("QA/after/h_Counts"), 2);
    hist.fill(HIST("QA/after/h_NTPC"), NTPC);
    hist.fill(HIST("QA/after/h_Cent"), Cent_FT0M);
    hist.fill(HIST("QA/after/h_NFT0M"), N_FT0M);
    hist.fill(HIST("QA/after/h2_NTPC_NFT0M"), N_FT0M, NTPC);
    hist.fill(HIST("QA/after/h2_NTPC_Cent"), Cent_FT0M, NTPC);
    hist.fill(HIST("QA/after/p_NTPC_Cent"), Cent_FT0M, NTPC);
    hist.fill(HIST("QA/after/p_NTPC_NFT0M"), N_FT0M, NTPC);
    hist.fill(HIST("QA/after/p_NFT0M_NTPC"), NTPC, N_FT0M);
    hist.fill(HIST("QA/after/h2_NTPC_Nch"), NTPC, Nch);

    hist.fill(HIST("Analysis/Charged/h_Mult_ch"), Nch);
    hist.fill(HIST("Analysis/Pion/h_Mult_pi"), N_Pi);
    hist.fill(HIST("Analysis/Kaon/h_Mult_ka"), N_Ka);
    hist.fill(HIST("Analysis/Proton/h_Mult_pr"), N_Pr);

    hist.fill(HIST("QA/after/Charged/h_Mult_ch"), Nch);
    hist.fill(HIST("QA/after/Pion/h_Mult_pi"), N_Pi);
    hist.fill(HIST("QA/after/Kaon/h_Mult_ka"), N_Ka);
    hist.fill(HIST("QA/after/Proton/h_Mult_pr"), N_Pr);

    // Charged Particles:
    if (Nch > 1) {
      auto Nch2 = static_cast<double>(Nch) * (static_cast<double>(Nch) - 1);
      auto mean_Q1 = Q1_ch / static_cast<double>(Nch);
      auto twopart = ((Q1_ch * Q1_ch) - Q2_ch);
      auto twopart1 = (twopart) / (Nch2);
      hist.fill(HIST("Analysis/Charged/h_mean_Q1_ch"), mean_Q1);
      hist.fill(HIST("Analysis/Charged/p_mean_Q1_ch"), NTPC, mean_Q1);
      hist.fill(HIST("Analysis/Charged/h_mean_Q1_Mult_ch"), NTPC, mean_Q1, N_FT0M);
      hist.fill(HIST("Analysis/Charged/h_twopart_Mult_ch"), NTPC, twopart1, N_FT0M);
    }

    if (Nch > 2) {
      auto Nch3 = static_cast<double>(Nch) * (static_cast<double>(Nch) - 1) * (static_cast<double>(Nch) - 2);
      auto threepart = ((Q1_ch * Q1_ch * Q1_ch) - (3 * Q2_ch * Q1_ch) + 2 * Q3_ch);
      auto threepart1 = threepart / Nch3;
      hist.fill(HIST("Analysis/Charged/h_threepart_Mult_ch"), NTPC, threepart1, N_FT0M);
    }

    if (Nch > 3) {
      auto Nch4 = static_cast<double>(Nch) * (static_cast<double>(Nch) - 1) * (static_cast<double>(Nch) - 2) * (static_cast<double>(Nch) - 3);
      auto fourpart = ((Q1_ch * Q1_ch * Q1_ch * Q1_ch) - (6 * Q2_ch * Q1_ch * Q1_ch) + (3 * Q2_ch * Q2_ch) + (8 * Q3_ch * Q1_ch) - 6 * Q4_ch);
      auto fourpart1 = fourpart / Nch4;
      hist.fill(HIST("Analysis/Charged/h_fourpart_Mult_ch"), NTPC, fourpart1, N_FT0M);
    }

    // Pions:
    if (N_Pi > 1) {
      auto Nch2_Pi = static_cast<double>(N_Pi) * (static_cast<double>(N_Pi) - 1);
      auto mean_Q1_Pi = Q1_Pi / static_cast<double>(N_Pi);
      auto twopart_Pi = ((Q1_Pi * Q1_Pi) - Q2_Pi);
      auto twopart1_Pi = (twopart_Pi) / (Nch2_Pi);
      hist.fill(HIST("Analysis/Pion/h_mean_Q1_pi"), mean_Q1_Pi);
      hist.fill(HIST("Analysis/Pion/p_mean_Q1_pi"), NTPC, mean_Q1_Pi);
      hist.fill(HIST("Analysis/Pion/h_mean_Q1_Mult_pi"), NTPC, mean_Q1_Pi, N_FT0M);
      hist.fill(HIST("Analysis/Pion/h_twopart_Mult_pi"), NTPC, twopart1_Pi, N_FT0M);
    }

    if (N_Pi > 2) {
      auto Nch3_Pi = static_cast<double>(N_Pi) * (static_cast<double>(N_Pi) - 1) * (static_cast<double>(N_Pi) - 2);
      auto threepart_Pi = ((Q1_Pi * Q1_Pi * Q1_Pi) - (3 * Q2_Pi * Q1_Pi) + 2 * Q3_Pi);
      auto threepart1_Pi = threepart_Pi / Nch3_Pi;
      hist.fill(HIST("Analysis/Pion/h_threepart_Mult_pi"), NTPC, threepart1_Pi, N_FT0M);
    }

    if (N_Pi > 3) {
      auto Nch4_Pi = static_cast<double>(N_Pi) * (static_cast<double>(N_Pi) - 1) * (static_cast<double>(N_Pi) - 2) * (static_cast<double>(N_Pi) - 3);
      auto fourpart_Pi = ((Q1_Pi * Q1_Pi * Q1_Pi * Q1_Pi) - (6 * Q2_Pi * Q1_Pi * Q1_Pi) + (3 * Q2_Pi * Q2_Pi) + (8 * Q3_Pi * Q1_Pi) - 6 * Q4_Pi);
      auto fourpart1_Pi = fourpart_Pi / Nch4_Pi;
      hist.fill(HIST("Analysis/Pion/h_fourpart_Mult_pi"), NTPC, fourpart1_Pi, N_FT0M);
    }

    // Kaons:
    if (N_Ka > 1) {
      auto Nch2_Ka = static_cast<double>(N_Ka) * (static_cast<double>(N_Ka) - 1);
      auto mean_Q1_Ka = Q1_Ka / static_cast<double>(N_Ka);
      auto twopart_Ka = ((Q1_Ka * Q1_Ka) - Q2_Ka);
      auto twopart1_Ka = (twopart_Ka) / (Nch2_Ka);
      hist.fill(HIST("Analysis/Kaon/h_mean_Q1_ka"), mean_Q1_Ka);
      hist.fill(HIST("Analysis/Kaon/p_mean_Q1_ka"), NTPC, mean_Q1_Ka);
      hist.fill(HIST("Analysis/Kaon/h_mean_Q1_Mult_ka"), NTPC, mean_Q1_Ka, N_FT0M);
      hist.fill(HIST("Analysis/Kaon/h_twopart_Mult_ka"), NTPC, twopart1_Ka, N_FT0M);
    }

    if (N_Ka > 2) {
      auto Nch3_Ka = static_cast<double>(N_Ka) * (static_cast<double>(N_Ka) - 1) * (static_cast<double>(N_Ka) - 2);
      auto threepart_Ka = ((Q1_Ka * Q1_Ka * Q1_Ka) - (3 * Q2_Ka * Q1_Ka) + 2 * Q3_Ka);
      auto threepart1_Ka = threepart_Ka / Nch3_Ka;
      hist.fill(HIST("Analysis/Kaon/h_threepart_Mult_ka"), NTPC, threepart1_Ka, N_FT0M);
    }

    if (N_Ka > 3) {
      auto Nch4_Ka = static_cast<double>(N_Ka) * (static_cast<double>(N_Ka) - 1) * (static_cast<double>(N_Ka) - 2) * (static_cast<double>(N_Ka) - 3);
      auto fourpart_Ka = ((Q1_Ka * Q1_Ka * Q1_Ka * Q1_Ka) - (6 * Q2_Ka * Q1_Ka * Q1_Ka) + (3 * Q2_Ka * Q2_Ka) + (8 * Q3_Ka * Q1_Ka) - 6 * Q4_Ka);
      auto fourpart1_Ka = fourpart_Ka / Nch4_Ka;
      hist.fill(HIST("Analysis/Kaon/h_fourpart_Mult_ka"), NTPC, fourpart1_Ka, N_FT0M);
    }

    // Protons:
    if (N_Pr > 1) {
      auto Nch2_Pr = static_cast<double>(N_Pr) * (static_cast<double>(N_Pr) - 1);
      auto mean_Q1_Pr = Q1_Pr / static_cast<double>(N_Pr);
      auto twopart_Pr = ((Q1_Pr * Q1_Pr) - Q2_Pr);
      auto twopart1_Pr = (twopart_Pr) / (Nch2_Pr);
      hist.fill(HIST("Analysis/Proton/h_mean_Q1_pr"), mean_Q1_Pr);
      hist.fill(HIST("Analysis/Proton/p_mean_Q1_pr"), NTPC, mean_Q1_Pr);
      hist.fill(HIST("Analysis/Proton/h_mean_Q1_Mult_pr"), NTPC, mean_Q1_Pr, N_FT0M);
      hist.fill(HIST("Analysis/Proton/h_twopart_Mult_pr"), NTPC, twopart1_Pr, N_FT0M);
    }

    if (N_Pr > 2) {
      auto Nch3_Pr = static_cast<double>(N_Pr) * (static_cast<double>(N_Pr) - 1) * (static_cast<double>(N_Pr) - 2);
      auto threepart_Pr = ((Q1_Pr * Q1_Pr * Q1_Pr) - (3 * Q2_Pr * Q1_Pr) + 2 * Q3_Pr);
      auto threepart1_Pr = threepart_Pr / Nch3_Pr;
      hist.fill(HIST("Analysis/Proton/h_threepart_Mult_pr"), NTPC, threepart1_Pr, N_FT0M);
    }

    if (N_Pr > 3) {
      auto Nch4_Pr = static_cast<double>(N_Pr) * (static_cast<double>(N_Pr) - 1) * (static_cast<double>(N_Pr) - 2) * (static_cast<double>(N_Pr) - 3);
      auto fourpart_Pr = ((Q1_Pr * Q1_Pr * Q1_Pr * Q1_Pr) - (6 * Q2_Pr * Q1_Pr * Q1_Pr) + (3 * Q2_Pr * Q2_Pr) + (8 * Q3_Pr * Q1_Pr) - 6 * Q4_Pr);
      auto fourpart1_Pr = fourpart_Pr / Nch4_Pr;
      hist.fill(HIST("Analysis/Proton/h_fourpart_Mult_pr"), NTPC, fourpart1_Pr, Cent_FT0M);
    }

    //----------------------------- TPC (No p cuts)---------------------------//
    if (N_Pi_tpc > 1) {
      double mean_Q1_Pi_tpc = Q1_Pi_tpc / static_cast<double>(N_Pi_tpc);
      double twopart_Pi_tpc = ((Q1_Pi_tpc * Q1_Pi_tpc) - Q2_Pi_tpc) / (static_cast<double>(N_Pi_tpc) * (static_cast<double>(N_Pi_tpc) - 1));
      hist.fill(HIST("Analysis/Pion/h_mean_Q1_Mult_pi_tpc"), NTPC, mean_Q1_Pi_tpc, N_FT0M);
      hist.fill(HIST("Analysis/Pion/h_twopart_Mult_pi_tpc"), NTPC, twopart_Pi_tpc, N_FT0M);
    }
    if (N_Ka_tpc > 1) {
      double mean_Q1_Ka_tpc = Q1_Ka_tpc / static_cast<double>(N_Ka_tpc);
      double twopart_Ka_tpc = ((Q1_Ka_tpc * Q1_Ka_tpc) - Q2_Ka_tpc) / (static_cast<double>(N_Ka_tpc) * (static_cast<double>(N_Ka_tpc) - 1));
      hist.fill(HIST("Analysis/Kaon/h_mean_Q1_Mult_ka_tpc"), NTPC, mean_Q1_Ka_tpc, N_FT0M);
      hist.fill(HIST("Analysis/Kaon/h_twopart_Mult_ka_tpc"), NTPC, twopart_Ka_tpc, N_FT0M);
    }
    if (N_Pr_tpc > 1) {
      double mean_Q1_Pr_tpc = Q1_Pr_tpc / static_cast<double>(N_Pr_tpc);
      double twopart_Pr_tpc = ((Q1_Pr_tpc * Q1_Pr_tpc) - Q2_Pr_tpc) / (static_cast<double>(N_Pr_tpc) * (static_cast<double>(N_Pr_tpc) - 1));
      hist.fill(HIST("Analysis/Proton/h_mean_Q1_Mult_pr_tpc"), NTPC, mean_Q1_Pr_tpc, N_FT0M);
      hist.fill(HIST("Analysis/Proton/h_twopart_Mult_pr_tpc"), NTPC, twopart_Pr_tpc, N_FT0M);
    }

    //-----------------------TPC + TOF (No p cuts)--------------------------//
    if (Nch_tof > 1) {
      double mean_Q1_tof = Q1_tof / static_cast<double>(Nch_tof);
      double twopart_tof = ((Q1_tof * Q1_tof) - Q2_tof) / (static_cast<double>(Nch_tof) * (static_cast<double>(Nch_tof) - 1));
      hist.fill(HIST("Analysis/Charged/h_mean_Q1_Mult_ch_tof"), NTPC, mean_Q1_tof, N_FT0M);
      hist.fill(HIST("Analysis/Charged/h_twopart_Mult_ch_tof"), NTPC, twopart_tof, N_FT0M);
    }
    if (N_Pi_tof > 1) {
      double mean_Q1_Pi_tof = Q1_Pi_tof / static_cast<double>(N_Pi_tof);
      double twopart_Pi_tof = ((Q1_Pi_tof * Q1_Pi_tof) - Q2_Pi_tof) / (static_cast<double>(N_Pi_tof) * (static_cast<double>(N_Pi_tof) - 1));
      hist.fill(HIST("Analysis/Pion/h_mean_Q1_Mult_pi_tof"), NTPC, mean_Q1_Pi_tof, N_FT0M);
      hist.fill(HIST("Analysis/Pion/h_twopart_Mult_pi_tof"), NTPC, twopart_Pi_tof, N_FT0M);
    }
    if (N_Ka_tof > 1) {
      double mean_Q1_Ka_tof = Q1_Ka_tof / static_cast<double>(N_Ka_tof);
      double twopart_Ka_tof = ((Q1_Ka_tof * Q1_Ka_tof) - Q2_Ka_tof) / (static_cast<double>(N_Ka_tof) * (static_cast<double>(N_Ka_tof) - 1));
      hist.fill(HIST("Analysis/Kaon/h_mean_Q1_Mult_ka_tof"), NTPC, mean_Q1_Ka_tof, N_FT0M);
      hist.fill(HIST("Analysis/Kaon/h_twopart_Mult_ka_tof"), NTPC, twopart_Ka_tof, N_FT0M);
    }
    if (N_Pr_tof > 1) {
      double mean_Q1_Pr_tof = Q1_Pr_tof / static_cast<double>(N_Pr_tof);
      double twopart_Pr_tof = ((Q1_Pr_tof * Q1_Pr_tof) - Q2_Pr_tof) / (static_cast<double>(N_Pr_tof) * (static_cast<double>(N_Pr_tof) - 1));
      hist.fill(HIST("Analysis/Proton/h_mean_Q1_Mult_pr_tof"), NTPC, mean_Q1_Pr_tof, N_FT0M);
      hist.fill(HIST("Analysis/Proton/h_twopart_Mult_pr_tof"), NTPC, twopart_Pr_tof, N_FT0M);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<meanPtFluc>(cfgc)};
}
