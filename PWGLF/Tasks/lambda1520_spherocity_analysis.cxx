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

/// \file lambda1520_spherocity_analysis.cxx
/// \brief Produce Spherocity table.
///        Invariant Mass Reconstruction of Lambda(1520) Resonance.
///
/// \author Yash Patley <yash.patley@cern.ch>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

double massKa = TDatabasePDG::Instance()->GetParticle(321)->Mass();
double massPr = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

namespace o2::aod
{
namespace mycols
{
DECLARE_SOA_COLUMN(Mult, mult, int);
DECLARE_SOA_COLUMN(Sp, sp, float);
} // namespace mycols
DECLARE_SOA_TABLE(MyCols, "AOD", "MYCOLS", mycols::Mult, mycols::Sp);
} // namespace o2::aod

struct myTable {

  Produces<aod::MyCols> rowMyCols;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axSp(120, -0.1, 1.1, "s_{0}");
    const AxisSpec axMult(10, 0, 10, "Mult");
    const AxisSpec axCtr(1, 0, 1, "CTR");
    histos.add("hCtr", "CTR", kTH1F, {axCtr});
    histos.add("hMult", "MULT", kTH1F, {axMult});
    histos.add("hSp", "Transverse Spherocity", kTH1F, {axSp});
  }

  void process(aod::ResoCollision const& col, aod::ResoTracks const& tracks)
  {
    int size = tracks.size();
    float Sp = 1;
    for (auto const& trk1 : tracks) {
      float sum1 = 0;
      float phi1 = trk1.phi();
      int ctr = 0;
      for (auto const& trk2 : tracks) {
        ++ctr;
        if (trk1.index() == trk2.index())
          continue;
        float phi2 = trk2.phi();
        sum1 += abs(sin(phi1 - phi2));
      }
      if (ctr == 0)
        histos.fill(HIST("hCtr"), 0.5);
      if (sum1 == 0)
        histos.fill(HIST("hMult"), size);
      float sph = pow(sum1 / static_cast<float>(size), 2);
      if (sph < Sp) {
        Sp = sph;
      }
    }
    float spher = pow(M_PI_2, 2) * Sp;
    histos.fill(HIST("hSp"), spher);
    rowMyCols(size, spher);
  }
};

struct lambdaTask {

  // Configurables for filters.
  Configurable<float> cEta{"cEta", 0.8, "Pseudorapidity bound."};
  Configurable<float> cMinPt{"cMinPt", 0.15, "Minimum transvesrse momentum."};
  Configurable<float> cMaxPt{"cMaxPt", 20., "Maximum transvesrse momentum."};
  Configurable<float> cDCAz{"cDCAz", 2., "DCA in z direction."};
  Configurable<float> cDCAxy{"cDCAxy", 0.12, "DCA in xy plane."};

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histogram."};
  Configurable<int> nBinsMult{"nBinsMult", 100,
                              "N bins in Multiplicity histograms."};
  Configurable<int> nBinsInvM{"nBinsInvM", 1600,
                              "N bins in InvMass histograms."};

  void init(InitContext const&)
  {
    // Define Axis.
    const AxisSpec axEv(1, 0, 1, "N_{Ev}");
    const AxisSpec axEta(200, -1, 1, "#eta");
    const AxisSpec axMult(nBinsMult, 0, 1000, "Multiplicity");
    const AxisSpec axSp(120, -0.1, 1.1, "S_{0}");
    const AxisSpec axPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axDCAz(500, -2.5, 2.5, {"DCA_{z}"});
    const AxisSpec axDCAxy(50, -0.25, 0.25, {"DCA_{xy}"});
    const AxisSpec axTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axTPCNSigma(28, -7, 7, {"TPC N^{Sigma}"});
    const AxisSpec axTOFNSigma(28, -7, 7, {"TOF N^{Sigma}"});
    const AxisSpec axInvM(nBinsInvM, 1.4, 3., {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    histos.add("Events/hEv", "Number of Events", kTH1F, {axEv});
    histos.add("Events/hMult", "TPCTemporal Multiplicity of Events", kTH1F,
               {axMult});
    histos.add("Events/hSp", "Transverse Spherocity of Events", kTH1F, {axSp});
    histos.add("Tracks/hPt", "Transverse Momentum", kTH1F, {axPt});
    histos.add("Tracks/hEta", "Pseudorapidity", kTH1F, {axEta});
    histos.add("QA_beforeCuts/Protons/hTPCNSigma", "TPC N^{Sigma} Protons",
               kTH2F, {{axPt}, {axTPCNSigma}});
    histos.add("QA_beforeCuts/Protons/hTOFNSigma", "TOF N^{Sigma} Protons",
               kTH2F, {{axPt}, {axTOFNSigma}});
    histos.add("QA_beforeCuts/Protons/hTpcTofMap", "TPC TOF N^{Sigma} Protons",
               kTH2F, {{axTOFNSigma}, {axTPCNSigma}});
    histos.add("QA_beforeCuts/Kaons/hTPCNSigma", "TPC N^{Sigma} Kaons", kTH2F,
               {{axPt}, {axTPCNSigma}});
    histos.add("QA_beforeCuts/Kaons/hTOFNSigma", "TOF N^{Sigma} Kaons", kTH2F,
               {{axPt}, {axTOFNSigma}});
    histos.add("QA_beforeCuts/Kaons/hTpcTofMap", "TPC TOF N^{Sigma} Kaons",
               kTH2F, {{axTOFNSigma}, {axTPCNSigma}});
    histos.add("QA_afterCuts/Protons/hDCAz", "DCAz Protons", kTH2F,
               {{axPt}, {axDCAz}});
    histos.add("QA_afterCuts/Protons/hDCAxy", "DCAxy Protons", kTH2F,
               {{axPt}, {axDCAxy}});
    histos.add("QA_afterCuts/Protons/hTPCNCls", "TPCNCls Protons", kTH2F,
               {{axPt}, {axTPCNCls}});
    histos.add("QA_afterCuts/Protons/hTPCNSigma", "TPC N^{Sigma} Protons",
               kTH2F, {{axPt}, {axTPCNSigma}});
    histos.add("QA_afterCuts/Protons/hTOFNSigma", "TOF N^{Sigma} Protons",
               kTH2F, {{axPt}, {axTOFNSigma}});
    histos.add("QA_afterCuts/Protons/hTpcTofMap", "TPC TOF N^{Sigma} Protons",
               kTH2F, {{axTOFNSigma}, {axTPCNSigma}});
    histos.add("QA_afterCuts/Kaons/hDCAz", "DCAz Kaons", kTH2F,
               {{axPt}, {axDCAz}});
    histos.add("QA_afterCuts/Kaons/hDCAxy", "DCAxy Kaons", kTH2F,
               {{axPt}, {axDCAxy}});
    histos.add("QA_afterCuts/Kaons/hTPCNCls", "TPCNCls Kaons", kTH2F,
               {{axPt}, {axTPCNCls}});
    histos.add("QA_afterCuts/Kaons/hTPCNSigma", "TPC N^{Sigma} Kaons", kTH2F,
               {{axPt}, {axTPCNSigma}});
    histos.add("QA_afterCuts/Kaons/hTOFNSigma", "TOF N^{Sigma} Kaons", kTH2F,
               {{axPt}, {axTOFNSigma}});
    histos.add("QA_afterCuts/Kaons/hTpcTofMap", "TPC TOF N^{Sigma} Kaons",
               kTH2F, {{axTOFNSigma}, {axTPCNSigma}});
    histos.add("Analysis/hPtProton", "Protons p_{T}", kTH1F, {axPt});
    histos.add("Analysis/hEtaProton", "Protons #eta", kTH1F, {axEta});
    histos.add("Analysis/hPtKaon", "Kaons p_{T}", kTH1F, {axPt});
    histos.add("Analysis/hEtaKaon", "Kaons #eta", kTH1F, {axEta});
    histos.add("Analysis/h4Lambda1", "THn #Lambda to p K^{-}", kTHnSparseF,
               {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4Lambda2", "THn #bar{#Lambda} to #bar{p} K^{+}",
               kTHnSparseF, {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4LikeSigns1", "THn Like Signs p K^{+}", kTHnSparseF,
               {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4LikeSigns2", "THn Like Signs #bar{p} K^{-}",
               kTHnSparseF, {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4Lambda",
               "THn for #Lambda(1520) and #bar{#Lambda}(1520)", kTHnSparseF,
               {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4LikeSigns", "THn Like Signs", kTHnSparseF,
               {axInvM, axPt, axSp, axMult});
    histos.add("Analysis/h4Mixed",
               "THn for Mixed Events #Lambda(1520) and #bar{#Lambda}(1520)",
               kTHnSparseF, {axInvM, axPt, axSp, axMult});
  }

  template <bool isMix, typename trackType>
  void fillHistos(int const& mult, float const& sp, trackType const& trkPr,
                  trackType const& trkKa)
  {
    for (auto const& [trk1, trk2] :
         soa::combinations(soa::CombinationsFullIndexPolicy(trkPr, trkKa))) {

      if (trk1.index() == trk2.index())
        continue;

      bool isPrSel = false;
      bool isKaSel = false;

      // Preselection Criterion.
      // TPC Only.
      if (abs(trk1.tpcNSigmaPr()) < 6) {
        histos.fill(HIST("QA_beforeCuts/Protons/hTPCNSigma"), trk1.pt(),
                    trk1.tpcNSigmaPr());
      }

      if (abs(trk2.tpcNSigmaKa()) < 6) {
        histos.fill(HIST("QA_beforeCuts/Kaons/hTPCNSigma"), trk2.pt(),
                    trk2.tpcNSigmaKa());
      }

      // Check TOF, TPC + TOF.
      if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) ==
          aod::resodaughter::kHasTOF) {
        // Proton.
        if (abs(trk1.tofNSigmaPr()) < 6 && abs(trk1.tpcNSigmaPr()) < 6) {
          histos.fill(HIST("QA_beforeCuts/Protons/hTOFNSigma"), trk1.pt(),
                      trk1.tofNSigmaPr());
          histos.fill(HIST("QA_beforeCuts/Protons/hTpcTofMap"),
                      trk1.tofNSigmaPr(), trk1.tpcNSigmaPr());
        }
      }
      if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) ==
          aod::resodaughter::kHasTOF) {
        // Kaon.
        if (abs(trk2.tofNSigmaKa()) < 6 && abs(trk2.tpcNSigmaKa()) < 6) {
          histos.fill(HIST("QA_beforeCuts/Kaons/hTOFNSigma"), trk2.pt(),
                      trk2.tofNSigmaKa());
          histos.fill(HIST("QA_beforeCuts/Kaons/hTpcTofMap"),
                      trk2.tofNSigmaKa(), trk2.tpcNSigmaKa());
        }
      }

      // TPC/TOF Cuts Selection Protons and Kaons.
      // TPC Only.
      if (abs(trk1.tpcNSigmaPr()) < 2) {
        histos.fill(HIST("QA_afterCuts/Protons/hDCAz"), trk1.pt(), trk1.dcaZ());
        histos.fill(HIST("QA_afterCuts/Protons/hDCAxy"), trk1.pt(),
                    trk1.dcaXY());
        histos.fill(HIST("QA_afterCuts/Protons/hTPCNCls"), trk1.pt(),
                    trk1.tpcNClsCrossedRows());
        histos.fill(HIST("QA_afterCuts/Protons/hTPCNSigma"), trk1.pt(),
                    trk1.tpcNSigmaPr());
        isPrSel = true;
      }

      if (abs(trk2.tpcNSigmaKa()) < 2) {
        histos.fill(HIST("QA_afterCuts/Kaons/hDCAz"), trk2.pt(), trk2.dcaZ());
        histos.fill(HIST("QA_afterCuts/Kaons/hDCAxy"), trk2.pt(), trk2.dcaXY());
        histos.fill(HIST("QA_afterCuts/Kaons/hTPCNCls"), trk2.pt(),
                    trk2.tpcNClsCrossedRows());
        histos.fill(HIST("QA_afterCuts/Kaons/hTPCNSigma"), trk2.pt(),
                    trk2.tpcNSigmaKa());
        isKaSel = true;
      }

      // Check TOF.
      if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) ==
          aod::resodaughter::kHasTOF) {
        // TPC + TOF.
        // Proton.
        if (abs(trk1.tofNSigmaPr()) < 2 && abs(trk1.tpcNSigmaPr()) < 3) {
          histos.fill(HIST("QA_afterCuts/Protons/hDCAz"), trk1.pt(),
                      trk1.dcaZ());
          histos.fill(HIST("QA_afterCuts/Protons/hDCAxy"), trk1.pt(),
                      trk1.dcaXY());
          histos.fill(HIST("QA_afterCuts/Protons/hTPCNCls"), trk1.pt(),
                      trk1.tpcNClsCrossedRows());
          histos.fill(HIST("QA_afterCuts/Protons/hTOFNSigma"), trk1.pt(),
                      trk1.tofNSigmaPr());
          histos.fill(HIST("QA_afterCuts/Protons/hTpcTofMap"),
                      trk1.tofNSigmaPr(), trk1.tpcNSigmaPr());
          isPrSel = true;
        }
      }

      if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) ==
          aod::resodaughter::kHasTOF) {
        // TPC + TOF.
        // Kaon.
        if (abs(trk2.tofNSigmaKa()) < 2 && abs(trk2.tpcNSigmaKa()) < 3) {
          histos.fill(HIST("QA_afterCuts/Kaons/hDCAz"), trk2.pt(), trk2.dcaZ());
          histos.fill(HIST("QA_afterCuts/Kaons/hDCAxy"), trk2.pt(),
                      trk2.dcaXY());
          histos.fill(HIST("QA_afterCuts/Kaons/hTPCNCls"), trk2.pt(),
                      trk2.tpcNClsCrossedRows());
          histos.fill(HIST("QA_afterCuts/Kaons/hTOFNSigma"), trk2.pt(),
                      trk2.tofNSigmaKa());
          histos.fill(HIST("QA_afterCuts/Kaons/hTpcTofMap"), trk2.tofNSigmaKa(),
                      trk2.tpcNSigmaKa());
          isKaSel = true;
        }
      }

      if (!isPrSel || !isKaSel)
        continue;

      // pT spectra of protons and kaons.
      histos.fill(HIST("Analysis/hPtProton"), trk1.pt());
      histos.fill(HIST("Analysis/hPtKaon"), trk2.pt());

      // #eta of protons and kaons.
      histos.fill(HIST("Analysis/hEtaProton"), trk1.eta());
      histos.fill(HIST("Analysis/hEtaKaon"), trk2.eta());

      TLorentzVector p, p1, p2;
      p1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPr);
      p2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      p = p1 + p2;
      if (abs(p.Rapidity()) > 0.5)
        continue;

      if (trk1.sign() * trk2.sign() < 0) {
        if (!isMix) {
          histos.fill(HIST("Analysis/h4Lambda"), p.M(), p.Pt(), sp, mult);
          if (trk1.sign() == +1 && trk2.sign() == -1) {
            histos.fill(HIST("Analysis/h4Lambda1"), p.M(), p.Pt(), sp, mult);
          } else {
            histos.fill(HIST("Analysis/h4Lambda2"), p.M(), p.Pt(), sp, mult);
          }
        } else {
          histos.fill(HIST("Analysis/h4Mixed"), p.M(), p.Pt(), sp, mult);
        }
      }

      if (trk1.sign() * trk2.sign() > 0 && !isMix) {
        histos.fill(HIST("Analysis/h4LikeSigns"), p.M(), p.Pt(), sp, mult);
        if (trk1.sign() == +1 && trk2.sign() == +1) {
          histos.fill(HIST("Analysis/h4LikeSigns1"), p.M(), p.Pt(), sp, mult);
        } else {
          histos.fill(HIST("Analysis/h4LikeSigns2"), p.M(), p.Pt(), sp, mult);
        }
      }
    }
  }

  // Collision and Track Table.
  using myCollisions =
    soa::Filtered<soa::Join<aod::ResoCollisions, aod::MyCols>>;
  using myTracks = soa::Filtered<aod::ResoTracks>;

  // Track Filters.
  Filter fTracks = aod::resodaughter::pt > cMinPt&& aod::resodaughter::pt <
                   cMaxPt&& nabs(aod::resodaughter::eta) < cEta;
  Filter fDCA =
    nabs(aod::track::dcaZ) < cDCAz && nabs(aod::track::dcaXY) < cDCAxy;
  Filter fNCLs = aod::resodaughter::tpcNClsCrossedRows >
                 static_cast<uint8_t>(70);

  // Collision Filters.
  Filter fCol = aod::mycols::mult > 2;

  void processData(myCollisions::iterator const& collision,
                   myTracks const& tracks)
  {
    histos.fill(HIST("Events/hEv"), 0.5);
    histos.fill(HIST("Events/hMult"), collision.mult());
    histos.fill(HIST("Events/hSp"), collision.sp());

    fillHistos<false>(collision.mult(), collision.sp(), tracks, tracks);

    for (auto const& track : tracks) {
      histos.fill(HIST("Tracks/hPt"), track.pt());
      histos.fill(HIST("Tracks/hEta"), track.eta());
    }
  }

  PROCESS_SWITCH(lambdaTask, processData, "Process Lambda(1520) Data", true);

  SliceCache cache;
  vector<double> vZBins{
    VARIABLE_WIDTH, -10., -8., -6., -4., -2., 2., 4., 6., 8., 10.};
  vector<double> vSBins{VARIABLE_WIDTH, 0.2, 0.4, 0.6, 0.8, 1.0};
  using BinningType =
    ColumnBinningPolicy<aod::ResoCollision::PosZ, aod::mycols::Sp>;
  BinningType binningOnPositions{{vZBins, vSBins}, true};

  void processMix(myCollisions& collisions, myTracks const& tracks)
  {

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<myCollisions, myTracks, BinningType> pair{
      binningOnPositions, 5, -1, collisions, tracksTuple, &cache};

    for (auto& [c1, t1, c2, t2] : pair) {
      fillHistos<true>(c1.mult(), c1.sp(), t1, t2);
    }
  }

  PROCESS_SWITCH(lambdaTask, processMix, "Process Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<myTable>(cfgc),
                      adaptAnalysisTask<lambdaTask>(cfgc)};
}
