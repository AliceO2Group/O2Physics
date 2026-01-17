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
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SGPIDSpectra {
  SGSelector sgSelector;
  ConfigurableAxis ptAxis{"ptAxis", {1000, 0.0, 20.0}, "p_{T}"};
  ConfigurableAxis sigmaAxis{"sigmaAxis", {2000, -20.0, 180.0}, "#sigma"};

  // configurables
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const AxisSpec axispt{ptAxis, "p_{T}"};
    const AxisSpec axistpc{sigmaAxis, "N_{#sigma}^{TPC}"};
    const AxisSpec axistof{sigmaAxis, "N_{#sigma}^{TOF}"};
    // Collision histograms
    registry.add("collisions/GapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("collisions/TrueGapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{4, -1.5, 2.5}}});
    registry.add("ttracks/pPion_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_Ka", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_Pr", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_Pi", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_Pr", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_Ka", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_Pi", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_Ka", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_Pr", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_Pi", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_Pr", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_El", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_Ka", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_Pi", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pPion_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nPion_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pKaon_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nKaon_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_p_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_n_0", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_p_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_n_1", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_p_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_n_2", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_p_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/pProton_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("ttracks/nProton_Pt_TPC_n_3", "", {HistType::kTH2F, {axispt, axistpc}});
  }

  // define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    // fill collision histograms
    registry.get<TH1>(HIST("collisions/GapSide"))->Fill(coll.gapSide(), 1.);
    // int truegapSide = sgSelector.trueGap(dgcand, FV0_cut, ZDC_cut);
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(coll, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.get<TH1>(HIST("collisions/TrueGapSide"))->Fill(truegapSide, 1.);
    // select PV contributors
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // check rho0 signals
    for (auto t : tracks) {
      if (trackselector(t, parameters) && t.hasTPC()) {
        if (truegapSide == 0) {
          if (t.sign() > 0) {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          } else {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_p_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_n_0"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          }
        } else if (truegapSide == 1) {
          if (t.sign() > 0) {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          } else {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_p_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_n_1"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          }
        } else if (truegapSide == 2) {
          if (t.sign() > 0) {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/pProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          } else {
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nPion_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaPr()) < 1.)
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_Pr"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
            if (std::abs(t.tpcNSigmaEl()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_El"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaKa()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Ka"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (std::abs(t.tpcNSigmaPi()) < 1.)
              registry.fill(HIST("ttracks/nProton_Pt_TPC_Pi"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_p_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_n_2"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          }
        } else {
          if (t.sign() > 0) {
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/pPion_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/pKaon_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/pProton_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          } else {
            if (t.pz() > 0) {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_p_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            } else {
              registry.fill(HIST("ttracks/nPion_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPi());
              registry.fill(HIST("ttracks/nKaon_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaKa());
              registry.fill(HIST("ttracks/nProton_Pt_TPC_n_3"), TMath::Abs(t.px() * t.px() + t.py() * t.py()), t.tpcNSigmaPr());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGPIDSpectra>(cfgc, TaskName{"sgpidspectra"}),
  };
}
