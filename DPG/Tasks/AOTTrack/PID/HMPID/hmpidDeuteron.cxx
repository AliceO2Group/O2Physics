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

// \brief Analysis task for the (anti-)deuteron absorption cross section using HMPID
// \authors rocco.liotino@cern.ch

#include "tableHMPID.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;

struct HmpidDeuteron {
  HistogramRegistry registryDA{"registryDA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.0, "nsigmaTOFMax"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};

  Configurable<float> qMipThreshold{"qMipThreshold", 120.0f, "minimum MIP charge required for a real signal"};
  Configurable<float> cutDistanceMipTrack{"cutDistanceMipTrack", 6.0f, "cut distance between MIP and track"};

  Configurable<float> offsetX{"offsetX", 0.696712f, "Residuals offset - X direction"};
  Configurable<float> offsetY{"offsetY", 2.32849f, "Residuals offset - Y direction"};

  Configurable<bool> collisionsPbPb{"collisionsPbPb", false, "Require PbPb collisions"};

  Configurable<float> centralityCut{"centralityCut", 20.0f, "centrality cut for PbPb collisions"};

  // Configurables for plotting
  //  Momentum
  Configurable<int> momNBins{"momNBins", 490, "Momentum bins"};
  Configurable<float> momMin{"momMin", 0.1f, "Momentum min"};
  Configurable<float> momMax{"momMax", 5.0f, "Momentum max"};

  // DeltaR
  Configurable<int> deltaRNBins{"deltaRNBins", 50, "DeltaR bins"};
  Configurable<float> deltaRMin{"deltaRMin", 0.0f, "DeltaR min"};
  Configurable<float> deltaRMax{"deltaRMax", 50.0f, "DeltaR max"};

  // nSigma TPC
  Configurable<int> nSigmaTPCNBins{"nSigmaTPCNBins", 100, "nSigmaTPC bins"};
  Configurable<float> nSigmaTPCMinPlot{"nSigmaTPCMinPlot", -5.0f, "nSigmaTPC min"};
  Configurable<float> nSigmaTPCMaxPlot{"nSigmaTPCMaxPlot", 5.0f, "nSigmaTPC max"};

  // nSigma TOF
  Configurable<int> nSigmaTOFNBins{"nSigmaTOFNBins", 100, "nSigmaTOF bins"};
  Configurable<float> nSigmaTOFMinPlot{"nSigmaTOFMinPlot", -5.0f, "nSigmaTOF min"};
  Configurable<float> nSigmaTOFMaxPlot{"nSigmaTOFMaxPlot", 5.0f, "nSigmaTOF max"};

  // Charge
  Configurable<int> chargeNBins{"chargeNBins", 1000, "Charge bins"};
  Configurable<float> chargeMin{"chargeMin", 200.0f, "Charge min"};
  Configurable<float> chargeMax{"chargeMax", 2200.0f, "Charge max"};

  // Eta
  Configurable<int> etaNBins{"etaNBins", 100, "Eta bins"};
  Configurable<float> etaMin{"etaMin", -5.0f, "Eta min"};
  Configurable<float> etaMax{"etaMax", 5.0f, "Eta max"};

  // Phi
  Configurable<int> phiNBins{"phiNBins", 100, "Phi bins"};
  Configurable<float> phiMin{"phiMin", -10.0f, "Phi min"};
  Configurable<float> phiMax{"phiMax", 10.0f, "Phi max"};

  // Mass
  Configurable<int> massNBins{"massNBins", 1000, "Mass bins"};
  Configurable<float> massMin{"massMin", -10.0f, "Mass min"};
  Configurable<float> massMax{"massMax", 10.0f, "Mass max"};

  // variables for chamber_number and HVs/PCs
  static constexpr int Rich0 = 0, Rich1 = 1, Rich2 = 2, Rich3 = 3, Rich4 = 4, Rich5 = 5, Rich6 = 6;
  static constexpr int Nch = 7;
  static constexpr float Nmean = 1.28947; // meanIdxRad(); method from param not working, using value from HMPIDBase/Param.cxx

  void init(o2::framework::InitContext const&)
  {

    // Axes for histograms
    AxisSpec axisMom{momNBins, momMin, momMax, "#it{p} (GeV/#it{c})", "axisMom"};
    AxisSpec axisDeltaR{deltaRNBins, deltaRMin, deltaRMax, "#Delta R (cm)", "axisDeltaR"};
    AxisSpec axisNsigmaTOF{nSigmaTOFNBins, nSigmaTOFMinPlot, nSigmaTOFMaxPlot, "n#sigma_{TOF}", "axisNsigmaTOF"};
    AxisSpec axisNsigmaTPC{nSigmaTPCNBins, nSigmaTPCMinPlot, nSigmaTPCMaxPlot, "n#sigma_{TPC}", "axisNsigmaTPC"};
    AxisSpec axisCharge{chargeNBins, chargeMin, chargeMax, "Charge (ADC)", "axisCharge"};
    AxisSpec axisEta{etaNBins, etaMin, etaMax, "#eta", "axisEta"};
    AxisSpec axisPhi{phiNBins, phiMin, phiMax, "#phi", "axisPhi"};
    AxisSpec axisMass{massNBins, massMin, massMax, "mass^{2} (GeV/#it{c}^{2})", "axisMass"};
    AxisSpec axisChamber{10, -1.5, 8.5, "Chamber number", "axisChamber"};

    // Deuteron Pos
    registryDA.add("incomingDe_Pos_8cm", "incomingDe_Pos_8cm", HistType::kTH1F, {axisMom});
    registryDA.add("incomingDe_Pos_4cm", "incomingDe_Pos_4cm", HistType::kTH1F, {axisMom});

    registryDA.add("De_Pos_deltaR_8cm", "De_Pos_deltaR_8cm", HistType::kTH1F, {axisDeltaR});
    registryDA.add("De_Pos_deltaR_4cm", "De_Pos_deltaR_4cm", HistType::kTH1F, {axisDeltaR});

    registryDA.add("nSigmaTOF_vs_momHMPID_DePos_8cm", "nSigmaTOF_vs_momHMPID_DePos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_DePos_4cm", "nSigmaTOF_vs_momHMPID_DePos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("survivingDe_Pos_8cm", "survivingDe_Pos_8cm", HistType::kTH2F, {axisMom, axisDeltaR});
    registryDA.add("survivingDe_Pos_4cm", "survivingDe_Pos_4cm", HistType::kTH2F, {axisMom, axisDeltaR});

    registryDA.add("De_Pos_Q_8cm", "De_Pos_Q_8cm", HistType::kTH2F, {axisMom, axisCharge});
    registryDA.add("De_Pos_Q_4cm", "De_Pos_Q_4cm", HistType::kTH2F, {axisMom, axisCharge});

    // nsigma plots
    registryDA.add("nSigmaTPC_vs_momHMPID_Cut_DePos", "nSigmaTPC_vs_momHMPID_Cut_DePos", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_Cut_DePos", "nSigmaTOF_vs_momHMPID_Cut_DePos", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    // check residuals distributions
    registryDA.add("residualsX_Vs_momentum_Neg_8cm", "residualsX_Vs_momentum_Neg_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Neg_8cm", "residualsY_Vs_momentum_Neg_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});
    registryDA.add("residualsX_Vs_momentum_Pos_8cm", "residualsX_Vs_momentum_Pos_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Pos_8cm", "residualsY_Vs_momentum_Pos_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});

    // Deuteron Neg
    registryDA.add("incomingDe_Neg_8cm", "incomingDe_Neg_8cm", HistType::kTH1F, {axisMom});
    registryDA.add("incomingDe_Neg_4cm", "incomingDe_Neg_4cm", HistType::kTH1F, {axisMom});

    registryDA.add("De_Neg_deltaR_8cm", "De_Neg_deltaR_8cm", HistType::kTH1F, {axisDeltaR});
    registryDA.add("De_Neg_deltaR_4cm", "De_Neg_deltaR_4cm", HistType::kTH1F, {axisDeltaR});

    registryDA.add("nSigmaTOF_vs_momHMPID_DeNeg_8cm", "nSigmaTOF_vs_momHMPID_DeNeg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_DeNeg_4cm", "nSigmaTOF_vs_momHMPID_DeNeg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("survivingDe_Neg_8cm", "survivingDe_Neg_8cm", HistType::kTH2F, {axisMom, axisDeltaR});
    registryDA.add("survivingDe_Neg_4cm", "survivingDe_Neg_4cm", HistType::kTH2F, {axisMom, axisDeltaR});
    registryDA.add("De_Neg_Q_8cm", "De_Neg_Q_8cm", HistType::kTH2F, {axisMom, axisCharge});
    registryDA.add("De_Neg_Q_4cm", "De_Neg_Q_4cm", HistType::kTH2F, {axisMom, axisCharge});
    registryDA.add("De_Neg_momentum", "De_Neg_momentum", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 5.0, "#it{p}_{hmpid} (GeV/#it{c})"}});

    // check residuals distributions
    registryDA.add("residualsX_Vs_momentum_Neg_4cm", "residualsX_Vs_momentum_Neg_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Neg_4cm", "residualsY_Vs_momentum_Neg_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});
    registryDA.add("residualsX_Vs_momentum_Pos_4cm", "residualsX_Vs_momentum_Pos_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Pos_4cm", "residualsY_Vs_momentum_Pos_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});

    // nsigma plots
    registryDA.add("nSigmaTPC_vs_momHMPID_Cut_DeNeg", "nSigmaTPC_vs_momHMPID_Cut_DeNeg", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_Cut_DeNeg", "nSigmaTOF_vs_momHMPID_Cut_DeNeg", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DePos", "nSigmaTOF_vs_momHMPID_TPCOnly_DePos", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg", "nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DePos_4cm", "nSigmaTOF_vs_momHMPID_TPCOnly_DePos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DePos_8cm", "nSigmaTOF_vs_momHMPID_TPCOnly_DePos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_4cm", "nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_8cm", "nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("pidTPC_TOF_p_ch_DePos", "pidTPC_TOF_p_ch_DePos", HistType::kTHnSparseF, {axisMom, axisNsigmaTPC, axisNsigmaTOF, axisChamber});
    registryDA.add("pidTPC_TOF_p_ch_DeNeg", "pidTPC_TOF_p_ch_DeNeg", HistType::kTHnSparseF, {axisMom, axisNsigmaTPC, axisNsigmaTOF, axisChamber});

    // general plots
    registryDA.add("hEta", "hEta", kTH1F, {axisEta});
    registryDA.add("hPhi", "hPhi", kTH1F, {axisPhi});
    registryDA.add("hMomentumHmpid", "hMomentumHmpid", kTH1F, {{1000, -10.f, 10.f, "#it{p}_{hmpid} (GeV/#it{c})"}});

    registryDA.add("hMass", "hMass", kTH1F, {axisMass});
    registryDA.add("hMass_postDeuteron", "hMass_postDeuteron", kTH1F, {axisMass});
    registryDA.add("hMass_Vs_centrality", "hMass_Vs_centrality", kTH2F, {axisMass, {100, 0.0f, 100.0f, "centrality (%)"}});

    for (int iCh = 0; iCh < Nch; iCh++) {
      registryDA.add(Form("De_Pos_deltaR_%d", iCh), Form("De_Pos_deltaR_%d", iCh), HistType::kTH1F, {axisDeltaR});
      registryDA.add(Form("De_Neg_deltaR_%d", iCh), Form("De_Neg_deltaR_%d", iCh), HistType::kTH1F, {axisDeltaR});
      registryDA.add(Form("hEta_%d", iCh), Form("hEta_%d", iCh), HistType::kTH1F, {axisEta});
      registryDA.add(Form("hPhi_%d", iCh), Form("hPhi_%d", iCh), HistType::kTH1F, {axisPhi});

      // tracking efficiency hmpid
      registryDA.add(Form("hAllTracks_%d", iCh), Form("hAllTracks_%d", iCh), HistType::kTH1F, {axisMom});
      registryDA.add(Form("hCutQandD_%d", iCh), Form("hCutQandD_%d", iCh), HistType::kTH1F, {axisMom});
      registryDA.add(Form("hAllTracks_Pos_%d", iCh), Form("hAllTracks_Pos_%d", iCh), HistType::kTH1F, {axisMom});
      registryDA.add(Form("hAllTracks_Neg_%d", iCh), Form("hAllTracks_Neg_%d", iCh), HistType::kTH1F, {axisMom});
      registryDA.add(Form("hCutQandD_Pos_%d", iCh), Form("hCutQandD_Pos_%d", iCh), HistType::kTH1F, {axisMom});
      registryDA.add(Form("hCutQandD_Neg_%d", iCh), Form("hCutQandD_Neg_%d", iCh), HistType::kTH1F, {axisMom});
    }
  }

  void process(aod::HmpidAnalysis const& hmpidTable)
  {
    for (const auto& hmpid : hmpidTable) {

      // -------------------------
      // track filters
      // -------------------------
      if (hmpid.itsNCluster() < minReqClusterITS) {
        continue;
      }
      if (hmpid.tpcNCluster() < minTPCnClsFound) {
        continue;
      }
      if (hmpid.tpcNClsCrossedRows() < minNCrossedRowsTPC) {
        continue;
      }
      if (hmpid.tpcChi2() > maxChi2TPC) {
        continue;
      }
      if (hmpid.itsChi2() > maxChi2ITS) {
        continue;
      }
      if (std::abs(hmpid.dcaXY()) > maxDCAxy) {
        continue;
      }
      if (std::abs(hmpid.dcaZ()) > maxDCAz) {
        continue;
      }

      // -------------------------
      // derived quantities
      // -------------------------
      const float momHmpid = hmpid.momentumHmpid();
      const float momAbs = std::abs(momHmpid);

      const float dx = hmpid.xMip() - hmpid.xTrack();
      const float dy = hmpid.yMip() - hmpid.yTrack();
      const float dr = std::hypot(dx - offsetX, dy - offsetY);
      const float dMipTrack = std::hypot(dx, dy);

      const float mass =
        std::pow(Nmean * momAbs * std::cos(hmpid.chAngle()), 2) -
        std::pow(momAbs, 2);

      const int chamber = static_cast<int>(hmpid.chamber());
      const bool isPos = momHmpid > 0;
      const bool isNeg = momHmpid < 0;

      const bool passQandD = hmpid.chargeMip() > qMipThreshold && dMipTrack < cutDistanceMipTrack;

      // -------------------------
      // global distributions
      // -------------------------
      registryDA.fill(HIST("hEta"), hmpid.etaTrack());
      registryDA.fill(HIST("hPhi"), hmpid.phiTrack());
      registryDA.fill(HIST("hMomentumHmpid"), momHmpid);

      if (hmpid.chAngle() > 0) {
        registryDA.fill(HIST("hMass"), mass);
        registryDA.fill(HIST("hMass_Vs_centrality"), mass, hmpid.centrality());
      }

      // -------------------------
      // Fill plots per chamber
      // -------------------------

      switch (chamber) {

        case Rich0:
          registryDA.fill(HIST("hEta_0"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_0"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_0"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_0"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_0"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_0"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_0"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_0"), momAbs);
            }
          }
          break;

        case Rich1:
          registryDA.fill(HIST("hEta_1"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_1"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_1"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_1"), momAbs);
          }

          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_1"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_1"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_1"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_1"), momAbs);
            }
          }
          break;

        case Rich2:
          registryDA.fill(HIST("hEta_2"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_2"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_2"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_2"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_2"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_2"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_2"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_2"), momAbs);
            }
          }
          break;

        case Rich3:
          registryDA.fill(HIST("hEta_3"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_3"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_3"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_3"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_3"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_3"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_3"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_3"), momAbs);
            }
          }
          break;

        case Rich4:
          registryDA.fill(HIST("hEta_4"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_4"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_4"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_4"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_4"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_4"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_4"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_4"), momAbs);
            }
          }
          break;

        case Rich5:
          registryDA.fill(HIST("hEta_5"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_5"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_5"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_5"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_5"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_5"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_5"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_5"), momAbs);
            }
          }
          break;

        case Rich6:
          registryDA.fill(HIST("hEta_6"), hmpid.etaTrack());
          registryDA.fill(HIST("hPhi_6"), hmpid.phiTrack());
          registryDA.fill(HIST("hAllTracks_6"), momAbs);
          if (passQandD) {
            registryDA.fill(HIST("hCutQandD_6"), momAbs);
          }
          if (isPos) {
            registryDA.fill(HIST("hAllTracks_Pos_6"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Pos_6"), momAbs);
            }
          }
          if (isNeg) {
            registryDA.fill(HIST("hAllTracks_Neg_6"), momAbs);
            if (passQandD) {
              registryDA.fill(HIST("hCutQandD_Neg_6"), momAbs);
            }
          }
          break;

        default:
          break;
      }

      // -------------------------
      // TPC + TOF PID distributions
      // -------------------------
      if (isPos) {
        registryDA.fill(HIST("pidTPC_TOF_p_ch_DePos"), momAbs, hmpid.tpcNSigmaDe(), hmpid.tofNSigmaDe(), chamber);
      }
      if (isNeg) {
        registryDA.fill(HIST("pidTPC_TOF_p_ch_DeNeg"), momAbs, hmpid.tpcNSigmaDe(), hmpid.tofNSigmaDe(), chamber);
      }

      // -------------------------
      // deuteron candidate cuts - TPC
      // -------------------------
      if (hmpid.tpcNSigmaDe() < nsigmaTPCMin ||
          hmpid.tpcNSigmaDe() > nsigmaTPCMax) {
        continue;
      }

      // -------------------------
      // post TPC cut
      // -------------------------
      if (isPos) {
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_Cut_DePos"), momAbs, hmpid.tpcNSigmaDe());
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DePos"), momAbs, hmpid.tofNSigmaDe());
      }
      if (isNeg) {
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_Cut_DeNeg"), momAbs, hmpid.tpcNSigmaDe());
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg"), momAbs, hmpid.tofNSigmaDe());
      }

      if (isPos && chamber == Rich2) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DePos_4cm"), momAbs, hmpid.tofNSigmaDe());
      }
      if (isPos && chamber == Rich4) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DePos_8cm"), momAbs, hmpid.tofNSigmaDe());
      }
      if (isNeg && chamber == Rich2) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_4cm"), momAbs, hmpid.tofNSigmaDe());
      }
      if (isNeg && chamber == Rich4) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_DeNeg_8cm"), momAbs, hmpid.tofNSigmaDe());
      }

      registryDA.fill(HIST("hMass_postDeuteron"), mass);

      // -------------------------
      // deuteron candidate cuts - TOF
      // -------------------------
      if (hmpid.tofNSigmaDe() < nsigmaTOFMin ||
          hmpid.tofNSigmaDe() > nsigmaTOFMax) {
        continue;
      }

      // -------------------------
      // TOF nsigma
      // -------------------------
      if (isPos) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Cut_DePos"),
                        momAbs, hmpid.tofNSigmaDe());
      }

      if (isNeg) {
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Cut_DeNeg"),
                        momAbs, hmpid.tofNSigmaDe());
      }

      // -------------------------
      // DEUTERON CANDIDATES AFTER ALL CUTS
      // -------------------------

      // Pb-Pb collisions centrality requirement
      if (collisionsPbPb) {
        if (hmpid.centrality() < centralityCut) {
          continue;
        }
      }

      // -------------------------
      // deltaR per chamber (Pos/Neg)
      // -------------------------
      if (chamber == Rich0) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_0"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_0"), dr);
        }
      }
      if (chamber == Rich1) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_1"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_1"), dr);
        }
      }
      if (chamber == Rich2) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_2"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_2"), dr);
        }
      }
      if (chamber == Rich3) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_3"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_3"), dr);
        }
      }
      if (chamber == Rich4) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_4"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_4"), dr);
        }
      }
      if (chamber == Rich5) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_5"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_5"), dr);
        }
      }
      if (chamber == Rich6) {
        if (isPos) {
          registryDA.fill(HIST("De_Pos_deltaR_6"), dr);
        }
        if (isNeg) {
          registryDA.fill(HIST("De_Neg_deltaR_6"), dr);
        }
      }

      // -------------------------
      // absorbers
      // -------------------------
      const int abs4cm = 2;
      const int abs8cm = 4;

      if (isPos && chamber == abs8cm) {
        registryDA.fill(HIST("incomingDe_Pos_8cm"), momAbs);
        registryDA.fill(HIST("survivingDe_Pos_8cm"), momAbs, dr);
        registryDA.fill(HIST("De_Pos_Q_8cm"), momAbs, hmpid.chargeMip());
        registryDA.fill(HIST("De_Pos_deltaR_8cm"), dr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_DePos_8cm"), momAbs, hmpid.tofNSigmaDe());

        registryDA.fill(HIST("residualsX_Vs_momentum_Pos_8cm"), momAbs, dx);
        registryDA.fill(HIST("residualsY_Vs_momentum_Pos_8cm"), momAbs, dy);
      }

      if (isNeg && chamber == abs8cm) {
        registryDA.fill(HIST("incomingDe_Neg_8cm"), momAbs);
        registryDA.fill(HIST("survivingDe_Neg_8cm"), momAbs, dr);
        registryDA.fill(HIST("De_Neg_Q_8cm"), momAbs, hmpid.chargeMip());
        registryDA.fill(HIST("De_Neg_deltaR_8cm"), dr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_DeNeg_8cm"), momAbs, hmpid.tofNSigmaDe());

        registryDA.fill(HIST("residualsX_Vs_momentum_Neg_8cm"), momAbs, dx);
        registryDA.fill(HIST("residualsY_Vs_momentum_Neg_8cm"), momAbs, dy);
      }

      if (isPos && chamber == abs4cm) {
        registryDA.fill(HIST("incomingDe_Pos_4cm"), momAbs);
        registryDA.fill(HIST("survivingDe_Pos_4cm"), momAbs, dr);
        registryDA.fill(HIST("De_Pos_Q_4cm"), momAbs, hmpid.chargeMip());
        registryDA.fill(HIST("De_Pos_deltaR_4cm"), dr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_DePos_4cm"), momAbs, hmpid.tofNSigmaDe());

        registryDA.fill(HIST("residualsX_Vs_momentum_Pos_4cm"), momAbs, dx);
        registryDA.fill(HIST("residualsY_Vs_momentum_Pos_4cm"), momAbs, dy);
      }

      if (isNeg && chamber == abs4cm) {
        registryDA.fill(HIST("incomingDe_Neg_4cm"), momAbs);
        registryDA.fill(HIST("survivingDe_Neg_4cm"), momAbs, dr);
        registryDA.fill(HIST("De_Neg_Q_4cm"), momAbs, hmpid.chargeMip());
        registryDA.fill(HIST("De_Neg_deltaR_4cm"), dr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_DeNeg_4cm"), momAbs, hmpid.tofNSigmaDe());

        registryDA.fill(HIST("residualsX_Vs_momentum_Neg_4cm"), momAbs, dx);
        registryDA.fill(HIST("residualsY_Vs_momentum_Neg_4cm"), momAbs, dy);
      }
    } // end of loop over hmpidTable
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HmpidDeuteron>(cfgc),
  };
}
