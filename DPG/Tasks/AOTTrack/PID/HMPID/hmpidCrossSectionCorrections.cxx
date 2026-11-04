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

struct HmpidCrossSectionCorrections {

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables for main track quality criteria
  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};

  // Configurables for QA selection
  Configurable<float> qMipMin{"qMipMin", 200.0f, "minimum MIP charge for QA selection"};
  Configurable<float> qMipMax{"qMipMax", 2000.0f, "maximum MIP charge for QA selection"};
  Configurable<float> cutDistanceMipTrack{"cutDistanceMipTrack", 10.0f, "maximum raw MIP-track distance for QA selection"};

  // Configurables for plotting
  Configurable<int> momNBins{"momNBins", 490, "Momentum bins"};
  Configurable<float> momMin{"momMin", 0.1f, "Momentum min"};
  Configurable<float> momMax{"momMax", 5.0f, "Momentum max"};

  Configurable<int> etaNBins{"etaNBins", 100, "Eta bins"};
  Configurable<float> etaMin{"etaMin", -5.0f, "Eta min"};
  Configurable<float> etaMax{"etaMax", 5.0f, "Eta max"};

  Configurable<int> phiNBins{"phiNBins", 100, "Phi bins"};
  Configurable<float> phiMin{"phiMin", -10.0f, "Phi min"};
  Configurable<float> phiMax{"phiMax", 10.0f, "Phi max"};

  Configurable<int> mapNBins{"mapNBins", 135, "Map bins (x or y)"};
  Configurable<float> mapMin{"mapMin", 0.0f, "Map min (x or y) in cm"};
  Configurable<float> mapMax{"mapMax", 135.0f, "Map max (x or y) in cm"};

  Configurable<int> residualsNBins{"residualsNBins", 400, "Residuals bins"};
  Configurable<float> residualsMin{"residualsMin", -20.0f, "Residuals min (x or y) in cm"};
  Configurable<float> residualsMax{"residualsMax", 20.0f, "Residuals max (x or y) in cm"};

  static constexpr int Rich0 = 0, Rich1 = 1, Rich2 = 2, Rich3 = 3, Rich4 = 4, Rich5 = 5, Rich6 = 6;
  static constexpr int Nch = 7;

  void init(o2::framework::InitContext const&)
  {
    // Axes for histograms
    AxisSpec axisHmpidMom{momNBins, momMin, momMax, "#it{p} (GeV/#it{c})", "axisHmpidMom"};
    AxisSpec axisGlobalMom{momNBins, momMin, momMax, "#it{p} (GeV/#it{c})", "axisGlobalMom"};
    AxisSpec axisMipX{mapNBins, mapMin, mapMax, "X_{MIP} (cm)", "axisMipX"};
    AxisSpec axisMipY{mapNBins, mapMin, mapMax, "Y_{MIP} (cm)", "axisMipY"};
    AxisSpec axisGlobalTrackX{mapNBins, mapMin, mapMax, "X_{track} (cm)", "axisGlobalTrackX"};
    AxisSpec axisGlobalTrackY{mapNBins, mapMin, mapMax, "Y_{track} (cm)", "axisGlobalTrackY"};
    AxisSpec axisEta{etaNBins, etaMin, etaMax, "#eta", "axisEta"};
    AxisSpec axisPhi{phiNBins, phiMin, phiMax, "#phi", "axisPhi"};
    AxisSpec axisResiduals{residualsNBins, residualsMin, residualsMax, "Residuals (cm)", "axisResiduals"};
    AxisSpec axisChamber{10, -1.5, 8.5, "Chamber number", "axisChamber"};

    for (int iCh = 0; iCh < Nch; iCh++) {
      registry.add(Form("hEtaRich%d", iCh), Form("hEtaRich%d", iCh), HistType::kTH1F, {axisEta});
      registry.add(Form("hPhiRich%d", iCh), Form("hPhiRich%d", iCh), HistType::kTH1F, {axisPhi});

      registry.add(Form("hmpidXYHmpidMipRich%dPos", iCh), "MIP XY - Q_{MIP} cut;X_{MIP} (cm);Y_{MIP} (cm)", HistType::kTHnSparseF, {axisHmpidMom, axisMipX, axisMipY});
      registry.add(Form("hmpidXYGlobalTrackRich%dPos", iCh), "Global track XY - all;X_{track} (cm);Y_{track} (cm)", HistType::kTHnSparseF, {axisGlobalMom, axisGlobalTrackX, axisGlobalTrackY});
      registry.add(Form("hmpidXYHmpidMipRich%dNeg", iCh), "MIP XY - Q_{MIP} cut;X_{MIP} (cm);Y_{MIP} (cm)", HistType::kTHnSparseF, {axisHmpidMom, axisMipX, axisMipY});
      registry.add(Form("hmpidXYGlobalTrackRich%dNeg", iCh), "Global track XY - all;X_{track} (cm);Y_{track} (cm)", HistType::kTHnSparseF, {axisGlobalMom, axisGlobalTrackX, axisGlobalTrackY});

      registry.add(Form("residualsXVsMomentumRich%dPos", iCh), Form("residualsXVsMomentumRich%dPos", iCh), HistType::kTH2F, {axisHmpidMom, axisResiduals});
      registry.add(Form("residualsYVsMomentumRich%dPos", iCh), Form("residualsYVsMomentumRich%dPos", iCh), HistType::kTH2F, {axisHmpidMom, axisResiduals});
      registry.add(Form("residualsXVsMomentumRich%dNeg", iCh), Form("residualsXVsMomentumRich%dNeg", iCh), HistType::kTH2F, {axisHmpidMom, axisResiduals});
      registry.add(Form("residualsYVsMomentumRich%dNeg", iCh), Form("residualsYVsMomentumRich%dNeg", iCh), HistType::kTH2F, {axisHmpidMom, axisResiduals});

      // Momentum and Position correlations
      registry.add(Form("hMomTrackVsHmpidRich%dPos", iCh), Form("hMomTrackVsHmpidRich%dPos", iCh), HistType::kTH2F, {axisGlobalMom, axisHmpidMom});
      registry.add(Form("hXTrackVsXMipRich%dPos", iCh), Form("hXTrackVsXMipRich%dPos", iCh), HistType::kTH2F, {axisGlobalTrackX, axisMipX});
      registry.add(Form("hYTrackVsYMipRich%dPos", iCh), Form("hYTrackVsYMipRich%dPos", iCh), HistType::kTH2F, {axisGlobalTrackY, axisMipY});
      registry.add(Form("hMomTrackVsHmpidRich%dNeg", iCh), Form("hMomTrackVsHmpidRich%dNeg", iCh), HistType::kTH2F, {axisGlobalMom, axisHmpidMom});
      registry.add(Form("hXTrackVsXMipRich%dNeg", iCh), Form("hXTrackVsXMipRich%dNeg", iCh), HistType::kTH2F, {axisGlobalTrackX, axisMipX});
      registry.add(Form("hYTrackVsYMipRich%dNeg", iCh), Form("hYTrackVsYMipRich%dNeg", iCh), HistType::kTH2F, {axisGlobalTrackY, axisMipY});
    } // end loop over chambers

  } // end init

  void process(aod::HmpidAnalysis const& hmpidTable)
  {
    for (const auto& hmpid : hmpidTable) {
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
      const int chamber = static_cast<int>(hmpid.chamber());

      const float momHmpid = hmpid.momentumHmpid();
      const float momTrackAbs = hmpid.momentumTrack();
      const float momHmpidAbs = std::abs(momHmpid);

      const bool isPos = momHmpid > 0;
      const bool isNeg = momHmpid < 0;
      const bool passQ = hmpid.chargeMip() > qMipMin && hmpid.chargeMip() < qMipMax;

      switch (chamber) {
        case Rich0:
          registry.fill(HIST("hEtaRich0"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich0"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich0Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich0Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich0Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich0Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich0Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich0Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich0Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich0Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich0Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich0Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich0Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich0Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich0Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich0Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich1:
          registry.fill(HIST("hEtaRich1"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich1"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich1Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich1Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich1Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich1Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich1Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich1Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich1Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich1Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich1Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich1Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich1Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich1Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich1Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich1Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich2:
          registry.fill(HIST("hEtaRich2"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich2"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich2Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich2Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich2Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich2Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich2Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich2Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich2Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich2Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich2Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich2Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich2Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich2Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich2Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich2Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich3:
          registry.fill(HIST("hEtaRich3"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich3"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich3Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich3Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich3Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich3Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich3Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich3Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich3Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich3Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich3Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich3Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich3Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich3Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich3Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich3Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich4:
          registry.fill(HIST("hEtaRich4"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich4"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich4Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich4Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich4Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich4Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich4Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich4Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich4Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich4Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich4Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich4Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich4Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich4Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich4Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich4Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich5:
          registry.fill(HIST("hEtaRich5"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich5"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich5Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich5Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich5Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich5Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich5Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich5Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich5Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich5Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich5Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich5Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich5Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich5Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich5Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich5Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        case Rich6:
          registry.fill(HIST("hEtaRich6"), hmpid.etaTrack());
          registry.fill(HIST("hPhiRich6"), hmpid.phiTrack());
          if (isPos) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich6Pos"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich6Pos"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich6Pos"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich6Pos"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich6Pos"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich6Pos"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich6Pos"), hmpid.yTrack(), hmpid.yMip());
          }
          if (isNeg) {
            if (passQ) {
              registry.fill(HIST("hmpidXYHmpidMipRich6Neg"), momHmpidAbs, hmpid.xMip(), hmpid.yMip());
            }
            registry.fill(HIST("hmpidXYGlobalTrackRich6Neg"), momHmpidAbs, hmpid.xTrack(), hmpid.yTrack());

            registry.fill(HIST("residualsXVsMomentumRich6Neg"), momHmpidAbs, hmpid.xTrack() - hmpid.xMip());
            registry.fill(HIST("residualsYVsMomentumRich6Neg"), momHmpidAbs, hmpid.yTrack() - hmpid.yMip());

            registry.fill(HIST("hMomTrackVsHmpidRich6Neg"), momTrackAbs, momHmpidAbs);
            registry.fill(HIST("hXTrackVsXMipRich6Neg"), hmpid.xTrack(), hmpid.xMip());
            registry.fill(HIST("hYTrackVsYMipRich6Neg"), hmpid.yTrack(), hmpid.yMip());
          }
          break;
        default:
          break;
      } // end switch on chamber
    } // end of loop over hmpidTable
  } // end process
}; // end struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HmpidCrossSectionCorrections>(cfgc),
  };
}
