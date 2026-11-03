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

// \brief Analysis task for the absorption cross section using HMPID (de,pi,ka,pr)
// \authors rocco.liotino@cern.ch

#include "tableHMPID.h"

#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct HmpidAbsCrossSection {

  HistogramRegistry registryDA{"registryDA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};

  // Configurables for PbPb collisions
  Configurable<bool> collisionsPbPb{"collisionsPbPb", false, "Require PbPb collisions"};
  Configurable<float> centralityCut{"centralityCut", 30.0f, "centrality cut for PbPb collisions"};

  // Configurables for PID cuts
  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.0, "nsigmaTOFMax"};

  // Particle species selection
  Configurable<std::string> particleSpecies{"particleSpecies", "de", "Particle species: de (deuteron), pr (proton), pi (pion), ka (kaon)"};

  // Configurable offset for residuals distributions
  Configurable<float> offsetX{"offsetX", 0.696712f, "Residuals offset - X direction"}; // mean value (old)
  Configurable<float> offsetY{"offsetY", 2.32849f, "Residuals offset - Y direction"};  // mean value (old)

  static constexpr float UnsetOffset = -999.f;
  static constexpr int MaxOffsetBins = 35;
  static constexpr int OffsetColumns = 8;

  static std::vector<std::string> makeOffsetRowLabels()
  {
    std::vector<std::string> labels;
    labels.reserve(MaxOffsetBins);
    for (int bin = 0; bin < MaxOffsetBins; ++bin) {
      labels.emplace_back("pBin" + std::to_string(bin));
    }
    return labels;
  }

  Configurable<float> offsetMomMin{"offsetMomMin", 0.5f, "Minimum momentum for residual offset correction"};
  Configurable<float> offsetMomMax{"offsetMomMax", 4.0f, "Maximum momentum for residual offset correction"};
  Configurable<float> offsetMomStep{"offsetMomStep", 0.5f, "Momentum bin width for residual offset correction"};

  Configurable<LabeledArray<float>> residualOffsets{"residualOffsets", LabeledArray<float>{Array2D<float>{std::vector<float>(MaxOffsetBins * OffsetColumns, UnsetOffset), MaxOffsetBins, OffsetColumns}, makeOffsetRowLabels(), {"xPosR2", "yPosR2", "xNegR2", "yNegR2", "xPosR4", "yPosR4", "xNegR4", "yNegR4"}}, "Residual offsets by momentum bin, chamber and charge"};

  // Configurables for plotting
  Configurable<int> momNBins{"momNBins", 450, "Momentum bins"};
  Configurable<float> momMin{"momMin", 0.5f, "Momentum min"};
  Configurable<float> momMax{"momMax", 5.0f, "Momentum max"};

  Configurable<int> deltaRNBins{"deltaRNBins", 100, "DeltaR bins"};
  Configurable<float> deltaRMin{"deltaRMin", 0.0f, "DeltaR min"};
  Configurable<float> deltaRMax{"deltaRMax", 50.0f, "DeltaR max"};

  Configurable<int> nSigmaTPCNBins{"nSigmaTPCNBins", 100, "nSigmaTPC bins"};
  Configurable<float> nSigmaTPCMinPlot{"nSigmaTPCMinPlot", -5.0f, "nSigmaTPC min"};
  Configurable<float> nSigmaTPCMaxPlot{"nSigmaTPCMaxPlot", 5.0f, "nSigmaTPC max"};

  Configurable<int> nSigmaTOFNBins{"nSigmaTOFNBins", 100, "nSigmaTOF bins"};
  Configurable<float> nSigmaTOFMinPlot{"nSigmaTOFMinPlot", -5.0f, "nSigmaTOF min"};
  Configurable<float> nSigmaTOFMaxPlot{"nSigmaTOFMaxPlot", 5.0f, "nSigmaTOF max"};

  static constexpr int Rich2 = 2, Rich4 = 4;
  // static constexpr Rich0 = 0, Rich1 = 1, Rich3 = 3, Rich5 = 5, Rich6 = 6;
  // static constexpr int Nch = 7;

  void init(o2::framework::InitContext const&)
  {
    // Axes for histograms
    AxisSpec axisMom{momNBins, momMin, momMax, "#it{p} (GeV/#it{c})", "axisMom"};
    AxisSpec axisDeltaR{deltaRNBins, deltaRMin, deltaRMax, "#Delta R (cm)", "axisDeltaR"};
    AxisSpec axisNsigmaTOF{nSigmaTOFNBins, nSigmaTOFMinPlot, nSigmaTOFMaxPlot, "n#sigma_{TOF}", "axisNsigmaTOF"};
    AxisSpec axisNsigmaTPC{nSigmaTPCNBins, nSigmaTPCMinPlot, nSigmaTPCMaxPlot, "n#sigma_{TPC}", "axisNsigmaTPC"};
    AxisSpec axisChamber{10, -1.5, 8.5, "Chamber number", "axisChamber"};

    // 1D integrated histograms (quick view)
    registryDA.add("incoming_Pos_8cm", "incoming_Pos_8cm", HistType::kTH1F, {axisMom});
    registryDA.add("incoming_Pos_4cm", "incoming_Pos_4cm", HistType::kTH1F, {axisMom});

    registryDA.add("DeltaR_Pos_8cm", "DeltaR_Pos_8cm", HistType::kTH1F, {axisDeltaR});
    registryDA.add("DeltaR_Pos_4cm", "DeltaR_Pos_4cm", HistType::kTH1F, {axisDeltaR});

    registryDA.add("incoming_Neg_8cm", "incoming_Neg_8cm", HistType::kTH1F, {axisMom});
    registryDA.add("incoming_Neg_4cm", "incoming_Neg_4cm", HistType::kTH1F, {axisMom});

    registryDA.add("DeltaR_Neg_8cm", "DeltaR_Neg_8cm", HistType::kTH1F, {axisDeltaR});
    registryDA.add("DeltaR_Neg_4cm", "DeltaR_Neg_4cm", HistType::kTH1F, {axisDeltaR});

    // Positive candidates histograms
    registryDA.add("nSigmaTPC_vs_momHMPID_PreCut_Pos_8cm", "nSigmaTPC_vs_momHMPID_PreCut_Pos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_PreCut_Pos_8cm", "nSigmaTOF_vs_momHMPID_PreCut_Pos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTPC_vs_momHMPID_PreCut_Pos_4cm", "nSigmaTPC_vs_momHMPID_PreCut_Pos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_PreCut_Pos_4cm", "nSigmaTOF_vs_momHMPID_PreCut_Pos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    // TPC only candidates histograms
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_Pos_8cm", "nSigmaTOF_vs_momHMPID_TPCOnly_Pos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_Pos_4cm", "nSigmaTOF_vs_momHMPID_TPCOnly_Pos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    // after TPC+TOF cuts candidates histograms
    registryDA.add("nSigmaTOF_vs_momHMPID_Pos_8cm", "nSigmaTOF_vs_momHMPID_Pos_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_Pos_4cm", "nSigmaTOF_vs_momHMPID_Pos_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("DeltaR_vs_momHMPID_Pos_8cm", "DeltaR_vs_momHMPID_Pos_8cm", HistType::kTH2F, {axisMom, axisDeltaR});
    registryDA.add("DeltaR_vs_momHMPID_Pos_4cm", "DeltaR_vs_momHMPID_Pos_4cm", HistType::kTH2F, {axisMom, axisDeltaR});

    // corrected residuals distributions
    registryDA.add("residualsX_Vs_momentum_Pos_8cm", "residualsX_Vs_momentum_Pos_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Pos_8cm", "residualsY_Vs_momentum_Pos_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});
    registryDA.add("residualsX_Vs_momentum_Pos_4cm", "residualsX_Vs_momentum_Pos_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Pos_4cm", "residualsY_Vs_momentum_Pos_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});

    // Negative candidates histograms
    registryDA.add("nSigmaTPC_vs_momHMPID_PreCut_Neg_8cm", "nSigmaTPC_vs_momHMPID_PreCut_Neg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_PreCut_Neg_8cm", "nSigmaTOF_vs_momHMPID_PreCut_Neg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTPC_vs_momHMPID_PreCut_Neg_4cm", "nSigmaTPC_vs_momHMPID_PreCut_Neg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTPC});
    registryDA.add("nSigmaTOF_vs_momHMPID_PreCut_Neg_4cm", "nSigmaTOF_vs_momHMPID_PreCut_Neg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    // TPC only candidates histograms
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_Neg_8cm", "nSigmaTOF_vs_momHMPID_TPCOnly_Neg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_TPCOnly_Neg_4cm", "nSigmaTOF_vs_momHMPID_TPCOnly_Neg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    // after TPC+TOF cuts candidates histograms
    registryDA.add("nSigmaTOF_vs_momHMPID_Neg_8cm", "nSigmaTOF_vs_momHMPID_Neg_8cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});
    registryDA.add("nSigmaTOF_vs_momHMPID_Neg_4cm", "nSigmaTOF_vs_momHMPID_Neg_4cm", HistType::kTH2F, {axisMom, axisNsigmaTOF});

    registryDA.add("DeltaR_vs_momHMPID_Neg_8cm", "DeltaR_vs_momHMPID_Neg_8cm", HistType::kTH2F, {axisMom, axisDeltaR});
    registryDA.add("DeltaR_vs_momHMPID_Neg_4cm", "DeltaR_vs_momHMPID_Neg_4cm", HistType::kTH2F, {axisMom, axisDeltaR});

    // corrected residuals distributions
    registryDA.add("residualsX_Vs_momentum_Neg_8cm", "residualsX_Vs_momentum_Neg_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Neg_8cm", "residualsY_Vs_momentum_Neg_8cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});
    registryDA.add("residualsX_Vs_momentum_Neg_4cm", "residualsX_Vs_momentum_Neg_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "x_{mip} - x_{track} (cm)"}});
    registryDA.add("residualsY_Vs_momentum_Neg_4cm", "residualsY_Vs_momentum_Neg_4cm", HistType::kTH2F, {axisMom, {400, -20.0f, 20.0f, "y_{mip} - y_{track} (cm)"}});

    // nsigma TPC vs nsigma TOF for candidates
    registryDA.add("pidTPC_TOF_p_ch_Pos", "pidTPC_TOF_p_ch_Pos", HistType::kTHnSparseF, {axisMom, axisNsigmaTPC, axisNsigmaTOF, axisChamber});
    registryDA.add("pidTPC_TOF_p_ch_Neg", "pidTPC_TOF_p_ch_Neg", HistType::kTHnSparseF, {axisMom, axisNsigmaTPC, axisNsigmaTOF, axisChamber});
  } // end init

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

      const int chamber = static_cast<int>(hmpid.chamber());
      const bool isPos = momHmpid > 0;
      const bool isNeg = momHmpid < 0;

      // Select chambers with absorber
      if (chamber != Rich2 && chamber != Rich4) {
        continue;
      }

      // Default residual offsets
      float usedOffsetX = offsetX.value;
      float usedOffsetY = offsetY.value;

      // Momentum-dependent offset correction
      if (momAbs >= offsetMomMin.value && momAbs < offsetMomMax.value && (isPos || isNeg)) {
        const int pBin = static_cast<int>((momAbs - offsetMomMin.value) / offsetMomStep.value);
        if (pBin >= 0 && pBin < static_cast<int>(residualOffsets->rows())) {
          const int col = (chamber == Rich2 ? 0 : 4) + (isPos ? 0 : 2);
          const float binnedOffsetX = residualOffsets->get(pBin, col);
          const float binnedOffsetY = residualOffsets->get(pBin, col + 1);

          if (binnedOffsetX != UnsetOffset) {
            usedOffsetX = binnedOffsetX;
          }
          if (binnedOffsetY != UnsetOffset) {
            usedOffsetY = binnedOffsetY;
          }
        }
      }

      // Corrected residuals
      const float dxCorr = (hmpid.xMip() - hmpid.xTrack()) - usedOffsetX;
      const float dyCorr = (hmpid.yMip() - hmpid.yTrack()) - usedOffsetY;
      const float drCorr = std::hypot(dxCorr, dyCorr);

      // check only HMPID chamber with an absorber (only RICH 2 (4cm abs) and RICH4 (8cm abs))
      if (chamber != Rich2 && chamber != Rich4) {
        continue;
      }

      // select particle species

      float nSigmaTPC = 0.f;
      float nSigmaTOF = 0.f;

      // Select particle species
      if (particleSpecies.value == "de") {
        nSigmaTPC = hmpid.tpcNSigmaDe();
        nSigmaTOF = hmpid.tofNSigmaDe();
      } else if (particleSpecies.value == "pr") {
        nSigmaTPC = hmpid.tpcNSigmaPr();
        nSigmaTOF = hmpid.tofNSigmaPr();
      } else if (particleSpecies.value == "pi") {
        nSigmaTPC = hmpid.tpcNSigmaPi();
        nSigmaTOF = hmpid.tofNSigmaPi();
      } else if (particleSpecies.value == "ka") {
        nSigmaTPC = hmpid.tpcNSigmaKa();
        nSigmaTOF = hmpid.tofNSigmaKa();
      } else {
        LOGF(fatal, "Invalid particle species: %s",
             particleSpecies.value.c_str());
      }

      // -------------------------
      // nsigma distributions pre-cuts
      // -------------------------
      if (isPos) {
        registryDA.fill(HIST("pidTPC_TOF_p_ch_Pos"), momAbs, nSigmaTPC, nSigmaTOF, chamber);
        if (chamber == Rich2) {
          registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_PreCut_Pos_4cm"), momAbs, nSigmaTPC);
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_PreCut_Pos_4cm"), momAbs, nSigmaTOF);
        }
        if (chamber == Rich4) {
          registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_PreCut_Pos_8cm"), momAbs, nSigmaTPC);
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_PreCut_Pos_8cm"), momAbs, nSigmaTOF);
        }
      }
      if (isNeg) {
        registryDA.fill(HIST("pidTPC_TOF_p_ch_Neg"), momAbs, nSigmaTPC, nSigmaTOF, chamber);
        if (chamber == Rich2) {
          registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_PreCut_Neg_4cm"), momAbs, nSigmaTPC);
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_PreCut_Neg_4cm"), momAbs, nSigmaTOF);
        }
        if (chamber == Rich4) {
          registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_PreCut_Neg_8cm"), momAbs, nSigmaTPC);
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_PreCut_Neg_8cm"), momAbs, nSigmaTOF);
        }
      }

      // -------------------------
      // deuteron candidate cuts - TPC
      // -------------------------
      if (nSigmaTPC < nsigmaTPCMin || nSigmaTPC > nsigmaTPCMax) {
        continue;
      }

      // -------------------------
      // post TPC cut
      // -------------------------
      if (isPos) {
        if (chamber == Rich2) {
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_Pos_8cm"), momAbs, nSigmaTOF);
        }
        if (chamber == Rich4) {
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_Pos_4cm"), momAbs, nSigmaTOF);
        }
      }
      if (isNeg) {
        if (chamber == Rich2) {
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_Neg_8cm"), momAbs, nSigmaTOF);
        }
        if (chamber == Rich4) {
          registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_TPCOnly_Neg_4cm"), momAbs, nSigmaTOF);
        }
      }

      // -------------------------
      // deuteron candidate cuts - TOF
      // -------------------------
      if (nSigmaTOF < nsigmaTOFMin || nSigmaTOF > nsigmaTOFMax) {
        continue;
      }

      // Pb-Pb collisions centrality requirement
      if (collisionsPbPb) {
        if (hmpid.centrality() < centralityCut) {
          continue;
        }
      }

      // Fill distributions for absorption cross section measurement
      if (isPos && chamber == Rich4) {
        registryDA.fill(HIST("incoming_Pos_8cm"), momAbs);
        registryDA.fill(HIST("DeltaR_Pos_8cm"), drCorr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Pos_8cm"), momAbs, nSigmaTOF);
        registryDA.fill(HIST("DeltaR_vs_momHMPID_Pos_8cm"), momAbs, drCorr);
        registryDA.fill(HIST("residualsX_Vs_momentum_Pos_8cm"), momAbs, dxCorr);
        registryDA.fill(HIST("residualsY_Vs_momentum_Pos_8cm"), momAbs, dyCorr);
      }

      if (isNeg && chamber == Rich4) {
        registryDA.fill(HIST("incoming_Neg_8cm"), momAbs);
        registryDA.fill(HIST("DeltaR_Neg_8cm"), drCorr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Neg_8cm"), momAbs, nSigmaTOF);
        registryDA.fill(HIST("DeltaR_vs_momHMPID_Neg_8cm"), momAbs, drCorr);

        registryDA.fill(HIST("residualsX_Vs_momentum_Neg_8cm"), momAbs, dxCorr);
        registryDA.fill(HIST("residualsY_Vs_momentum_Neg_8cm"), momAbs, dyCorr);
      }

      if (isPos && chamber == Rich2) {
        registryDA.fill(HIST("incoming_Pos_4cm"), momAbs);
        registryDA.fill(HIST("DeltaR_Pos_4cm"), drCorr);
        registryDA.fill(HIST("DeltaR_vs_momHMPID_Pos_4cm"), momAbs, drCorr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Pos_4cm"), momAbs, nSigmaTOF);

        registryDA.fill(HIST("residualsX_Vs_momentum_Pos_4cm"), momAbs, dxCorr);
        registryDA.fill(HIST("residualsY_Vs_momentum_Pos_4cm"), momAbs, dyCorr);
      }

      if (isNeg && chamber == Rich2) {
        registryDA.fill(HIST("incoming_Neg_4cm"), momAbs);
        registryDA.fill(HIST("DeltaR_Neg_4cm"), drCorr);
        registryDA.fill(HIST("DeltaR_vs_momHMPID_Neg_4cm"), momAbs, drCorr);
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Neg_4cm"), momAbs, nSigmaTOF);

        registryDA.fill(HIST("residualsX_Vs_momentum_Neg_4cm"), momAbs, dxCorr);
        registryDA.fill(HIST("residualsY_Vs_momentum_Neg_4cm"), momAbs, dyCorr);
      }
    } // end of loop over hmpidTable
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HmpidAbsCrossSection>(cfgc),
  };
}
