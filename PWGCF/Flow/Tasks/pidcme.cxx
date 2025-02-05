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

/// \author ZhengqingWang(zhengqing.wang@cern.ch)
/// \file   pidcme.cxx
/// \brief  task to calculate the pikp cme signal and bacground.
// C++/ROOT includes.
// o2-linter: disable=name/workflow-file
#include <chrono>
#include <string>
#include <vector>
#include <utility>
#include <TF1.h>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>

// o2Physics includes.
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"

#include "CommonConstants/PhysicsConstants.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace cme_track_pid_columns
{
DECLARE_SOA_COLUMN(NPidFlag, nPidFlag, int8_t); // Flag tracks without proper binning as -1, and indicate type of particle 0->un-Id, 1->pion, 2->kaon, 3->proton
DECLARE_SOA_COLUMN(AverClusterSizeCosl, averClusterSizeCosl, float);
DECLARE_SOA_COLUMN(NSigmaPiITS, nSigmaPiITS, float);
DECLARE_SOA_COLUMN(NSigmaKaITS, nSigmaKaITS, float);
DECLARE_SOA_COLUMN(NSigmaPrITS, nSigmaPrITS, float);
DECLARE_SOA_COLUMN(NSigmaPiTPC, nSigmaPiTPC, float);
DECLARE_SOA_COLUMN(NSigmaKaTPC, nSigmaKaTPC, float);
DECLARE_SOA_COLUMN(NSigmaPrTPC, nSigmaPrTPC, float);
} // namespace cme_track_pid_columns
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", cme_track_pid_columns::NPidFlag, cme_track_pid_columns::AverClusterSizeCosl, cme_track_pid_columns::NSigmaPiITS, cme_track_pid_columns::NSigmaKaITS, cme_track_pid_columns::NSigmaPrITS, cme_track_pid_columns::NSigmaPiTPC, cme_track_pid_columns::NSigmaKaTPC, cme_track_pid_columns::NSigmaPrTPC);
} // namespace o2::aod

using TracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
using CollisionPID = soa::Join<aod::Collisions, aod::CentFT0Cs>;
struct FillPIDcolums {
  Configurable<float> cfgPtMaxforTPCOnlyPID{"cfgPtMaxforTPCOnlyPID", 0.4, "Maxmium track pt for TPC only PID,only when onlyTOF and onlyTOFHIT closed"};
  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.8, "Maximum track #eta for PID"};
  Configurable<float> cfgAveClusSizeCoslMinPi{"cfgAveClusSizeCoslMinPi", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Pions"};
  Configurable<float> cfgAveClusSizeCoslMaxPi{"cfgAveClusSizeCoslMaxPi", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Pions"};
  Configurable<float> cfgAveClusSizeCoslMinKa{"cfgAveClusSizeCoslMinKa", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Kaons"};
  Configurable<float> cfgAveClusSizeCoslMaxKa{"cfgAveClusSizeCoslMaxKa", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Kaons"};
  Configurable<float> cfgAveClusSizeCoslMinPr{"cfgAveClusSizeCoslMinPr", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Protons"};
  Configurable<float> cfgAveClusSizeCoslMaxPr{"cfgAveClusSizeCoslMaxPr", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Protons"};

  ConfigurableAxis cfgrigidityBins{"cfgrigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis cfgdedxBins{"cfgdedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsTOF{"cfgnSigmaBinsTOF", {200, -5.f, 5.f}, "Binning for n sigma TOF"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma ITS"};
  ConfigurableAxis cfgnSigmaBinsCom{"cfgnSigmaBinsCom", {100, 0.f, 10.f}, "Combination Binning for TPC&TOF nsigma"};
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgaxispPID{"cfgaxispPID", {50, 0, 5}, "Binning for P PID"};
  ConfigurableAxis cfgaxisAverClusterCosl{"cfgaxisAverClusterCosl", {50, 0, 10}, "Binning for average cluster size x cos(#lambda)"};
  ConfigurableAxis cfgaxisAverClusterCoslnSigma{"cfgaxisAverClusterCoslnSigma", {50, 0, 5}, "Binning for average cluster size x cos(#lambda) vs nSigam"};
  ConfigurableAxis cfgaxisetaPID{"cfgaxisetaPID", {90, -0.9, 0.9}, "Binning for Pt QA"};

  Configurable<bool> cfgOpenITSCut{"cfgOpenITSCut", true, "open ITSnsigma cut"};
  Configurable<bool> cfgOpenDetailPlots{"cfgOpenDetailPlots", true, "open detail TH3D plots for nSigmaTPC-ITS Pt-eta-Phi nSigmaITS-clustersize"};

  Configurable<std::vector<float>> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", {3, 3, 3}, "TPC cut for pi k p respectively at low pt"};
  Configurable<std::vector<float>> cfgnSigmaCutRMS{"cfgnSigmaCutRMS", {3, 3, 3}, "TPC_TOF combined cut for pi k p respectively at high pt"};
  Configurable<std::vector<float>> cfgnSigmaCutITS{"cfgnSigmaCutITS", {3, 2.5, 2}, "TPC_TOF combined cut for pi k p respectively at high pt"};

  static float averageClusterSizeCosl(uint32_t itsClusterSizes, float eta)
  {
    float average = 0;
    int nclusters = 0;
    const float cosl = 1. / std::cosh(eta);

    for (int layer = 0; layer < 7; layer++) {
      if ((itsClusterSizes >> (layer * 4)) & 0xf) {
        nclusters++;
        average += (itsClusterSizes >> (layer * 4)) & 0xf;
      }
    }
    if (nclusters == 0) {
      return 0;
    }
    return average * cosl / nclusters;
  };

  template <typename TrackType>
  bool selTrackPid(const TrackType track)
  {
    if (!(track.pt() > cfgMinPtPID))
      return false;
    if (!(std::abs(track.eta()) < cfgMaxEtaPID))
      return false;
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;
    return true;
  }

  template <typename T>
  int selectionPidtpctof(const T& candidate)
  {
    float averClusSizeCosl = averageClusterSizeCosl(candidate.itsClusterSizes(), candidate.eta());
    std::array<float, 3> nSigmaTPC = {candidate.tpcNSigmaPi(), candidate.tpcNSigmaKa(), candidate.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(candidate.tpcNSigmaPi(), candidate.tofNSigmaPi()), std::hypot(candidate.tpcNSigmaKa(), candidate.tofNSigmaKa()), std::hypot(candidate.tpcNSigmaPr(), candidate.tofNSigmaPr())};
    // Choose which nSigma to use
    int pid = -1;
    std::vector<float> pidVector = (candidate.pt() > cfgPtMaxforTPCOnlyPID && candidate.hasTOF()) ? cfgnSigmaCutRMS.value : cfgnSigmaCutTPC.value;
    float nsigma = pidVector[0];
    std::array<float, 3> nSigmaToUse = (candidate.pt() > cfgPtMaxforTPCOnlyPID && candidate.hasTOF()) ? nSigmaCombined : nSigmaTPC;
    // Fill cross pid QA
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < pidVector[i]) {
        if (i == 0) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Pi"), candidate.tpcNSigmaPi(), candidate.tofNSigmaPi(), candidate.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pi"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Pi"), candidate.tpcNSigmaPi());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pi"), candidate.pt(), candidate.tpcNSigmaPi());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pi"), candidate.itsNSigmaPi());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pi"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), candidate.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), candidate.tpcNSigmaPi(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pi"), candidate.phi());
          histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pi"), candidate.tpcNSigmaPi(), candidate.itsNSigmaPi());
        }
        if (i == 1) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Ka"), candidate.tpcNSigmaKa(), candidate.tofNSigmaKa(), candidate.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Ka"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Ka"), candidate.tpcNSigmaKa());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Ka"), candidate.pt(), candidate.tpcNSigmaKa());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Ka"), candidate.itsNSigmaKa());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Ka"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), candidate.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), candidate.tpcNSigmaKa(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Ka"), candidate.phi());
          histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Ka"), candidate.tpcNSigmaKa(), candidate.itsNSigmaKa());
        }
        if (i == 2) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Pr"), candidate.tpcNSigmaPr(), candidate.tofNSigmaPr(), candidate.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pr"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Pr"), candidate.tpcNSigmaPr());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pr"), candidate.pt(), candidate.tpcNSigmaPr());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pr"), candidate.itsNSigmaPr());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pr"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), candidate.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), candidate.tpcNSigmaKa(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pr"), candidate.phi());
          histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pr"), candidate.tpcNSigmaPr(), candidate.itsNSigmaPr());
        }
      }
    }
    // Select particle with the lowest nsigma
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma && std::abs(nSigmaToUse[i]) < pidVector[i]) {
        pid = i;
        nsigma = std::abs(nSigmaToUse[i]);
      }
    }
    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  HistogramRegistry histosQA{"histosQAPID", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    AxisSpec axisRigidity{cfgrigidityBins, "#it{p}^{TPC}/#it{z}"};
    AxisSpec axisdEdx{cfgdedxBins, "d#it{E}/d#it{x}"};
    AxisSpec axisnSigmaTPC{cfgnSigmaBinsTPC, "n_{#sigma}TPC"};
    AxisSpec axisnSigmaTOF{cfgnSigmaBinsTOF, "n_{#sigma}TOF"};
    AxisSpec axisnSigmaITS{cfgnSigmaBinsITS, "n_{#sigma}ITS"};
    AxisSpec axisnSigmaCom{cfgnSigmaBinsCom, "hypot(n_{#sigma}TPC,TOF)"};
    AxisSpec axisPtPID{cfgaxisptPID, "#it{p}_{T}"};
    AxisSpec axisPPID{cfgaxispPID, "#it{p}"};
    AxisSpec axisClusterSize{cfgaxisAverClusterCosl, "<ITS Cluster Size> x <cos(#lambda)>"};
    AxisSpec axisClusterSizenSigma{cfgaxisAverClusterCoslnSigma, "<ITS Cluster Size> x <cos(#lambda)>"};
    AxisSpec axisPhi = {100, 0, 2.1 * constants::math::PI, "#phi"};
    // TH3D NSigmaTPC,NSigmaTOF,pt
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Pi"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Ka"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Pr"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Pi"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Ka"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Pr"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_AfterITS_Pi"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_AfterITS_Ka"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_AfterITS_Pr"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    // Hist for PID Averge Cluster Size ITS related
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisnSigmaITS}});
    // Hist for Nsigma TPC TOF
    histosQA.add(Form("QA/PID/histdEdxTPC_All"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Pi"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Pi"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Ka"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Ka"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Pr"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Pr"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_com_Pi"), "", {HistType::kTH1F, {axisnSigmaCom}});
    histosQA.add(Form("QA/PID/histnSigma_com_Ka"), "", {HistType::kTH1F, {axisnSigmaCom}});
    histosQA.add(Form("QA/PID/histnSigma_com_Pr"), "", {HistType::kTH1F, {axisnSigmaCom}});
    histosQA.add(Form("QA/PID/histnSigma_TOF_Pi"), "", {HistType::kTH1F, {axisnSigmaTOF}});
    histosQA.add(Form("QA/PID/histnSigma_TOF_Ka"), "", {HistType::kTH1F, {axisnSigmaTOF}});
    histosQA.add(Form("QA/PID/histnSigma_TOF_Pr"), "", {HistType::kTH1F, {axisnSigmaTOF}});
    histosQA.add(Form("QA/PID/histdEdxTPC_cross_Pi"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_cross_Pi"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histdEdxTPC_cross_Ka"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_cross_Ka"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histdEdxTPC_cross_Pr"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_cross_Pr"), "", {HistType::kTH1F, {axisnSigmaTPC}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
    // Hist for nSigma ITS
    histosQA.add(Form("QA/PID/histnSigma_ITS_Pi"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_Ka"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_Pr"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Pi"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Ka"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Pr"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_AfterITS_Pi"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_AfterITS_Ka"), "", {HistType::kTH1F, {axisnSigmaITS}});
    histosQA.add(Form("QA/PID/histnSigma_ITS_AfterITS_Pr"), "", {HistType::kTH1F, {axisnSigmaITS}});
    // Hist for checking the PID phi distribution
    histosQA.add(Form("QA/PID/histPhi_Dis_Pi"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_Ka"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_Pr"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pi"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Ka"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pr"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Pi"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Ka"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Pr"), "", {HistType::kTH1F, {axisPhi}});
    // Hist 3D for PID check
    if (cfgOpenDetailPlots) {
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPi"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPi"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosKa"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegKa"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPr"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPr"), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhi, axisPtPID, cfgaxisetaPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPi_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosKa_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPr_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPi_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosKa_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPr_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPi_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegKa_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPr_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPi_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegKa_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
      histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPr_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
    }
  }
  Produces<aod::Flags> pidCmeTable;
  void process(TracksPID const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<TracksPID,
                                        aod::pidits::ITSNSigmaPi,
                                        aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr>(tracks);
    int8_t pidFlag;
    for (const auto& track : tracksWithITSPid) {
      float averClusSizeCosl = averageClusterSizeCosl(track.itsClusterSizes(), track.eta());
      if (!selTrackPid(track)) {
        pidFlag = -1;
      } else {
        histosQA.fill(HIST("QA/PID/histdEdxTPC_All"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
        pidFlag = selectionPidtpctof(track);
        switch (pidFlag) { // no cross pid track QA plots
          case 1:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Pi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pi"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pi"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Pi"), track.tpcNSigmaPi(), track.itsNSigmaPi());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pi"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pi"), track.tpcNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pi"), track.pt(), track.tpcNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Pi"), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pi"), track.tofNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pi"), track.itsNSigmaPi());
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Pi"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPi_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPi_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
          case 2:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Ka"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Ka"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Ka"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Ka"), track.tpcNSigmaKa(), track.itsNSigmaKa());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Ka"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Ka"), track.tpcNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Ka"), track.pt(), track.tpcNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Ka"), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Ka"), track.tofNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Ka"), track.itsNSigmaKa());
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Ka"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosKa_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegKa_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
          case 3:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pr"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pr"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), track.tpcNSigmaPr(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Pr"), track.tpcNSigmaPr(), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pr"), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pr"), track.pt(), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Pr"), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pr"), track.tofNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pr"), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Pr"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPr_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPr_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
        }
        if (cfgOpenITSCut) {
          switch (pidFlag) {
            case 1:
              if (!(std::abs(track.itsNSigmaPi()) < cfgnSigmaCutITS.value[0] && averClusSizeCosl > cfgAveClusSizeCoslMinPi && averClusSizeCosl < cfgAveClusSizeCoslMaxPi)) {
                pidFlag = 4;
              }
              break;
            case 2:
              if (!(std::abs(track.itsNSigmaKa()) < cfgnSigmaCutITS.value[1] && averClusSizeCosl > cfgAveClusSizeCoslMinKa && averClusSizeCosl < cfgAveClusSizeCoslMaxKa)) {
                pidFlag = 5;
              }
              break;
            case 3:
              if (!(std::abs(track.itsNSigmaPr()) < cfgnSigmaCutITS.value[2] && averClusSizeCosl > cfgAveClusSizeCoslMinPr && averClusSizeCosl < cfgAveClusSizeCoslMaxPr)) {
                pidFlag = 6;
              }
              break;
          }
        }
        switch (pidFlag) { // after ITScut QA plots
          case 1:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_AfterITS_Pi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Pi"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pi"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pi"), track.tpcNSigmaPi(), track.itsNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_AfterITS_Pi"), track.itsNSigmaPi());
            histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Pi"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_PosPi"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPi_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_NegPi"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPi_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
          case 2:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_AfterITS_Ka"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_AfterITS_Ka"), track.itsNSigmaKa());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Ka"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Ka"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Ka"), track.tpcNSigmaKa(), track.itsNSigmaKa());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Ka"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_PosKa"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosKa_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_NegKa"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegKa_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
          case 3:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_AfterITS_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_AfterITS_Pr"), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Pr"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pr"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pr"), track.tpcNSigmaPr(), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pr"), track.tpcNSigmaPr(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Pr"), track.phi());
            if (cfgOpenDetailPlots) {
              if (track.sign() > 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_PosPr"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPr_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              } else if (track.sign() < 0) {
                histosQA.fill(HIST("QA/PID/histPhi_Pt_Eta_NegPr"), track.phi(), track.pt(), track.eta());
                histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPr_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
              }
            }
            break;
        }
      }
      pidCmeTable(pidFlag, averClusSizeCosl, track.itsNSigmaPi(), track.itsNSigmaKa(), track.itsNSigmaPr(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr());
    }
  }
};

struct QAProcessCent {
  HistogramRegistry histosQA{"histosQAwithcent", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::vector<int>> cfgCentralitybinsforQA{"cfgCentralitybinsforQA", {0, 30, 60}, "Centrality bins for track phi and TPC_ITS matching check"};
  Configurable<bool> cfgOpenDetailPlotscme{"cfgOpenDetailPlotscme", true, "open detail TH3D plots for nSigmaTPC-ITS Pt-eta-Phi nSigmaITS-clustersize"};
  Configurable<bool> cfgOpenPtEtaPhi{"cfgOpenPtEtaPhi", true, "open pt-#eta-#phi PID QA  (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenITSTPCnSigma{"cfgOpenITSTPCnSigma", true, "open ITS-TPC nSigma QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenClusSizenSigmaTPC{"cfgOpenClusSizenSigmaTPC", true, "open ITSClustersize-TPCnsigma QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenPi{"cfgOpenPi", true, "open Pion QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenKa{"cfgOpenKa", true, "open Kaon QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenPr{"cfgOpenPr", true, "open Proton QA (Optional for limited memory usage)"};
  ConfigurableAxis cfgaxisetaPIDcme{"cfgaxisetaPIDcme", {90, -0.9, 0.9}, "Binning for Pt QA"};
  ConfigurableAxis cfgaxisptPIDcme{"cfgaxisptPIDcme", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgnSigmaBinsTPCcme{"cfgnSigmaBinsTPCcme", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsITScme{"cfgnSigmaBinsITScme", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgaxisAverClusterCoslnSigmacme{"cfgaxisAverClusterCoslnSigmacme", {50, 0, 5}, "Binning for average cluster size x cos(#lambda) vs nSigam"};

  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaPosPiCen;
  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaNegPiCen;
  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaPosKaCen;
  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaNegKaCen;
  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaPosPrCen;
  std::vector<std::shared_ptr<TH3>> vhistPhiPtEtaNegPrCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtPosPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaITSTPCPtNegPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtPosPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistAverClusterSizeCoslnSigmaTPCPtNegPrAfterCen;
  Filter trackPIDfilter = aod::cme_track_pid_columns::nPidFlag > (int8_t)0;
  void init(InitContext const&)
  {
    AxisSpec axisPhicme = {100, 0, 2.1 * constants::math::PI, "#phi"};
    // Additional QA histograms for PID
    if (cfgOpenDetailPlotscme) {
      for (int i = 0; i < static_cast<int>(cfgCentralitybinsforQA.value.size()) - 1; ++i) {
        if (cfgOpenPtEtaPhi) {
          if (cfgOpenPi) {
            auto hPhiPtEtaPosPi = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPi_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            auto hPhiPtEtaNegPi = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPi_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            vhistPhiPtEtaPosPiCen.push_back(std::move(hPhiPtEtaPosPi));
            vhistPhiPtEtaNegPiCen.push_back(std::move(hPhiPtEtaNegPi));
          }
          if (cfgOpenKa) {
            auto hPhiPtEtaPosKa = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosKa_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            auto hPhiPtEtaNegKa = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegKa_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            vhistPhiPtEtaPosKaCen.push_back(std::move(hPhiPtEtaPosKa));
            vhistPhiPtEtaNegKaCen.push_back(std::move(hPhiPtEtaNegKa));
          }
          if (cfgOpenPr) {
            auto hPhiPtEtaPosPr = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPr_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            auto hPhiPtEtaNegPr = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPr_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPIDcme, cfgaxisetaPIDcme}});
            vhistPhiPtEtaPosPrCen.push_back(std::move(hPhiPtEtaPosPr));
            vhistPhiPtEtaNegPrCen.push_back(std::move(hPhiPtEtaNegPr));
          }
        }
        if (cfgOpenITSTPCnSigma) {
          if (cfgOpenPi) {
            auto hnSigmaITSTPCPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            vhistnSigmaITSTPCPtPosPiBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosPiBefore));
            vhistnSigmaITSTPCPtNegPiBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegPiBefore));
            vhistnSigmaITSTPCPtPosPiAfterCen.push_back(std::move(hnSigmaITSTPCPtPosPiAfter));
            vhistnSigmaITSTPCPtNegPiAfterCen.push_back(std::move(hnSigmaITSTPCPtNegPiAfter));
          }
          if (cfgOpenKa) {
            auto hnSigmaITSTPCPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            vhistnSigmaITSTPCPtPosKaBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosKaBefore));
            vhistnSigmaITSTPCPtNegKaBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegKaBefore));
            vhistnSigmaITSTPCPtPosKaAfterCen.push_back(std::move(hnSigmaITSTPCPtPosKaAfter));
            vhistnSigmaITSTPCPtNegKaAfterCen.push_back(std::move(hnSigmaITSTPCPtNegKaAfter));
          }
          if (cfgOpenPr) {
            auto hnSigmaITSTPCPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            auto hnSigmaITSTPCPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgnSigmaBinsITScme, cfgaxisptPIDcme}});
            vhistnSigmaITSTPCPtPosPrBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosPrBefore));
            vhistnSigmaITSTPCPtNegPrBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegPrBefore));
            vhistnSigmaITSTPCPtPosPrAfterCen.push_back(std::move(hnSigmaITSTPCPtPosPrAfter));
            vhistnSigmaITSTPCPtNegPrAfterCen.push_back(std::move(hnSigmaITSTPCPtNegPrAfter));
          }
        }
        if (cfgOpenClusSizenSigmaTPC) {
          if (cfgOpenPi) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPiBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPiBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtPosPiAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPiAfter));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPiAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPiAfter));
          }
          if (cfgOpenKa) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosKaBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegKaBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtPosKaAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosKaAfter));
            vhistAverClusterSizeCoslnSigmaTPCPtNegKaAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegKaAfter));
          }
          if (cfgOpenPr) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPCcme, cfgaxisAverClusterCoslnSigmacme, cfgaxisptPIDcme}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPrBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPrBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtPosPrAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPrAfter));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPrAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPrAfter));
          }
        }
      }
    }
  }
  void process(CollisionPID::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::Flags>> const& tracks)
  {
    if (cfgOpenDetailPlotscme) {
      const auto cent = collision.centFT0C();
      int currentBin = -1;
      for (int i = 0; i < static_cast<int>(cfgCentralitybinsforQA.value.size()) - 1; ++i) {
        if (cent >= cfgCentralitybinsforQA.value[i] && cent < cfgCentralitybinsforQA.value[i + 1]) {
          currentBin = i;
          break;
        }
      }
      if (currentBin >= 0) {
        for (const auto& trk : tracks) {
          int8_t pididx = trk.nPidFlag();
          switch (pididx) {
            case 1:
              if (cfgOpenPi) {
                if (trk.sign() > 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosPiCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                    vhistnSigmaITSTPCPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegPiCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                    vhistnSigmaITSTPCPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
            case 2:
              if (cfgOpenKa) {
                if (trk.sign() > 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosKaCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                    vhistnSigmaITSTPCPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegKaCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                    vhistnSigmaITSTPCPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
            case 3:
              if (cfgOpenPr) {
                if (trk.sign() > 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosPrCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                    vhistnSigmaITSTPCPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegPrCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                    vhistnSigmaITSTPCPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
            case 4:
              if (cfgOpenPi) {
                if (trk.sign() > 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
            case 5:
              if (cfgOpenKa) {
                if (trk.sign() > 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
            case 6:
              if (cfgOpenPr) {
                if (trk.sign() > 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                } else if (trk.sign() < 0) {
                  if (cfgOpenITSTPCnSigma) {
                    vhistnSigmaITSTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                  }
                  if (cfgOpenClusSizenSigmaTPC) {
                    vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                  }
                }
              }
              break;
          }
        }
      }
    }
  }
};

struct pidcme { // o2-linter: disable=name/struct
  HistogramRegistry histosQA{"histosmain", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low boundary cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High boundary cut on TPC occupancy"};

  Configurable<float> cfgVtzCut{"cfgVtzCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Centrality min"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Centrality max"};
  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.1, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0, "Maximum longitudinal DCA"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {100, 0, 100}, ""};

  ConfigurableAxis cfgaxiscos{"cfgaxiscos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, ""};
  ConfigurableAxis cfgaxisCentMerged{"cfgaxisCentMerged", {20, 0, 100}, ""};
  ConfigurableAxis cfgaxisCentForQA{"cfgaxisCentForQA", {100, 0, 100}, "centrality for event QA"};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis cfgaxisT0C{"cfgaxisT0C", {70, 0, 70000}, "N_{ch} (T0C)"};
  ConfigurableAxis cfgaxisT0A{"cfgaxisT0A", {200, 0, 200000}, "N_{ch} (T0A)"};
  ConfigurableAxis cfgaxisNchPV{"cfgaxisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};

  ConfigurableAxis cfgaxissumpt{"cfgaxissumpt", {7, 1, 8}, "Binning for #gamma and #delta pt(particle1 + particle2)"};
  ConfigurableAxis cfgaxisdeltaeta{"cfgaxisdeltaeta", {5, 0, 1}, "Binning for #gamma and #delta |#eta(particle1 - particle2)|"};

  Configurable<bool> cfgUseAdditionalEventCut{"cfgUseAdditionalEventCut", true, "Use additional event cut beyond sel8"};
  Configurable<bool> cfgOpenEvSelkIsGoodZvtxFT0vsPV{"cfgOpenEvSelkIsGoodZvtxFT0vsPV", true, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution"};
  Configurable<bool> cfgOpenEvSelkNoSameBunchPileup{"cfgOpenEvSelkNoSameBunchPileup", true, "rejects collisions which are associated with the same found-by-T0 bunch crossing"};
  Configurable<bool> cfgOpenEvSelkNoCollInTimeRangeStandard{"cfgOpenEvSelkNoCollInTimeRangeStandard", true, "no collisions in specified time range"};
  Configurable<bool> cfgOpenEvSelkIsGoodITSLayersAll{"cfgOpenEvSelkIsGoodITSLayersAll", true, "cut time intervals with dead ITS staves"};
  Configurable<bool> cfgOpenEvSelkNoCollInRofStandard{"cfgOpenEvSelkNoCollInRofStandard", true, "no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
  Configurable<bool> cfgOpenEvSelkNoHighMultCollInPrevRof{"cfgOpenEvSelkNoHighMultCollInPrevRof", true, "veto an event if FT0C amplitude in previous ITS ROF is above threshold"};
  Configurable<bool> cfgOpenEvSelOccupancy{"cfgOpenEvSelOccupancy", true, "Occupancy cut"};
  Configurable<bool> cfgOpenEvSelMultCorrelation{"cfgOpenEvSelMultCorrelation", true, "Multiplicity correlation cut"};
  Configurable<bool> cfgOpenEvSelV0AT0ACut{"cfgOpenEvSelV0AT0ACut", true, "V0A T0A 5 sigma cut"};
  Configurable<bool> cfgOpenFullEventQA{"cfgOpenFullEventQA", true, "Open full QA plots for event QA"};
  Configurable<bool> cfgkOpenCME{"cfgkOpenCME", true, "open PID CME"};
  Configurable<bool> cfgkOpenPiPi{"cfgkOpenPiPi", true, "open Pi-Pi"};
  Configurable<bool> cfgkOpenKaKa{"cfgkOpenKaKa", false, "open Ka-Ka"};
  Configurable<bool> cfgkOpenPrPr{"cfgkOpenPrPr", false, "open Pr-Pr"};
  Configurable<bool> cfgkOpenPiKa{"cfgkOpenPiKa", true, "open Pi-Ka"};
  Configurable<bool> cfgkOpenPiPr{"cfgkOpenPiPr", true, "open Pi-Pr"};
  Configurable<bool> cfgkOpenKaPr{"cfgkOpenKaPr", true, "open Ka-Pr"};
  Configurable<bool> cfgkOpenHaHa{"cfgkOpenHaHa", false, "open Ha-Ha"};
  Configurable<bool> cfgkOpenSsOsCrossCheck{"cfgkOpenSsOsCrossCheck", false, "open check for matter an antimatter #gamma#delta"};

  EventPlaneHelper helperEP;
  SliceCache cache;

  unsigned int mult1, mult2, mult3;
  int detId;
  int refAId;
  int refBId;
  // Additional Event selection cuts - Copy from flowGenericFramework.cxx
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  template <typename T>
  int getDetId(const T& name)
  {
    if (name.value == "BPos" || name.value == "BNeg" || name.value == "BTot") {
      LOGF(warning, "Using deprecated label: %s. Please use TPCpos, TPCneg, TPCall instead.", name.value);
    }
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos" || name.value == "BPos") {
      return 4;
    } else if (name.value == "TPCneg" || name.value == "BNeg") {
      return 5;
    } else if (name.value == "TPCall" || name.value == "BTot") {
      return 6;
    } else {
      return 0;
    }
  }

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVtzCut) && (aod::cent::centFT0C > cfgCentMin) && (aod::cent::centFT0C < cfgCentMax);
  Filter ptfilter = aod::track::pt > cfgMinPt;
  Filter etafilter = aod::track::eta < cfgMaxEta;
  Filter properPIDfilter = aod::cme_track_pid_columns::nPidFlag > (int8_t)0;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet1 = aod::cme_track_pid_columns::nPidFlag == 1;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet2 = aod::cme_track_pid_columns::nPidFlag == 2;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet3 = aod::cme_track_pid_columns::nPidFlag == 3;
  // Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet4 = aod::cme_track_pid_columns::nPidFlag == 4;
  // Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet5 = aod::cme_track_pid_columns::nPidFlag == 5;
  // Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet6 = aod::cme_track_pid_columns::nPidFlag == 6;
  void init(InitContext const&)
  {

    detId = getDetId(cfgDetName);
    refAId = getDetId(cfgRefAName);
    refBId = getDetId(cfgRefBName);

    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }

    if (cfgUseAdditionalEventCut) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x - 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x + 3.5*([5]+[6]*x+[7]*x*x+[8]*x*x*x+[9]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(3257.29, -121.848, 1.98492, -0.0172128, 6.47528e-05, 154.756, -1.86072, -0.0274713, 0.000633499, -3.37757e-06);

      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutLow->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultCutHigh->SetParameters(1654.46, -47.2379, 0.449833, -0.0014125, 150.773, -3.67334, 0.0530503, -0.000614061, 3.15956e-06);

      fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgaxiscos, "angle function"};
    AxisSpec axisPt{cfgaxispt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgaxisCentMerged, "merged centrality for cme and PID v2"};

    AxisSpec axissumpt{cfgaxissumpt, "#it{p}_{T}^{sum}}"};
    AxisSpec axisdeltaeta{cfgaxisdeltaeta, "#Delta#eta"};
    AxisSpec axisvertexz = {100, -15., 15., "vrtx_{Z} [cm]"};

    histosQA.add(Form("QA/histEventCount"), "", {HistType::kTH1F, {{3, 0.0, 3.0}}});
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(3, "after additional event cut");
    if (cfgUseAdditionalEventCut) {
      histosQA.add(Form("QA/histEventCountDetail"), "Number of Event;; Count", {HistType::kTH1F, {{10, 0, 10}}});
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(1, "after sel8");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(8, "occupancy");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(9, "MultCorrelation");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(10, "cfgEvSelV0AT0ACut");
    }
    if (cfgOpenFullEventQA) {
      histosQA.add("QA/hist_globalTracks_centT0C_before", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNch}});
      histosQA.add("QA/hist_PVTracks_centT0C_before", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNchPV}});
      histosQA.add("QA/hist_globalTracks_PVTracks_before", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histosQA.add("QA/hist_globalTracks_multT0A_before", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histosQA.add("QA/hist_globalTracks_multV0A_before", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histosQA.add("QA/hist_multV0A_multT0A_before", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histosQA.add("QA/hist_multT0C_centT0C_before", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisT0C}});
      histosQA.add("QA/hist_globalTracks_centT0C_after", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNch}});
      histosQA.add("QA/hist_PVTracks_centT0C_after", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNchPV}});
      histosQA.add("QA/hist_globalTracks_PVTracks_after", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histosQA.add("QA/hist_globalTracks_multT0A_after", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histosQA.add("QA/hist_globalTracks_multV0A_after", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histosQA.add("QA/hist_multV0A_multT0A_after", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histosQA.add("QA/hist_multT0C_centT0C_after", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisT0C}});
    }
    histosQA.add(Form("QA/histVertexZRec"), "", {HistType::kTH1F, {axisvertexz}});
    histosQA.add(Form("QA/histCentrality"), "", {HistType::kTH1F, {axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL0_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL1_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL2_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL3_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL0_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL1_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL2_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL3_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histQvecRes_SigRefAV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvecRes_SigRefBV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvecRes_RefARefBV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});

    histosQA.add(Form("V2/histCosDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/histSinDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});

    histosQA.add(Form("V2/PID/histCosDetV2_Pi"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Ka"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pr"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pi_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Ka_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pr_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});

    if (cfgkOpenCME) {
      if (cfgkOpenPiPi) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPi_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPi_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPi_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPi_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenKaKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPrPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PrPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PrPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PrPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PrPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPiKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPiPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenKaPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenHaHa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_HaHa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_HaHa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_HaHa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_HaHa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
    }
  }

  template <typename CollType>
  bool selEvent(const CollType& collision, const int multTrk, const float centrality)
  {
    histosQA.fill(HIST("QA/histEventCountDetail"), 0.5);
    if (cfgOpenEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (cfgOpenEvSelkIsGoodZvtxFT0vsPV) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 1.5);
    }
    if (cfgOpenEvSelkNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgOpenEvSelkNoSameBunchPileup) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 2.5);
    }
    if (cfgOpenEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (cfgOpenEvSelkNoCollInTimeRangeStandard) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 3.5);
    }
    if (cfgOpenEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return 0;
    }
    if (cfgOpenEvSelkIsGoodITSLayersAll) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 4.5);
    }
    if (cfgOpenEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return 0;
    }
    if (cfgOpenEvSelkNoCollInRofStandard) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 5.5);
    }
    if (cfgOpenEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return 0;
    }
    if (cfgOpenEvSelkNoHighMultCollInPrevRof) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 6.5);
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgOpenEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)) {
      return 0;
    }
    if (cfgOpenEvSelOccupancy) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 7.5);
    }
    if (cfgOpenEvSelMultCorrelation) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
    }
    if (cfgOpenEvSelMultCorrelation) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 8.5);
    }
    if (cfgOpenEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return 0;
    }
    if (cfgOpenEvSelV0AT0ACut) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 9.5);
    }
    return 1;
  }

  template <typename TrackType>
  bool selTrack(const TrackType track)
  {
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;
    return true;
  }

  template <typename CollType>
  void fillHistosQvec(const CollType& collision, int nmode)
  {
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int refAInd = refAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int refBInd = refBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if (nmode == 2) {
      if (collision.qvecAmp()[detId] > 1e-8) {
        histosQA.fill(HIST("QA/histQvec_CorrL0_V2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL1_V2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL2_V2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL3_V2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL0_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL1_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL2_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL3_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.centFT0C());
      }
      if (collision.qvecAmp()[detId] > 1e-8 && collision.qvecAmp()[refAId] > 1e-8 && collision.qvecAmp()[refBId] > 1e-8) {
        histosQA.fill(HIST("QA/histQvecRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histQvecRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histQvecRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode), collision.centFT0C());
      }
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlowGammaDelta(const CollType& collision, const TrackType& track1, const TrackType& track2, const TrackType& track3, int nmode)
  {
    if (collision.qvecAmp()[detId] < 1e-8) {
      return;
    }
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    float psiN = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
    for (const auto& trk : track1) {
      if (!selTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pi"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pi_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        }
      }
    }
    for (const auto& trk : track2) {
      if (!selTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Ka"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Ka_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        }
      }
    }
    for (const auto& trk : track3) {
      if (!selTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pr"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pr_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
        }
      }
    }
    if (cfgkOpenCME) {
      if (cfgkOpenPiPi) {
        for (const auto& trk1 : track1) {
          for (const auto& trk2 : track1) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_PiPi_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPi_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPi_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPi_PP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_PP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPi_NN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_NN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiPi_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPi_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPi_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPi_PN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_PN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPi_NP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_NP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenKaKa) {
        for (const auto& trk1 : track2) {
          for (const auto& trk2 : track2) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_KaKa_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaKa_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_KaKa_PP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_PP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_KaKa_NN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_NN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_KaKa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaKa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_KaKa_PN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_PN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_KaKa_NP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_NP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPrPr) {
        for (const auto& trk1 : track3) {
          for (const auto& trk2 : track3) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_PrPr_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PrPr_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PrPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PrPr_PP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_PP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PrPr_NN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_NN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PrPr_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PrPr_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PrPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PrPr_PN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_PN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PrPr_NP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_NP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPiKa) {
        for (const auto& trk1 : track1) {
          for (const auto& trk2 : track2) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_PiKa_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiKa_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiKa_PP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_PP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiKa_NN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_NN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiKa_PN"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_PN"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiKa_NP"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_NP"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPiPr) {
        for (const auto& trk1 : track1) {
          for (const auto& trk3 : track3) {
            if (trk1.globalIndex() == trk3.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk3.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_PiPr_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPr_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() - trk3.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk3.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPr_PP"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_PP"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPr_NN"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_NN"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() - trk3.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk3.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPr_PN"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_PN"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_PiPr_NP"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_NP"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenKaPr) {
        for (const auto& trk2 : track2) {
          for (const auto& trk3 : track3) {
            if (trk2.globalIndex() == trk3.globalIndex())
              continue;
            if (nmode == 2) {
              if (trk2.sign() == trk3.sign()) {
                histosQA.fill(HIST("PIDCME/histgamama_KaPr_ss"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaPr_ss"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaPr_ss_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_ss_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() - trk3.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk2.sign() > 0 && trk3.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_KaPr_PP"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_PP"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_KaPr_NN"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_NN"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() - trk3.phi())));
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk2.sign() > 0 && trk3.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamama_KaPr_PN"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_PN"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamama_KaPr_NP"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_NP"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::Qvectors>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    const auto cent = collision.centFT0C();
    histosQA.fill(HIST("QA/histEventCount"), 0.5);
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    histosQA.fill(HIST("QA/histEventCount"), 1.5);
    if (cfgOpenFullEventQA) {
      histosQA.fill(HIST("QA/hist_globalTracks_centT0C_before"), cent, tracks.size());
      histosQA.fill(HIST("QA/hist_PVTracks_centT0C_before"), cent, collision.multNTracksPV());
      histosQA.fill(HIST("QA/hist_globalTracks_PVTracks_before"), collision.multNTracksPV(), tracks.size());
      histosQA.fill(HIST("QA/hist_globalTracks_multT0A_before"), collision.multFT0A(), tracks.size());
      histosQA.fill(HIST("QA/hist_globalTracks_multV0A_before"), collision.multFV0A(), tracks.size());
      histosQA.fill(HIST("QA/hist_multV0A_multT0A_before"), collision.multFT0A(), collision.multFV0A());
      histosQA.fill(HIST("QA/hist_multT0C_centT0C_before"), cent, collision.multFT0C());
    }
    if (cfgUseAdditionalEventCut && !selEvent(collision, tracks.size(), cent)) {
      return;
    }
    histosQA.fill(HIST("QA/histEventCount"), 2.5);
    histosQA.fill(HIST("QA/histCentrality"), cent);
    histosQA.fill(HIST("QA/histVertexZRec"), collision.posZ());
    if (cfgOpenFullEventQA) {
      histosQA.fill(HIST("QA/hist_globalTracks_centT0C_after"), cent, tracks.size());
      histosQA.fill(HIST("QA/hist_PVTracks_centT0C_after"), cent, collision.multNTracksPV());
      histosQA.fill(HIST("QA/hist_globalTracks_PVTracks_after"), collision.multNTracksPV(), tracks.size());
      histosQA.fill(HIST("QA/hist_globalTracks_multT0A_after"), collision.multFT0A(), tracks.size());
      histosQA.fill(HIST("QA/hist_globalTracks_multV0A_after"), collision.multFV0A(), tracks.size());
      histosQA.fill(HIST("QA/hist_multV0A_multT0A_after"), collision.multFT0A(), collision.multFV0A());
      histosQA.fill(HIST("QA/hist_multT0C_centT0C_after"), cent, collision.multFT0C());
    }
    auto tracks1 = tracksSet1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracks2 = tracksSet2->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracks3 = tracksSet3->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult1 = tracks1.size();
    mult2 = tracks2.size();
    mult3 = tracks3.size();
    if (mult1 < 1 || mult2 < 1 || mult3 < 1) // Reject Collisions without sufficient particles
      return;
    for (auto i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
      int detIndGlobal = detId * 4 + cfgnTotalSystem * 4 * (cfgnMods->at(i) - 2);
      float psiNGlobal = helperEP.GetEventPlane(collision.qvecRe()[detIndGlobal + 3], collision.qvecIm()[detIndGlobal + 3], cfgnMods->at(i));
      for (const auto& trk : tracks) {
        if (!selTrack(trk))
          continue;
        histosQA.fill(HIST("V2/histSinDetV2"), cent, trk.pt(),
                      std::sin(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
        histosQA.fill(HIST("V2/histCosDetV2"), cent, trk.pt(),
                      std::cos(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
      }
      if (cfgkOpenCME && cfgkOpenHaHa && cfgnMods->at(i) == 2) {
        for (const auto& trk1 : tracks) {
          for (const auto& trk2 : tracks) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (trk1.sign() == trk2.sign()) {
              histosQA.fill(HIST("PIDCME/histgamama_HaHa_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_HaHa_ss_Dif"), cent, trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_ss_Dif"), cent, trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
              if (cfgkOpenSsOsCrossCheck) {
                if (trk1.sign() > 0 && trk2.sign() > 0) {
                  histosQA.fill(HIST("PIDCME/histgamama_HaHa_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                } else {
                  histosQA.fill(HIST("PIDCME/histgamama_HaHa_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                }
              }
            } else {
              histosQA.fill(HIST("PIDCME/histgamama_HaHa_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_os"), cent, std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_HaHa_os_Dif"), cent, trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_os_Dif"), cent, trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
              if (cfgkOpenSsOsCrossCheck) {
                if (trk1.sign() > 0 && trk2.sign() < 0) {
                  histosQA.fill(HIST("PIDCME/histgamama_HaHa_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                } else {
                  histosQA.fill(HIST("PIDCME/histgamama_HaHa_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_NP"), cent, std::cos((trk1.phi() - trk2.phi())));
                }
              }
            }
          }
        }
      }
      fillHistosQvec(collision, cfgnMods->at(i));
      fillHistosFlowGammaDelta(collision, tracks1, tracks2, tracks3, cfgnMods->at(i));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FillPIDcolums>(cfgc),
    adaptAnalysisTask<QAProcessCent>(cfgc),
    adaptAnalysisTask<pidcme>(cfgc),
  };
}
