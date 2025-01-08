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
} // namespace cme_track_pid_columns
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", cme_track_pid_columns::NPidFlag);
} // namespace o2::aod

using TracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
struct FillPIDcolums {
  Configurable<float> cfgnSigmaCutTPCPi{"cfgnSigmaCutTPCPi", 3.0, "Value of the TPC Nsigma cut for pions"};
  Configurable<float> cfgnSigmaCutTPCKa{"cfgnSigmaCutTPCKa", 3.0, "Value of the TPC Nsigma cut for kaons"};
  Configurable<float> cfgnSigmaCutTPCPr{"cfgnSigmaCutTPCPr", 3.0, "Value of the TPC Nsigma cut for protons"};
  Configurable<float> cfgnSigmaCutTOFPi{"cfgnSigmaCutTOFPi", 3.0, "Value of the TOF Nsigma cut for pions"};
  Configurable<float> cfgnSigmaCutTOFKa{"cfgnSigmaCutTOFKa", 3.0, "Value of the TOF Nsigma cut for kaons"};
  Configurable<float> cfgnSigmaCutTOFPr{"cfgnSigmaCutTOFPr", 3.0, "Value of the TOF Nsigma cut for protons"};
  Configurable<float> cfgnSigmaCutCombine{"cfgnSigmaCutCombine", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<float> cfgPtMaxforTPCOnlyPID{"cfgPtMaxforTPCOnlyPID", 0.4, "Maxmium track pt for TPC only PID,only when onlyTOF and onlyTOFHIT closed"};
  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.8, "Maximum track #eta for PID"};
  Configurable<float> cfgnSigmaCutITSPi{"cfgnSigmaCutITSPi", 3.0, "Value of the ITS Nsigma cut for Pions"};
  Configurable<float> cfgnSigmaCutITSKa{"cfgnSigmaCutITSKa", 2.5, "Value of the ITS Nsigma cut for Kaons"};
  Configurable<float> cfgnSigmaCutITSPr{"cfgnSigmaCutITSPr", 2.0, "Value of the ITS Nsigma cut for Protons"};
  Configurable<float> cfgAveClusSizeCoslMinPi{"cfgAveClusSizeCoslMinPi", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Pions"};
  Configurable<float> cfgAveClusSizeCoslMaxPi{"cfgAveClusSizeCoslMaxPi", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Pions"};
  Configurable<float> cfgAveClusSizeCoslMinKa{"cfgAveClusSizeCoslMinKa", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Kaons"};
  Configurable<float> cfgAveClusSizeCoslMaxKa{"cfgAveClusSizeCoslMaxKa", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Kaons"};
  Configurable<float> cfgAveClusSizeCoslMinPr{"cfgAveClusSizeCoslMinPr", 0, "Base line for minmum ITS cluster size x cos(#lambda) for Protons"};
  Configurable<float> cfgAveClusSizeCoslMaxPr{"cfgAveClusSizeCoslMaxPr", 1e9, "Base line for maxmum ITS cluster size x cos(#lambda) for Protons"};

  ConfigurableAxis cfgrigidityBins{"cfgrigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis cfgdedxBins{"cfgdedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis cfgnSigmaBins{"cfgnSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis cfgnSigmaBinsCom{"cfgnSigmaBinsCom", {100, 0.f, 10.f}, "Combination Binning for TPC&TOF nsigma"};
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {24, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgaxispPID{"cfgaxispPID", {50, 0, 5}, "Binning for P PID"};
  ConfigurableAxis cfgaxisAverClusterCosl{"cfgaxisAverClusterCosl", {50, 0, 10}, "Binning for average cluster size x cos(#lambda)"};

  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  Configurable<bool> openITSCut{"openITSCut", true, "open ITSnsigma cut"};
  bool onlyTPC = true;

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
  bool selectionPid(const T& candidate, int8_t PID, float clustersize)
  {
    if (candidate.pt() > cfgPtMaxforTPCOnlyPID) {
      onlyTPC = false;
    }

    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOFPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOFPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPi()) < cfgnSigmaCutITSPi && clustersize > cfgAveClusSizeCoslMinPi && clustersize < cfgAveClusSizeCoslMaxPi)
              return true;
          } else {
            return true;
          }
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOFKa) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOFKa) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaKa()) < cfgnSigmaCutITSKa && clustersize > cfgAveClusSizeCoslMinKa && clustersize < cfgAveClusSizeCoslMaxKa)
              return true;
          } else {
            return true;
          }
        }
      }
    } else if (PID == 2) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOFPr) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOFPr) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          if (openITSCut) {
            if (std::abs(candidate.itsNSigmaPr()) < cfgnSigmaCutITSPr && clustersize > cfgAveClusSizeCoslMinPr && clustersize < cfgAveClusSizeCoslMaxPr)
              return true;
          } else {
            return true;
          }
        }
      }
    }
    return false;
  }

  HistogramRegistry histosQA{"histosQAPID", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    AxisSpec axisRigidity{cfgrigidityBins, "#it{p}^{TPC}/#it{z}"};
    AxisSpec axisdEdx{cfgdedxBins, "d#it{E}/d#it{x}"};
    AxisSpec axisnSigmaTPC{cfgnSigmaBins, "n_{#sigma}TPC"};
    AxisSpec axisnSigmaTOF{cfgnSigmaBins, "n_{#sigma}TOF"};
    AxisSpec axisnSigmaITS{cfgnSigmaBins, "n_{#sigma}TOF"};
    AxisSpec axisnSigmaCom{cfgnSigmaBinsCom, "hypot(n_{#sigma}TPC,TOF)"};
    AxisSpec axisPtPID{cfgaxisptPID, "#it{p}_{T}"};
    AxisSpec axisPPID{cfgaxispPID, "#it{p}"};
    AxisSpec axisClusterSize{cfgaxisAverClusterCosl, "<ITS Cluster Size> x <cos(#lambda)>"};
    AxisSpec axisPhi = {100, 0, 2.1 * constants::math::PI, "#phi"};
    // TH3D NSigmaTPC,NSigmaTOF,pt
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Pi"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Ka"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_Pr"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Pi"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Ka"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    histosQA.add(Form("QA/PID/histnSigma_TPC_TOF_cross_Pr"), "", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaTOF, axisPtPID}});
    // Hist for PID Averge Cluster Size
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
    histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSize}});
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
    // Hist for checking the PID phi distribution
    histosQA.add(Form("QA/PID/histPhi_Dis_Pi"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_Ka"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_Pr"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pi"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Ka"), "", {HistType::kTH1F, {axisPhi}});
    histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pr"), "", {HistType::kTH1F, {axisPhi}});
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
        if (selectionPid(track, 0, averClusSizeCosl)) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Pi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pi"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Pi"), track.tpcNSigmaPi());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pi"), track.pt(), track.tpcNSigmaPi());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pi"), track.itsNSigmaPi());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pi"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), track.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pi"), track.phi());
        }
        if (selectionPid(track, 1, averClusSizeCosl)) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Ka"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Ka"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Ka"), track.tpcNSigmaKa());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Ka"), track.pt(), track.tpcNSigmaKa());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Ka"), track.itsNSigmaKa());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Ka"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), track.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Ka"), track.phi());
        }
        if (selectionPid(track, 2, averClusSizeCosl)) {
          histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_cross_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
          histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
          histosQA.fill(HIST("QA/PID/histnSigma_cross_Pr"), track.tpcNSigmaPr());
          histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pr"), track.pt(), track.tpcNSigmaPr());
          histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pr"), track.itsNSigmaPr());
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pr"), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), track.p(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), track.tpcNSigmaKa(), averClusSizeCosl);
          histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pr"), track.phi());
        }
        histosQA.fill(HIST("QA/PID/histdEdxTPC_All"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
        float nSigmaArray[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
        pidFlag = 0;
        for (int8_t i = 0; i < 3; i++) {
          if (selectionPid(track, i, averClusSizeCosl))
            pidFlag = pidFlag * 10 + i + 1;
          if (pidFlag > 10) {                                                                          // If a track is identified as two different tracks.
            if (std::abs(nSigmaArray[(pidFlag / 10) - 1]) < std::abs(nSigmaArray[(pidFlag % 10) - 1])) // The track is identified as the particle whose |nsigma| is the least.
              pidFlag /= 10;
            else
              pidFlag %= 10;
          }
        }

        switch (pidFlag) {
          case 1:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Pi"), track.tpcNSigmaPi(), track.tofNSigmaPi(), track.pt());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pi"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pi"), track.tpcNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pi"), track.pt(), track.tpcNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Pi"), (track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pi"), track.tofNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pi"), track.itsNSigmaPi());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pi"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pi"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Pi"), track.phi());
            break;
          case 2:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Ka"), track.tpcNSigmaKa(), track.tofNSigmaKa(), track.pt());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Ka"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Ka"), track.tpcNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Ka"), track.pt(), track.tpcNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Ka"), (track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Ka"), track.tofNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Ka"), track.itsNSigmaKa());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Ka"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Ka"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Ka"), track.phi());
            break;
          case 3:
            histosQA.fill(HIST("QA/PID/histnSigma_TPC_TOF_Pr"), track.tpcNSigmaPr(), track.tofNSigmaPr(), track.pt());
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pr"), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pr"), track.pt(), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_com_Pr"), (track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr()));
            histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pr"), track.tofNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pr"), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pr"), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pr"), track.p(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), track.tpcNSigmaPr(), averClusSizeCosl);
            histosQA.fill(HIST("QA/PID/histPhi_Dis_Pr"), track.phi());
            break;
        }
      }
      pidCmeTable(pidFlag);
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

  ConfigurableAxis cfgaxissumpt{"cfgaxissumpt", {7, 1, 8}, "Binning for #gamma and #delta pt(particle1 + particle2)"};
  ConfigurableAxis cfgaxisdeltaeta{"cfgaxisdeltaeta", {5, 0, 1}, "Binning for #gamma and #delta |#eta(particle1 - particle2)|"};

  Configurable<bool> cfgkOpeanCME{"cfgkOpeanCME", true, "open PID CME"};
  Configurable<bool> cfgkOpeanPiPi{"cfgkOpeanPiPi", true, "open Pi-Pi"};
  Configurable<bool> cfgkOpeanKaKa{"cfgkOpeanKaKa", false, "open Ka-Ka"};
  Configurable<bool> cfgkOpeanPrPr{"cfgkOpeanPrPr", false, "open Pr-Pr"};
  Configurable<bool> cfgkOpeanPiKa{"cfgkOpeanPiKa", true, "open Pi-Ka"};
  Configurable<bool> cfgkOpeanPiPr{"cfgkOpeanPiPr", true, "open Pi-Pr"};
  Configurable<bool> cfgkOpeanKaPr{"cfgkOpeanKaPr", true, "open Ka-Pr"};
  Configurable<bool> cfgkOpeanHaHa{"cfgkOpeanHaHa", false, "open Ha-Ha"};

  EventPlaneHelper helperEP;
  SliceCache cache;

  unsigned int mult1, mult2, mult3;
  int detId;
  int refAId;
  int refBId;

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

  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > cfgMinPt;
  Filter etafilter = aod::track::eta < cfgMaxEta;
  Filter properPIDfilter = aod::cme_track_pid_columns::nPidFlag != -1;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet1 = aod::cme_track_pid_columns::nPidFlag == 1;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet2 = aod::cme_track_pid_columns::nPidFlag == 2;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> tracksSet3 = aod::cme_track_pid_columns::nPidFlag == 3;
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

    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgaxiscos, "angle function"};
    AxisSpec axisPt{cfgaxispt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgaxisCentMerged, "merged centrality"};

    AxisSpec axissumpt{cfgaxissumpt, "#it{p}_{T}^{sum}}"};
    AxisSpec axisdeltaeta{cfgaxisdeltaeta, "#Delta#eta"};
    AxisSpec axisvertexz = {100, -15., 15., "vrtx_{Z} [cm]"};

    histosQA.add(Form("QA/histEventCount"), "", {HistType::kTH1F, {{2, 0.0, 2.0}}});
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(1, "Not selected events");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(2, "Selected events");
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

    if (cfgkOpeanCME) {
      if (cfgkOpeanPiPi) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPi_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPi_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPi_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPi_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanKaKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanPrPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PrPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PrPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PrPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PrPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanPiKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanPiPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanKaPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
      if (cfgkOpeanHaHa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamama_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_HaHa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_HaHa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_HaHa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
        histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_HaHa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      }
    }
  }

  template <typename CollType>
  bool selEvent(const CollType& collision)
  {
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
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
    if (cfgkOpeanCME) {
      if (cfgkOpeanPiPi) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiPi_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPi_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPi_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
              }
            }
          }
        }
      }
      if (cfgkOpeanKaKa) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_KaKa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaKa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
              }
            }
          }
        }
      }
      if (cfgkOpeanPrPr) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PrPr_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PrPr_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PrPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
              }
            }
          }
        }
      }
      if (cfgkOpeanPiKa) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                              std::cos((trk1.phi() - trk2.phi())));
              }
            }
          }
        }
      }
      if (cfgkOpeanPiPr) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                              std::cos((trk1.phi() - trk3.phi())));
              }
            }
          }
        }
      }
      if (cfgkOpeanKaPr) {
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
              } else {
                histosQA.fill(HIST("PIDCME/histgamama_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
                histosQA.fill(HIST("PIDCME/Differential/histgamama_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                              std::cos((trk2.phi() - trk3.phi())));
              }
            }
          }
        }
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Qvectors>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    histosQA.fill(HIST("QA/histEventCount"), 0.5);
    if (!selEvent(collision)) {
      return;
    }
    histosQA.fill(HIST("QA/histEventCount"), 1.5);
    histosQA.fill(HIST("QA/histCentrality"), collision.centFT0C());
    histosQA.fill(HIST("QA/histVertexZRec"), collision.posZ());
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
        histosQA.fill(HIST("V2/histSinDetV2"), collision.centFT0C(), trk.pt(),
                      std::sin(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
        histosQA.fill(HIST("V2/histCosDetV2"), collision.centFT0C(), trk.pt(),
                      std::cos(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
      }
      if (cfgkOpeanCME && cfgkOpeanHaHa && cfgnMods->at(i) == 2) {
        for (const auto& trk1 : tracks) {
          for (const auto& trk2 : tracks) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (trk1.sign() == trk2.sign()) {
              histosQA.fill(HIST("PIDCME/histgamama_HaHa_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_HaHa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
            } else {
              histosQA.fill(HIST("PIDCME/histgamama_HaHa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_HaHa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
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
    adaptAnalysisTask<pidcme>(cfgc),
  };
}
