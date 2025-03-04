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
#include <CCDB/BasicCCDBManager.h>
#include <chrono>
#include <string>
#include <vector>
#include <utility>
#include <memory>
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
DECLARE_SOA_COLUMN(NPidFlag, nPidFlag, int8_t); // Flag tracks without proper binning as -1, and indicate type of particle [0->(un_Id) 1->(pi_only), 2->(ka_only), 3->(pr_only), 4->(pi_ITSleft), 5->(ka_ITSleft), 6->(pr_ITSleft), 7->(pi_ka), 8->(pi_Pr), 9->(ka_pr), 10->(pi_ka_pr), 11->(pi_ka_ITSleft), 12->(pi_pr_ITSleft), 13->(ka_pr_ITSleft), 14->(pi_ka_pr_ITSleft)]
DECLARE_SOA_COLUMN(AverClusterSizeCosl, averClusterSizeCosl, float);
DECLARE_SOA_COLUMN(NSigmaPiITS, nSigmaPiITS, float);
DECLARE_SOA_COLUMN(NSigmaKaITS, nSigmaKaITS, float);
DECLARE_SOA_COLUMN(NSigmaPrITS, nSigmaPrITS, float);
DECLARE_SOA_COLUMN(NSigmaPiTPC, nSigmaPiTPC, float);
DECLARE_SOA_COLUMN(NSigmaKaTPC, nSigmaKaTPC, float);
DECLARE_SOA_COLUMN(NSigmaPrTPC, nSigmaPrTPC, float);
DECLARE_SOA_COLUMN(NSigmaPiTOF, nSigmaPiTOF, float);
DECLARE_SOA_COLUMN(NSigmaKaTOF, nSigmaKaTOF, float);
DECLARE_SOA_COLUMN(NSigmaPrTOF, nSigmaPrTOF, float);
} // namespace cme_track_pid_columns
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", cme_track_pid_columns::NPidFlag);
DECLARE_SOA_TABLE(PidInfo, "AOD", "PidInfo", cme_track_pid_columns::AverClusterSizeCosl, cme_track_pid_columns::NSigmaPiITS, cme_track_pid_columns::NSigmaKaITS, cme_track_pid_columns::NSigmaPrITS, cme_track_pid_columns::NSigmaPiTPC, cme_track_pid_columns::NSigmaKaTPC, cme_track_pid_columns::NSigmaPrTPC, cme_track_pid_columns::NSigmaPiTOF, cme_track_pid_columns::NSigmaKaTOF, cme_track_pid_columns::NSigmaPrTOF);
} // namespace o2::aod

using TracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
using CollisionPID = soa::Join<aod::Collisions, aod::CentFT0Cs>;
struct FillPIDcolums {
  Configurable<float> cfgPtMaxforTPCOnlyPID{"cfgPtMaxforTPCOnlyPID", 0.4, "Maxmium track pt for TPC only PID,only when onlyTOF and onlyTOFHIT closed"};
  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxPtPID{"cfgMaxPtPID", 99.9, "Maximum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.8, "Maximum track #eta for PID"};
  Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 2.5, "Maximum chi2 per cluster TPC for PID if not use costom track cuts"};
  Configurable<float> cfgMaxChi2NClITS{"cfgMaxChi2NClITS", 2.5, "Maximum chi2 per cluster ITS for PID if not use costom track cuts"};
  Configurable<float> cfgMinTPCCls{"cfgMinTPCCls", 70, "Minimum TPC clusters for PID if not use costom track cuts"};
  Configurable<float> cfgMinITSCls{"cfgMinITSCls", 5, "Minimum ITS clusters for PID if not use costom track cuts"};
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

  Configurable<bool> cfgQuietMode{"cfgQuietMode", false, "open quiet mode for saving cpu cost and only do some basic QA plots"};
  Configurable<bool> cfgOpenPlotnSigmaTOFITSPt{"cfgOpenPlotnSigmaTOFITSPt", true, "plot nSigmaTOF vs nSigmaITS vs Pt"};
  Configurable<bool> cfgOpenPlotnSigmaITSTPCPt{"cfgOpenPlotnSigmaITSTPCPt", true, "plot nSigmaITS vs nSigmaTOF vs Pt"};
  Configurable<bool> cfgOpenPlotnSigmaTOFTPCPt{"cfgOpenPlotnSigmaTOFTPCPt", true, "plot nSigmaTOF vs nSigmaTPC vs Pt"};
  Configurable<bool> cfgOpenPlotAverClus{"cfgOpenPlotAverClus", true, "plot average cluster size x cos(#lambda)"};
  Configurable<bool> cfgOpenPlotAverClusP{"cfgOpenPlotAverClusP", true, "plot average cluster size x cos(#lambda) vs p"};
  Configurable<bool> cfgOpenPlotAverClusnSigmaTPC{"cfgOpenPlotAverClusnSigmaTPC", true, "plot average cluster size x cos(#lambda) vs nSigmaTPC"};
  Configurable<bool> cfgOpenPlotPhiDis{"cfgOpenPlotPhiDis", true, "plot phi distribution QA"};
  Configurable<bool> cfgOpenPlotPhiDisPtEta{"cfgOpenPlotPhiDisPtEta", true, "plot phi pt eta distribution QA"};
  Configurable<bool> cfgOpenITSCut{"cfgOpenITSCut", true, "open ITSnsigma cut"};
  Configurable<bool> cfgOpenDetailPlotsTPCITSContaimination{"cfgOpenDetailPlotsTPCITSContaimination", false, "open detail TH3D plots for nSigmaTPC-ITS Pt-eta-Phi nSigmaITS-clustersize"};
  Configurable<bool> cfgOpenAllowCrossTrack{"cfgOpenAllowCrossTrack", false, "Allow one track to be identified as different kind of PID particles"};
  Configurable<bool> cfgOpenCrossTrackQAPlots{"cfgOpenCrossTrackQAPlots", true, "open cross pid track QA plots"};
  Configurable<bool> cfgOpenTOFOnlyPID{"cfgOpenTOFOnlyPID", true, "only accept tracks who has TOF infomation and use TOFnsigma for PID(priority greater than TPConly and combined)"};
  Configurable<bool> cfgOpenTPCOnlyPID{"cfgOpenTPCOnlyPID", false, "only use TPCnsigma for PID(priority grater than combined less than TOFOnly)"};
  Configurable<bool> cfgUseCostomTrackCuts{"cfgUseCostomTrackCuts", true, "use track cuts from default track selection table producer"};
  Configurable<bool> cfgOpenPtRangedTOFnSigmacutPi{"cfgOpenPtRangedTOFnSigmacutPi", false, "use nSigma TOF cut for different pt Pion"};
  Configurable<bool> cfgOpenPtRangedTPCnSigmacutPi{"cfgOpenPtRangedTPCnSigmacutPi", false, "use nSigma TPC cut for different pt Pion"};
  Configurable<bool> cfgOpenPtRangedITSnSigmacutPi{"cfgOpenPtRangedITSnSigmacutPi", false, "use nSigma ITS cut for different pt Pion"};
  Configurable<bool> cfgOpenPtRangedTOFnSigmacutKa{"cfgOpenPtRangedTOFnSigmacutKa", false, "use nSigma TOF cut for different pt Kaon"};
  Configurable<bool> cfgOpenPtRangedTPCnSigmacutKa{"cfgOpenPtRangedTPCnSigmacutKa", false, "use nSigma TPC cut for different pt Kaon"};
  Configurable<bool> cfgOpenPtRangedITSnSigmacutKa{"cfgOpenPtRangedITSnSigmacutKa", false, "use nSigma ITS cut for different pt Kaon"};
  Configurable<bool> cfgOpenPtRangedTOFnSigmacutPr{"cfgOpenPtRangedTOFnSigmacutPr", false, "use nSigma TOF cut for different pt Proton"};
  Configurable<bool> cfgOpenPtRangedTPCnSigmacutPr{"cfgOpenPtRangedTPCnSigmacutPr", false, "use nSigma TPC cut for different pt Proton"};
  Configurable<bool> cfgOpenPtRangedITSnSigmacutPr{"cfgOpenPtRangedITSnSigmacutPr", false, "use nSigma ITS cut for different pt Proton"};

  Configurable<std::vector<float>> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", {3, 3, 3}, "TPC nsigma cut for pi k p respectively at low pt and for the TPCOnly case"};
  Configurable<std::vector<float>> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", {1.5, 1.5, 1.5}, "TOF nsigma cut for pi k p respectively for the TOFonly case"};
  Configurable<std::vector<float>> cfgnSigmaCutRMS{"cfgnSigmaCutRMS", {3, 3, 3}, "TPC_TOF combined cut for pi k p respectively at high pt"};
  Configurable<std::vector<float>> cfgnSigmaCutITS{"cfgnSigmaCutITS", {3, 2.5, 2}, "ITS nSigma cut for pi k p"};
  Configurable<std::vector<float>> cfgPtBinPionPID{"cfgPtBinPionPID", {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0}, "pt bin for pion PIDnsigma"};
  Configurable<std::vector<float>> cfgPtBinKaonPID{"cfgPtBinKaonPID", {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0}, "pt bin for pion PIDnsigma"};
  Configurable<std::vector<float>> cfgPtBinProtonPID{"cfgPtBinProtonPID", {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0}, "pt bin for pion PIDnsigma"};
  Configurable<std::vector<float>> cfgnSigmaTPCPionPt{"cfgnSigmaTPCPionPt", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaTPC cut anchored to pion pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFPionPt{"cfgnSigmaTOFPionPt", {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}, "nSigmaTOF cut anchored to pion pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSPionPt{"cfgnSigmaITSPionPt", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaITS cut anchored to pion pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTPCKaonPt{"cfgnSigmaTPCKaonPt", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaTPC cut anchored to kaon pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFKaonPt{"cfgnSigmaTOFKaonPt", {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}, "nSigmaTOF cut anchored to kaon pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSKaonPt{"cfgnSigmaITSKaonPt", {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5}, "nSigmaITS cut anchored to kaon pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTPCProtonPt{"cfgnSigmaTPCProtonPt", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaTPC cut anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFProtonPt{"cfgnSigmaTOFProtonPt", {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}, "nSigmaTOF cut anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSProtonPt{"cfgnSigmaITSProtonPt", {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, "nSigmaITS cut anchored to proton pt bins"};

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
    if ((track.pt() < cfgMinPtPID) || (track.pt() > cfgMaxPtPID))
      return false;
    if (std::abs(track.eta()) > cfgMaxEtaPID)
      return false;
    if (cfgUseCostomTrackCuts) {
      if (!track.passedITSNCls())
        return false;
      if (!track.passedITSChi2NDF())
        return false;
      if (!track.passedITSHits())
        return false;
      if (!track.passedTPCChi2NDF())
        return false;
    } else {
      if (track.tpcChi2NCl() > cfgMaxTPCChi2NCl)
        return false;
      if (track.tpcNClsFound() < cfgMinTPCCls)
        return false;
      if (track.itsChi2NCl() < cfgMaxChi2NClITS)
        return false;
      if (track.itsNCls() < cfgMinITSCls)
        return false;
    }
    if (!track.passedTPCCrossedRowsOverNCls())
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
    // initialization for basic parameter
    float averClusSizeCosl = averageClusterSizeCosl(candidate.itsClusterSizes(), candidate.eta());
    std::array<float, 3> nSigmaTPC = {candidate.tpcNSigmaPi(), candidate.tpcNSigmaKa(), candidate.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {candidate.tofNSigmaPi(), candidate.tofNSigmaKa(), candidate.tofNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(candidate.tpcNSigmaPi(), candidate.tofNSigmaPi()), std::hypot(candidate.tpcNSigmaKa(), candidate.tofNSigmaKa()), std::hypot(candidate.tpcNSigmaPr(), candidate.tofNSigmaPr())};
    std::array<float, 3> nSigmaToUse;
    std::vector<float> pidVector;
    std::vector<float> pidVectorTOFPt;
    std::vector<float> pidVectorTPCPt;
    int pid = -1;
    bool kIsPi = false, kIsKa = false, kIsPr = false;
    int currentPtBinPi = -1, currentPtBinKa = -1, currentPtBinPr = -1;
    if (cfgOpenPtRangedTOFnSigmacutPi || cfgOpenPtRangedTPCnSigmacutPi) {
      for (int i = 0; i < static_cast<int>(cfgPtBinPionPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinPionPID.value[i] && candidate.pt() < cfgPtBinPionPID.value[i + 1]) {
          currentPtBinPi = i;
          break;
        }
      }
    }
    if (cfgOpenPtRangedTOFnSigmacutKa || cfgOpenPtRangedTPCnSigmacutKa) {
      for (int i = 0; i < static_cast<int>(cfgPtBinKaonPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinKaonPID.value[i] && candidate.pt() < cfgPtBinKaonPID.value[i + 1]) {
          currentPtBinKa = i;
          break;
        }
      }
    }
    if (cfgOpenPtRangedTOFnSigmacutPr || cfgOpenPtRangedTPCnSigmacutPr) {
      for (int i = 0; i < static_cast<int>(cfgPtBinProtonPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinProtonPID.value[i] && candidate.pt() < cfgPtBinProtonPID.value[i + 1]) {
          currentPtBinPr = i;
          break;
        }
      }
    }
    float nSigmaTOFCutPiPt = (currentPtBinPi == -1) ? cfgnSigmaCutTOF.value[0] : cfgnSigmaTOFPionPt.value[currentPtBinPi];
    float nSigmaTOFCutKaPt = (currentPtBinKa == -1) ? cfgnSigmaCutTOF.value[1] : cfgnSigmaTOFKaonPt.value[currentPtBinKa];
    float nSigmaTOFCutPrPt = (currentPtBinPr == -1) ? cfgnSigmaCutTOF.value[2] : cfgnSigmaTOFProtonPt.value[currentPtBinPr];
    float nSigmaTPCCutPiPt = (currentPtBinPi == -1) ? cfgnSigmaCutTPC.value[0] : cfgnSigmaTPCPionPt.value[currentPtBinPi];
    float nSigmaTPCCutKaPt = (currentPtBinKa == -1) ? cfgnSigmaCutTPC.value[1] : cfgnSigmaTPCKaonPt.value[currentPtBinKa];
    float nSigmaTPCCutPrPt = (currentPtBinPr == -1) ? cfgnSigmaCutTPC.value[2] : cfgnSigmaTPCProtonPt.value[currentPtBinPr];
    pidVectorTOFPt.push_back(nSigmaTOFCutPiPt);
    pidVectorTOFPt.push_back(nSigmaTOFCutKaPt);
    pidVectorTOFPt.push_back(nSigmaTOFCutPrPt);
    pidVectorTPCPt.push_back(nSigmaTPCCutPiPt);
    pidVectorTPCPt.push_back(nSigmaTPCCutKaPt);
    pidVectorTPCPt.push_back(nSigmaTPCCutPrPt);
    // Choose which nSigma array and PIDcut array to use
    if (cfgOpenTOFOnlyPID) {
      if (!candidate.hasTOF())
        return -1;
      nSigmaToUse = nSigmaTOF;
      pidVector = pidVectorTOFPt;
    } else if (cfgOpenTPCOnlyPID) {
      nSigmaToUse = nSigmaTPC;
      pidVector = pidVectorTPCPt;
    } else {
      nSigmaToUse = (candidate.pt() > cfgPtMaxforTPCOnlyPID && candidate.hasTOF()) ? nSigmaCombined : nSigmaTPC;
      pidVector = (candidate.pt() > cfgPtMaxforTPCOnlyPID && candidate.hasTOF()) ? cfgnSigmaCutRMS.value : cfgnSigmaCutTPC.value;
    }
    float nsigma = pidVector[0];
    // Fill cross pid QA
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < pidVector[i]) {
        if (i == 0) {
          kIsPi = true;
          if (!cfgQuietMode) {
            if (cfgOpenCrossTrackQAPlots) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pi"), candidate.itsNSigmaPi(), candidate.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_cross_Pi"), candidate.tofNSigmaPi(), candidate.itsNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Pi"), candidate.tofNSigmaPi(), candidate.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pi"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_cross_Pi"), candidate.tofNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_cross_Pi"), candidate.pt(), candidate.tofNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_cross_Pi"), candidate.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pi"), candidate.pt(), candidate.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pi"), candidate.itsNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_cross_Pi"), candidate.pt(), candidate.itsNSigmaPi());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_cross_Pi"), candidate.pt(), candidate.tofNSigmaPi(), candidate.itsNSigmaPi());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_cross_Pi"), candidate.pt(), candidate.tofNSigmaPi(), candidate.tpcNSigmaPi());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_cross_Pi"), candidate.pt(), candidate.itsNSigmaPi(), candidate.tpcNSigmaPi());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pi"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), candidate.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), candidate.tpcNSigmaPi(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pi"), candidate.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_cross_Pi"), candidate.phi(), candidate.pt(), candidate.eta());
              }
            }
          }
        }
        if (i == 1) {
          kIsKa = true;
          if (!cfgQuietMode) {
            if (cfgOpenCrossTrackQAPlots) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Ka"), candidate.itsNSigmaKa(), candidate.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_cross_Ka"), candidate.tofNSigmaKa(), candidate.itsNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Ka"), candidate.tofNSigmaKa(), candidate.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Ka"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_cross_Ka"), candidate.tofNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_cross_Ka"), candidate.pt(), candidate.tofNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_cross_Ka"), candidate.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Ka"), candidate.pt(), candidate.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Ka"), candidate.itsNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_cross_Ka"), candidate.pt(), candidate.itsNSigmaKa());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_cross_Ka"), candidate.pt(), candidate.tofNSigmaKa(), candidate.itsNSigmaKa());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_cross_Ka"), candidate.pt(), candidate.tofNSigmaKa(), candidate.tpcNSigmaKa());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_cross_Ka"), candidate.pt(), candidate.itsNSigmaKa(), candidate.tpcNSigmaKa());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Ka"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), candidate.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), candidate.tpcNSigmaKa(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Ka"), candidate.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_cross_Ka"), candidate.phi(), candidate.pt(), candidate.eta());
              }
            }
          }
        }
        if (i == 2) {
          kIsPr = true;
          if (!cfgQuietMode) {
            if (cfgOpenCrossTrackQAPlots) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pr"), candidate.itsNSigmaPr(), candidate.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_cross_Pr"), candidate.tofNSigmaPr(), candidate.itsNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Pr"), candidate.tofNSigmaPr(), candidate.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_cross_Pr"), candidate.sign() * candidate.tpcInnerParam(), candidate.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_cross_Pr"), candidate.tofNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_cross_Pr"), candidate.pt(), candidate.tofNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_cross_Pr"), candidate.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_Pt_cross_Pr"), candidate.pt(), candidate.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_cross_Pr"), candidate.itsNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_cross_Pr"), candidate.pt(), candidate.itsNSigmaPr());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_cross_Pr"), candidate.pt(), candidate.tofNSigmaPr(), candidate.itsNSigmaPr());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_cross_Pr"), candidate.pt(), candidate.tofNSigmaPr(), candidate.tpcNSigmaPr());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_cross_Pr"), candidate.pt(), candidate.itsNSigmaPr(), candidate.tpcNSigmaPr());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_cross_Pr"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), candidate.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), candidate.tpcNSigmaPr(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_cross_Pr"), candidate.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_cross_Pr"), candidate.phi(), candidate.pt(), candidate.eta());
              }
            }
          }
        }
      }
    }
    if (cfgOpenAllowCrossTrack) {
      // one track can be recognized as different PID particles
      int index = (kIsPr << 2) | (kIsKa << 1) | kIsPi;
      const int map[] = {0, 1, 2, 7, 3, 8, 9, 10};
      return map[index];
    } else {
      // Select particle with the lowest nsigma (If not allow cross track)
      for (int i = 0; i < 3; ++i) {
        if (std::abs(nSigmaToUse[i]) < nsigma && std::abs(nSigmaToUse[i]) < pidVector[i]) {
          pid = i;
          nsigma = std::abs(nSigmaToUse[i]);
        }
      }
      return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
    }
    // Clear the vectors
    std::vector<float>().swap(pidVector);
    std::vector<float>().swap(pidVectorTOFPt);
    std::vector<float>().swap(pidVectorTPCPt);
  }

  template <typename T>
  bool selectionITS(const T& candidate, int mode, float avgclssize)
  {
    std::array<float, 3> nSigmaITSToUse;
    int currentPtBinPi = -1, currentPtBinKa = -1, currentPtBinPr = -1;
    if (cfgOpenPtRangedITSnSigmacutPi) {
      for (int i = 0; i < static_cast<int>(cfgPtBinPionPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinPionPID.value[i] && candidate.pt() < cfgPtBinPionPID.value[i + 1]) {
          currentPtBinPi = i;
          break;
        }
      }
    }
    if (cfgOpenPtRangedITSnSigmacutKa) {
      for (int i = 0; i < static_cast<int>(cfgPtBinKaonPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinKaonPID.value[i] && candidate.pt() < cfgPtBinKaonPID.value[i + 1]) {
          currentPtBinKa = i;
          break;
        }
      }
    }
    if (cfgOpenPtRangedITSnSigmacutPr) {
      for (int i = 0; i < static_cast<int>(cfgPtBinProtonPID.value.size()) - 1; ++i) {
        if (candidate.pt() >= cfgPtBinProtonPID.value[i] && candidate.pt() < cfgPtBinProtonPID.value[i + 1]) {
          currentPtBinPr = i;
          break;
        }
      }
    }
    nSigmaITSToUse[0] = (currentPtBinPi == -1) ? cfgnSigmaCutITS.value[0] : cfgnSigmaITSPionPt.value[currentPtBinPi];
    nSigmaITSToUse[1] = (currentPtBinKa == -1) ? cfgnSigmaCutITS.value[1] : cfgnSigmaITSKaonPt.value[currentPtBinKa];
    nSigmaITSToUse[2] = (currentPtBinPr == -1) ? cfgnSigmaCutITS.value[2] : cfgnSigmaITSProtonPt.value[currentPtBinPr];
    switch (mode) {
      case 1: // For Pion
        if (!(std::abs(candidate.itsNSigmaPi()) < nSigmaITSToUse[0] && avgclssize > cfgAveClusSizeCoslMinPi && avgclssize < cfgAveClusSizeCoslMaxPi)) {
          return false;
        } else {
          return true;
        }
        break;

      case 2: // For Kaon
        if (!(std::abs(candidate.itsNSigmaKa()) < nSigmaITSToUse[1] && avgclssize > cfgAveClusSizeCoslMinKa && avgclssize < cfgAveClusSizeCoslMaxKa)) {
          return false;
        } else {
          return true;
        }
        break;

      case 3: // For Proton
        if (!(std::abs(candidate.itsNSigmaPr()) < nSigmaITSToUse[2] && avgclssize > cfgAveClusSizeCoslMinPr && avgclssize < cfgAveClusSizeCoslMaxPr)) {
          return false;
        } else {
          return true;
        }
        break;
    }
    return false;
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
    AxisSpec axisEtaPID{cfgaxisetaPID, "#it{#eta}"};
    AxisSpec axisClusterSize{cfgaxisAverClusterCosl, "<ITS Cluster Size> x <cos(#lambda)>"};
    AxisSpec axisClusterSizenSigma{cfgaxisAverClusterCoslnSigma, "<ITS Cluster Size> x <cos(#lambda)>"};
    AxisSpec axisPhi = {100, 0, 2.1 * constants::math::PI, "#phi"};
    if (!cfgQuietMode) {
      // TH3D NSigmaTPC,NSigmaTOF,NSigmaITS combo vs pt(if necessary for whole centrality)
      if (cfgOpenPlotnSigmaTOFITSPt) {
        histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_cross_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_cross_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_cross_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaITS, axisPtPID}});
        }
      }
      if (cfgOpenPlotnSigmaTOFTPCPt) {
        histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_cross_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_cross_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_cross_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Pi"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Ka"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Pr"), "", {HistType::kTH3F, {axisnSigmaTOF, axisnSigmaTPC, axisPtPID}});
        }
      }
      if (cfgOpenPlotnSigmaITSTPCPt) {
        histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_Pi"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_Ka"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_Pr"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_cross_Pi"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_cross_Ka"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_cross_Pr"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Pi"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Ka"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
          histosQA.add(Form("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Pr"), "", {HistType::kTH3F, {axisnSigmaITS, axisnSigmaTPC, axisPtPID}});
        }
      }
      if (cfgOpenPlotAverClus) {
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_cross_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Pi"), "", {HistType::kTH1F, {axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Ka"), "", {HistType::kTH1F, {axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_AfterITS_Pr"), "", {HistType::kTH1F, {axisClusterSize}});
        }
      }
      if (cfgOpenPlotAverClusP) {
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_cross_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pi"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Ka"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pr"), "", {HistType::kTH2F, {axisPPID, axisClusterSize}});
        }
      }
      if (cfgOpenPlotAverClusnSigmaTPC) {
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
        histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
          histosQA.add(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTPC, axisClusterSizenSigma}});
        }
      }
      if (cfgOpenPlotPhiDis) {
        histosQA.add(Form("QA/PID/histPhi_Dis_Pi"), "", {HistType::kTH1F, {axisPhi}});
        histosQA.add(Form("QA/PID/histPhi_Dis_Ka"), "", {HistType::kTH1F, {axisPhi}});
        histosQA.add(Form("QA/PID/histPhi_Dis_Pr"), "", {HistType::kTH1F, {axisPhi}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pi"), "", {HistType::kTH1F, {axisPhi}});
          histosQA.add(Form("QA/PID/histPhi_Dis_cross_Ka"), "", {HistType::kTH1F, {axisPhi}});
          histosQA.add(Form("QA/PID/histPhi_Dis_cross_Pr"), "", {HistType::kTH1F, {axisPhi}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Pi"), "", {HistType::kTH1F, {axisPhi}});
          histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Ka"), "", {HistType::kTH1F, {axisPhi}});
          histosQA.add(Form("QA/PID/histPhi_Dis_AfterITS_Pr"), "", {HistType::kTH1F, {axisPhi}});
        }
      }
      if (cfgOpenPlotPhiDisPtEta) {
        histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_Pi"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
        histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_Ka"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
        histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_Pr"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
        if (cfgOpenCrossTrackQAPlots) {
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_cross_Pi"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_cross_Ka"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_cross_Pr"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
        }
        if (cfgOpenITSCut) {
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Pi"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Ka"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
          histosQA.add(Form("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Pr"), "", {HistType::kTH3F, {axisPhi, axisPtPID, axisEtaPID}});
        }
      }
      // some basic plots should be ploted (except for the quite mode)
      // nSigma TPC TOF ITS combo plots
      histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
      if (cfgOpenCrossTrackQAPlots) {
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_cross_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
      }
      if (cfgOpenITSCut) {
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaITS, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Pi"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Ka"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Pr"), "", {HistType::kTH2F, {axisnSigmaTOF, axisnSigmaTPC}});
      }
      // nSigma TPC TOF ITS signle and some simple QA plots
      histosQA.add(Form("QA/PID/histdEdxTPC_All"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
      histosQA.add(Form("QA/PID/histdEdxTPC_Pi"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
      histosQA.add(Form("QA/PID/histdEdxTPC_Ka"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
      histosQA.add(Form("QA/PID/histdEdxTPC_Pr"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Pi"), "", {HistType::kTH1F, {axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Ka"), "", {HistType::kTH1F, {axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Pr"), "", {HistType::kTH1F, {axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
      histosQA.add(Form("QA/PID/histnSigma_com_Pi"), "", {HistType::kTH1F, {axisnSigmaCom}});
      histosQA.add(Form("QA/PID/histnSigma_com_Ka"), "", {HistType::kTH1F, {axisnSigmaCom}});
      histosQA.add(Form("QA/PID/histnSigma_com_Pr"), "", {HistType::kTH1F, {axisnSigmaCom}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Pi"), "", {HistType::kTH1F, {axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Ka"), "", {HistType::kTH1F, {axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Pr"), "", {HistType::kTH1F, {axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Pi"), "", {HistType::kTH1F, {axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Ka"), "", {HistType::kTH1F, {axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Pr"), "", {HistType::kTH1F, {axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
      histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
      if (cfgOpenCrossTrackQAPlots) {
        histosQA.add(Form("QA/PID/histdEdxTPC_cross_Pi"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
        histosQA.add(Form("QA/PID/histdEdxTPC_cross_Ka"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
        histosQA.add(Form("QA/PID/histdEdxTPC_cross_Pr"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_cross_Pi"), "", {HistType::kTH1F, {axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_cross_Ka"), "", {HistType::kTH1F, {axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_cross_Pr"), "", {HistType::kTH1F, {axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_cross_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_cross_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_cross_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_cross_Pi"), "", {HistType::kTH1F, {axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_cross_Ka"), "", {HistType::kTH1F, {axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_cross_Pr"), "", {HistType::kTH1F, {axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_Pt_cross_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Pi"), "", {HistType::kTH1F, {axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Ka"), "", {HistType::kTH1F, {axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_cross_Pr"), "", {HistType::kTH1F, {axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_cross_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_cross_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_cross_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
      }
      if (cfgOpenITSCut) {
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_AfterITS_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_AfterITS_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TOF_Pt_AfterITS_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTOF}});
        histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_AfterITS_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_AfterITS_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_TPC_Pt_AfterITS_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaTPC}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_AfterITS_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_AfterITS_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
        histosQA.add(Form("QA/PID/histnSigma_ITS_Pt_AfterITS_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigmaITS}});
      }
      // plots for TPC-ITS contamination (whole centrality)
      if (cfgOpenDetailPlotsTPCITSContaimination) {
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_Before"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        if (cfgOpenITSCut) {
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisnSigmaITS, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPi_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosKa_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_PosPr_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPi_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegKa_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
          histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSIgmaTPC_Pt_NegPr_After"), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {axisnSigmaTPC, axisClusterSizenSigma, axisPtPID}});
        }
      }
    }
  }
  Produces<aod::Flags> pidCmeTable;
  Produces<aod::PidInfo> pidInfoTable;
  void process(TracksPID const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<TracksPID,
                                        aod::pidits::ITSNSigmaPi,
                                        aod::pidits::ITSNSigmaKa,
                                        aod::pidits::ITSNSigmaPr>(tracks);
    int8_t pidFlag;
    for (const auto& track : tracksWithITSPid) {
      const float averClusSizeCosl = averageClusterSizeCosl(track.itsClusterSizes(), track.eta());
      if (!selTrackPid(track)) {
        pidFlag = -1;
      } else {
        if (!cfgQuietMode) {
          histosQA.fill(HIST("QA/PID/histdEdxTPC_All"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
        }
        pidFlag = selectionPidtpctof(track);
        if (!(pidFlag == 0 || pidFlag == -1)) {
          // First fill ITS uncut plots
          if (!cfgQuietMode) {
            if ((pidFlag == 1) || (pidFlag == 7) || (pidFlag == 8) || (pidFlag == 10)) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Pi"), track.itsNSigmaPi(), track.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_Pi"), track.tofNSigmaPi(), track.itsNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_Pi"), track.tofNSigmaPi(), track.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_Pi"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pi"), track.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_Pi"), track.pt(), track.tpcNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_com_Pi"), std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()));
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pi"), track.tofNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_Pi"), track.pt(), track.tofNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pi"), track.itsNSigmaPi());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_Pi"), track.pt(), track.itsNSigmaPi());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_Pi"), track.tofNSigmaPi(), track.itsNSigmaPi(), track.pt());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_Pi"), track.tofNSigmaPi(), track.tpcNSigmaPi(), track.pt());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_Pi"), track.itsNSigmaPi(), track.tpcNSigmaPi(), track.pt());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pi"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pi"), track.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pi"), track.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_Pi"), track.phi(), track.pt(), track.eta());
              }
              if (cfgOpenDetailPlotsTPCITSContaimination) {
                if (track.sign() > 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
                } else if (track.sign() < 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_Before"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
                }
              }
            }
            if ((pidFlag == 2) || (pidFlag == 7) || (pidFlag == 9) || (pidFlag == 10)) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Ka"), track.itsNSigmaKa(), track.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_Ka"), track.tofNSigmaKa(), track.itsNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_Ka"), track.tofNSigmaKa(), track.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_Ka"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Ka"), track.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_Ka"), track.pt(), track.tpcNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_com_Ka"), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()));
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Ka"), track.tofNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_Ka"), track.pt(), track.tofNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Ka"), track.itsNSigmaKa());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_Ka"), track.pt(), track.itsNSigmaKa());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_Ka"), track.tofNSigmaKa(), track.itsNSigmaKa(), track.pt());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_Ka"), track.tofNSigmaKa(), track.tpcNSigmaKa(), track.pt());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_Ka"), track.itsNSigmaKa(), track.tpcNSigmaKa(), track.pt());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Ka"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Ka"), track.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Ka"), track.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_Ka"), track.phi(), track.pt(), track.eta());
              }
              if (cfgOpenDetailPlotsTPCITSContaimination) {
                if (track.sign() > 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before"), track.tpcNSigmaKa(), track.itsNSigmaKa(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_Before"), track.tpcNSigmaKa(), averClusSizeCosl, track.pt());
                } else if (track.sign() < 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before"), track.tpcNSigmaKa(), track.itsNSigmaKa(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_Before"), track.tpcNSigmaKa(), averClusSizeCosl, track.pt());
                }
              }
            }
            if ((pidFlag == 3) || (pidFlag == 8) || (pidFlag == 9) || (pidFlag == 10)) {
              histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_Pr"), track.itsNSigmaPr(), track.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_Pr"), track.tofNSigmaPr(), track.itsNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_Pr"), track.tofNSigmaPr(), track.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histdEdxTPC_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pr"), track.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_Pr"), track.pt(), track.tpcNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_com_Pr"), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pr"), track.tofNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_Pr"), track.pt(), track.tofNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pr"), track.itsNSigmaPr());
              histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_Pr"), track.pt(), track.itsNSigmaPr());
              if (cfgOpenPlotnSigmaTOFITSPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_Pr"), track.tofNSigmaPr(), track.itsNSigmaPr(), track.pt());
              }
              if (cfgOpenPlotnSigmaTOFTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_Pr"), track.tofNSigmaPr(), track.tpcNSigmaPr(), track.pt());
              }
              if (cfgOpenPlotnSigmaITSTPCPt) {
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_Pr"), track.itsNSigmaPr(), track.tpcNSigmaPr(), track.pt());
              }
              if (cfgOpenPlotAverClus) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_Pr"), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusP) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_Pr"), track.p(), averClusSizeCosl);
              }
              if (cfgOpenPlotAverClusnSigmaTPC) {
                histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pr"), track.tpcNSigmaPr(), averClusSizeCosl);
              }
              if (cfgOpenPlotPhiDis) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pr"), track.phi());
              }
              if (cfgOpenPlotPhiDisPtEta) {
                histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_Pr"), track.phi(), track.pt(), track.eta());
              }
              if (cfgOpenDetailPlotsTPCITSContaimination) {
                if (track.sign() > 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before"), track.tpcNSigmaPr(), track.itsNSigmaPr(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_Before"), track.tpcNSigmaPr(), averClusSizeCosl, track.pt());
                } else if (track.sign() < 0) {
                  histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before"), track.tpcNSigmaPr(), track.itsNSigmaPr(), track.pt());
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_Before"), track.tpcNSigmaPr(), averClusSizeCosl, track.pt());
                }
              }
            }
          }
          // Second proform ITS cut
          if (cfgOpenITSCut) {
            int idx = -1;
            switch (pidFlag) {
              case 1:
                if (!selectionITS(track, 1, averClusSizeCosl)) {
                  pidFlag = 4;
                }
                break;
              case 2:
                if (!selectionITS(track, 2, averClusSizeCosl)) {
                  pidFlag = 5;
                }
                break;
              case 3:
                if (!selectionITS(track, 3, averClusSizeCosl)) {
                  pidFlag = 6;
                }
                break;
              case 7:
                idx = (selectionITS(track, 1, averClusSizeCosl) << 1) | selectionITS(track, 2, averClusSizeCosl);
                switch (idx) {
                  case 0:
                    pidFlag = 11;
                    break;

                  case 1:
                    pidFlag = 2;
                    break;

                  case 2:
                    pidFlag = 1;
                    break;
                }
                break;
              case 8:
                idx = (selectionITS(track, 1, averClusSizeCosl) << 1) | selectionITS(track, 3, averClusSizeCosl);
                switch (idx) {
                  case 0:
                    pidFlag = 12;
                    break;

                  case 1:
                    pidFlag = 3;
                    break;

                  case 2:
                    pidFlag = 1;
                    break;
                }
                break;
              case 9:
                idx = (selectionITS(track, 2, averClusSizeCosl) << 1) | selectionITS(track, 3, averClusSizeCosl);
                switch (idx) {
                  case 0:
                    pidFlag = 13;
                    break;

                  case 1:
                    pidFlag = 3;
                    break;

                  case 2:
                    pidFlag = 2;
                    break;
                }
                break;
              case 10:
                idx = (selectionITS(track, 1, averClusSizeCosl) << 2) | (selectionITS(track, 2, averClusSizeCosl) << 1) | selectionITS(track, 3, averClusSizeCosl);
                switch (idx) {
                  case 0:
                    pidFlag = 14;
                    break;

                  case 1:
                    pidFlag = 3;
                    break;

                  case 2:
                    pidFlag = 2;
                    break;

                  case 3:
                    pidFlag = 9;
                    break;

                  case 4:
                    pidFlag = 1;
                    break;

                  case 5:
                    pidFlag = 8;
                    break;

                  case 6:
                    pidFlag = 7;
                    break;
                }
                break;
            }
          }
          // Third Fill ITS cut plots
          if (cfgOpenITSCut) {
            if (!cfgQuietMode) {
              if ((pidFlag == 1) || (pidFlag == 7) || (pidFlag == 8) || (pidFlag == 10)) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pi"), track.itsNSigmaPi(), track.tpcNSigmaPi());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Pi"), track.tofNSigmaPi(), track.itsNSigmaPi());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Pi"), track.tofNSigmaPi(), track.tpcNSigmaPi());
                histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_AfterITS_Pi"), track.pt(), track.tpcNSigmaPi());
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_AfterITS_Pi"), track.pt(), track.tofNSigmaPi());
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_AfterITS_Pi"), track.pt(), track.itsNSigmaPi());
                if (cfgOpenPlotnSigmaTOFITSPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Pi"), track.tofNSigmaPi(), track.itsNSigmaPi(), track.pt());
                }
                if (cfgOpenPlotnSigmaTOFTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Pi"), track.tofNSigmaPi(), track.tpcNSigmaPi(), track.pt());
                }
                if (cfgOpenPlotnSigmaITSTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Pi"), track.itsNSigmaPi(), track.tpcNSigmaPi(), track.pt());
                }
                if (cfgOpenPlotAverClus) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Pi"), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusP) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pi"), track.p(), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusnSigmaTPC) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pi"), track.tpcNSigmaPi(), averClusSizeCosl);
                }
                if (cfgOpenPlotPhiDis) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Pi"), track.phi());
                }
                if (cfgOpenPlotPhiDisPtEta) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Pi"), track.phi(), track.pt(), track.eta());
                }
                if (cfgOpenDetailPlotsTPCITSContaimination) {
                  if (track.sign() > 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
                  } else if (track.sign() < 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After"), track.tpcNSigmaPi(), track.itsNSigmaPi(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_After"), track.tpcNSigmaPi(), averClusSizeCosl, track.pt());
                  }
                }
              }
              if ((pidFlag == 2) || (pidFlag == 7) || (pidFlag == 9) || (pidFlag == 10)) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Ka"), track.itsNSigmaKa(), track.tpcNSigmaKa());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Ka"), track.tofNSigmaKa(), track.itsNSigmaKa());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Ka"), track.tofNSigmaKa(), track.tpcNSigmaKa());
                histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_AfterITS_Ka"), track.pt(), track.tpcNSigmaKa());
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_AfterITS_Ka"), track.pt(), track.tofNSigmaKa());
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_AfterITS_Ka"), track.pt(), track.itsNSigmaKa());
                if (cfgOpenPlotnSigmaTOFITSPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Ka"), track.tofNSigmaKa(), track.itsNSigmaKa(), track.pt());
                }
                if (cfgOpenPlotnSigmaTOFTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Ka"), track.tofNSigmaKa(), track.tpcNSigmaKa(), track.pt());
                }
                if (cfgOpenPlotnSigmaITSTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Ka"), track.itsNSigmaKa(), track.tpcNSigmaKa(), track.pt());
                }
                if (cfgOpenPlotAverClus) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Ka"), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusP) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Ka"), track.p(), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusnSigmaTPC) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Ka"), track.tpcNSigmaKa(), averClusSizeCosl);
                }
                if (cfgOpenPlotPhiDis) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Ka"), track.phi());
                }
                if (cfgOpenPlotPhiDisPtEta) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Ka"), track.phi(), track.pt(), track.eta());
                }
                if (cfgOpenDetailPlotsTPCITSContaimination) {
                  if (track.sign() > 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After"), track.tpcNSigmaKa(), track.itsNSigmaKa(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_After"), track.tpcNSigmaKa(), averClusSizeCosl, track.pt());
                  } else if (track.sign() < 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After"), track.tpcNSigmaKa(), track.itsNSigmaKa(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_After"), track.tpcNSigmaKa(), averClusSizeCosl, track.pt());
                  }
                }
              }
              if ((pidFlag == 3) || (pidFlag == 8) || (pidFlag == 9) || (pidFlag == 10)) {
                histosQA.fill(HIST("QA/PID/histnSigmaITS_nSigmaTPC_AfterITS_Pr"), track.itsNSigmaPr(), track.tpcNSigmaPr());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaITS_AfterITS_Pr"), track.tofNSigmaPr(), track.itsNSigmaPr());
                histosQA.fill(HIST("QA/PID/histnSigmaTOF_nSigmaTPC_AfterITS_Pr"), track.tofNSigmaPr(), track.tpcNSigmaPr());
                histosQA.fill(HIST("QA/PID/histnSigma_TPC_Pt_AfterITS_Pr"), track.pt(), track.tpcNSigmaPr());
                histosQA.fill(HIST("QA/PID/histnSigma_TOF_Pt_AfterITS_Pr"), track.pt(), track.tofNSigmaPr());
                histosQA.fill(HIST("QA/PID/histnSigma_ITS_Pt_AfterITS_Pr"), track.pt(), track.itsNSigmaPr());
                if (cfgOpenPlotnSigmaTOFITSPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_ITS_Pt_AfterITS_Pr"), track.tofNSigmaPr(), track.itsNSigmaPr(), track.pt());
                }
                if (cfgOpenPlotnSigmaTOFTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_TOF_TPC_Pt_AfterITS_Pr"), track.tofNSigmaPr(), track.tpcNSigmaPr(), track.pt());
                }
                if (cfgOpenPlotnSigmaITSTPCPt) {
                  histosQA.fill(HIST("QA/PID/histnSigma_ITS_TPC_Pt_AfterITS_Pr"), track.itsNSigmaPr(), track.tpcNSigmaPr(), track.pt());
                }
                if (cfgOpenPlotAverClus) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_AfterITS_Pr"), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusP) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_P_AfterITS_Pr"), track.p(), averClusSizeCosl);
                }
                if (cfgOpenPlotAverClusnSigmaTPC) {
                  histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_AfterITS_Pr"), track.tpcNSigmaPr(), averClusSizeCosl);
                }
                if (cfgOpenPlotPhiDis) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_AfterITS_Pr"), track.phi());
                }
                if (cfgOpenPlotPhiDisPtEta) {
                  histosQA.fill(HIST("QA/PID/histPhi_Dis_Pt_Eta_AfterITS_Pr"), track.phi(), track.pt(), track.eta());
                }
                if (cfgOpenDetailPlotsTPCITSContaimination) {
                  if (track.sign() > 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After"), track.tpcNSigmaPr(), track.itsNSigmaPr(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_After"), track.tpcNSigmaPr(), averClusSizeCosl, track.pt());
                  } else if (track.sign() < 0) {
                    histosQA.fill(HIST("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After"), track.tpcNSigmaPr(), track.itsNSigmaPr(), track.pt());
                    histosQA.fill(HIST("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_After"), track.tpcNSigmaPr(), averClusSizeCosl, track.pt());
                  }
                }
              }
            }
          }
        }
      }
      pidCmeTable(pidFlag);
      pidInfoTable(averClusSizeCosl, track.itsNSigmaPi(), track.itsNSigmaKa(), track.itsNSigmaPr(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
    }
  }
};

struct QAProcessCent {
  HistogramRegistry histosQA{"histosQAwithcent", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<std::vector<int>> cfgCentralitybinsforQA{"cfgCentralitybinsforQA", {0, 30, 60}, "Centrality bins for track phi and TPC_ITS matching check"};
  Configurable<bool> cfgOpenDetailPlots{"cfgOpenDetailPlots", true, "open detail TH3D plots for nSigmaTPC-ITS Pt-eta-Phi nSigmaITS-clustersize"};
  Configurable<bool> cfgOpenITSafter{"cfgOpenITSafter", true, "open check for after ITS cut check(if used ITScut in table producer it need else close to save cpu usage)"};
  Configurable<bool> cfgOpenPtEtaPhi{"cfgOpenPtEtaPhi", true, "open pt-#eta-#phi PID QA  (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenITSTPCnSigma{"cfgOpenITSTPCnSigma", true, "open ITS-TPC nSigma QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenClusSizenSigmaTPC{"cfgOpenClusSizenSigmaTPC", true, "open ITSClustersize-TPCnsigma QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenPi{"cfgOpenPi", true, "open Pion QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenKa{"cfgOpenKa", true, "open Kaon QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenPr{"cfgOpenPr", true, "open Proton QA (Optional for limited memory usage)"};
  Configurable<bool> cfgOpenTOFITSnSigma{"cfgOpenTOFITSnSigma", false, "open TOF-ITS nsigma 2D plots vs pt at a certain centrality"};
  Configurable<bool> cfgOpenTOFTPCnSigma{"cfgOpenTOFTPCnSigma", false, "open TOF-TPC nsigma 2D plots vs pt at a certain centrality"};
  ConfigurableAxis cfgaxisetaPID{"cfgaxisetaPID", {90, -0.9, 0.9}, "Binning for Pt QA"};
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma ITS"};
  ConfigurableAxis cfgnSigmaBinsTOF{"cfgnSigmaBinsTOF", {200, -5.f, 5.f}, "Binning for n sigma TOF"};
  ConfigurableAxis cfgaxisAverClusterCoslnSigma{"cfgaxisAverClusterCoslnSigma", {50, 0, 5}, "Binning for average cluster size x cos(#lambda) vs nSigam"};

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
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtPosPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFITSPtNegPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegPiBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegKaBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegPrBeforeCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegPiAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegKaAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtPosPrAfterCen;
  std::vector<std::shared_ptr<TH3>> vhistnSigmaTOFTPCPtNegPrAfterCen;
  Filter trackPIDfilter = aod::cme_track_pid_columns::nPidFlag > (int8_t)0;
  void init(InitContext const&)
  {
    AxisSpec axisPhicme = {100, 0, 2.1 * constants::math::PI, "#phi"};
    // Additional QA histograms for PID
    if (cfgOpenDetailPlots) {
      for (int i = 0; i < static_cast<int>(cfgCentralitybinsforQA.value.size()) - 1; ++i) {
        if (cfgOpenPtEtaPhi) {
          if (cfgOpenPi) {
            auto hPhiPtEtaPosPi = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPi_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            auto hPhiPtEtaNegPi = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPi_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            vhistPhiPtEtaPosPiCen.push_back(std::move(hPhiPtEtaPosPi));
            vhistPhiPtEtaNegPiCen.push_back(std::move(hPhiPtEtaNegPi));
          }
          if (cfgOpenKa) {
            auto hPhiPtEtaPosKa = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosKa_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            auto hPhiPtEtaNegKa = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegKa_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            vhistPhiPtEtaPosKaCen.push_back(std::move(hPhiPtEtaPosKa));
            vhistPhiPtEtaNegKaCen.push_back(std::move(hPhiPtEtaNegKa));
          }
          if (cfgOpenPr) {
            auto hPhiPtEtaPosPr = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_PosPr_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            auto hPhiPtEtaNegPr = histosQA.add<TH3>(Form("QA/PID/histPhi_Pt_Eta_NegPr_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";#phi;#p_{t};#eta", {HistType::kTH3F, {axisPhicme, cfgaxisptPID, cfgaxisetaPID}});
            vhistPhiPtEtaPosPrCen.push_back(std::move(hPhiPtEtaPosPr));
            vhistPhiPtEtaNegPrCen.push_back(std::move(hPhiPtEtaNegPr));
          }
        }
        if (cfgOpenITSTPCnSigma) {
          if (cfgOpenPi) {
            auto hnSigmaITSTPCPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaITSTPCPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaITSTPCPtPosPiBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosPiBefore));
            vhistnSigmaITSTPCPtNegPiBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegPiBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaITSTPCPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaITSTPCPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaITSTPCPtPosPiAfterCen.push_back(std::move(hnSigmaITSTPCPtPosPiAfter));
              vhistnSigmaITSTPCPtNegPiAfterCen.push_back(std::move(hnSigmaITSTPCPtNegPiAfter));
            }
          }
          if (cfgOpenKa) {
            auto hnSigmaITSTPCPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaITSTPCPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaITSTPCPtPosKaBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosKaBefore));
            vhistnSigmaITSTPCPtNegKaBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegKaBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaITSTPCPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaITSTPCPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaITSTPCPtPosKaAfterCen.push_back(std::move(hnSigmaITSTPCPtPosKaAfter));
              vhistnSigmaITSTPCPtNegKaAfterCen.push_back(std::move(hnSigmaITSTPCPtNegKaAfter));
            }
          }
          if (cfgOpenPr) {
            auto hnSigmaITSTPCPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaITSTPCPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaITSTPCPtPosPrBeforeCen.push_back(std::move(hnSigmaITSTPCPtPosPrBefore));
            vhistnSigmaITSTPCPtNegPrBeforeCen.push_back(std::move(hnSigmaITSTPCPtNegPrBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaITSTPCPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaITSTPCPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaITS_TPC_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaITSTPCPtPosPrAfterCen.push_back(std::move(hnSigmaITSTPCPtPosPrAfter));
              vhistnSigmaITSTPCPtNegPrAfterCen.push_back(std::move(hnSigmaITSTPCPtNegPrAfter));
            }
          }
        }
        if (cfgOpenClusSizenSigmaTPC) {
          if (cfgOpenPi) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPiBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPiBefore));
            if (cfgOpenITSafter) {
              auto hAverClusterSizeCoslnSigmaTPCPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              auto hAverClusterSizeCoslnSigmaTPCPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              vhistAverClusterSizeCoslnSigmaTPCPtPosPiAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPiAfter));
              vhistAverClusterSizeCoslnSigmaTPCPtNegPiAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPiAfter));
            }
          }
          if (cfgOpenKa) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosKaBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegKaBefore));
            if (cfgOpenITSafter) {
              auto hAverClusterSizeCoslnSigmaTPCPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              auto hAverClusterSizeCoslnSigmaTPCPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              vhistAverClusterSizeCoslnSigmaTPCPtPosKaAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosKaAfter));
              vhistAverClusterSizeCoslnSigmaTPCPtNegKaAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegKaAfter));
            }
          }
          if (cfgOpenPr) {
            auto hAverClusterSizeCoslnSigmaTPCPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            auto hAverClusterSizeCoslnSigmaTPCPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
            vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPrBefore));
            vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPrBefore));
            if (cfgOpenITSafter) {
              auto hAverClusterSizeCoslnSigmaTPCPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              auto hAverClusterSizeCoslnSigmaTPCPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histAverClusterSizeCosl_nSigmaTPC_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TPC};<ITS Cluster Size> x <cos(#lambda)>;#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgaxisAverClusterCoslnSigma, cfgaxisptPID}});
              vhistAverClusterSizeCoslnSigmaTPCPtPosPrAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtPosPrAfter));
              vhistAverClusterSizeCoslnSigmaTPCPtNegPrAfterCen.push_back(std::move(hAverClusterSizeCoslnSigmaTPCPtNegPrAfter));
            }
          }
        }
        if (cfgOpenTOFITSnSigma) {
          if (cfgOpenPi) {
            auto hnSigmaTOFITSPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaTOFITSPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaTOFITSPtPosPiBeforeCen.push_back(std::move(hnSigmaTOFITSPtPosPiBefore));
            vhistnSigmaTOFITSPtNegPiBeforeCen.push_back(std::move(hnSigmaTOFITSPtNegPiBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFITSPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaTOFITSPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaTOFITSPtPosPiAfterCen.push_back(std::move(hnSigmaTOFITSPtPosPiAfter));
              vhistnSigmaTOFITSPtNegPiAfterCen.push_back(std::move(hnSigmaTOFITSPtNegPiAfter));
            }
          }
          if (cfgOpenKa) {
            auto hnSigmaTOFITSPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaTOFITSPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaTOFITSPtPosKaBeforeCen.push_back(std::move(hnSigmaTOFITSPtPosKaBefore));
            vhistnSigmaTOFITSPtNegKaBeforeCen.push_back(std::move(hnSigmaTOFITSPtNegKaBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFITSPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaTOFITSPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaTOFITSPtPosKaAfterCen.push_back(std::move(hnSigmaTOFITSPtPosKaAfter));
              vhistnSigmaTOFITSPtNegKaAfterCen.push_back(std::move(hnSigmaTOFITSPtNegKaAfter));
            }
          }
          if (cfgOpenPr) {
            auto hnSigmaTOFITSPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            auto hnSigmaTOFITSPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
            vhistnSigmaTOFITSPtPosPrBeforeCen.push_back(std::move(hnSigmaTOFITSPtPosPrBefore));
            vhistnSigmaTOFITSPtNegPrBeforeCen.push_back(std::move(hnSigmaTOFITSPtNegPrBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFITSPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              auto hnSigmaTOFITSPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_ITS_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
              vhistnSigmaTOFITSPtPosPrAfterCen.push_back(std::move(hnSigmaTOFITSPtPosPrAfter));
              vhistnSigmaTOFITSPtNegPrAfterCen.push_back(std::move(hnSigmaTOFITSPtNegPrAfter));
            }
          }
        }
        if (cfgOpenTOFTPCnSigma) {
          if (cfgOpenPi) {
            auto hnSigmaTOFTPCPtPosPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            auto hnSigmaTOFTPCPtNegPiBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegPi_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            vhistnSigmaTOFTPCPtPosPiBeforeCen.push_back(std::move(hnSigmaTOFTPCPtPosPiBefore));
            vhistnSigmaTOFTPCPtNegPiBeforeCen.push_back(std::move(hnSigmaTOFTPCPtNegPiBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFTPCPtPosPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              auto hnSigmaTOFTPCPtNegPiAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegPi_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              vhistnSigmaTOFTPCPtPosPiAfterCen.push_back(std::move(hnSigmaTOFTPCPtPosPiAfter));
              vhistnSigmaTOFTPCPtNegPiAfterCen.push_back(std::move(hnSigmaTOFTPCPtNegPiAfter));
            }
          }
          if (cfgOpenKa) {
            auto hnSigmaTOFTPCPtPosKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            auto hnSigmaTOFTPCPtNegKaBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegKa_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            vhistnSigmaTOFTPCPtPosKaBeforeCen.push_back(std::move(hnSigmaTOFTPCPtPosKaBefore));
            vhistnSigmaTOFTPCPtNegKaBeforeCen.push_back(std::move(hnSigmaTOFTPCPtNegKaBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFTPCPtPosKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              auto hnSigmaTOFTPCPtNegKaAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegKa_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              vhistnSigmaTOFTPCPtPosKaAfterCen.push_back(std::move(hnSigmaTOFTPCPtPosKaAfter));
              vhistnSigmaTOFTPCPtNegKaAfterCen.push_back(std::move(hnSigmaTOFTPCPtNegKaAfter));
            }
          }
          if (cfgOpenPr) {
            auto hnSigmaTOFTPCPtPosPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            auto hnSigmaTOFTPCPtNegPrBefore = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegPr_Before_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
            vhistnSigmaTOFTPCPtPosPrBeforeCen.push_back(std::move(hnSigmaTOFTPCPtPosPrBefore));
            vhistnSigmaTOFTPCPtNegPrBeforeCen.push_back(std::move(hnSigmaTOFTPCPtNegPrBefore));
            if (cfgOpenITSafter) {
              auto hnSigmaTOFTPCPtPosPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_PosPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              auto hnSigmaTOFTPCPtNegPrAfter = histosQA.add<TH3>(Form("QA/PID/histnSigmaTOF_TPC_Pt_NegPr_After_Cen_%d_%d", cfgCentralitybinsforQA.value[i], cfgCentralitybinsforQA.value[i + 1]), ";n#sigma_{TOF};n#sigma_{TPC};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
              vhistnSigmaTOFTPCPtPosPrAfterCen.push_back(std::move(hnSigmaTOFTPCPtPosPrAfter));
              vhistnSigmaTOFTPCPtNegPrAfterCen.push_back(std::move(hnSigmaTOFTPCPtNegPrAfter));
            }
          }
        }
      }
    }
  }
  void process(CollisionPID::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::Flags, aod::PidInfo>> const& tracks)
  {
    if (cfgOpenDetailPlots) {
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
          int8_t pidFlag = trk.nPidFlag();
          if (cfgOpenPi) {
            if ((pidFlag == 1) || (pidFlag == 4) || (pidFlag == 7) || (pidFlag == 8) || (pidFlag == 10) || (pidFlag == 11) || (pidFlag == 12) || (pidFlag == 14)) {
              if (trk.sign() > 0) {
                if (!((pidFlag == 4) || (pidFlag == 11) || (pidFlag == 12) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosPiCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtPosPiAfterCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtPosPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiTPC(), trk.pt());
                }
              } else if (trk.sign() < 0) {
                if (!((pidFlag == 4) || (pidFlag == 11) || (pidFlag == 12) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegPiCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtNegPiAfterCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.nSigmaPiITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtNegPiBeforeCen[currentBin]->Fill(trk.nSigmaPiTOF(), trk.nSigmaPiTPC(), trk.pt());
                }
              }
            }
          }
          if (cfgOpenKa) {
            if ((pidFlag == 2) || (pidFlag == 5) || (pidFlag == 7) || (pidFlag == 9) || (pidFlag == 10) || (pidFlag == 11) || (pidFlag == 13) || (pidFlag == 14)) {
              if (trk.sign() > 0) {
                if (!((pidFlag == 5) || (pidFlag == 11) || (pidFlag == 13) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosKaCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtPosKaAfterCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtPosKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaTPC(), trk.pt());
                }
              } else if (trk.sign() < 0) {
                if (!((pidFlag == 5) || (pidFlag == 11) || (pidFlag == 13) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegKaCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtNegKaAfterCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.nSigmaKaITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtNegKaBeforeCen[currentBin]->Fill(trk.nSigmaKaTOF(), trk.nSigmaKaTPC(), trk.pt());
                }
              }
            }
          }
          if (cfgOpenPr) {
            if ((pidFlag == 3) || (pidFlag == 6) || (pidFlag == 8) || (pidFlag == 9) || (pidFlag == 10) || (pidFlag == 12) || (pidFlag == 13) || (pidFlag == 14)) {
              if (trk.sign() > 0) {
                if (!((pidFlag == 6) || (pidFlag == 12) || (pidFlag == 13) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaPosPrCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtPosPrAfterCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtPosPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrTPC(), trk.pt());
                }
              } else if (trk.sign() < 0) {
                if (!((pidFlag == 6) || (pidFlag == 12) || (pidFlag == 13) || (pidFlag == 14))) {
                  if (cfgOpenPtEtaPhi) {
                    vhistPhiPtEtaNegPrCen[currentBin]->Fill(trk.phi(), trk.pt(), trk.eta());
                  }
                  if (cfgOpenITSafter) {
                    if (cfgOpenITSTPCnSigma) {
                      vhistnSigmaITSTPCPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                    }
                    if (cfgOpenClusSizenSigmaTPC) {
                      vhistAverClusterSizeCoslnSigmaTPCPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                    }
                    if (cfgOpenTOFITSnSigma) {
                      vhistnSigmaTOFITSPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrITS(), trk.pt());
                    }
                    if (cfgOpenTOFTPCnSigma) {
                      vhistnSigmaTOFTPCPtNegPrAfterCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrTPC(), trk.pt());
                    }
                  }
                }
                if (cfgOpenITSTPCnSigma) {
                  vhistnSigmaITSTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
                }
                if (cfgOpenClusSizenSigmaTPC) {
                  vhistAverClusterSizeCoslnSigmaTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTPC(), trk.averClusterSizeCosl(), trk.pt());
                }
                if (cfgOpenTOFITSnSigma) {
                  vhistnSigmaTOFITSPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrITS(), trk.pt());
                }
                if (cfgOpenTOFTPCnSigma) {
                  vhistnSigmaTOFTPCPtNegPrBeforeCen[currentBin]->Fill(trk.nSigmaPrTOF(), trk.nSigmaPrTPC(), trk.pt());
                }
              }
            }
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
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma TPC"};

  ConfigurableAxis cfgaxissumpt{"cfgaxissumpt", {16, 0, 16}, "Binning for #gamma and #delta sum p_{t}(particle1 + particle2)"};
  ConfigurableAxis cfgaxisdeltaeta{"cfgaxisdeltaeta", {16, -1.6, 1.6}, "Binning for #gamma and #delta #eta(particle1 - particle2)"};
  ConfigurableAxis cfgaxisdeltapt{"cfgaxisdeltapt", {16, -8, 8}, "Binning for #gamma and #delta p_{t}(particle1 - particle2)"};

  Configurable<bool> cfgUseAdditionalEventCut{"cfgUseAdditionalEventCut", true, "Use additional event cut beyond sel8"};
  Configurable<bool> cfgOpenEvSelkIsGoodZvtxFT0vsPV{"cfgOpenEvSelkIsGoodZvtxFT0vsPV", true, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution"};
  Configurable<bool> cfgOpenEvSelkNoSameBunchPileup{"cfgOpenEvSelkNoSameBunchPileup", true, "rejects collisions which are associated with the same found-by-T0 bunch crossing"};
  Configurable<bool> cfgOpenEvSelkNoCollInTimeRangeStandard{"cfgOpenEvSelkNoCollInTimeRangeStandard", true, "no collisions in specified time range"};
  Configurable<bool> cfgOpenEvSelkIsGoodITSLayersAll{"cfgOpenEvSelkIsGoodITSLayersAll", true, "cut time intervals with dead ITS staves"};
  Configurable<bool> cfgOpenEvSelkNoCollInRofStandard{"cfgOpenEvSelkNoCollInRofStandard", true, "no other collisions in this Readout Frame with per-collision multiplicity above threshold"};
  Configurable<bool> cfgOpenEvSelkNoHighMultCollInPrevRof{"cfgOpenEvSelkNoHighMultCollInPrevRof", true, "veto an event if FT0C amplitude in previous ITS ROF is above threshold"};
  Configurable<bool> cfgOpenEvSelOccupancy{"cfgOpenEvSelOccupancy", true, "Occupancy cut"};
  Configurable<bool> cfgOpenEvSelMultCorrelationPVTracks{"cfgOpenEvSelMultCorrelationPVTracks", true, "Multiplicity correlation cut for PVtracks vs centrality(FT0C)"};
  Configurable<bool> cfgOpenEvSelMultCorrelationGlobalTracks{"cfgOpenEvSelMultCorrelationGlobalTracks", false, "Multiplicity correlation cut for Globaltracks vs centrality(FT0C)"};
  Configurable<bool> cfgOpenEvSelV0AT0ACut{"cfgOpenEvSelV0AT0ACut", true, "V0A T0A 5 sigma cut"};
  Configurable<bool> cfgOpenFullEventQA{"cfgOpenFullEventQA", true, "Open full QA plots for event QA"};
  Configurable<bool> cfgkOpenV2{"cfgkOpenV2", true, "open V2 plots"};
  Configurable<bool> cfgkOpenCME{"cfgkOpenCME", true, "open PID CME"};
  Configurable<bool> cfgkOpenCMEDifferential{"cfgkOpenCMEDifferential", true, "open Differential plot(#delta pt, #delta eta, #average pt) for cme #delta and #gamma"};
  Configurable<bool> cfgkOpenDeltaPt{"cfgkOpenDeltaPt", true, "open CME Differential #Delta Pt"};
  Configurable<bool> cfgkOpenDeltaEta{"cfgkOpenDeltaEta", true, "open CME Differential #Delta Eta"};
  Configurable<bool> cfgkOpenAveragePt{"cfgkOpenAveragePt", true, "open CME Differential #Average Pt"};
  Configurable<bool> cfgkOpenPiPi{"cfgkOpenPiPi", true, "open Pi-Pi"};
  Configurable<bool> cfgkOpenKaKa{"cfgkOpenKaKa", true, "open Ka-Ka"};
  Configurable<bool> cfgkOpenPrPr{"cfgkOpenPrPr", true, "open Pr-Pr"};
  Configurable<bool> cfgkOpenPiKa{"cfgkOpenPiKa", true, "open Pi-Ka"};
  Configurable<bool> cfgkOpenPiPr{"cfgkOpenPiPr", true, "open Pi-Pr"};
  Configurable<bool> cfgkOpenKaPr{"cfgkOpenKaPr", true, "open Ka-Pr"};
  Configurable<bool> cfgkOpenHaHa{"cfgkOpenHaHa", false, "open Ha-Ha"};
  Configurable<bool> cfgkOpenSsOsCrossCheck{"cfgkOpenSsOsCrossCheck", false, "open check for matter an antimatter #gamma#delta"};
  Configurable<bool> cfgkOpenTPCITSPurityCut{"cfgkOpenTPCITSPurityCut", true, "open ITS-TPC purity cut"};
  Configurable<bool> cfgkOpenTPCITSPurityCutQA{"cfgkOpenTPCITSPurityCutQA", true, "open ITS-TPC purity cut QA plots"};
  Configurable<bool> cfgkOpenDebugPIDCME{"cfgkOpenDebugPIDCME", false, "open pidcme workflow debug mode"};
  Configurable<bool> cfgOpenPlotITSNcls{"cfgOpenPlotITSNcls", true, "open QA for overall ITSNcls distribution"};
  Configurable<bool> cfgOpenPlotITSNclsPtCent{"cfgOpenPlotITSNclsPtCent", false, "open QA for ITSNcls distribution vs centality and pt"};
  Configurable<bool> cfgOpenTPCNcls{"cfgOpenTPCNcls", true, "open QA for overall TPCNcls distribution"};
  Configurable<bool> cfgOpenTPCNclsPtCent{"cfgOpenTPCNclsPtCent", false, "open QA for TPCNcls distribution vs centality and pt"};
  Configurable<bool> cfgOpenCustomTrackCutAssurance{"cfgOpenCustomTrackCutAssurance", true, "Assure track using for v2 and cme pass the custom track cuts"};

  Configurable<std::vector<int>> cfgITSPurityCen{"cfgITSPurityCen", {20, 30}, "ITS purity cut centrality"};
  Configurable<std::vector<float>> cfgPtPrCut{"cfgPtPrCut", {0.5, 0.6, 0.7, 0.8, 0.9}, "pt binings for proton ITS purity cut"};
  Configurable<std::string> cfgCCDBPurityPath{"cfgCCDBPurityPath", "Users/z/zhengqiw/PurityCut", "CCDB path for nsigmaITS - nSigmaTPC purity cut"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Service<ccdb::BasicCCDBManager> ccdb;
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
  // vectors for ITS-TPC contanmination cut
  std::vector<std::vector<std::shared_ptr<TH1>>> hPosPrCut;
  std::vector<std::vector<std::shared_ptr<TH1>>> hNegPrCut;

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
  Filter properPIDfilter = aod::cme_track_pid_columns::nPidFlag >= (int8_t)0;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags, aod::PidInfo>>> tracksSet1 = ((aod::cme_track_pid_columns::nPidFlag == 1) || (aod::cme_track_pid_columns::nPidFlag == 7) || (aod::cme_track_pid_columns::nPidFlag == 8) || (aod::cme_track_pid_columns::nPidFlag == 10));
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags, aod::PidInfo>>> tracksSet2 = ((aod::cme_track_pid_columns::nPidFlag == 2) || (aod::cme_track_pid_columns::nPidFlag == 7) || (aod::cme_track_pid_columns::nPidFlag == 9) || (aod::cme_track_pid_columns::nPidFlag == 10));
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags, aod::PidInfo>>> tracksSet3 = ((aod::cme_track_pid_columns::nPidFlag == 3) || (aod::cme_track_pid_columns::nPidFlag == 8) || (aod::cme_track_pid_columns::nPidFlag == 9) || (aod::cme_track_pid_columns::nPidFlag == 10));
  void init(InitContext const&)
  {

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

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

    AxisSpec axissumpt{cfgaxissumpt, "#it{p}_{T}^{sum}"};
    AxisSpec axisdeltaeta{cfgaxisdeltaeta, "#Delta#eta"};
    AxisSpec axisdeltapt{cfgaxisdeltapt, "#Delta#it{p}_{T}"};
    AxisSpec axisvertexz = {100, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec axisITSNcls = {10, -0.5, 9.5, "ITSNcls"};
    AxisSpec axisTPCNcls = {160, 0, 160, "TPCNcls"};

    histosQA.add(Form("QA/histEventCount"), "", {HistType::kTH1F, {{3, 0.0, 3.0}}});
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(3, "after additional event cut");
    if (cfgUseAdditionalEventCut) {
      histosQA.add(Form("QA/histEventCountDetail"), "Number of Event;; Count", {HistType::kTH1F, {{11, 0, 11}}});
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(1, "after sel8");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(8, "occupancy");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(9, "MultCorrelationPVTracks");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(10, "MultCorrelationGlobalTracks");
      histosQA.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(11, "cfgEvSelV0AT0ACut");
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

    if (cfgkOpenV2) {
      histosQA.add(Form("V2/histCosDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/histSinDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Pi"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Ka"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Pr"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Pi_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Ka_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
      histosQA.add(Form("V2/PID/histCosDetV2_Pr_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    }

    if (cfgkOpenTPCITSPurityCut && cfgkOpenTPCITSPurityCutQA) {
      histosQA.add(Form("QA/histITSPuritycheck_Pr_Pos_Cen_20_30"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
      histosQA.add(Form("QA/histITSPuritycheck_Pr_Neg_Cen_20_30"), ";n#sigma_{TPC};n#sigma_{ITS};#p_{t}", {HistType::kTH3F, {cfgnSigmaBinsTPC, cfgnSigmaBinsITS, cfgaxisptPID}});
    }
    if (cfgOpenPlotITSNcls) {
      histosQA.add(Form("QA/histITSNcls_PosPi"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add(Form("QA/histITSNcls_NegPi"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add(Form("QA/histITSNcls_PosKa"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add(Form("QA/histITSNcls_NegKa"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add(Form("QA/histITSNcls_PosPr"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add(Form("QA/histITSNcls_NegPr"), ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
    }
    if (cfgOpenPlotITSNclsPtCent) {
      histosQA.add(Form("QA/histITSNclsPtCent_PosPi"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histITSNclsPtCent_NegPi"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histITSNclsPtCent_PosKa"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histITSNclsPtCent_NegKa"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histITSNclsPtCent_PosPr"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histITSNclsPtCent_NegPr"), ";ITSNcls;Pt;Centrality", {HistType::kTH3F, {axisITSNcls, axisPt, axisCentMerged}});
    }
    if (cfgOpenTPCNcls) {
      histosQA.add(Form("QA/histTPCNcls_PosPi"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add(Form("QA/histTPCNcls_NegPi"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add(Form("QA/histTPCNcls_PosKa"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add(Form("QA/histTPCNcls_NegKa"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add(Form("QA/histTPCNcls_PosPr"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add(Form("QA/histTPCNcls_NegPr"), ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
    }
    if (cfgOpenTPCNclsPtCent) {
      histosQA.add(Form("QA/histTPCNclsPtCent_PosPi"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histTPCNclsPtCent_NegPi"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histTPCNclsPtCent_PosKa"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histTPCNclsPtCent_NegKa"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histTPCNclsPtCent_PosPr"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
      histosQA.add(Form("QA/histTPCNclsPtCent_NegPr"), ";TPCNcls;Pt;Centrality", {HistType::kTH3F, {axisTPCNcls, axisPt, axisCentMerged}});
    }

    if (cfgkOpenCME) {
      if (cfgkOpenPiPi) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_PiPi_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_PiPi_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_PiPi_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_PiPi_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_PiPi_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_PiPi_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPi_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPi_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPi_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenKaKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_KaKa_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_KaKa_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_KaKa_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_KaKa_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_KaKa_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_KaKa_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_KaKa_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPrPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_PrPr_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_PrPr_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_PrPr_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_PrPr_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_PrPr_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_PrPr_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PrPr_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PrPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PrPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPiKa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_PiKa_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_PiKa_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_PiKa_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_PiKa_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_PiKa_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_PiKa_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiKa_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenPiPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_PiPr_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_PiPr_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_PiPr_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_PiPr_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_PiPr_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_PiPr_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_PiPr_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_PiPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenKaPr) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_KaPr_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_KaPr_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_KaPr_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_KaPr_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_KaPr_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_KaPr_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_KaPr_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_KaPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
      if (cfgkOpenHaHa) {
        histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
        histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_os"), "", {HistType::kTProfile, {axisCentMerged}});
        if (cfgkOpenCMEDifferential) {
          if (cfgkOpenDeltaPt) {
            histosQA.add("PIDCME/Differential/histgamma_HaHa_ss_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histgamma_HaHa_os_DPt", ";centrality;#Delta #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_ss_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_os_DPt", ";centrality;#Delta #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltapt}});
          }
          if (cfgkOpenDeltaEta) {
            histosQA.add("PIDCME/Differential/histgamma_HaHa_ss_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histgamma_HaHa_os_DEt", ";centrality;#Delta #eta;#gamma", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_ss_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_os_DEt", ";centrality;#Delta #eta;#delta", {HistType::kTProfile2D, {axisCentMerged, axisdeltaeta}});
          }
          if (cfgkOpenAveragePt) {
            histosQA.add("PIDCME/Differential/histgamma_HaHa_ss_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histgamma_HaHa_os_SPt", ";centrality;Sum #p_{T};#gamma", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_ss_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
            histosQA.add("PIDCME/Differential/histdelta_HaHa_os_SPt", ";centrality;Sum #p_{T};#delta", {HistType::kTProfile2D, {axisCentMerged, axissumpt}});
          }
        }
        if (cfgkOpenSsOsCrossCheck) {
          histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histgamma_HaHa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_PP"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_NN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_PN"), "", {HistType::kTProfile, {axisCentMerged}});
          histosQA.add<TProfile>(Form("PIDCME/histdelta_HaHa_NP"), "", {HistType::kTProfile, {axisCentMerged}});
        }
      }
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    std::string fullPath;
    auto timestamp = bc.timestamp();
    // hPosPrCut.clear();
    // hNegPrCut.clear();
    for (auto i = 0; i < static_cast<int>(cfgITSPurityCen->size()) - 1; i++) {
      std::vector<std::shared_ptr<TH1>> hPosPrCutCen;
      std::vector<std::shared_ptr<TH1>> hNegPrCutCen;
      for (auto j = 0; j < static_cast<int>(cfgPtPrCut->size()) - 1; j++) {
        fullPath = cfgCCDBPurityPath;
        fullPath += Form("/ProtonPos/Cen_%d_%d/Pt_%d_%d", cfgITSPurityCen->at(i), cfgITSPurityCen->at(i + 1), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j))), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j + 1))));
        auto posPrHist = ccdb->getForTimeStamp<TH1F>(fullPath, timestamp);
        fullPath = cfgCCDBPurityPath;
        fullPath += Form("/ProtonNeg/Cen_%d_%d/Pt_%d_%d", cfgITSPurityCen->at(i), cfgITSPurityCen->at(i + 1), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j))), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j + 1))));
        auto negPrHist = ccdb->getForTimeStamp<TH1F>(fullPath, timestamp);
        if (!posPrHist) {
          LOGF(fatal, Form("could not load Pos Proton ITS TPC purity hist for Cent_%d_%d Pt_%d_%d(MeV)", cfgITSPurityCen->at(i), cfgITSPurityCen->at(i + 1), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j))), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j + 1)))));
        }
        if (!negPrHist) {
          LOGF(fatal, Form("could not load Neg Proton ITS TPC purity hist for Cent_%d_%d Pt_%d_%d(MeV)", cfgITSPurityCen->at(i), cfgITSPurityCen->at(i + 1), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j))), static_cast<int>(std::round(1e3 * cfgPtPrCut->at(j + 1)))));
        }
        std::shared_ptr<TH1> sharedPosPrHist(posPrHist);
        std::shared_ptr<TH1> sharedNegPrHist(negPrHist);
        hPosPrCutCen.push_back(std::move(sharedPosPrHist));
        hNegPrCutCen.push_back(std::move(sharedNegPrHist));
      }
      hPosPrCut.push_back(std::move(hPosPrCutCen));
      hNegPrCut.push_back(std::move(hNegPrCutCen));
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
    if (cfgOpenEvSelMultCorrelationPVTracks) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return 0;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return 0;
    }
    if (cfgOpenEvSelMultCorrelationPVTracks) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 8.5);
    }
    if (cfgOpenEvSelMultCorrelationGlobalTracks) {
      if (multTrk < fMultCutLow->Eval(centrality))
        return 0;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return 0;
    }
    if (cfgOpenEvSelMultCorrelationGlobalTracks) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 9.5);
    }
    if (cfgOpenEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > 5 * fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return 0;
    }
    if (cfgOpenEvSelV0AT0ACut) {
      histosQA.fill(HIST("QA/histEventCountDetail"), 10.5);
    }
    return 1;
  }

  template <typename TrackType>
  bool selTrack(const TrackType track, float centrality)
  {
    if (cfgkOpenDebugPIDCME) {
      LOGF(info, "====================Entering track selection=============================");
    }
    if (cfgOpenCustomTrackCutAssurance) {
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
    }
    if (cfgkOpenTPCITSPurityCut) {
      int cenBin = -1;
      int ptBin = -1;
      for (int i = 0; i < static_cast<int>(cfgITSPurityCen.value.size()) - 1; ++i) {
        if (centrality >= cfgITSPurityCen.value[i] && centrality < cfgITSPurityCen.value[i + 1]) {
          cenBin = i;
          break;
        }
      }
      for (int i = 0; i < static_cast<int>(cfgPtPrCut.value.size()) - 1; ++i) {
        if (track.pt() >= cfgPtPrCut.value[i] && track.pt() < cfgPtPrCut.value[i + 1]) {
          ptBin = i;
          break;
        }
      }
      if ((cenBin >= 0) && (ptBin >= 0)) {
        if ((track.nPidFlag() == 3) || (track.nPidFlag() == 8) || (track.nPidFlag() == 9) || (track.nPidFlag() == 10)) {
          if (cfgkOpenDebugPIDCME) {
            LOGF(info, Form("=========cen_bin: %d pt_bin: %d=========", cenBin, ptBin));
          }
          float nSigmaITSPr = track.nSigmaPrITS();
          int xBin = hPosPrCut[cenBin][ptBin]->GetXaxis()->FindBin(nSigmaITSPr);
          float binContentPosPr = hPosPrCut[cenBin][ptBin]->GetBinContent(xBin);
          float binContentNegPr = hNegPrCut[cenBin][ptBin]->GetBinContent(xBin);
          if (track.sign() > 0) {
            if ((binContentPosPr != 0) && (track.nSigmaPrTPC() < binContentPosPr)) {
              if (cfgkOpenDebugPIDCME) {
                LOGF(info, "====================Track selection Finished with cut=============================");
              }
              return false;
            }
          } else {
            if ((binContentNegPr != 0) && (track.nSigmaPrTPC() < binContentNegPr)) {
              if (cfgkOpenDebugPIDCME) {
                LOGF(info, "====================Track selection Finished with cut=============================");
              }
              return false;
            }
          }
        }
      }
    }
    if (cfgkOpenDebugPIDCME) {
      LOGF(info, "====================Track selection Finished without cut=============================");
    }
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
    auto cent = collision.centFT0C();
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    float psiN = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode);
    if (cfgkOpenV2) {
      for (const auto& trk : track1) {
        if (!selTrack(trk, cent))
          continue;
        if (nmode == 2) {
          if (trk.sign() > 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Pi"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_PosPi"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_PosPi"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_PosPi"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_PosPi"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          } else if (trk.sign() < 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Pi_Neg"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_NegPi"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_NegPi"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_NegPi"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_NegPi"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          }
        }
      }
      for (const auto& trk : track2) {
        if (!selTrack(trk, cent))
          continue;
        if (nmode == 2) {
          if (trk.sign() > 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Ka"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_PosKa"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_PosKa"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_PosKa"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_PosKa"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          } else if (trk.sign() < 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Ka_Neg"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_NegKa"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_NegKa"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_NegKa"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_NegKa"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          }
        }
      }
      for (const auto& trk : track3) {
        if (!selTrack(trk, cent))
          continue;
        if (nmode == 2) {
          if (trk.sign() > 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Pr"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_PosPr"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_PosPr"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_PosPr"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_PosPr"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          } else if (trk.sign() < 0) {
            histosQA.fill(HIST("V2/PID/histCosDetV2_Pr_Neg"), cent, trk.pt(),
                          std::cos(static_cast<float>(nmode) * (trk.phi() - psiN)));
            if (cfgOpenPlotITSNcls) {
              histosQA.fill(HIST("QA/histITSNcls_NegPr"), trk.itsNCls());
            }
            if (cfgOpenPlotITSNclsPtCent) {
              histosQA.fill(HIST("QA/histITSNclsPtCent_NegPr"), trk.itsNCls(), trk.pt(), cent);
            }
            if (cfgOpenTPCNcls) {
              histosQA.fill(HIST("QA/histTPCNcls_NegPr"), trk.tpcNClsFound());
            }
            if (cfgOpenTPCNclsPtCent) {
              histosQA.fill(HIST("QA/histTPCNclsPtCent_NegPr"), trk.tpcNClsFound(), trk.pt(), cent);
            }
          }
        }
      }
    }
    if (cfgkOpenCME) {
      if (cfgkOpenPiPi) {
        for (const auto& trk1 : track1) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk2 : track1) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (!selTrack(trk2, cent))
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_PiPi_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPi_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPi_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPi_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_PiPi_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPi_os"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPi_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPi_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPi_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPi_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPi_NP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenKaKa) {
        for (const auto& trk1 : track2) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk2 : track2) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (!selTrack(trk2, cent))
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_KaKa_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaKa_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_KaKa_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_KaKa_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_KaKa_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaKa_os"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaKa_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaKa_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_KaKa_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_KaKa_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaKa_NP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPrPr) {
        for (const auto& trk1 : track3) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk2 : track3) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (!selTrack(trk2, cent))
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_PrPr_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PrPr_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PrPr_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PrPr_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_PrPr_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PrPr_os"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PrPr_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PrPr_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PrPr_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PrPr_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PrPr_NP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPiKa) {
        for (const auto& trk1 : track1) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk2 : track2) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (!selTrack(trk2, cent))
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk2.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_PiKa_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiKa_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_ss_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_ss_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_ss_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiKa_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiKa_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_PiKa_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiKa_os"), cent, std::cos((trk1.phi() - trk2.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_DPt"), cent, trk1.pt() - trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_DEt"), cent, trk1.eta() - trk2.eta(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiKa_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_SPt"), cent, trk1.pt() + trk2.pt(),
                                  std::cos((trk1.phi() - trk2.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk2.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiKa_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiKa_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiKa_NP"), cent, std::cos((trk1.phi() - trk2.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenPiPr) {
        for (const auto& trk1 : track1) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk3 : track3) {
            if (trk1.globalIndex() == trk3.globalIndex())
              continue;
            if (!selTrack(trk3, cent))
              continue;
            if (nmode == 2) {
              if (trk1.sign() == trk3.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_PiPr_ss"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPr_ss"), cent, std::cos((trk1.phi() - trk3.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_ss_DPt"), cent, trk1.pt() - trk3.pt(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_ss_DPt"), cent, trk1.pt() - trk3.pt(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_ss_DEt"), cent, trk1.eta() - trk3.eta(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_ss_DEt"), cent, trk1.eta() - trk3.eta(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_ss_SPt"), cent, trk1.pt() + trk3.pt(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_ss_SPt"), cent, trk1.pt() + trk3.pt(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk3.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPr_PP"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_PP"), cent, std::cos((trk1.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPr_NN"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_NN"), cent, std::cos((trk1.phi() - trk3.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_PiPr_os"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_PiPr_os"), cent, std::cos((trk1.phi() - trk3.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_os_DPt"), cent, trk1.pt() - trk3.pt(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_DPt"), cent, trk1.pt() - trk3.pt(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_os_DEt"), cent, trk1.eta() - trk3.eta(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_DEt"), cent, trk1.eta() - trk3.eta(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_PiPr_os_SPt"), cent, trk1.pt() + trk3.pt(),
                                  std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_SPt"), cent, trk1.pt() + trk3.pt(),
                                  std::cos((trk1.phi() - trk3.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk1.sign() > 0 && trk3.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPr_PN"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_PN"), cent, std::cos((trk1.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_PiPr_NP"), cent, std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_PiPr_NP"), cent, std::cos((trk1.phi() - trk3.phi())));
                  }
                }
              }
            }
          }
        }
      }
      if (cfgkOpenKaPr) {
        for (const auto& trk2 : track2) {
          if (!selTrack(trk2, cent))
            continue;
          for (const auto& trk3 : track3) {
            if (trk2.globalIndex() == trk3.globalIndex())
              continue;
            if (!selTrack(trk3, cent))
              continue;
            if (nmode == 2) {
              if (trk2.sign() == trk3.sign()) {
                histosQA.fill(HIST("PIDCME/histgamma_KaPr_ss"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaPr_ss"), cent, std::cos((trk2.phi() - trk3.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_ss_DPt"), cent, trk2.pt() - trk3.pt(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_ss_DPt"), cent, trk2.pt() - trk3.pt(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_ss_DEt"), cent, trk2.eta() - trk3.eta(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_ss_DEt"), cent, trk2.eta() - trk3.eta(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_ss_SPt"), cent, trk2.pt() + trk3.pt(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_ss_SPt"), cent, trk2.pt() + trk3.pt(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk2.sign() > 0 && trk3.sign() > 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_KaPr_PP"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_PP"), cent, std::cos((trk2.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_KaPr_NN"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_NN"), cent, std::cos((trk2.phi() - trk3.phi())));
                  }
                }
              } else {
                histosQA.fill(HIST("PIDCME/histgamma_KaPr_os"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                histosQA.fill(HIST("PIDCME/histdelta_KaPr_os"), cent, std::cos((trk2.phi() - trk3.phi())));
                if (cfgkOpenCMEDifferential) {
                  if (cfgkOpenDeltaPt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_os_DPt"), cent, trk2.pt() - trk3.pt(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_DPt"), cent, trk2.pt() - trk3.pt(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                  if (cfgkOpenDeltaEta) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_os_DEt"), cent, trk2.eta() - trk3.eta(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_DEt"), cent, trk2.eta() - trk3.eta(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                  if (cfgkOpenAveragePt) {
                    histosQA.fill(HIST("PIDCME/Differential/histgamma_KaPr_os_SPt"), cent, trk2.pt() + trk3.pt(),
                                  std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_SPt"), cent, trk2.pt() + trk3.pt(),
                                  std::cos((trk2.phi() - trk3.phi())));
                  }
                }
                if (cfgkOpenSsOsCrossCheck) {
                  if (trk2.sign() > 0 && trk3.sign() < 0) {
                    histosQA.fill(HIST("PIDCME/histgamma_KaPr_PN"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_PN"), cent, std::cos((trk2.phi() - trk3.phi())));
                  } else {
                    histosQA.fill(HIST("PIDCME/histgamma_KaPr_NP"), cent, std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * psiN)));
                    histosQA.fill(HIST("PIDCME/histdelta_KaPr_NP"), cent, std::cos((trk2.phi() - trk3.phi())));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::Qvectors>>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags, aod::PidInfo>> const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if ((hPosPrCut.empty()) || (hNegPrCut.empty())) {
      initCCDB(bc);
      LOGF(info, "==================CCDB file successfuly applied====================");
      LOGF(info, Form("size of hPosPrCut is %lu x %lu", hPosPrCut.size(), hPosPrCut[0].size()));
      LOGF(info, Form("size of hNegPrCut is %lu x %lu", hNegPrCut.size(), hNegPrCut[0].size()));
      LOGF(info, "===================================================================");
    }
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
    if (cfgkOpenDebugPIDCME) {
      LOGF(info, "==================Event Cut Finished====================");
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
        if (!selTrack(trk, cent))
          continue;
        if (cfgkOpenTPCITSPurityCut && cfgkOpenTPCITSPurityCutQA) {
          if (cent >= 20 && cent < 30) {
            if ((trk.nPidFlag() == 3) || (trk.nPidFlag() == 8) || (trk.nPidFlag() == 9) || (trk.nPidFlag() == 10)) {
              if (trk.sign() > 0) {
                histosQA.fill(HIST("QA/histITSPuritycheck_Pr_Pos_Cen_20_30"), trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
              } else {
                histosQA.fill(HIST("QA/histITSPuritycheck_Pr_Neg_Cen_20_30"), trk.nSigmaPrTPC(), trk.nSigmaPrITS(), trk.pt());
              }
            }
          }
        }
        if (cfgkOpenV2) {
          histosQA.fill(HIST("V2/histSinDetV2"), cent, trk.pt(),
                        std::sin(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
          histosQA.fill(HIST("V2/histCosDetV2"), cent, trk.pt(),
                        std::cos(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - psiNGlobal)));
        }
      }
      if (cfgkOpenCME && cfgkOpenHaHa && cfgnMods->at(i) == 2) {
        for (const auto& trk1 : tracks) {
          if (!selTrack(trk1, cent))
            continue;
          for (const auto& trk2 : tracks) {
            if (trk1.globalIndex() == trk2.globalIndex())
              continue;
            if (!selTrack(trk2, cent))
              continue;
            if (trk1.sign() == trk2.sign()) {
              histosQA.fill(HIST("PIDCME/histgamma_HaHa_ss"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_ss"), cent, std::cos((trk1.phi() - trk2.phi())));
              if (cfgkOpenCMEDifferential) {
                if (cfgkOpenDeltaPt) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_ss_DPt"), cent, trk2.pt() - trk2.pt(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_ss_DPt"), cent, trk2.pt() - trk2.pt(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
                if (cfgkOpenDeltaEta) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_ss_DEt"), cent, trk2.eta() - trk2.eta(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_ss_DEt"), cent, trk2.eta() - trk2.eta(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
                if (cfgkOpenAveragePt) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_ss_SPt"), cent, trk2.pt() + trk2.pt(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_ss_SPt"), cent, trk2.pt() + trk2.pt(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
              }
              if (cfgkOpenSsOsCrossCheck) {
                if (trk1.sign() > 0 && trk2.sign() > 0) {
                  histosQA.fill(HIST("PIDCME/histgamma_HaHa_PP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_PP"), cent, std::cos((trk1.phi() - trk2.phi())));
                } else {
                  histosQA.fill(HIST("PIDCME/histgamma_HaHa_NN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_NN"), cent, std::cos((trk1.phi() - trk2.phi())));
                }
              }
            } else {
              histosQA.fill(HIST("PIDCME/histgamma_HaHa_os"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
              histosQA.fill(HIST("PIDCME/histdelta_HaHa_os"), cent, std::cos((trk1.phi() - trk2.phi())));
              if (cfgkOpenCMEDifferential) {
                if (cfgkOpenDeltaPt) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_os_DPt"), cent, trk2.pt() - trk2.pt(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_os_DPt"), cent, trk2.pt() - trk2.pt(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
                if (cfgkOpenDeltaEta) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_os_DEt"), cent, trk2.eta() - trk2.eta(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_os_DEt"), cent, trk2.eta() - trk2.eta(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
                if (cfgkOpenAveragePt) {
                  histosQA.fill(HIST("PIDCME/Differential/histgamma_HaHa_os_SPt"), cent, trk2.pt() + trk2.pt(),
                                std::cos((trk2.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/Differential/histdelta_HaHa_os_SPt"), cent, trk2.pt() + trk2.pt(),
                                std::cos((trk2.phi() - trk2.phi())));
                }
              }
              if (cfgkOpenSsOsCrossCheck) {
                if (trk1.sign() > 0 && trk2.sign() < 0) {
                  histosQA.fill(HIST("PIDCME/histgamma_HaHa_PN"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
                  histosQA.fill(HIST("PIDCME/histdelta_HaHa_PN"), cent, std::cos((trk1.phi() - trk2.phi())));
                } else {
                  histosQA.fill(HIST("PIDCME/histgamma_HaHa_NP"), cent, std::cos((trk1.phi() + trk2.phi() - static_cast<float>(cfgnMods->at(i)) * psiNGlobal)));
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
