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
/// \file   flowEsePHe3.cxx
/// \brief  task to calculate the P He3 flow correlation.
// C++/ROOT includes.
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace flow_ese_p_he3
{
DECLARE_SOA_COLUMN(NPidFlag, nPidFlag, int8_t); // unqualified -1, hadron 0, proton 1, he3 2, proton+he3 3
} // namespace flow_ese_p_he3
DECLARE_SOA_TABLE(PHe3ESEFlags, "AOD", "PHe3ESEFlags", flow_ese_p_he3::NPidFlag);
} // namespace o2::aod

namespace pid_flags
{
// constexpr int8_t kUnqualified = -1;
// constexpr int8_t kUnPOIHadron = 0;
constexpr int8_t kProton = 1;
constexpr int8_t kHe3 = 2;
constexpr int8_t kProtonHe3 = 3;
} // namespace pid_flags

namespace event_selection
{
constexpr int kFT0AV0ASigma = 5;
}

namespace fourier_mode
{
// constexpr int kMode1 = 1;
constexpr int kMode2 = 2;
// constexpr int kMode3 = 3;
} // namespace fourier_mode

using TracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullHe, aod::pidTOFFullHe>;
struct FillPIDcolums {

  HistogramRegistry histosQA{"histosQAPID", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxPtPID{"cfgMaxPtPID", 99.9, "Maximum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.8, "Maximum track #eta for PID"};
  Configurable<float> cfgMinTPCChi2NCl{"cfgMinTPCChi2NCl", 0, "Minimum chi2 per cluster TPC for PID if not use costom track cuts"};
  Configurable<float> cfgMinChi2NClITS{"cfgMinChi2NClITS", 0, "Minimum chi2 per cluster ITS for PID if not use costom track cuts"};
  Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4, "Maximum chi2 per cluster TPC for PID if not use costom track cuts"};
  Configurable<float> cfgMaxChi2NClITS{"cfgMaxChi2NClITS", 36, "Maximum chi2 per cluster ITS for PID if not use costom track cuts"};
  Configurable<float> cfgMinTPCCls{"cfgMinTPCCls", 70, "Minimum TPC clusters for PID if not use costom track cuts"};
  Configurable<float> cfgMinITSCls{"cfgMinITSCls", 1, "Minimum ITS clusters for PID if not use costom track cuts"};
  Configurable<float> cfgMaxTPCCls{"cfgMaxTPCCls", 999, "Max TPC clusters for PID if not use costom track cuts"};
  Configurable<float> cfgMaxITSCls{"cfgMaxITSCls", 999, "Max ITS clusters for PID if not use costom track cuts"};
  Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 99, "Maxium DCAxy for standard PID tracking"};
  Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 2, "Maxium DCAz for standard PID tracking"};
  Configurable<float> cfgPtMaxforTPCOnlyPIDPrton{"cfgPtMaxforTPCOnlyPIDPrton", 0.4, "Maxmium track pt for TPC only PID, at RMS PID mode for proton"};
  Configurable<float> cfgPtMaxforTPCOnlyPIDHe3{"cfgPtMaxforTPCOnlyPIDHe3", 0.5, "Maxmium track pt for TPC only PID, at RMS PID mode for he3"};

  Configurable<int> cfgProtonPIDMode{"cfgProtonPIDMode", 2, "Proton PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};
  Configurable<int> cfgHe3PIDMode{"cfgHe3PIDMode", 1, "He3 PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};

  Configurable<bool> cfgOpenpassedITSNCls{"cfgOpenpassedITSNCls", false, "useTrackSelectionTables passedITSNCls for basic track selection"};
  Configurable<bool> cfgOpenpassedITSChi2NDF{"cfgOpenpassedITSChi2NDF", false, "useTrackSelectionTables passedITSChi2NDF for basic track selection"};
  Configurable<bool> cfgOpenpassedITSHits{"cfgOpenpassedITSHits", false, "useTrackSelectionTables passedITSHits for basic track selection"};
  Configurable<bool> cfgOpenpassedTPCChi2NDF{"cfgOpenpassedTPCChi2NDF", false, "useTrackSelectionTables passedTPCChi2NDF for basic track selection"};
  Configurable<bool> cfgOpenpassedTPCCrossedRowsOverNCls{"cfgOpenpassedTPCCrossedRowsOverNCls", false, "useTrackSelectionTables passedTPCCrossedRowsOverNCls for basic track selection"};
  Configurable<bool> cfgOpenpassedDCAxy{"cfgOpenpassedDCAxy", false, "useTrackSelectionTables passedDCAxy for basic track selection"};
  Configurable<bool> cfgOpenpassedDCAz{"cfgOpenpassedDCAz", false, "useTrackSelectionTables passedDCAz for basic track selection"};

  Configurable<bool> cfgQuietMode{"cfgQuietMode", false, "open quiet mode for saving cpu cost and only do some basic QA plots"};
  Configurable<bool> cfgOpenPIDITSProton{"cfgOpenPIDITSProton", true, "open ITS assistance cut for proton PID"};
  Configurable<bool> cfgOpenPIDITSHe3{"cfgOpenPIDITSHe3", false, "open ITS assistance cut for He3 PID"};
  Configurable<bool> cfgOpenPIDByPtProtonMain{"cfgOpenPIDByPtProtonMain", false, "Selection Proton by pt its pt binnings for main selection"};
  Configurable<bool> cfgOpenPIDByPtHe3Main{"cfgOpenPIDByPtHe3Main", false, "Selection He3 by pt its pt binnings for main selection"};
  Configurable<bool> cfgOpenPIDByPtProtonITS{"cfgOpenPIDByPtProtonITS", false, "Selection Proton by pt its pt binnings for ITS selection"};
  Configurable<bool> cfgOpenPIDByPtHe3ITS{"cfgOpenPIDByPtHe3ITS", false, "Selection He3 by pt its pt binnings for ITS selection"};
  Configurable<bool> cfgOpenHe3ITSPtCut{"cfgOpenHe3ITSPtCut", true, "Do He3 ITS contamination cut"};
  Configurable<bool> cfgOpenAllowCrossTrack{"cfgOpenAllowCrossTrack", false, "Allow one track to be identified as different kind of PID particles"};

  Configurable<bool> cfgOpenPlotnSigmaTOFITSPt{"cfgOpenPlotnSigmaTOFITSPt", true, "plot nSigmaTOF vs nSigmaITS vs Pt"};
  Configurable<bool> cfgOpenPlotnSigmaITSTPCPt{"cfgOpenPlotnSigmaITSTPCPt", true, "plot nSigmaITS vs nSigmaTOF vs Pt"};
  Configurable<bool> cfgOpenPlotnSigmaTOFTPCPt{"cfgOpenPlotnSigmaTOFTPCPt", true, "plot nSigmaTOF vs nSigmaTPC vs Pt"};

  Configurable<std::vector<float>> cfgPtCutProton{"cfgPtCutProton", {0.15, 99.}, "Pt limit for Proton"};
  Configurable<std::vector<float>> cfgPtCutHe3{"cfgPtCutHe3", {0.15, 99.}, "Pt limit for He3"};
  Configurable<std::vector<float>> cfgnSigmaCutTPCProton{"cfgnSigmaCutTPCProton", {-3, 3}, "TPC nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutTPCHe3{"cfgnSigmaCutTPCHe3", {-2, 2}, "TPC nsigma cut limit for He3"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFProton{"cfgnSigmaCutTOFProton", {-1.5, 1.5}, "TOF nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFHe3{"cfgnSigmaCutTOFHe3", {-1.5, 1.5}, "TOF nsigma cut limit for He3"};
  Configurable<std::vector<float>> cfgnSigmaCutITSProton{"cfgnSigmaCutITSProton", {-3, 3}, "ITS nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutITSHe3{"cfgnSigmaCutITSHe3", {-3, 3}, "ITS nsigma cut limit for He3"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSProton{"cfgnSigmaCutRMSProton", {-3, 3}, "RMS nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSHe3{"cfgnSigmaCutRMSHe3", {-3, 3}, "RMS nsigma cut limit for He3"};

  Configurable<std::vector<float>> cfgPtBinProtonPID{"cfgPtBinProtonPID", {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0}, "pt bin for pion PIDnsigma"};
  Configurable<std::vector<float>> cfgPtBinHe3PID{"cfgPtBinHe3PID", {2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.6, 4, 4.4, 4.8, 5.2, 5.6, 6, 6.4, 7.2, 8, 10}, "pt bin for pion PIDnsigma"};

  Configurable<std::vector<float>> cfgnSigmaTPCProtonPtUpper{"cfgnSigmaTPCProtonPtUpper", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaTPC cut upper limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFProtonPtUpper{"cfgnSigmaTOFProtonPtUpper", {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}, "nSigmaTOF cut upper limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSProtonPtUpper{"cfgnSigmaITSProtonPtUpper", {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, "nSigmaITS cut upper limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaRMSProtonPtUpper{"cfgnSigmaRMSProtonPtUpper", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaRMS cut upper limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTPCProtonPtLower{"cfgnSigmaTPCProtonPtLower", {-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3}, "nSigmaTPC cut lower limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFProtonPtLower{"cfgnSigmaTOFProtonPtLower", {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5}, "nSigmaTOF cut lower limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSProtonPtLower{"cfgnSigmaITSProtonPtLower", {-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2}, "nSigmaITS cut lower limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaRMSProtonPtLower{"cfgnSigmaRMSProtonPtLower", {-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3}, "nSigmaRMS cut lower limit anchored to proton pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTPCHe3PtUpper{"cfgnSigmaTPCHe3PtUpper", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaTPC cut upper limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFHe3PtUpper{"cfgnSigmaTOFHe3PtUpper", {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}, "nSigmaTOF cut upper limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSHe3PtUpper{"cfgnSigmaITSHe3PtUpper", {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, "nSigmaITS cut upper limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaRMSHe3PtUpper{"cfgnSigmaRMSHe3PtUpper", {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}, "nSigmaRMS cut upper limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTPCHe3PtLower{"cfgnSigmaTPCHe3PtLower", {-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3}, "nSigmaTPC cut lower limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaTOFHe3PtLower{"cfgnSigmaTOFHe3PtLower", {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5}, "nSigmaTOF cut lower limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaITSHe3PtLower{"cfgnSigmaITSHe3PtLower", {-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2}, "nSigmaITS cut lower limit anchored to He3 pt bins"};
  Configurable<std::vector<float>> cfgnSigmaRMSHe3PtLower{"cfgnSigmaRMSHe3PtLower", {-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3}, "nSigmaRMS cut lower limit anchored to He3 pt bins"};

  ConfigurableAxis cfgrigidityBins{"cfgrigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis cfgdedxBins{"cfgdedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsTOF{"cfgnSigmaBinsTOF", {200, -5.f, 5.f}, "Binning for n sigma TOF"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma ITS"};
  ConfigurableAxis cfgnSigmaBinsRMS{"cfgnSigmaBinsRMS", {100, 0.f, 10.f}, "Combination Binning for TPC&TOF nsigma"};
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgaxispPID{"cfgaxispPID", {50, 0, 5}, "Binning for P PID"};
  ConfigurableAxis cfgaxisetaPID{"cfgaxisetaPID", {90, -0.9, 0.9}, "Binning for Pt QA"};
  ConfigurableAxis cfgaxisDCAz{"cfgaxisDCAz", {200, -1, 1}, "Binning for DCAz"};
  ConfigurableAxis cfgaxisDCAxy{"cfgaxisDCAxy", {100, -0.5, 0.5}, "Binning for DCAxy"};
  ConfigurableAxis cfgaxisChi2Ncls{"cfgaxisChi2Ncls", {100, 0, 100}, "Binning for Chi2Ncls TPC/ITS"};

  // Function for He3 TPC-ITS mismatching cuts referd to by chiara's slides
  TF1* fNSigmaITSPt = nullptr;

  template <typename TrackType>
  bool trackSelBasic(const TrackType track)
  {
    if ((track.pt() < cfgMinPtPID) || (track.pt() > cfgMaxPtPID))
      return false;
    if (std::abs(track.eta()) > cfgMaxEtaPID)
      return false;
    if (cfgOpenpassedITSNCls) {
      if (!track.passedITSNCls())
        return false;
    } else {
      if (track.itsNCls() < cfgMinITSCls || track.itsNCls() > cfgMaxITSCls)
        return false;
    }
    if (cfgOpenpassedITSChi2NDF) {
      if (!track.passedITSChi2NDF())
        return false;
    } else {
      if (track.itsChi2NCl() < cfgMinChi2NClITS || track.itsChi2NCl() > cfgMaxChi2NClITS)
        return false;
    }
    if (cfgOpenpassedITSHits) {
      if (!track.passedITSHits())
        return false;
    }
    if (cfgOpenpassedTPCChi2NDF) {
      if (!track.passedTPCChi2NDF())
        return false;
    } else {
      if (track.tpcChi2NCl() < cfgMinTPCChi2NCl || track.tpcChi2NCl() > cfgMaxTPCChi2NCl)
        return false;
    }
    if (cfgOpenpassedTPCCrossedRowsOverNCls) {
      if (!track.passedTPCCrossedRowsOverNCls())
        return false;
    }
    if (cfgOpenpassedDCAxy) {
      if (!track.passedDCAxy())
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cfgMaxDCAxy)
        return false;
    }
    if (cfgOpenpassedDCAz) {
      if (!track.passedDCAz())
        return false;
    } else {
      if (std::abs(track.dcaZ()) > cfgMaxDCAz)
        return false;
    }
    if (track.tpcNClsFound() < cfgMinTPCCls || track.tpcNClsFound() > cfgMaxTPCCls)
      return false;
    return true;
  }

  template <typename TrackType>
  bool pidProtonSel(const TrackType track, float nSigmaLower, float nSigmaUpper)
  { // proton == 1 , He3 == 2
    if (track.pt() < cfgPtCutProton.value[0] || track.pt() > cfgPtCutProton.value[1])
      return false;
    float nSigmaUse = -999;
    switch (cfgProtonPIDMode) {
      case 0: // RMS
        nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDPrton) ? std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()) : track.tpcNSigmaPr();
        break;
      case 1: // TPC only
        nSigmaUse = track.tpcNSigmaPr();
        break;
      case 2: // TOF only
        nSigmaUse = track.tofNSigmaPr();
        break;
    }
    if (nSigmaUse < nSigmaLower || nSigmaUse > nSigmaUpper) {
      return false;
    } else {
      return true;
    }
  }

  template <typename TrackType>
  bool pidHe3Sel(const TrackType track, float nSigmaLower, float nSigmaUpper)
  { // proton == 1 , He3 == 2
    if (track.pt() < cfgPtCutHe3.value[0] || track.pt() > cfgPtCutHe3.value[1])
      return false;
    float nSigmaUse = -999;
    switch (cfgHe3PIDMode) {
      case 0: // RMS
        nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDHe3) ? std::hypot(track.tpcNSigmaHe(), track.tofNSigmaHe()) : track.tpcNSigmaHe();
        break;
      case 1: // TPC only
        nSigmaUse = track.tpcNSigmaHe();
        break;
      case 2: // TOF only
        nSigmaUse = track.tofNSigmaHe();
        break;
    }
    if (nSigmaUse < nSigmaLower || nSigmaUse > nSigmaUpper) {
      return false;
    } else {
      return true;
    }
  }

  template <typename TrackType>
  int crossTrackID(const TrackType track)
  {
    if (track.tpcNSigmaPr() < track.tpcNSigmaHe()) {
      return 0;
    } else {
      return 1;
    }
  }

  void init(InitContext const&)
  {
    if (cfgOpenHe3ITSPtCut) {
      fNSigmaITSPt = new TF1("fNSigmaITSPt", "[0]/pow(x,0.5) - [2]", 0.02, 1000);
      fNSigmaITSPt->SetParameters(4.6, 0.5, 4.5);
    }
    AxisSpec axisITSNcls = {10, -1.5, 8.5, "ITSNcls"};
    AxisSpec axisTPCNcls = {160, 0, 160, "TPCNcls"};
    if (!cfgQuietMode) {
      histosQA.add("QA/hist_dEdxTPC_All", ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", {HistType::kTH2F, {cfgrigidityBins, cfgdedxBins}});
      histosQA.add("QA/Proton/hist_dEdxTPC_Pr", ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", {HistType::kTH2F, {cfgrigidityBins, cfgdedxBins}});
      histosQA.add("QA/He3/hist_dEdxTPC_He3", ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", {HistType::kTH2F, {cfgrigidityBins, cfgdedxBins}});
      histosQA.add("QA/hist_pt_All", ";#it{p}_{T};counts", {HistType::kTH1F, {cfgaxisptPID}});
      histosQA.add("QA/Proton/hist_pt_Pr", ";#it{p}_{T};counts", {HistType::kTH1F, {cfgaxisptPID}});
      histosQA.add("QA/He3/hist_pt_He3", ";#it{p}_{T};counts", {HistType::kTH1F, {cfgaxisptPID}});
      histosQA.add("QA/hist_eta_All", ";#it{#eta};counts", {HistType::kTH1F, {cfgaxisetaPID}});
      histosQA.add("QA/Proton/hist_eta_Pr", ";#it{#eta};counts", {HistType::kTH1F, {cfgaxisetaPID}});
      histosQA.add("QA/He3/hist_eta_He3", ";#it{#eta};counts", {HistType::kTH1F, {cfgaxisetaPID}});
      histosQA.add("QA/hist_ITSNcls_All", ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add("QA/Proton/hist_ITSNcls_Pr", ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add("QA/He3/hist_ITSNcls_He3", ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
      histosQA.add("QA/hist_TPCNcls_All", ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add("QA/Proton/hist_TPCNcls_Pr", ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add("QA/He3/hist_TPCNcls_He3", ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
      histosQA.add("QA/hist_ITSChi2NDF_All", ";ITS#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/Proton/hist_ITSChi2NDF_Pr", ";ITS#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/He3/hist_ITSChi2NDF_He3", ";ITS#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/hist_TPCChi2NDF_All", ";TPC#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/Proton/hist_TPCChi2NDF_Pr", ";TPC#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/He3/hist_TPCChi2NDF_He3", ";TPC#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
      histosQA.add("QA/hist_DCAxy_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAxy}});
      histosQA.add("QA/Proton/hist_DCAxy_Pr", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAxy}});
      histosQA.add("QA/He3/hist_DCAxy_He3", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAxy}});
      histosQA.add("QA/hist_DCAz_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAz}});
      histosQA.add("QA/Proton/hist_DCAz_Pr", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAz}});
      histosQA.add("QA/He3/hist_DCAz_He3", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAz}});
      histosQA.add("QA/Proton/hist_nSigmaTPC_Pr", ";n_{#sigma}TPC", {HistType::kTH1F, {cfgnSigmaBinsTPC}});
      histosQA.add("QA/Proton/hist_nSigmaTPCPt_Pr", ";#it{p}_{T};n_{#sigma}TPC", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsTPC}});
      histosQA.add("QA/Proton/hist_nSigmaTOF_Pr", ";n_{#sigma}TOF", {HistType::kTH1F, {cfgnSigmaBinsTOF}});
      histosQA.add("QA/Proton/hist_nSigmaTOFPt_Pr", ";#it{p}_{T};n_{#sigma}TOF", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsTOF}});
      histosQA.add("QA/Proton/hist_nSigmaITS_Pr", ";n_{#sigma}ITS", {HistType::kTH1F, {cfgnSigmaBinsITS}});
      histosQA.add("QA/Proton/hist_nSigmaITSPt_Pr", ";#it{p}_{T};n_{#sigma}ITS", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsITS}});
      histosQA.add("QA/Proton/hist_nSigmaRMS_Pr", ";n_{#sigma}RMS", {HistType::kTH1F, {cfgnSigmaBinsRMS}});
      histosQA.add("QA/Proton/hist_nSigmaRMSPt_Pr", ";#it{p}_{T};n_{#sigma}RMS", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsRMS}});
      histosQA.add("QA/He3/hist_nSigmaTPC_He3", ";n_{#sigma}TPC", {HistType::kTH1F, {cfgnSigmaBinsTPC}});
      histosQA.add("QA/He3/hist_nSigmaTPCPt_He3", ";#it{p}_{T};n_{#sigma}TPC", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsTPC}});
      histosQA.add("QA/He3/hist_nSigmaTOF_He3", ";n_{#sigma}TOF", {HistType::kTH1F, {cfgnSigmaBinsTOF}});
      histosQA.add("QA/He3/hist_nSigmaTOFPt_He3", ";#it{p}_{T};n_{#sigma}TOF", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsTOF}});
      histosQA.add("QA/He3/hist_nSigmaITS_He3", ";n_{#sigma}ITS", {HistType::kTH1F, {cfgnSigmaBinsITS}});
      histosQA.add("QA/He3/hist_nSigmaITSPt_He3", ";#it{p}_{T};n_{#sigma}ITS", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsITS}});
      histosQA.add("QA/He3/hist_nSigmaRMS_He3", ";n_{#sigma}RMS", {HistType::kTH1F, {cfgnSigmaBinsRMS}});
      histosQA.add("QA/He3/hist_nSigmaRMSPt_He3", ";#it{p}_{T};n_{#sigma}RMS", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsRMS}});
      if (cfgOpenHe3ITSPtCut) {
        histosQA.add("QA/He3/hist_nSigmaITSPt_He3_unCuted", ";#it{p}_{T};n_{#sigma}ITS", {HistType::kTH2F, {cfgaxisptPID, cfgnSigmaBinsITS}});
      }
      if (cfgOpenPlotnSigmaTOFITSPt) {
        histosQA.add("QA/Proton/hist_nSigmaTOFITSPt_Pr", ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
        histosQA.add("QA/He3/hist_nSigmaTOFITSPt_He3", ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxisptPID}});
      }
      if (cfgOpenPlotnSigmaITSTPCPt) {
        histosQA.add("QA/Proton/hist_nSigmaITSTPCPt_Pr", ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxisptPID}});
        histosQA.add("QA/He3/hist_nSigmaITSTPCPt_He3", ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxisptPID}});
      }
      if (cfgOpenPlotnSigmaTOFTPCPt) {
        histosQA.add("QA/Proton/hist_nSigmaTOFTPCPt_Pr", ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
        histosQA.add("QA/He3/hist_nSigmaTOFTPCPt_He3", ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", {HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxisptPID}});
      }
    }
  }
  Produces<aod::PHe3ESEFlags> pidEsePHe3Table;
  void process(TracksPID const& tracks)
  {
    auto tracksWithITSPid = soa::Attach<TracksPID,
                                        aod::pidits::ITSNSigmaPr,
                                        aod::pidits::ITSNSigmaHe>(tracks);
    int8_t pidFlag;
    for (const auto& track : tracksWithITSPid) {
      histosQA.fill(HIST("QA/hist_dEdxTPC_All"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
      histosQA.fill(HIST("QA/hist_pt_All"), track.pt());
      histosQA.fill(HIST("QA/hist_eta_All"), track.eta());
      histosQA.fill(HIST("QA/hist_ITSNcls_All"), track.itsNCls());
      histosQA.fill(HIST("QA/hist_TPCNcls_All"), track.tpcNClsFound());
      histosQA.fill(HIST("QA/hist_ITSChi2NDF_All"), track.itsChi2NCl());
      histosQA.fill(HIST("QA/hist_TPCChi2NDF_All"), track.tpcChi2NCl());
      histosQA.fill(HIST("QA/hist_DCAxy_All"), track.dcaXY());
      histosQA.fill(HIST("QA/hist_DCAz_All"), track.dcaZ());
      if (!trackSelBasic(track)) {
        pidFlag = -1;
      } else {
        int currentPtBinPr = -1, currentPtBinHe3 = -1;
        if (cfgOpenPIDByPtProtonMain || (cfgOpenPIDByPtProtonITS && cfgOpenPIDITSProton)) {
          for (int i = 0; i < static_cast<int>(cfgPtBinProtonPID.value.size()) - 1; ++i) {
            if (track.pt() >= cfgPtBinProtonPID.value[i] && track.pt() < cfgPtBinProtonPID.value[i + 1]) {
              currentPtBinPr = i;
              break;
            }
          }
        }
        if (cfgOpenPIDByPtHe3Main || (cfgOpenPIDByPtHe3ITS && cfgOpenPIDITSHe3)) {
          for (int i = 0; i < static_cast<int>(cfgPtBinHe3PID.value.size()) - 1; ++i) {
            if (track.pt() >= cfgPtBinHe3PID.value[i] && track.pt() < cfgPtBinHe3PID.value[i + 1]) {
              currentPtBinHe3 = i;
              break;
            }
          }
        }
        float nSigmaTPCCutPrPtLower = (currentPtBinPr == -1) ? cfgnSigmaCutTPCProton.value[0] : cfgnSigmaTPCProtonPtLower.value[currentPtBinPr];
        float nSigmaTPCCutPrPtUpper = (currentPtBinPr == -1) ? cfgnSigmaCutTPCProton.value[1] : cfgnSigmaTPCProtonPtUpper.value[currentPtBinPr];
        float nSigmaTOFCutPrPtLower = (currentPtBinPr == -1) ? cfgnSigmaCutTOFProton.value[0] : cfgnSigmaTOFProtonPtLower.value[currentPtBinPr];
        float nSigmaTOFCutPrPtUpper = (currentPtBinPr == -1) ? cfgnSigmaCutTOFProton.value[1] : cfgnSigmaTOFProtonPtUpper.value[currentPtBinPr];
        float nSigmaRMSCutPrPtLower = (currentPtBinPr == -1) ? cfgnSigmaCutRMSProton.value[0] : cfgnSigmaRMSProtonPtLower.value[currentPtBinPr];
        float nSigmaRMSCutPrPtUpper = (currentPtBinPr == -1) ? cfgnSigmaCutRMSProton.value[1] : cfgnSigmaRMSProtonPtUpper.value[currentPtBinPr];
        float nSigmaITSCutPrPtLower = (currentPtBinPr == -1) ? cfgnSigmaCutITSProton.value[0] : cfgnSigmaITSProtonPtLower.value[currentPtBinPr];
        float nSigmaITSCutPrPtUpper = (currentPtBinPr == -1) ? cfgnSigmaCutITSProton.value[1] : cfgnSigmaITSProtonPtUpper.value[currentPtBinPr];
        float nSigmaTPCCutHe3PtLower = (currentPtBinHe3 == -1) ? cfgnSigmaCutTPCHe3.value[0] : cfgnSigmaTPCHe3PtLower.value[currentPtBinHe3];
        float nSigmaTPCCutHe3PtUpper = (currentPtBinHe3 == -1) ? cfgnSigmaCutTPCHe3.value[1] : cfgnSigmaTPCHe3PtUpper.value[currentPtBinHe3];
        float nSigmaTOFCutHe3PtLower = (currentPtBinHe3 == -1) ? cfgnSigmaCutTOFHe3.value[0] : cfgnSigmaTOFHe3PtLower.value[currentPtBinHe3];
        float nSigmaTOFCutHe3PtUpper = (currentPtBinHe3 == -1) ? cfgnSigmaCutTOFHe3.value[1] : cfgnSigmaTOFHe3PtUpper.value[currentPtBinHe3];
        float nSigmaRMSCutHe3PtLower = (currentPtBinHe3 == -1) ? cfgnSigmaCutRMSHe3.value[0] : cfgnSigmaRMSHe3PtLower.value[currentPtBinHe3];
        float nSigmaRMSCutHe3PtUpper = (currentPtBinHe3 == -1) ? cfgnSigmaCutRMSHe3.value[1] : cfgnSigmaRMSHe3PtUpper.value[currentPtBinHe3];
        float nSigmaITSCutHe3PtLower = (currentPtBinHe3 == -1) ? cfgnSigmaCutITSHe3.value[0] : cfgnSigmaITSHe3PtLower.value[currentPtBinHe3];
        float nSigmaITSCutHe3PtUpper = (currentPtBinHe3 == -1) ? cfgnSigmaCutITSHe3.value[1] : cfgnSigmaITSHe3PtUpper.value[currentPtBinHe3];
        float nSigmaMainLowerPr = -999, nSigmaMainUpperPr = -999;
        float nSigmaMainLowerHe3 = -999, nSigmaMainUpperHe3 = -999;
        switch (cfgProtonPIDMode) {
          case 0:
            nSigmaMainLowerPr = nSigmaRMSCutPrPtLower;
            nSigmaMainUpperPr = nSigmaRMSCutPrPtUpper;
            break;
          case 1:
            nSigmaMainLowerPr = nSigmaTPCCutPrPtLower;
            nSigmaMainUpperPr = nSigmaTPCCutPrPtUpper;
            break;
          case 2:
            nSigmaMainLowerPr = nSigmaTOFCutPrPtLower;
            nSigmaMainUpperPr = nSigmaTOFCutPrPtUpper;
            break;
        }
        switch (cfgHe3PIDMode) {
          case 0:
            nSigmaMainLowerHe3 = nSigmaRMSCutHe3PtLower;
            nSigmaMainUpperHe3 = nSigmaRMSCutHe3PtUpper;
            break;
          case 1:
            nSigmaMainLowerHe3 = nSigmaTPCCutHe3PtLower;
            nSigmaMainUpperHe3 = nSigmaTPCCutHe3PtUpper;
            break;
          case 2:
            nSigmaMainLowerHe3 = nSigmaTOFCutHe3PtLower;
            nSigmaMainUpperHe3 = nSigmaTOFCutHe3PtUpper;
            break;
        }
        bool kIsPr = false, kIsHe3 = false;
        // Identify Proton
        if (pidProtonSel(track, nSigmaMainLowerPr, nSigmaMainUpperPr)) {
          kIsPr = true;
          if (cfgOpenPIDITSProton) {
            if (track.itsNSigmaPr() < nSigmaITSCutPrPtLower || track.itsNSigmaPr() > nSigmaITSCutPrPtUpper) {
              kIsPr = false;
            }
          }
        }
        // Identify He3
        if (pidHe3Sel(track, nSigmaMainLowerHe3, nSigmaMainUpperHe3)) {
          kIsHe3 = true;
          if (cfgOpenPIDITSHe3) {
            if (track.itsNSigmaHe() < nSigmaITSCutHe3PtLower || track.itsNSigmaHe() > nSigmaITSCutHe3PtUpper) {
              kIsHe3 = false;
            }
          }
        }
        // Cross track rejection
        if (!cfgOpenAllowCrossTrack) {
          if (kIsPr && kIsHe3) {
            switch (crossTrackID(track)) {
              case 0:
                kIsPr = true;
                kIsHe3 = false;
                break;
              case 1:
                kIsPr = false;
                kIsHe3 = true;
                break;
            }
          }
        }
        // Filter He3 contaimination
        if (cfgOpenHe3ITSPtCut && kIsHe3) {
          if (!cfgQuietMode) {
            histosQA.fill(HIST("QA/He3/hist_nSigmaITSPt_He3_unCuted"), track.pt(), track.itsNSigmaHe());
          }
          if (track.itsNSigmaHe() < fNSigmaITSPt->Eval(track.pt())) {
            kIsHe3 = false;
          }
        }
        pidFlag = (kIsHe3 << 1) | kIsPr;
        // Fill QA histograms
        if (!cfgQuietMode) {
          if (kIsPr) {
            histosQA.fill(HIST("QA/Proton/hist_dEdxTPC_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/Proton/hist_pt_Pr"), track.pt());
            histosQA.fill(HIST("QA/Proton/hist_eta_Pr"), track.eta());
            histosQA.fill(HIST("QA/Proton/hist_ITSNcls_Pr"), track.itsNCls());
            histosQA.fill(HIST("QA/Proton/hist_TPCNcls_Pr"), track.tpcNClsFound());
            histosQA.fill(HIST("QA/Proton/hist_ITSChi2NDF_Pr"), track.itsChi2NCl());
            histosQA.fill(HIST("QA/Proton/hist_TPCChi2NDF_Pr"), track.tpcChi2NCl());
            histosQA.fill(HIST("QA/Proton/hist_DCAxy_Pr"), track.dcaXY());
            histosQA.fill(HIST("QA/Proton/hist_DCAz_Pr"), track.dcaZ());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaTPC_Pr"), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaTPCPt_Pr"), track.pt(), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaTOF_Pr"), track.tofNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaTOFPt_Pr"), track.pt(), track.tofNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaITS_Pr"), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaITSPt_Pr"), track.pt(), track.itsNSigmaPr());
            histosQA.fill(HIST("QA/Proton/hist_nSigmaRMS_Pr"), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));
            histosQA.fill(HIST("QA/Proton/hist_nSigmaRMSPt_Pr"), track.pt(), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()));
            if (cfgOpenPlotnSigmaTOFITSPt) {
              histosQA.fill(HIST("QA/Proton/hist_nSigmaTOFITSPt_Pr"), track.tofNSigmaPr(), track.itsNSigmaPr(), track.pt());
            }
            if (cfgOpenPlotnSigmaITSTPCPt) {
              histosQA.fill(HIST("QA/Proton/hist_nSigmaITSTPCPt_Pr"), track.itsNSigmaPr(), track.tpcNSigmaPr(), track.pt());
            }
            if (cfgOpenPlotnSigmaTOFTPCPt) {
              histosQA.fill(HIST("QA/Proton/hist_nSigmaTOFTPCPt_Pr"), track.tofNSigmaPr(), track.tpcNSigmaPr(), track.pt());
            }
          }
          if (kIsHe3) {
            histosQA.fill(HIST("QA/He3/hist_dEdxTPC_He3"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/He3/hist_pt_He3"), track.pt());
            histosQA.fill(HIST("QA/He3/hist_eta_He3"), track.eta());
            histosQA.fill(HIST("QA/He3/hist_ITSNcls_He3"), track.itsNCls());
            histosQA.fill(HIST("QA/He3/hist_TPCNcls_He3"), track.tpcNClsFound());
            histosQA.fill(HIST("QA/He3/hist_ITSChi2NDF_He3"), track.itsChi2NCl());
            histosQA.fill(HIST("QA/He3/hist_TPCChi2NDF_He3"), track.tpcChi2NCl());
            histosQA.fill(HIST("QA/He3/hist_DCAxy_He3"), track.dcaXY());
            histosQA.fill(HIST("QA/He3/hist_DCAz_He3"), track.dcaZ());
            histosQA.fill(HIST("QA/He3/hist_nSigmaTPC_He3"), track.tpcNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaTPCPt_He3"), track.pt(), track.tpcNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaTOF_He3"), track.tofNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaTOFPt_He3"), track.pt(), track.tofNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaITS_He3"), track.itsNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaITSPt_He3"), track.pt(), track.itsNSigmaHe());
            histosQA.fill(HIST("QA/He3/hist_nSigmaRMS_He3"), std::hypot(track.tpcNSigmaHe(), track.tofNSigmaHe()));
            histosQA.fill(HIST("QA/He3/hist_nSigmaRMSPt_He3"), track.pt(), std::hypot(track.tpcNSigmaHe(), track.tofNSigmaHe()));
            if (cfgOpenPlotnSigmaTOFITSPt) {
              histosQA.fill(HIST("QA/He3/hist_nSigmaTOFITSPt_He3"), track.tofNSigmaHe(), track.itsNSigmaHe(), track.pt());
            }
            if (cfgOpenPlotnSigmaITSTPCPt) {
              histosQA.fill(HIST("QA/He3/hist_nSigmaITSTPCPt_He3"), track.itsNSigmaHe(), track.tpcNSigmaHe(), track.pt());
            }
            if (cfgOpenPlotnSigmaTOFTPCPt) {
              histosQA.fill(HIST("QA/He3/hist_nSigmaTOFTPCPt_He3"), track.tofNSigmaHe(), track.tpcNSigmaHe(), track.pt());
            }
          }
        }
      }
      pidEsePHe3Table(pidFlag);
    }
  }
};

struct FlowEsePHe3 {
  HistogramRegistry histos{"histosmain", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};

  Configurable<float> cfgVtzCut{"cfgVtzCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Centrality min"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Centrality max"};

  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low boundary cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High boundary cut on TPC occupancy"};
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
  Configurable<bool> cfgOpenv2q{"cfgOpenv2q", true, "Open v2(EP)and q calculation for Proton and He3"};
  Configurable<bool> cfgOpenESE{"cfgOpenESE", true, "Open ESE process"};
  Configurable<bool> cfgOpenESEChargeSeperation{"cfgOpenESEChargeSeperation", true, "Open ESE for postive and negative charge repectivily"};
  Configurable<bool> cfgOpenESEProton{"cfgOpenESEProton", true, "Open ESE Proton process"};
  Configurable<bool> cfgOpenESEHe3{"cfgOpenESEHe3", true, "Open ESE He3 process"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, ""};
  ConfigurableAxis cfgaxisCentForQA{"cfgaxisCentForQA", {100, 0, 100}, "centrality for event QA"};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis cfgaxisT0C{"cfgaxisT0C", {70, 0, 70000}, "N_{ch} (T0C)"};
  ConfigurableAxis cfgaxisT0A{"cfgaxisT0A", {200, 0, 200000}, "N_{ch} (T0A)"};
  ConfigurableAxis cfgaxisNchPV{"cfgaxisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};
  ConfigurableAxis cfgaxisq2{"cfgaxisq2", {120, 0, 12}, "Binning for P_{t} PID"};

  EventPlaneHelper helperEP;
  SliceCache cache;

  int detId;
  int refAId;
  int refBId;
  // Additional Event selection cuts - Copy from flowGenericFramework.cxx
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVtzCut) && (aod::cent::centFT0C > cfgCentMin) && (aod::cent::centFT0C < cfgCentMax);
  Filter properPIDfilter = aod::flow_ese_p_he3::nPidFlag >= (int8_t)0; // Only POI

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::PHe3ESEFlags>>> protonTrackSet = ((aod::flow_ese_p_he3::nPidFlag == pid_flags::kProton) || (aod::flow_ese_p_he3::nPidFlag == pid_flags::kProtonHe3));
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::PHe3ESEFlags>>> he3TrackSet = ((aod::flow_ese_p_he3::nPidFlag == pid_flags::kHe3) || (aod::flow_ese_p_he3::nPidFlag == pid_flags::kProtonHe3));

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

  template <typename CollType>
  bool selEvent(const CollType& collision, const int multTrk, const float centrality)
  {
    histos.fill(HIST("QA/histEventCountDetail"), 0.5);
    if (cfgOpenEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (cfgOpenEvSelkIsGoodZvtxFT0vsPV) {
      histos.fill(HIST("QA/histEventCountDetail"), 1.5);
    }
    if (cfgOpenEvSelkNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (cfgOpenEvSelkNoSameBunchPileup) {
      histos.fill(HIST("QA/histEventCountDetail"), 2.5);
    }
    if (cfgOpenEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (cfgOpenEvSelkNoCollInTimeRangeStandard) {
      histos.fill(HIST("QA/histEventCountDetail"), 3.5);
    }
    if (cfgOpenEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (cfgOpenEvSelkIsGoodITSLayersAll) {
      histos.fill(HIST("QA/histEventCountDetail"), 4.5);
    }
    if (cfgOpenEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (cfgOpenEvSelkNoCollInRofStandard) {
      histos.fill(HIST("QA/histEventCountDetail"), 5.5);
    }
    if (cfgOpenEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (cfgOpenEvSelkNoHighMultCollInPrevRof) {
      histos.fill(HIST("QA/histEventCountDetail"), 6.5);
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgOpenEvSelOccupancy && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)) {
      return false;
    }
    if (cfgOpenEvSelOccupancy) {
      histos.fill(HIST("QA/histEventCountDetail"), 7.5);
    }
    if (cfgOpenEvSelMultCorrelationPVTracks) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return false;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return false;
    }
    if (cfgOpenEvSelMultCorrelationPVTracks) {
      histos.fill(HIST("QA/histEventCountDetail"), 8.5);
    }
    if (cfgOpenEvSelMultCorrelationGlobalTracks) {
      if (multTrk < fMultCutLow->Eval(centrality))
        return false;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return false;
    }
    if (cfgOpenEvSelMultCorrelationGlobalTracks) {
      histos.fill(HIST("QA/histEventCountDetail"), 9.5);
    }
    if (cfgOpenEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > event_selection::kFT0AV0ASigma * fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (cfgOpenEvSelV0AT0ACut) {
      histos.fill(HIST("QA/histEventCountDetail"), 10.5);
    }
    return true;
  }

  template <typename CollType>
  void fillHistosQvec(const CollType& collision, int nmode)
  {
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int refAInd = refAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int refBInd = refBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if (nmode == fourier_mode::kMode2) {
      if (collision.qvecAmp()[detId] > 1e-8) {
        histos.fill(HIST("QA/histQvec_CorrL0_V2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.centFT0C());
        histos.fill(HIST("QA/histQvec_CorrL1_V2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.centFT0C());
        histos.fill(HIST("QA/histQvec_CorrL2_V2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.centFT0C());
        histos.fill(HIST("QA/histQvec_CorrL3_V2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.centFT0C());
        histos.fill(HIST("QA/histEvtPl_CorrL0_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], nmode), collision.centFT0C());
        histos.fill(HIST("QA/histEvtPl_CorrL1_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], nmode), collision.centFT0C());
        histos.fill(HIST("QA/histEvtPl_CorrL2_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], nmode), collision.centFT0C());
        histos.fill(HIST("QA/histEvtPl_CorrL3_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), collision.centFT0C());
      }
      if (collision.qvecAmp()[detId] > 1e-8 && collision.qvecAmp()[refAId] > 1e-8 && collision.qvecAmp()[refBId] > 1e-8) {
        histos.fill(HIST("QA/histQvecRes_SigRefAV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), nmode));
        histos.fill(HIST("QA/histQvecRes_SigRefBV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode));
        histos.fill(HIST("QA/histQvecRes_RefARefBV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], nmode), nmode));
      }
    }
  }

  template <typename TrackType>
  float calculateq2(const TrackType tracks, float psi2, float cent, int pidmode) // pidmode 1 for proton , 2 for he3
  {
    int multi = tracks.size();
    if (multi > 0) {
      float q2x = 0, q2y = 0;
      for (const auto& track : tracks) {
        q2x += std::cos(2 * track.phi());
        q2y += std::sin(2 * track.phi());
        if (pidmode == pid_flags::kProton) {
          if (track.sign() > 0) {
            histos.fill(HIST("V2/histCosV2EP_Pr_Pos"), track.pt(), cent, std::cos(2 * (track.phi() - psi2)));
          } else {
            histos.fill(HIST("V2/histCosV2EP_Pr_Neg"), track.pt(), cent, std::cos(2 * (track.phi() - psi2)));
          }
        }
        if (pidmode == pid_flags::kHe3) {
          if (track.sign() > 0) {
            histos.fill(HIST("V2/histCosV2EP_He3_Pos"), track.pt(), cent, std::cos(2 * (track.phi() - psi2)));
          } else {
            histos.fill(HIST("V2/histCosV2EP_He3_Neg"), track.pt(), cent, std::cos(2 * (track.phi() - psi2)));
          }
        }
      }
      return std::hypot(q2x, q2y) / std::sqrt(multi);
    } else {
      return 0;
    }
  }

  template <typename TrackType>
  void processESE(const TrackType tracks, float psi2, float q2, float cent, int pidmode, bool spcharge) // pidmode 1 for proton , 2 for he3
  {
    for (const auto& track : tracks) {
      if (pidmode == pid_flags::kProton) {
        if (spcharge) {
          if (track.sign() > 0) {
            histos.fill(HIST("ESE/hist_v2PosPr_Cent_Pt_q2He3"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
          } else {
            histos.fill(HIST("ESE/hist_v2NegPr_Cent_Pt_q2He3"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
          }
        } else {
          histos.fill(HIST("ESE/hist_v2Pr_Cent_Pt_q2He3"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
        }
      }
      if (pidmode == pid_flags::kHe3) {
        if (spcharge) {
          if (track.sign() > 0) {
            histos.fill(HIST("ESE/hist_v2PosHe3_Cent_Pt_q2Pr"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
          } else {
            histos.fill(HIST("ESE/hist_v2NegHe3_Cent_Pt_q2Pr"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
          }
        } else {
          histos.fill(HIST("ESE/hist_v2He3_Cent_Pt_q2Pr"), track.pt(), cent, q2, std::cos(2 * (track.phi() - psi2)));
        }
      }
    }
  }

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
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec axisvertexz = {100, -15., 15., "vrtx_{Z} [cm]"};
    histos.add("QA/histEventCount", "", {HistType::kTH1F, {{3, 0.0, 3.0}}});
    histos.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(1, "Filtered event");
    histos.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    histos.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(3, "after additional event cut");
    if (cfgUseAdditionalEventCut) {
      histos.add("QA/histEventCountDetail", "Number of Event;; Count", {HistType::kTH1F, {{11, 0, 11}}});
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(1, "after sel8");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(8, "occupancy");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(9, "MultCorrelationPVTracks");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(10, "MultCorrelationGlobalTracks");
      histos.get<TH1>(HIST("QA/histEventCountDetail"))->GetXaxis()->SetBinLabel(11, "cfgEvSelV0AT0ACut");
    }
    if (cfgOpenFullEventQA) {
      histos.add("QA/hist_globalTracks_centT0C_before", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNch}});
      histos.add("QA/hist_PVTracks_centT0C_before", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNchPV}});
      histos.add("QA/hist_globalTracks_PVTracks_before", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histos.add("QA/hist_globalTracks_multT0A_before", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histos.add("QA/hist_globalTracks_multV0A_before", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histos.add("QA/hist_multV0A_multT0A_before", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histos.add("QA/hist_multT0C_centT0C_before", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisT0C}});
      histos.add("QA/hist_globalTracks_centT0C_after", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNch}});
      histos.add("QA/hist_PVTracks_centT0C_after", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisNchPV}});
      histos.add("QA/hist_globalTracks_PVTracks_after", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histos.add("QA/hist_globalTracks_multT0A_after", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histos.add("QA/hist_globalTracks_multV0A_after", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histos.add("QA/hist_multV0A_multT0A_after", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histos.add("QA/hist_multT0C_centT0C_after", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {cfgaxisCentForQA, cfgaxisT0C}});
    }
    histos.add("QA/histVertexZRec", ";vrtx_{Z} [cm];counts", {HistType::kTH1F, {axisvertexz}});
    histos.add("QA/histCentrality", ";Centrality;counts", {HistType::kTH1F, {cfgaxisCentForQA}});
    histos.add("QA/histProtonNum", "ProtonNum;counts", {HistType::kTH1F, {{100, 0, 100}}});
    histos.add("QA/histHe3Num", "He3Num;counts", {HistType::kTH1F, {{20, 0, 20}}});
    histos.add("QA/histQvec_CorrL0_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histos.add("QA/histQvec_CorrL1_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histos.add("QA/histQvec_CorrL2_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histos.add("QA/histQvec_CorrL3_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histos.add("QA/histEvtPl_CorrL0_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histos.add("QA/histEvtPl_CorrL1_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histos.add("QA/histEvtPl_CorrL2_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histos.add("QA/histEvtPl_CorrL3_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histos.add("QA/histQvecRes_SigRefAV2", ";Centrality;Cos(Sig-RefA)", {HistType::kTProfile, {cfgaxisCent}});
    histos.add("QA/histQvecRes_SigRefBV2", ";Centrality;Cos(Sig-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    histos.add("QA/histQvecRes_RefARefBV2", ";Centrality;Cos(RefA-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    if (cfgOpenv2q) {
      histos.add("V2/histCosV2EP_Pr_Pos", ";#it{p}_{T};Centrality", {HistType::kTProfile2D, {cfgaxispt, cfgaxisCent}});
      histos.add("V2/histCosV2EP_Pr_Neg", ";#it{p}_{T};Centrality", {HistType::kTProfile2D, {cfgaxispt, cfgaxisCent}});
      histos.add("V2/histCosV2EP_He3_Pos", ";#it{p}_{T};Centrality", {HistType::kTProfile2D, {cfgaxispt, cfgaxisCent}});
      histos.add("V2/histCosV2EP_He3_Neg", ";#it{p}_{T};Centrality", {HistType::kTProfile2D, {cfgaxispt, cfgaxisCent}});
      histos.add("q2/hist_q2_Cen_Pr", ";q_{2} (TPC);Centrality", {HistType::kTH2F, {cfgaxisq2, cfgaxisCent}});
      histos.add("q2/hist_q2_Cen_He3", ";q_{2} (TPC);Centrality", {HistType::kTH2F, {cfgaxisq2, cfgaxisCent}});
      histos.add("q2/hist_q2_Pr", ";q_{2} (TPC);counts", {HistType::kTH1F, {cfgaxisq2}});
      histos.add("q2/hist_q2_He3", ";q_{2} (TPC);counts", {HistType::kTH1F, {cfgaxisq2}});
    }
    if (cfgOpenESE) {
      if (cfgOpenESEChargeSeperation) {
        if (cfgOpenESEProton) {
          histos.add("ESE/hist_v2PosPr_Cent_Pt_q2He3", ";#it{p}_{T};q_{2}(He3);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
          histos.add("ESE/hist_v2NegPr_Cent_Pt_q2He3", ";#it{p}_{T};q_{2}(He3);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
        }
        if (cfgOpenESEHe3) {
          histos.add("ESE/hist_v2PosHe3_Cent_Pt_q2Pr", ";#it{p}_{T};q_{2}(Proton);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
          histos.add("ESE/hist_v2NegHe3_Cent_Pt_q2Pr", ";#it{p}_{T};q_{2}(Proton);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
        }
      } else {
        if (cfgOpenESEProton) {
          histos.add("ESE/hist_v2Pr_Cent_Pt_q2He3", ";#it{p}_{T};q_{2}(He3);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
        }
        if (cfgOpenESEHe3) {
          histos.add("ESE/hist_v2He3_Cent_Pt_q2Pr", ";#it{p}_{T};q_{2}(Proton);Centrality", HistType::kTProfile3D, {cfgaxispt, cfgaxisq2, cfgaxisCent});
        }
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::Qvectors>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::PHe3ESEFlags>> const& tracks)
  {
    const auto cent = collision.centFT0C();
    histos.fill(HIST("QA/histEventCount"), 0.5);
    if (!collision.sel8())
      return;
    if (tracks.size() < 1)
      return;
    histos.fill(HIST("QA/histEventCount"), 1.5);
    if (cfgOpenFullEventQA) {
      histos.fill(HIST("QA/hist_globalTracks_centT0C_before"), cent, tracks.size());
      histos.fill(HIST("QA/hist_PVTracks_centT0C_before"), cent, collision.multNTracksPV());
      histos.fill(HIST("QA/hist_globalTracks_PVTracks_before"), collision.multNTracksPV(), tracks.size());
      histos.fill(HIST("QA/hist_globalTracks_multT0A_before"), collision.multFT0A(), tracks.size());
      histos.fill(HIST("QA/hist_globalTracks_multV0A_before"), collision.multFV0A(), tracks.size());
      histos.fill(HIST("QA/hist_multV0A_multT0A_before"), collision.multFT0A(), collision.multFV0A());
      histos.fill(HIST("QA/hist_multT0C_centT0C_before"), cent, collision.multFT0C());
    }
    if (cfgUseAdditionalEventCut && !selEvent(collision, tracks.size(), cent)) {
      return;
    }
    histos.fill(HIST("QA/histEventCount"), 2.5);
    histos.fill(HIST("QA/histCentrality"), cent);
    histos.fill(HIST("QA/histVertexZRec"), collision.posZ());
    if (cfgOpenFullEventQA) {
      histos.fill(HIST("QA/hist_globalTracks_centT0C_after"), cent, tracks.size());
      histos.fill(HIST("QA/hist_PVTracks_centT0C_after"), cent, collision.multNTracksPV());
      histos.fill(HIST("QA/hist_globalTracks_PVTracks_after"), collision.multNTracksPV(), tracks.size());
      histos.fill(HIST("QA/hist_globalTracks_multT0A_after"), collision.multFT0A(), tracks.size());
      histos.fill(HIST("QA/hist_globalTracks_multV0A_after"), collision.multFV0A(), tracks.size());
      histos.fill(HIST("QA/hist_multV0A_multT0A_after"), collision.multFT0A(), collision.multFV0A());
      histos.fill(HIST("QA/hist_multT0C_centT0C_after"), cent, collision.multFT0C());
    }
    auto tracksPr = protonTrackSet->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksHe3 = he3TrackSet->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int multiPr = tracksPr.size();
    int multiHe3 = tracksHe3.size();
    // LOGF(info, Form("Collison ID + 1; Proton Num:%d; He3 Num:%d;\n", multiPr, multiHe3));
    histos.fill(HIST("QA/histProtonNum"), multiPr);
    histos.fill(HIST("QA/histHe3Num"), multiHe3);
    if (multiPr < 1 && multiHe3 < 1)
      return; // Reject Collisions without enough POI
    for (auto i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
      int detIndGlobal = detId * 4 + cfgnTotalSystem * 4 * (cfgnMods->at(i) - 2);
      float psiNGlobal = helperEP.GetEventPlane(collision.qvecRe()[detIndGlobal + 3], collision.qvecIm()[detIndGlobal + 3], cfgnMods->at(i));
      if (cfgnMods->at(i) == fourier_mode::kMode2) {
        // LOGF(info, "Process q2\n");
        float q2Proton = calculateq2(tracksPr, psiNGlobal, cent, 1);
        float q2He3 = calculateq2(tracksHe3, psiNGlobal, cent, 2);
        histos.fill(HIST("q2/hist_q2_Pr"), q2Proton);
        histos.fill(HIST("q2/hist_q2_He3"), q2He3);
        histos.fill(HIST("q2/hist_q2_Cen_Pr"), q2Proton, cent);
        histos.fill(HIST("q2/hist_q2_Cen_He3"), q2He3, cent);
        if (cfgOpenESE && multiPr > 0 && multiHe3 > 0) {
          // LOGF(info, "Process ESE\n");
          if (cfgOpenESEProton) {
            processESE(tracksPr, psiNGlobal, q2Proton, cent, 1, cfgOpenESEChargeSeperation);
          }
          if (cfgOpenESEHe3) {
            processESE(tracksHe3, psiNGlobal, q2He3, cent, 2, cfgOpenESEChargeSeperation);
          }
        }
        // LOGF(info, "Process for this event over\n");
      }
      fillHistosQvec(collision, cfgnMods->at(i));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FillPIDcolums>(cfgc),
    adaptAnalysisTask<FlowEsePHe3>(cfgc),
  };
}
