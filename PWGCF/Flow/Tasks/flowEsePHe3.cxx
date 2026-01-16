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

#include <TComplex.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>

#include <chrono>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// o2Physics includes.
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "MathUtils/BetheBlochAleph.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace ese_var_table
{
DECLARE_SOA_COLUMN(EseVtz, eseVtz, float);
DECLARE_SOA_COLUMN(EseCentFT0C, eseCentFT0C, float);
DECLARE_SOA_COLUMN(EsePsi2FT0C, esePsi2FT0C, float);
DECLARE_SOA_COLUMN(Eseq2Tar, eseq2Tar, float);
DECLARE_SOA_COLUMN(Eseq2Ref, eseq2Ref, float);
DECLARE_SOA_COLUMN(Eseq2FT0C, eseq2FT0C, float);
DECLARE_SOA_COLUMN(EseTarSign, eseTarSign, int8_t);
DECLARE_SOA_COLUMN(EseTarTPCInnerParam, eseTarTPCInnerParam, float);
DECLARE_SOA_COLUMN(EseTarTPCSignal, eseTarTPCSignal, float);
DECLARE_SOA_COLUMN(EseTarPt, eseTarPt, float);
DECLARE_SOA_COLUMN(EseTarEta, eseTarEta, float);
DECLARE_SOA_COLUMN(EseTarPhi, eseTarPhi, float);
DECLARE_SOA_COLUMN(EseTarDCAxy, eseTarDCAxy, float);
DECLARE_SOA_COLUMN(EseTarDCAz, eseTarDCAz, float);
DECLARE_SOA_COLUMN(EseTarTPCNcls, eseTarTPCNcls, uint8_t);
DECLARE_SOA_COLUMN(EseTarITSNcls, eseTarITSNcls, uint8_t);
DECLARE_SOA_COLUMN(EseTarTPCChi2NDF, eseTarTPCChi2NDF, float);
DECLARE_SOA_COLUMN(EseTarITSChi2NDF, eseTarITSChi2NDF, float);
DECLARE_SOA_COLUMN(EseTarTPCNSigma, eseTarTPCNSigma, float);
DECLARE_SOA_COLUMN(EseTarTOFNSigma, eseTarTOFNSigma, float);
DECLARE_SOA_COLUMN(EseTarITSNSigma, eseTarITSNSigma, float);
DECLARE_SOA_COLUMN(EseTarITSClusSize, eseTarITSClusSize, uint32_t);
} // namespace ese_var_table

DECLARE_SOA_TABLE(ESETable, "AOD", "ESETable",
                  ese_var_table::EseVtz,
                  ese_var_table::EseCentFT0C,
                  ese_var_table::EsePsi2FT0C,
                  ese_var_table::Eseq2Tar,
                  ese_var_table::Eseq2Ref,
                  ese_var_table::Eseq2FT0C,
                  ese_var_table::EseTarSign,
                  ese_var_table::EseTarTPCInnerParam,
                  ese_var_table::EseTarTPCSignal,
                  ese_var_table::EseTarPt,
                  ese_var_table::EseTarEta,
                  ese_var_table::EseTarPhi,
                  ese_var_table::EseTarDCAxy,
                  ese_var_table::EseTarDCAz,
                  ese_var_table::EseTarTPCNcls,
                  ese_var_table::EseTarITSNcls,
                  ese_var_table::EseTarTPCChi2NDF,
                  ese_var_table::EseTarITSChi2NDF,
                  ese_var_table::EseTarTPCNSigma,
                  ese_var_table::EseTarTOFNSigma,
                  ese_var_table::EseTarITSNSigma,
                  ese_var_table::EseTarITSClusSize);
} // namespace o2::aod

struct ESECandidate {
  float vtz;
  float centFT0C;
  float psi2FT0C;
  float q2Tar;
  float q2Ref;
  float q2FT0C;
  int8_t signTar;
  float tpcInnerParamTar;
  float tpcSignalTar;
  float ptTar;
  float etaTar;
  float phiTar;
  float dcaXYTar;
  float dcaZTar;
  uint8_t tpcNclsTar;
  uint8_t itsNclsTar;
  float tpcChi2NDFTar;
  float itsChi2NDFTar;
  float tpcNSigmaTar;
  float tofNSigmaTar;
  float itsNSigmaTar;
  uint32_t itsClusSizeTar;
};

struct ESEReference {
  int8_t signRef;
  float ptRef;
  float v2Ref;
};

namespace ese_parameters
{
constexpr uint8_t kProton = 0;
constexpr uint8_t kDeuteron = 1;
constexpr uint8_t kTriton = 2;
constexpr uint8_t kHe3 = 3;
constexpr uint8_t kAlpha = 4;
constexpr uint8_t kPion = 5;
constexpr uint8_t kKaon = 6;
constexpr uint8_t kHadron = 7;
constexpr int kFT0AV0ASigma = 5;
constexpr int kRMSMode = 0;
constexpr int kTPCMode = 1;
constexpr int kTOFOnlyMode = 2;
constexpr float Amplitudelow = 1e-8;
constexpr float Charges[8]{1.f, 1.f, 1.f, 2.f, 2.f, 1.f, 1.f, 1.f};
constexpr float Masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};
constexpr double BetheBlochDefault[5][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-239.99, 1.155, 1.099, 1.137, 1.006, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09},
  {-586.66, 1.859, 4.435, 0.282, 3.201, 0.09}};
constexpr double BbMomScalingDefault[5][2]{// 0:poscharged 1:negcharged
                                           {1., 1.},
                                           {1., 1.},
                                           {1., 1.},
                                           {1., 1.},
                                           {1., 1.}};
constexpr int Open3DPIDPlots[3][1]{
  {0},
  {0},
  {0}};
constexpr int OpenEvSel[10][1]{
  {1},
  {1},
  {1},
  {1},
  {1},
  {1},
  {1},
  {1},
  {0},
  {1}};
constexpr int OpenTrackSel[7][1]{
  {1},
  {1},
  {1},
  {1},
  {1},
  {1},
  {1}};
constexpr double TPCnSigmaCutDefault[7][2]{
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.}};
constexpr double ITSnSigmaCutDefault[7][2]{
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.}};
constexpr double POverZPreselection[7][2]{
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.}};
constexpr double EtaPreselection[7][2]{
  {0.9},
  {0.9},
  {0.9},
  {0.8},
  {0.9},
  {0.9},
  {0.9}};
constexpr double TPCNclsPreselection[7][2]{
  {50, 160},
  {50, 160},
  {50, 160},
  {100, 160},
  {50, 160},
  {50, 160},
  {50, 160}};
constexpr double ITSNclsPreselection[7][2]{
  {5, 7},
  {5, 7},
  {5, 7},
  {5, 7},
  {5, 7},
  {5, 7},
  {5, 7}};
constexpr double TPCChi2Preselection[7][2]{
  {0, 10},
  {0, 10},
  {0, 10},
  {0.5, 4},
  {0, 10},
  {0, 10},
  {0, 10}};
constexpr double ITSChi2Preselection[7][2]{
  {0, 36},
  {0, 36},
  {0, 36},
  {0, 36},
  {0, 36},
  {0, 36},
  {0, 36}};
constexpr double DCAxyPreselection[7][2]{
  {1},
  {1},
  {1},
  {0.1},
  {1},
  {1},
  {1},
};
constexpr double DCAzPreselection[7][2]{
  {5},
  {5},
  {5},
  {1},
  {5},
  {5},
  {5},
};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> namesFull{"proton", "deuteron", "triton", "He3", "alpha", "pion", "kaon"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> plot3DPIDNames{"TOF vs ITS", "ITS vs TPC", "TOF vs TPC"};
static const std::vector<std::string> openEventSelNames{"EvSelkIsGoodZvtxFT0vsPV", "EvSelkNoSameBunchPileup", "EvSelkNoCollInTimeRangeStandard", "EvSelkIsGoodITSLayersAll", "EvSelkNoCollInRofStandard", "EvSelkNoHighMultCollInPrevRof", "EvSelOccupancy", "EvSelMultCorrelationPVTracks", "EvSelMultCorrelationGlobalTracks", "EvSelV0AT0ACut"};
static const std::vector<std::string> openTrackSelNames{"passedITSNCls", "passedITSChi2NDF", "passedITSHits", "passedTPCChi2NDF", "passedTPCCrossedRowsOverNCls", "passedDCAxy", "passedDCAz"};
static const std::vector<std::string> plot3DConfigNames{"Open related 3D nSigma plots"};
static const std::vector<std::string> openEventSelConfigNames{"Open related event selection options"};
static const std::vector<std::string> openTrackSelConfigNames{"Open track selection from TrackSelection table"};
static const std::vector<std::string> pidTPCnSigmaNames{"n#sigma_{TPC} Low", "n#sigma_{TPC} High"};
static const std::vector<std::string> pidITSnSigmaNames{"n#sigma_{ITS} Low", "n#sigma_{ITS} High"};
static const std::vector<std::string> pidPOverZNames{"p/z Low", "p/z High"};
static const std::vector<std::string> pidEtaNames{"Abs Eta Max"};
static const std::vector<std::string> pidTPCNclsNames{"TPCNcls Low", "TPCNcls High"};
static const std::vector<std::string> pidITSNclsNames{"ITSNcls Low", "ITSNcls High"};
static const std::vector<std::string> pidTPCChi2Names{"TPCChi2 Low", "TPCChi2 High"};
static const std::vector<std::string> pidITSChi2Names{"ITSChi2 Low", "ITSChi2 High"};
static const std::vector<std::string> pidDCAxyNames{"Abs DCAxy Max"};
static const std::vector<std::string> pidDCAzNames{"Abs DCAz Max"};
std::vector<ESECandidate> eseCandidates;
std::vector<ESEReference> eseReferences;
// Tar ptr
std::shared_ptr<TH1> hPIDQATar1D[12];
std::shared_ptr<TH2> hPIDQATar2D[4];
std::shared_ptr<TH3> hPIDQATar3D[3];
std::shared_ptr<TProfile3D> hv2Tar[2];
std::shared_ptr<TH1> hESEQATar1D[8];
std::shared_ptr<TH2> hESEQATar2D[2];
std::shared_ptr<THnSparse> hESETar;
// Ref ptr
std::shared_ptr<TH1> hPIDQARef1D[12];
std::shared_ptr<TH2> hPIDQARef2D[4];
std::shared_ptr<TH3> hPIDQARef3D[3];
std::shared_ptr<TProfile3D> hv2Ref[2];
std::shared_ptr<TH1> hESEQARef1D[8];
std::shared_ptr<TH2> hESEQARef2D[2];
// FT0C general q2 QA
std::shared_ptr<TH1> hESEQAFT0C1D[2];
std::shared_ptr<TH2> hESEQAFT0C2D[2];
} // namespace ese_parameters

using TracksPIDFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>;

struct FlowEsePHe3 {

  EventPlaneHelper helperEP;
  o2::aod::ITSResponse itsResponse;
  HistogramRegistry histsESE{"histsESE", {}, OutputObjHandlingPolicy::AnalysisObject};
  // process POI control
  Configurable<std::string> cfgTarName{"cfgTarName", "kHe3", "Name of the v2 particle: kProton, kDeuteron, kTriton, kHe3, kAlpha"};
  Configurable<std::string> cfgRefName{"cfgRefName", "kProton", "Name of the q2 reference particle: kProton, kDeuteron, kTriton, kHe3, kAlpha"};
  // total control config
  Configurable<bool> cfgOpenAllowCrossTrack{"cfgOpenAllowCrossTrack", true, "Allow one track to be identified as different kind of PID particles"};
  Configurable<bool> cfgOpenFullEventQA{"cfgOpenFullEventQA", true, "Open full QA plots for event QA"};
  Configurable<bool> cfgOpenPIDQA{"cfgOpenPIDQA", true, "Open PID QA plots"};
  Configurable<bool> cfgOpenv2Tar{"cfgOpenv2Tar", true, "Open v2(EP) for Tar patricle"};
  Configurable<bool> cfgOpenv2Ref{"cfgOpenv2Ref", true, "Open v2(EP) for Ref patricle"};
  Configurable<bool> cfgOpenESE{"cfgOpenESE", true, "Open ESE plots"};
  Configurable<bool> cfgOpenESEQA{"cfgOpenESEQA", true, "Open ESE QA plots"};
  Configurable<LabeledArray<int>> cfgOpen3DPIDPlots{"cfgOpen3DPIDPlots", {ese_parameters::Open3DPIDPlots[0], 3, 1, ese_parameters::plot3DPIDNames, ese_parameters::plot3DConfigNames}, "3D PID QA Plots switch configuration"};
  // Qvec configs
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};
  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};
  // pre event selection(filter)
  Configurable<float> cfgVtzCut{"cfgVtzCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Centrality min"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.0f, "Centrality max"};
  // event selection configs
  Configurable<int> cfgCutOccupancyLow{"cfgCutOccupancyLow", 0, "Low boundary cut on TPC occupancy"};
  Configurable<int> cfgCutOccupancyHigh{"cfgCutOccupancyHigh", 3000, "High boundary cut on TPC occupancy"};
  Configurable<LabeledArray<int>> cfgOpenEvSel{"cfgOpenEvSel", {ese_parameters::OpenEvSel[0], 10, 1, ese_parameters::openEventSelNames, ese_parameters::openEventSelConfigNames}, "Event selection switch configuration"};
  // track selection configs
  Configurable<LabeledArray<int>> cfgOpenTrackSel{"cfgOpenTrackSel", {ese_parameters::OpenTrackSel[0], 7, 1, ese_parameters::openTrackSelNames, ese_parameters::openTrackSelConfigNames}, "Track selection switch configuration"};
  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxPtPID{"cfgMaxPtPID", 99.9, "Maximum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.9, "Maximum track #eta for PID"};
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
  Configurable<float> cfgPtMaxforTPCOnlyPIDProton{"cfgPtMaxforTPCOnlyPIDProton", 0.4, "Maxmium track pt for TPC only PID, at RMS PID mode for proton"};
  Configurable<float> cfgPtMaxforTPCOnlyPIDPion{"cfgPtMaxforTPCOnlyPIDPion", 0.4, "Maxmium track pt for TPC only PID, at RMS PID mode for pion"};
  Configurable<float> cfgPtMaxforTPCOnlyPIDKaon{"cfgPtMaxforTPCOnlyPIDKaon", 0.4, "Maxmium track pt for TPC only PID, at RMS PID mode for kaon"};
  // PID configs
  Configurable<bool> cfgOpenITSPreselection{"cfgOpenITSPreselection", false, "Use nSigma ITS preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgPOverZPreselection{"cfgPOverZPreselection", {ese_parameters::POverZPreselection[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidPOverZNames}, "P/Z preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgEtaPreselection{"cfgEtaPreselection", {ese_parameters::EtaPreselection[0], 7, 1, ese_parameters::namesFull, ese_parameters::pidEtaNames}, "Eta preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgTPCNclsPreselection{"cfgTPCNclsPreselection", {ese_parameters::TPCNclsPreselection[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidTPCNclsNames}, "TPCNcls preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgITSNclsPreselection{"cfgITSNclsPreselection", {ese_parameters::ITSNclsPreselection[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidITSNclsNames}, "ITSNcls preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgTPCChi2Preselection{"cfgTPCChi2Preselection", {ese_parameters::TPCChi2Preselection[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidTPCChi2Names}, "TPCChi2 preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgITSChi2Preselection{"cfgITSChi2Preselection", {ese_parameters::ITSChi2Preselection[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidITSChi2Names}, "ITSChi2 preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgDCAxyPreselection{"cfgDCAxyPreselection", {ese_parameters::DCAxyPreselection[0], 7, 1, ese_parameters::namesFull, ese_parameters::pidDCAxyNames}, "DCAxy preselection for light nuclei"};
  Configurable<LabeledArray<double>> cfgDCAzPreselection{"cfgDCAzPreselection", {ese_parameters::DCAzPreselection[0], 7, 1, ese_parameters::namesFull, ese_parameters::pidDCAzNames}, "DCAz preselection for light nuclei"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFProton{"cfgnSigmaCutTOFProton", {-1.5, 1.5}, "TOF nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSProton{"cfgnSigmaCutRMSProton", {-3, 3}, "RMS nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFPion{"cfgnSigmaCutTOFPion", {-1.5, 1.5}, "TOF nsigma cut limit for Pion"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSPion{"cfgnSigmaCutRMSPion", {-3, 3}, "RMS nsigma cut limit for Pion"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFKaon{"cfgnSigmaCutTOFKaon", {-1.5, 1.5}, "TOF nsigma cut limit for Kaon"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSKaon{"cfgnSigmaCutRMSKaon", {-3, 3}, "RMS nsigma cut limit for Kaon"};
  Configurable<bool> cfgUseSelfnSigmaTPCProton{"cfgUseSelfnSigmaTPCProton", true, "Use self nSigma TPC for Proton PID"};
  Configurable<int> cfgProtonPIDMode{"cfgProtonPIDMode", 2, "Proton PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};
  Configurable<int> cfgPionPIDMode{"cfgPionPIDMode", 2, "Pion PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};
  Configurable<int> cfgKaonPIDMode{"cfgKaonPIDMode", 2, "Kaon PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};
  Configurable<LabeledArray<double>> cfgnSigmaTPC{"cfgnSigmaTPC", {ese_parameters::TPCnSigmaCutDefault[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidTPCnSigmaNames}, "TPC nSigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgnSigmaITS{"cfgnSigmaITS", {ese_parameters::ITSnSigmaCutDefault[0], 7, 2, ese_parameters::namesFull, ese_parameters::pidITSnSigmaNames}, "ITS nSigma selection for light nuclei"};
  // PID BBself paras config
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {ese_parameters::BetheBlochDefault[0], 5, 6, ese_parameters::names, ese_parameters::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {ese_parameters::BbMomScalingDefault[0], 5, 2, ese_parameters::names, ese_parameters::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", true, "If true, divide tpcInnerParam by the electric charge"};
  // Axias configs
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsTOF{"cfgnSigmaBinsTOF", {200, -5.f, 5.f}, "Binning for n sigma TOF"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma ITS"};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, "Binning for P_{t}"};
  ConfigurableAxis cfgaxisDCAz{"cfgaxisDCAz", {200, -1, 1}, "Binning for DCAz"};
  ConfigurableAxis cfgaxisDCAxy{"cfgaxisDCAxy", {100, -0.5, 0.5}, "Binning for DCAxy"};
  ConfigurableAxis cfgaxisChi2Ncls{"cfgaxisChi2Ncls", {100, 0, 30}, "Binning for Chi2Ncls TPC/ITS"};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};
  ConfigurableAxis cfgaxisq2Tar{"cfgaxisq2Tar", {100, 0, 2}, "Binning for q_{2} traget particle"};
  ConfigurableAxis cfgaxisq2Ref{"cfgaxisq2Ref", {120, 0, 12}, "Binning for q_{2} reference particle"};
  ConfigurableAxis cfgaxisq2FT0C{"cfgaxisq2FT0C", {900, 0, 900}, "Binning for q_{2} FT0C"};
  ConfigurableAxis cfgq2NumeratorTar{"cfgq2NumeratorTar", {20, 0, 2}, "q2 Numerator bin for tar particle"};
  ConfigurableAxis cfgq2NumeratorRef{"cfgq2NumeratorRef", {100, 0, 10}, "q2 Numerator bin for ref particle"};
  ConfigurableAxis cfgq2DenominatorTar{"cfgq2DenominatorTar", {20, 0, 2}, "q2 Denominator bin for tar particle"};
  ConfigurableAxis cfgq2DenominatorRef{"cfgq2DenominatorRef", {50, 0, 5}, "q2 Denominator bin for ref particle"};

  uint8_t poiTar;
  uint8_t poiRef;
  int detId;
  int refAId;
  int refBId;
  int detInd;
  int refAInd;
  int refBInd;
  // Function for He3 purity cut refered from luca's slides
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fT0AV0AMean = nullptr;
  TF1* fT0AV0ASigma = nullptr;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgVtzCut) && (aod::cent::centFT0C > cfgCentMin) && (aod::cent::centFT0C < cfgCentMax);

  Produces<o2::aod::ESETable> eseTable;

  template <typename TrackType>
  float getNSigmaTPCSelfBB(const TrackType track, uint8_t POI)
  {
    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
    const int iC{track.sign() < 0};
    switch (POI) {
      case ese_parameters::kProton: {
        const double bgScaling[2]{ese_parameters::Charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / ese_parameters::Masses[0], ese_parameters::Charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / ese_parameters::Masses[0]};
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(0u, 0u), cfgBetheBlochParams->get(0u, 1u), cfgBetheBlochParams->get(0u, 2u), cfgBetheBlochParams->get(0u, 3u), cfgBetheBlochParams->get(0u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(0u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kDeuteron: {
        const double bgScaling[2]{ese_parameters::Charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / ese_parameters::Masses[1], ese_parameters::Charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / ese_parameters::Masses[1]};
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(1u, 0u), cfgBetheBlochParams->get(1u, 1u), cfgBetheBlochParams->get(1u, 2u), cfgBetheBlochParams->get(1u, 3u), cfgBetheBlochParams->get(1u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(1u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kTriton: {
        const double bgScaling[2]{ese_parameters::Charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / ese_parameters::Masses[2], ese_parameters::Charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / ese_parameters::Masses[2]};
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(2u, 0u), cfgBetheBlochParams->get(2u, 1u), cfgBetheBlochParams->get(2u, 2u), cfgBetheBlochParams->get(2u, 3u), cfgBetheBlochParams->get(2u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(2u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kHe3: {
        const double bgScaling[2]{ese_parameters::Charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / ese_parameters::Masses[3], ese_parameters::Charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / ese_parameters::Masses[3]};
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(3u, 0u), cfgBetheBlochParams->get(3u, 1u), cfgBetheBlochParams->get(3u, 2u), cfgBetheBlochParams->get(3u, 3u), cfgBetheBlochParams->get(3u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(3u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kAlpha: {
        const double bgScaling[2]{ese_parameters::Charges[4] * cfgMomentumScalingBetheBloch->get(4u, 0u) / ese_parameters::Masses[4], ese_parameters::Charges[4] * cfgMomentumScalingBetheBloch->get(4u, 1u) / ese_parameters::Masses[4]};
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(4u, 0u), cfgBetheBlochParams->get(4u, 1u), cfgBetheBlochParams->get(4u, 2u), cfgBetheBlochParams->get(4u, 3u), cfgBetheBlochParams->get(4u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(4u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      default:
        return -99.f;
    }
  }

  template <typename TrackType>
  float getNSigmaTPC(const TrackType track, uint8_t POI)
  {
    switch (POI) {
      case ese_parameters::kProton: {
        float nSigmaUse = (cfgUseSelfnSigmaTPCProton ? getNSigmaTPCSelfBB(track, ese_parameters::kProton) : track.tpcNSigmaPr());
        return nSigmaUse;
      }

      case ese_parameters::kDeuteron: {
        return getNSigmaTPCSelfBB(track, ese_parameters::kDeuteron);
      }

      case ese_parameters::kTriton: {
        return getNSigmaTPCSelfBB(track, ese_parameters::kTriton);
      }

      case ese_parameters::kHe3: {
        return getNSigmaTPCSelfBB(track, ese_parameters::kHe3);
      }

      case ese_parameters::kAlpha: {
        return getNSigmaTPCSelfBB(track, ese_parameters::kAlpha);
      }

      case ese_parameters::kPion: {
        return track.tpcNSigmaPi();
      }

      case ese_parameters::kKaon: {
        return track.tpcNSigmaKa();
      }

      default:
        LOGF(error, "Unknown POI: %d", POI);
        return 0;
    }
  }

  template <typename TrackType>
  float getNSigmaTOF(const TrackType track, uint8_t POI)
  {
    switch (POI) {
      case ese_parameters::kProton: {
        return track.tofNSigmaPr();
      }
      case ese_parameters::kDeuteron: {
        return track.tofNSigmaDe();
      }
      case ese_parameters::kTriton: {
        return track.tofNSigmaTr();
      }
      case ese_parameters::kHe3: {
        return track.tofNSigmaHe();
      }
      case ese_parameters::kAlpha: {
        return track.tofNSigmaAl();
      }
      case ese_parameters::kPion: {
        return track.tofNSigmaPi();
      }
      case ese_parameters::kKaon: {
        return track.tofNSigmaKa();
      }
      default:
        return -99.f;
    }
  }

  template <typename TrackType>
  float getNSigmaITS(const TrackType track, uint8_t POI)
  {
    switch (POI) {
      case ese_parameters::kProton: {
        return itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
      }
      case ese_parameters::kDeuteron: {
        return itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track);
      }
      case ese_parameters::kTriton: {
        return itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
      }
      case ese_parameters::kHe3: {
        return itsResponse.nSigmaITS<o2::track::PID::Helium3>(track);
      }
      case ese_parameters::kAlpha: {
        return itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
      }
      case ese_parameters::kPion: {
        return itsResponse.nSigmaITS<o2::track::PID::Pion>(track);
      }
      case ese_parameters::kKaon: {
        return itsResponse.nSigmaITS<o2::track::PID::Kaon>(track);
      }
      default:
        return -99.f;
    }
  }

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

  template <typename T>
  uint8_t getPOI(const T& POI)
  {
    if (POI.value == "kProton") {
      return ese_parameters::kProton;
    } else if (POI.value == "kDeuteron") {
      return ese_parameters::kDeuteron;
    } else if (POI.value == "kTriton") {
      return ese_parameters::kTriton;
    } else if (POI.value == "kHe3") {
      return ese_parameters::kHe3;
    } else if (POI.value == "kAlpha") {
      return ese_parameters::kAlpha;
    } else if (POI.value == "kPion") {
      return ese_parameters::kPion;
    } else if (POI.value == "kKaon") {
      return ese_parameters::kKaon;
    } else if (POI.value == "kHadron") {
      return ese_parameters::kHadron;
    } else {
      LOGF(error, "Unknown POIstr: %s", POI.value.c_str());
      return 0;
    }
  }

  float calculateq2(const float qx, const float qy, const int multi)
  {
    if (multi <= 0) {
      return 0.f;
    } else {
      return std::hypot(qx, qy) / std::sqrt(static_cast<float>(multi));
    }
  }

  template <typename CollType>
  bool eventSelBasic(const CollType& collision, const int64_t multTrk, const float centrality, bool fillQA)
  {
    if (!collision.sel8())
      return false;
    if (fillQA) {
      histsESE.fill(HIST("EventQA/histEventCount"), 0.5);
    }
    if (cfgOpenEvSel->get(0u) && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(0u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 1.5);
    }
    if (cfgOpenEvSel->get(1u) && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(1u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 2.5);
    }
    if (cfgOpenEvSel->get(2u) && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(2u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 3.5);
    }
    if (cfgOpenEvSel->get(3u) && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(3u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 4.5);
    }
    if (cfgOpenEvSel->get(4u) && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(4u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 5.5);
    }
    if (cfgOpenEvSel->get(5u) && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(5u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 6.5);
    }
    auto multNTracksPV = collision.multNTracksPV();
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgOpenEvSel->get(6u) && (occupancy < cfgCutOccupancyLow || occupancy > cfgCutOccupancyHigh)) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(6u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 7.5);
    }
    if (cfgOpenEvSel->get(7u)) {
      if (multNTracksPV < fMultPVCutLow->Eval(centrality))
        return false;
      if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
        return false;
    }
    if (fillQA && cfgOpenEvSel->get(7u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 8.5);
    }
    if (cfgOpenEvSel->get(8u)) {
      if (multTrk < fMultCutLow->Eval(centrality))
        return false;
      if (multTrk > fMultCutHigh->Eval(centrality))
        return false;
    }
    if (fillQA && cfgOpenEvSel->get(8u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 9.5);
    }
    if (cfgOpenEvSel->get(9u) && (std::fabs(collision.multFV0A() - fT0AV0AMean->Eval(collision.multFT0A())) > ese_parameters::kFT0AV0ASigma * fT0AV0ASigma->Eval(collision.multFT0A()))) {
      return false;
    }
    if (fillQA && cfgOpenEvSel->get(9u)) {
      histsESE.fill(HIST("EventQA/histEventCount"), 10.5);
    }
    if (fillQA) {
      histsESE.fill(HIST("EventQA/histVtz"), collision.posZ());
      histsESE.fill(HIST("EventQA/histCent"), centrality);
    }
    return true;
  }

  template <typename TrackType>
  bool trackSelBasic(const TrackType track)
  {
    if ((track.pt() < cfgMinPtPID) || (track.pt() > cfgMaxPtPID))
      return false;
    if (std::abs(track.eta()) > cfgMaxEtaPID)
      return false;
    if (cfgOpenTrackSel->get(0u)) {
      if (!track.passedITSNCls())
        return false;
    } else {
      if (track.itsNCls() < cfgMinITSCls || track.itsNCls() > cfgMaxITSCls)
        return false;
    }
    if (cfgOpenTrackSel->get(1u)) {
      if (!track.passedITSChi2NDF())
        return false;
    } else {
      if (track.itsChi2NCl() < cfgMinChi2NClITS || track.itsChi2NCl() > cfgMaxChi2NClITS)
        return false;
    }
    if (cfgOpenTrackSel->get(2u)) {
      if (!track.passedITSHits())
        return false;
    }
    if (cfgOpenTrackSel->get(3u)) {
      if (!track.passedTPCChi2NDF())
        return false;
    } else {
      if (track.tpcChi2NCl() < cfgMinTPCChi2NCl || track.tpcChi2NCl() > cfgMaxTPCChi2NCl)
        return false;
    }
    if (cfgOpenTrackSel->get(4u)) {
      if (!track.passedTPCCrossedRowsOverNCls())
        return false;
    }
    if (cfgOpenTrackSel->get(5u)) {
      if (!track.passedDCAxy())
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cfgMaxDCAxy)
        return false;
    }
    if (cfgOpenTrackSel->get(6u)) {
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

  template <typename CollType>
  void fillEventQAhistBe(const CollType collision, const int multTrk, const float centrality)
  {
    histsESE.fill(HIST("EventQA/hist_globalTracks_centT0C_before"), centrality, multTrk);
    histsESE.fill(HIST("EventQA/hist_PVTracks_centT0C_before"), centrality, collision.multNTracksPV());
    histsESE.fill(HIST("EventQA/hist_globalTracks_PVTracks_before"), collision.multNTracksPV(), multTrk);
    histsESE.fill(HIST("EventQA/hist_globalTracks_multT0A_before"), collision.multFT0A(), multTrk);
    histsESE.fill(HIST("EventQA/hist_globalTracks_multV0A_before"), collision.multFV0A(), multTrk);
    histsESE.fill(HIST("EventQA/hist_multV0A_multT0A_before"), collision.multFT0A(), collision.multFV0A());
    histsESE.fill(HIST("EventQA/hist_multT0C_centT0C_before"), centrality, collision.multFT0C());
  }

  template <typename CollType>
  void fillEventQAhistAf(const CollType collision, const int multTrk, const float centrality)
  {
    histsESE.fill(HIST("EventQA/hist_globalTracks_centT0C_after"), centrality, multTrk);
    histsESE.fill(HIST("EventQA/hist_PVTracks_centT0C_after"), centrality, collision.multNTracksPV());
    histsESE.fill(HIST("EventQA/hist_globalTracks_PVTracks_after"), collision.multNTracksPV(), multTrk);
    histsESE.fill(HIST("EventQA/hist_globalTracks_multT0A_after"), collision.multFT0A(), multTrk);
    histsESE.fill(HIST("EventQA/hist_globalTracks_multV0A_after"), collision.multFV0A(), multTrk);
    histsESE.fill(HIST("EventQA/hist_multV0A_multT0A_after"), collision.multFT0A(), collision.multFV0A());
    histsESE.fill(HIST("EventQA/hist_multT0C_centT0C_after"), centrality, collision.multFT0C());
  }

  template <typename TrackType>
  void fillTrackQAhist(const TrackType track)
  {
    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
    histsESE.fill(HIST("TrackQA/hist_dEdxTPC_All"), track.sign() * correctedTpcInnerParam, track.tpcSignal());
    histsESE.fill(HIST("TrackQA/hist_pt_All"), track.pt());
    histsESE.fill(HIST("TrackQA/hist_eta_All"), track.eta());
    histsESE.fill(HIST("TrackQA/hist_phi_All"), track.phi());
    histsESE.fill(HIST("TrackQA/hist_DCAxy_All"), track.dcaXY());
    histsESE.fill(HIST("TrackQA/hist_DCAz_All"), track.dcaZ());
    histsESE.fill(HIST("TrackQA/hist_ITSNcls_All"), track.itsNCls());
    histsESE.fill(HIST("TrackQA/hist_TPCNcls_All"), track.tpcNClsFound());
    histsESE.fill(HIST("TrackQA/hist_ITSChi2NDF_All"), track.itsChi2NCl());
    histsESE.fill(HIST("TrackQA/hist_TPCChi2NDF_All"), track.tpcChi2NCl());
    histsESE.fill(HIST("TrackQA/hist_MomRes_All"), track.p(), 1 - (correctedTpcInnerParam / track.p()));
    if (heliumPID) {
      histsESE.fill(HIST("TrackQA/hist_He3AlphaTrackcounts_All"), 1.5);
    } else {
      histsESE.fill(HIST("TrackQA/hist_He3AlphaTrackcounts_All"), 0.5);
    }
  }

  template <typename CollType>
  void fillHistosQvec(const CollType& collision)
  {
    if (collision.qvecAmp()[detId] > ese_parameters::Amplitudelow) {
      histsESE.fill(HIST("PlanQA/histQvec_CorrL0_V2"), collision.qvecRe()[detInd], collision.qvecIm()[detInd], collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histQvec_CorrL1_V2"), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histQvec_CorrL2_V2"), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histQvec_CorrL3_V2"), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histEvtPl_CorrL0_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd], collision.qvecIm()[detInd], 2), collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histEvtPl_CorrL1_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1], 2), collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histEvtPl_CorrL2_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2], 2), collision.centFT0C());
      histsESE.fill(HIST("PlanQA/histEvtPl_CorrL3_V2"), helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], 2), collision.centFT0C());
    }
    if (collision.qvecAmp()[detId] > ese_parameters::Amplitudelow && collision.qvecAmp()[refAId] > ese_parameters::Amplitudelow && collision.qvecAmp()[refBId] > ese_parameters::Amplitudelow) {
      histsESE.fill(HIST("PlanQA/histQvecRes_SigRefAV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], 2), helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], 2), 2));
      histsESE.fill(HIST("PlanQA/histQvecRes_SigRefBV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], 2), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], 2), 2));
      histsESE.fill(HIST("PlanQA/histQvecRes_RefARefBV2"), collision.centFT0C(), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[refAInd + 3], collision.qvecIm()[refAInd + 3], 2), helperEP.GetEventPlane(collision.qvecRe()[refBInd + 3], collision.qvecIm()[refBInd + 3], 2), 2));
    }
  }

  template <typename TrackType>
  bool pidSel(const TrackType& track, uint8_t POI)
  {
    // Hadron
    if (POI == ese_parameters::kHadron)
      return true;
    // PID particles
    bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
    if (correctedTpcInnerParam < cfgPOverZPreselection->get(POI, 0u) || correctedTpcInnerParam > cfgPOverZPreselection->get(POI, 1u)) {
      return false;
    }
    if (std::abs(track.eta()) > cfgEtaPreselection->get(POI)) {
      return false;
    }
    if (track.tpcNClsFound() < cfgTPCNclsPreselection->get(POI, 0u) || track.tpcNClsFound() > cfgTPCNclsPreselection->get(POI, 1u)) {
      return false;
    }
    if (track.itsNCls() < cfgITSNclsPreselection->get(POI, 0u) || track.itsNCls() > cfgITSNclsPreselection->get(POI, 1u)) {
      return false;
    }
    if (track.tpcChi2NCl() < cfgTPCChi2Preselection->get(POI, 0u) || track.tpcChi2NCl() > cfgTPCChi2Preselection->get(POI, 1u)) {
      return false;
    }
    if (track.itsChi2NCl() < cfgITSChi2Preselection->get(POI, 0u) || track.itsChi2NCl() > cfgITSChi2Preselection->get(POI, 1u)) {
      return false;
    }
    if (std::abs(track.dcaXY()) > cfgDCAxyPreselection->get(POI)) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cfgDCAzPreselection->get(POI)) {
      return false;
    }
    float nSigmaTPC{0.f};
    switch (POI) {
      case ese_parameters::kProton: {
        if (cfgProtonPIDMode == ese_parameters::kRMSMode) { // RMS mode
          float nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDProton) ? std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()) : track.tpcNSigmaPr();
          if (nSigmaUse < cfgnSigmaCutRMSProton.value[0] || nSigmaUse > cfgnSigmaCutRMSProton.value[1]) {
            return false;
          }
        } else if (cfgProtonPIDMode == ese_parameters::kTPCMode) { // TPC mode
          nSigmaTPC = getNSigmaTPC(track, ese_parameters::kProton);
        } else if (cfgProtonPIDMode == ese_parameters::kTOFOnlyMode) { // TOF only mode
          if (!track.hasTOF())
            return false;
          if (track.tofNSigmaPr() < cfgnSigmaCutTOFProton.value[0] || track.tofNSigmaPr() > cfgnSigmaCutTOFProton.value[1]) {
            return false;
          }
        }
        break;
      }

      case ese_parameters::kDeuteron: {
        nSigmaTPC = getNSigmaTPC(track, ese_parameters::kDeuteron);
        break;
      }

      case ese_parameters::kTriton: {
        nSigmaTPC = getNSigmaTPC(track, ese_parameters::kTriton);
        break;
      }

      case ese_parameters::kHe3: {
        nSigmaTPC = getNSigmaTPC(track, ese_parameters::kHe3);
        break;
      }

      case ese_parameters::kAlpha: {
        nSigmaTPC = getNSigmaTPC(track, ese_parameters::kAlpha);
        break;
      }

      case ese_parameters::kPion: {
        if (cfgPionPIDMode == ese_parameters::kRMSMode) { // RMS mode
          float nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDPion) ? std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()) : track.tpcNSigmaPi();
          if (nSigmaUse < cfgnSigmaCutRMSPion.value[0] || nSigmaUse > cfgnSigmaCutRMSPion.value[1]) {
            return false;
          }
        } else if (cfgPionPIDMode == ese_parameters::kTPCMode) { // TPC mode
          nSigmaTPC = getNSigmaTPC(track, ese_parameters::kPion);
        } else if (cfgPionPIDMode == ese_parameters::kTOFOnlyMode) { // TOF only mode
          if (!track.hasTOF())
            return false;
          if (track.tofNSigmaPi() < cfgnSigmaCutTOFPion.value[0] || track.tofNSigmaPi() > cfgnSigmaCutTOFPion.value[1]) {
            return false;
          }
        }
        break;
      }

      case ese_parameters::kKaon: {
        if (cfgKaonPIDMode == ese_parameters::kRMSMode) { // RMS mode
          float nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDKaon) ? std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()) : track.tpcNSigmaKa();
          if (nSigmaUse < cfgnSigmaCutRMSKaon.value[0] || nSigmaUse > cfgnSigmaCutRMSKaon.value[1]) {
            return false;
          }
        } else if (cfgKaonPIDMode == ese_parameters::kTPCMode) { // TPC mode
          nSigmaTPC = getNSigmaTPC(track, ese_parameters::kKaon);
        } else if (cfgKaonPIDMode == ese_parameters::kTOFOnlyMode) { // TOF only mode
          if (!track.hasTOF())
            return false;
          if (track.tofNSigmaKa() < cfgnSigmaCutTOFKaon.value[0] || track.tofNSigmaKa() > cfgnSigmaCutTOFKaon.value[1]) {
            return false;
          }
        }
        break;
      }

      default:
        LOGF(error, "Unknown POI: %d", POI);
        return false;
    }
    if (nSigmaTPC < cfgnSigmaTPC->get(POI, 0u) || nSigmaTPC > cfgnSigmaTPC->get(POI, 1u)) {
      return false;
    }
    if (cfgOpenITSPreselection) {
      float nSigmaITS{getNSigmaITS(track, POI)};
      if (nSigmaITS < cfgnSigmaITS->get(POI, 0u) || nSigmaITS > cfgnSigmaITS->get(POI, 1u)) {
        return false;
      }
    }
    return true;
  }

  template <typename Tcoll, typename Ttrks>
  void fillESECandidates(Tcoll const& collision, Ttrks const& tracks, float& q2Tarx, float& q2Tary, int& multiTar, float& q2Refx, float& q2Refy, int& multiRef)
  {
    float psi2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], 2);
    float q2Tarinit{0.f};
    float q2Refinit{0.f};
    float q2FT0Cinit{0.f};
    for (const auto& track : tracks) { // loop on tracks
      if (!trackSelBasic(track)) {
        continue;
      }
      // we fill the track info QA in the main process
      fillTrackQAhist(track);
      bool kIsTar{false};
      bool kIsRef{false};
      kIsTar = pidSel(track, poiTar);
      kIsRef = pidSel(track, poiRef);
      if (kIsTar && kIsRef && !cfgOpenAllowCrossTrack && poiRef != ese_parameters::kHadron) {
        if (getNSigmaTPC(track, poiTar) < getNSigmaTPC(track, poiRef)) {
          kIsRef = false;
        } else {
          kIsTar = false;
        }
      }
      if (kIsTar) {
        multiTar++;
        q2Tarx += std::cos(2 * track.phi());
        q2Tary += std::sin(2 * track.phi());
        bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
        float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
        float nSigmaTPCTar{getNSigmaTPC(track, poiTar)};
        float nSigmaTOFTar{getNSigmaTOF(track, poiTar)};
        float nSigmaITSTar{getNSigmaITS(track, poiTar)};
        if (cfgOpenPIDQA) {
          ese_parameters::hPIDQATar1D[0]->Fill(ese_parameters::Charges[poiTar] * track.pt());
          ese_parameters::hPIDQATar1D[1]->Fill(track.eta());
          ese_parameters::hPIDQATar1D[2]->Fill(track.phi());
          ese_parameters::hPIDQATar1D[3]->Fill(track.itsNCls());
          ese_parameters::hPIDQATar1D[4]->Fill(track.tpcNClsFound());
          ese_parameters::hPIDQATar1D[5]->Fill(track.itsChi2NCl());
          ese_parameters::hPIDQATar1D[6]->Fill(track.tpcChi2NCl());
          ese_parameters::hPIDQATar1D[7]->Fill(track.dcaXY());
          ese_parameters::hPIDQATar1D[8]->Fill(track.dcaZ());
          ese_parameters::hPIDQATar1D[9]->Fill(nSigmaTPCTar);
          ese_parameters::hPIDQATar1D[10]->Fill(nSigmaTOFTar);
          ese_parameters::hPIDQATar1D[11]->Fill(nSigmaITSTar);
          ese_parameters::hPIDQATar2D[0]->Fill(track.sign() * correctedTpcInnerParam, track.tpcSignal());
          ese_parameters::hPIDQATar2D[1]->Fill(ese_parameters::Charges[poiTar] * track.pt(), nSigmaTPCTar);
          ese_parameters::hPIDQATar2D[2]->Fill(ese_parameters::Charges[poiTar] * track.pt(), nSigmaTOFTar);
          ese_parameters::hPIDQATar2D[3]->Fill(ese_parameters::Charges[poiTar] * track.pt(), nSigmaITSTar);
          if (cfgOpen3DPIDPlots->get(0u)) {
            ese_parameters::hPIDQATar3D[0]->Fill(nSigmaTOFTar, nSigmaITSTar, ese_parameters::Charges[poiTar] * track.pt());
          }
          if (cfgOpen3DPIDPlots->get(1u)) {
            ese_parameters::hPIDQATar3D[1]->Fill(nSigmaITSTar, nSigmaTPCTar, ese_parameters::Charges[poiTar] * track.pt());
          }
          if (cfgOpen3DPIDPlots->get(2u)) {
            ese_parameters::hPIDQATar3D[2]->Fill(nSigmaTOFTar, nSigmaTPCTar, ese_parameters::Charges[poiTar] * track.pt());
          }
        }
        ese_parameters::eseCandidates.emplace_back(ESECandidate{
          collision.posZ(), collision.centFT0C(), psi2, q2Tarinit, q2Refinit, q2FT0Cinit, static_cast<int8_t>(track.sign()), correctedTpcInnerParam, track.tpcSignal(), ese_parameters::Charges[poiTar] * track.pt(), track.eta(), track.phi(),
          track.dcaXY(), track.dcaZ(), static_cast<uint8_t>(track.tpcNClsFound()), track.itsNCls(), track.tpcChi2NCl(), track.itsChi2NCl(),
          nSigmaTPCTar, nSigmaTOFTar, nSigmaITSTar, track.itsClusterSizes()});
      }
      if (kIsRef) {
        multiRef++;
        q2Refx += std::cos(2 * track.phi());
        q2Refy += std::sin(2 * track.phi());
        if (cfgOpenPIDQA && poiRef != ese_parameters::kHadron) {
          bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
          float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
          float nSigmaTPCRef{getNSigmaTPC(track, poiRef)};
          float nSigmaTOFRef{getNSigmaTOF(track, poiRef)};
          float nSigmaITSRef{getNSigmaITS(track, poiRef)};
          ese_parameters::hPIDQARef1D[0]->Fill(ese_parameters::Charges[poiRef] * track.pt());
          ese_parameters::hPIDQARef1D[1]->Fill(track.eta());
          ese_parameters::hPIDQARef1D[2]->Fill(track.phi());
          ese_parameters::hPIDQARef1D[3]->Fill(track.itsNCls());
          ese_parameters::hPIDQARef1D[4]->Fill(track.tpcNClsFound());
          ese_parameters::hPIDQARef1D[5]->Fill(track.itsChi2NCl());
          ese_parameters::hPIDQARef1D[6]->Fill(track.tpcChi2NCl());
          ese_parameters::hPIDQARef1D[7]->Fill(track.dcaXY());
          ese_parameters::hPIDQARef1D[8]->Fill(track.dcaZ());
          ese_parameters::hPIDQARef1D[9]->Fill(nSigmaTPCRef);
          ese_parameters::hPIDQARef1D[10]->Fill(nSigmaTOFRef);
          ese_parameters::hPIDQARef1D[11]->Fill(nSigmaITSRef);
          ese_parameters::hPIDQARef2D[0]->Fill(track.sign() * correctedTpcInnerParam, track.tpcSignal());
          ese_parameters::hPIDQARef2D[1]->Fill(ese_parameters::Charges[poiRef] * track.pt(), nSigmaTPCRef);
          ese_parameters::hPIDQARef2D[2]->Fill(ese_parameters::Charges[poiRef] * track.pt(), nSigmaTOFRef);
          ese_parameters::hPIDQARef2D[3]->Fill(ese_parameters::Charges[poiRef] * track.pt(), nSigmaITSRef);
          if (cfgOpen3DPIDPlots->get(0u)) {
            ese_parameters::hPIDQARef3D[0]->Fill(nSigmaTOFRef, nSigmaITSRef, ese_parameters::Charges[poiRef] * track.pt());
          }
          if (cfgOpen3DPIDPlots->get(1u)) {
            ese_parameters::hPIDQARef3D[1]->Fill(nSigmaITSRef, nSigmaTPCRef, ese_parameters::Charges[poiRef] * track.pt());
          }
          if (cfgOpen3DPIDPlots->get(2u)) {
            ese_parameters::hPIDQARef3D[2]->Fill(nSigmaTOFRef, nSigmaTPCRef, ese_parameters::Charges[poiRef] * track.pt());
          }
        }
        if (cfgOpenv2Ref) {
          ese_parameters::eseReferences.emplace_back(ESEReference{static_cast<int8_t>(track.sign()), ese_parameters::Charges[poiRef] * track.pt(), std::cos(2 * (track.phi() - psi2))});
        }
      }
    }
  }

  void init(InitContext const&)
  {
    poiTar = getPOI(cfgTarName);
    if (poiTar == ese_parameters::kHadron) {
      LOGF(info, "This work flow is not designed to run over inclusive hadrons v2 ESE, you can only set the reference to be hadron, Reset Tar to He3");
      poiTar = ese_parameters::kHe3;
    }
    if (poiTar == poiRef) {
      LOGF(error, "The POI and the reference cannot be the same, please change the configuration, Reset Tar to He3 Ref to Proton");
      poiTar = ese_parameters::kHe3;
      poiRef = ese_parameters::kProton;
    }
    poiRef = getPOI(cfgRefName);
    detId = getDetId(cfgDetName);
    refAId = getDetId(cfgRefAName);
    refBId = getDetId(cfgRefBName);
    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }
    detInd = detId * 4 + cfgnTotalSystem * 4 * (2 - 2);
    refAInd = refAId * 4 + cfgnTotalSystem * 4 * (2 - 2);
    refBInd = refBId * 4 + cfgnTotalSystem * 4 * (2 - 2);

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

    AxisSpec axisITSNcls = {10, -1.5, 8.5, "ITSNcls"};
    AxisSpec axisTPCNcls = {160, 0, 160, "TPCNcls"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec axisPhi = {200, -2.1 * constants::math::PI, 2.1 * constants::math::PI};
    AxisSpec axisCentForQA = {100, 0, 100};
    AxisSpec axisCharge = {4, -2, 2, "Charge"};
    AxisSpec axisRigidity = {200, -10, 10, "#it{p}^{TPC}/#it{z}"};
    AxisSpec axisdEdx = {1400, 0, 1400, "dE/dx [arb. units]"};
    AxisSpec axisEta = {200, -1.0, 1.0, "#eta"};
    AxisSpec axisQvecF = {300, -1, 1, "Qvec"};
    AxisSpec axisNch = {4000, 0, 4000, "N_{ch}"};
    AxisSpec axisT0C = {70, 0, 70000, "N_{ch} (T0C)"};
    AxisSpec axisT0A = {70, 0, 70000, "N_{ch} (T0A)"};
    AxisSpec axisNchPV = {4000, 0, 4000, "N_{ch} (PV)"};
    AxisSpec axisCos = {102, -1.02, 1.02, "Cos"};
    // hists for event level QA
    histsESE.add("EventQA/histEventCount", ";Event Count;Counts", {HistType::kTH1F, {{11, 0, 11}}});
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(1, "after sel8");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(2, "kIsGoodZvtxFT0vsPV");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(3, "kNoSameBunchPileup");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(4, "kNoCollInTimeRangeStandard");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(5, "kIsGoodITSLayersAll");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(6, "kNoCollInRofStandard");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(7, "kNoHighMultCollInPrevRof");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(8, "occupancy");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(9, "MultCorrelationPVTracks");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(10, "MultCorrelationGlobalTracks");
    histsESE.get<TH1>(HIST("EventQA/histEventCount"))->GetXaxis()->SetBinLabel(11, "cfgEvSelV0AT0ACut");
    histsESE.add("EventQA/histVtz", ";#it{Vtz} (cm);Counts", {HistType::kTH1F, {{200, -20., +20.}}});
    histsESE.add("EventQA/histCent", ";Centrality (%);Counts", {HistType::kTH1F, {{100, 0., 100.}}});
    if (cfgOpenFullEventQA) {
      histsESE.add("EventQA/hist_globalTracks_centT0C_before", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      histsESE.add("EventQA/hist_PVTracks_centT0C_before", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
      histsESE.add("EventQA/hist_globalTracks_PVTracks_before", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
      histsESE.add("EventQA/hist_globalTracks_multT0A_before", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      histsESE.add("EventQA/hist_globalTracks_multV0A_before", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      histsESE.add("EventQA/hist_multV0A_multT0A_before", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      histsESE.add("EventQA/hist_multT0C_centT0C_before", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
      histsESE.add("EventQA/hist_globalTracks_centT0C_after", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, axisNch}});
      histsESE.add("EventQA/hist_PVTracks_centT0C_after", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, axisNchPV}});
      histsESE.add("EventQA/hist_globalTracks_PVTracks_after", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {axisNchPV, axisNch}});
      histsESE.add("EventQA/hist_globalTracks_multT0A_after", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      histsESE.add("EventQA/hist_globalTracks_multV0A_after", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {axisT0A, axisNch}});
      histsESE.add("EventQA/hist_multV0A_multT0A_after", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {axisT0A, axisT0A}});
      histsESE.add("EventQA/hist_multT0C_centT0C_after", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, axisT0C}});
    }
    histsESE.add("PlanQA/histQvec_CorrL0_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {axisQvecF, axisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL1_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {axisQvecF, axisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL2_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {axisQvecF, axisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL3_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {axisQvecF, axisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL0_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL1_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL2_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL3_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_SigRefAV2", ";Centrality;Cos(Sig-RefA)", {HistType::kTProfile, {cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_SigRefBV2", ";Centrality;Cos(Sig-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_RefARefBV2", ";Centrality;Cos(RefA-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    // hists for track level QA
    histsESE.add("TrackQA/hist_dEdxTPC_All", ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x} [arb. units]", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histsESE.add("TrackQA/hist_pt_All", ";#it{p}_{T};counts", {HistType::kTH1F, {cfgaxispt}});
    histsESE.add("TrackQA/hist_eta_All", ";#it{#eta};counts", {HistType::kTH1F, {axisEta}});
    histsESE.add("TrackQA/hist_phi_All", ";#it{#phi};counts", {HistType::kTH1F, {axisPhi}});
    histsESE.add("TrackQA/hist_ITSNcls_All", ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
    histsESE.add("TrackQA/hist_TPCNcls_All", ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
    histsESE.add("TrackQA/hist_ITSChi2NDF_All", ";ITS#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
    histsESE.add("TrackQA/hist_TPCChi2NDF_All", ";TPC#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
    histsESE.add("TrackQA/hist_DCAxy_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAxy}});
    histsESE.add("TrackQA/hist_DCAz_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAz}});
    histsESE.add("TrackQA/hist_MomRes_All", ";#it{p};Res(1 - corrted_p/track_p)", {HistType::kTH2F, {{100, 0, 10}, {100, -1, 1}}});
    histsESE.add("TrackQA/hist_He3AlphaTrackcounts_All", ";Track counts;Counts", {HistType::kTH1F, {{2, 0, 2}}});
    histsESE.get<TH1>(HIST("TrackQA/hist_He3AlphaTrackcounts_All"))->GetXaxis()->SetBinLabel(1, "All Tracks");
    histsESE.get<TH1>(HIST("TrackQA/hist_He3AlphaTrackcounts_All"))->GetXaxis()->SetBinLabel(2, "He3+Alpha");
    // v2 and ESEPlots
    /// QA plots
    if (cfgOpenPIDQA) {
      ese_parameters::hPIDQATar2D[0] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_dEdxTPC_%s", cfgTarName.value.c_str()), ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x} [arb. units]", HistType::kTH2F, {axisRigidity, axisdEdx});
      ese_parameters::hPIDQATar1D[0] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_pt_%s", cfgTarName.value.c_str()), ";#it{p}_{T};counts", HistType::kTH1F, {cfgaxispt});
      ese_parameters::hPIDQATar1D[1] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_eta_%s", cfgTarName.value.c_str()), ";#it{#eta};counts", HistType::kTH1F, {axisEta});
      ese_parameters::hPIDQATar1D[2] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_phi_%s", cfgTarName.value.c_str()), ";#it{#phi};counts", HistType::kTH1F, {axisPhi});
      ese_parameters::hPIDQATar1D[3] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_ITSNcls_%s", cfgTarName.value.c_str()), ";ITSNcls;counts", HistType::kTH1F, {axisITSNcls});
      ese_parameters::hPIDQATar1D[4] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_TPCNcls_%s", cfgTarName.value.c_str()), ";TPCNcls;counts", HistType::kTH1F, {axisTPCNcls});
      ese_parameters::hPIDQATar1D[5] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_ITSChi2NDF_%s", cfgTarName.value.c_str()), ";ITS#it{#chi^{2}}/NDF;counts", HistType::kTH1F, {cfgaxisChi2Ncls});
      ese_parameters::hPIDQATar1D[6] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_TPCChi2NDF_%s", cfgTarName.value.c_str()), ";TPC#it{#chi^{2}}/NDF;counts", HistType::kTH1F, {cfgaxisChi2Ncls});
      ese_parameters::hPIDQATar1D[7] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_DCAxy_%s", cfgTarName.value.c_str()), ";#it{DCA_{xy}};counts", HistType::kTH1F, {cfgaxisDCAxy});
      ese_parameters::hPIDQATar1D[8] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_DCAz_%s", cfgTarName.value.c_str()), ";#it{DCA_{xy}};counts", HistType::kTH1F, {cfgaxisDCAz});
      ese_parameters::hPIDQATar1D[9] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaTPC_%s", cfgTarName.value.c_str()), ";n#sigmaTPC;counts", HistType::kTH1F, {cfgnSigmaBinsTPC});
      ese_parameters::hPIDQATar2D[1] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaTPC_pt_%s", cfgTarName.value.c_str()), ";#it{p}_{T};n#sigmaTPC", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsTPC});
      ese_parameters::hPIDQATar1D[10] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaTOF_%s", cfgTarName.value.c_str()), ";n#sigmaTOF;counts", HistType::kTH1F, {cfgnSigmaBinsTOF});
      ese_parameters::hPIDQATar2D[2] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaTOF_pt_%s", cfgTarName.value.c_str()), ";#it{p}_{T};n#sigmaTOF", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsTOF});
      ese_parameters::hPIDQATar1D[11] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaITS_%s", cfgTarName.value.c_str()), ";n#sigmaITS;counts", HistType::kTH1F, {cfgnSigmaBinsITS});
      ese_parameters::hPIDQATar2D[3] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaITS_pt_%s", cfgTarName.value.c_str()), ";#it{p}_{T};n#sigmaITS", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsITS});
      if (cfgRefName.value != "kHadron") {
        ese_parameters::hPIDQARef2D[0] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_dEdxTPC_%s", cfgRefName.value.c_str()), ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x} [arb. units]", HistType::kTH2F, {axisRigidity, axisdEdx});
        ese_parameters::hPIDQARef1D[0] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_pt_%s", cfgRefName.value.c_str()), ";#it{p}_{T};counts", HistType::kTH1F, {cfgaxispt});
        ese_parameters::hPIDQARef1D[1] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_eta_%s", cfgRefName.value.c_str()), ";#it{#eta};counts", HistType::kTH1F, {axisEta});
        ese_parameters::hPIDQARef1D[2] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_phi_%s", cfgRefName.value.c_str()), ";#it{#phi};counts", HistType::kTH1F, {axisPhi});
        ese_parameters::hPIDQARef1D[3] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_ITSNcls_%s", cfgRefName.value.c_str()), ";ITSNcls;counts", HistType::kTH1F, {axisITSNcls});
        ese_parameters::hPIDQARef1D[4] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_TPCNcls_%s", cfgRefName.value.c_str()), ";TPCNcls;counts", HistType::kTH1F, {axisTPCNcls});
        ese_parameters::hPIDQARef1D[5] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_ITSChi2NDF_%s", cfgRefName.value.c_str()), ";ITS#it{#chi^{2}}/NDF;counts", HistType::kTH1F, {cfgaxisChi2Ncls});
        ese_parameters::hPIDQARef1D[6] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_TPCChi2NDF_%s", cfgRefName.value.c_str()), ";TPC#it{#chi^{2}}/NDF;counts", HistType::kTH1F, {cfgaxisChi2Ncls});
        ese_parameters::hPIDQARef1D[7] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_DCAxy_%s", cfgRefName.value.c_str()), ";#it{DCA_{xy}};counts", HistType::kTH1F, {cfgaxisDCAxy});
        ese_parameters::hPIDQARef1D[8] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_DCAz_%s", cfgRefName.value.c_str()), ";#it{DCA_{xy}};counts", HistType::kTH1F, {cfgaxisDCAz});
        ese_parameters::hPIDQARef1D[9] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaTPC_%s", cfgRefName.value.c_str()), ";n#sigmaTPC;counts", HistType::kTH1F, {cfgnSigmaBinsTPC});
        ese_parameters::hPIDQARef2D[1] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaTPC_pt_%s", cfgRefName.value.c_str()), ";#it{p}_{T};n#sigmaTPC", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsTPC});
        ese_parameters::hPIDQARef1D[10] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaTOF_%s", cfgRefName.value.c_str()), ";n#sigmaTOF;counts", HistType::kTH1F, {cfgnSigmaBinsTOF});
        ese_parameters::hPIDQARef2D[2] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaTOF_pt_%s", cfgRefName.value.c_str()), ";#it{p}_{T};n#sigmaTOF", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsTOF});
        ese_parameters::hPIDQARef1D[11] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_nSigmaITS_%s", cfgRefName.value.c_str()), ";n#sigmaITS;counts", HistType::kTH1F, {cfgnSigmaBinsITS});
        ese_parameters::hPIDQARef2D[3] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_nSigmaITS_pt_%s", cfgRefName.value.c_str()), ";#it{p}_{T};n#sigmaITS", HistType::kTH2F, {cfgaxispt, cfgnSigmaBinsITS});
      }
      if (cfgOpen3DPIDPlots->get(0u)) {
        ese_parameters::hPIDQATar3D[0] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFITSPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxispt});
        if (cfgRefName.value != "kHadron") {
          ese_parameters::hPIDQARef3D[0] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFITSPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxispt});
        }
      }
      if (cfgOpen3DPIDPlots->get(1u)) {
        ese_parameters::hPIDQATar3D[1] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaITSTPCPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxispt});
        if (cfgRefName.value != "kHadron") {
          ese_parameters::hPIDQARef3D[1] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaITSTPCPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxispt});
        }
      }
      if (cfgOpen3DPIDPlots->get(2u)) {
        ese_parameters::hPIDQATar3D[2] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFTPCPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxispt});
        if (cfgRefName.value != "kHadron") {
          ese_parameters::hPIDQARef3D[2] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFTPCPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxispt});
        }
      }
    }
    // v2 plots
    if (cfgOpenv2Tar) {
      ese_parameters::hv2Tar[0] = histsESE.add<TProfile3D>(Form("ESE/V2/hist_%sPosV2", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2}", {HistType::kTProfile3D, {cfgaxispt, cfgaxisCent, cfgaxisq2Tar}});
      ese_parameters::hv2Tar[1] = histsESE.add<TProfile3D>(Form("ESE/V2/hist_%sNegV2", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2}", {HistType::kTProfile3D, {cfgaxispt, cfgaxisCent, cfgaxisq2Tar}});
    }
    if (cfgOpenv2Ref) {
      ese_parameters::hv2Ref[0] = histsESE.add<TProfile3D>(Form("ESE/V2/hist_%sPosV2", cfgRefName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2}", {HistType::kTProfile3D, {cfgaxispt, cfgaxisCent, cfgaxisq2Ref}});
      ese_parameters::hv2Ref[1] = histsESE.add<TProfile3D>(Form("ESE/V2/hist_%sNegV2", cfgRefName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2}", {HistType::kTProfile3D, {cfgaxispt, cfgaxisCent, cfgaxisq2Ref}});
    }
    // ESE plots
    if (cfgOpenESEQA) {
      ese_parameters::hESEQATar1D[0] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum", cfgTarName.value.c_str()), ";Num_{Proton}/Event;counts", HistType::kTH1F, {{100, 0, 100}});
      ese_parameters::hESEQATar1D[1] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2", cfgTarName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2Tar});
      ese_parameters::hESEQATar1D[2] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Numerator", cfgTarName.value.c_str()), ";#it{q}_{2} numerator;counts", HistType::kTH1F, {cfgq2NumeratorTar});
      ese_parameters::hESEQATar1D[3] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Denominator", cfgTarName.value.c_str()), ";#it{q}_{2} denominator;counts", HistType::kTH1F, {cfgq2DenominatorTar});
      ese_parameters::hESEQATar1D[4] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum_Survived", cfgTarName.value.c_str()), ";Num_{Proton}/Event;counts", HistType::kTH1F, {{100, 0, 100}});
      ese_parameters::hESEQATar1D[5] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2_Survived", cfgTarName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2Tar});
      ese_parameters::hESEQATar1D[6] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Numerator_Survived", cfgTarName.value.c_str()), ";#it{q}_{2} numerator;counts", HistType::kTH1F, {cfgq2NumeratorTar});
      ese_parameters::hESEQATar1D[7] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Denominator_Survived", cfgTarName.value.c_str()), ";#it{q}_{2} denominator;counts", HistType::kTH1F, {cfgq2DenominatorTar});
      ese_parameters::hESEQATar2D[0] = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent", cfgTarName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2Tar, cfgaxisCent});
      ese_parameters::hESEQATar2D[1] = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent_Survived", cfgTarName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2Tar, cfgaxisCent});
      ese_parameters::hESEQARef1D[0] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum", cfgRefName.value.c_str()), ";Num_{He3}/Event;counts", HistType::kTH1F, {{10, 0, 10}});
      ese_parameters::hESEQARef1D[1] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2", cfgRefName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2Ref});
      ese_parameters::hESEQARef1D[2] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Numerator", cfgRefName.value.c_str()), ";#it{q}_{2} numerator;counts", HistType::kTH1F, {cfgq2NumeratorRef});
      ese_parameters::hESEQARef1D[3] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Denominator", cfgRefName.value.c_str()), ";#it{q}_{2} denominator;counts", HistType::kTH1F, {cfgq2DenominatorRef});
      ese_parameters::hESEQARef1D[4] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum_Survived", cfgRefName.value.c_str()), ";Num_{He3}/Event;counts", HistType::kTH1F, {{10, 0, 10}});
      ese_parameters::hESEQARef1D[5] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2_Survived", cfgRefName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2Ref});
      ese_parameters::hESEQARef1D[6] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Numerator_Survived", cfgRefName.value.c_str()), ";#it{q}_{2} numerator;counts", HistType::kTH1F, {cfgq2NumeratorRef});
      ese_parameters::hESEQARef1D[7] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2Denominator_Survived", cfgRefName.value.c_str()), ";#it{q}_{2} denominator;counts", HistType::kTH1F, {cfgq2DenominatorRef});
      ese_parameters::hESEQARef2D[0] = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent", cfgRefName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2Ref, cfgaxisCent});
      ese_parameters::hESEQARef2D[1] = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent_Survived", cfgRefName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2Ref, cfgaxisCent});
      ese_parameters::hESEQAFT0C1D[0] = histsESE.add<TH1>("ESE/ESEQA/hist_q2FT0C", ";#it{q}_{2} (FT0C)", HistType::kTH1F, {cfgaxisq2FT0C});
      ese_parameters::hESEQAFT0C1D[1] = histsESE.add<TH1>("ESE/ESEQA/hist_q2FT0C_Survived", ";#it{q}_{2} (FT0C)", HistType::kTH1F, {cfgaxisq2FT0C});
      ese_parameters::hESEQAFT0C2D[0] = histsESE.add<TH2>("ESE/ESEQA/hist_q2FT0C_Cent", ";#it{q}_{2} (FT0C);Centrality (%)", HistType::kTH2F, {cfgaxisq2FT0C, cfgaxisCent});
      ese_parameters::hESEQAFT0C2D[1] = histsESE.add<TH2>("ESE/ESEQA/hist_q2FT0C_Cent_Survived", ";#it{q}_{2} (FT0C);Centrality (%)", HistType::kTH2F, {cfgaxisq2FT0C, cfgaxisCent});
    }
    if (cfgOpenESE) {
      ese_parameters::hESETar = histsESE.add<THnSparse>(Form("ESE/ESE/histESE_%s", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2};cos(#phi-#Psi_{2});Charge", HistType::kTHnSparseF, {cfgaxispt, cfgaxisCent, cfgaxisq2Ref, axisCos, axisCharge});
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::Qvectors, aod::QvectorFT0CVecs>>::iterator const& collision, TracksPIDFull const& tracks)
  {
    ese_parameters::eseCandidates.clear();
    ese_parameters::eseReferences.clear();
    const float centrality{collision.centFT0C()};
    const int64_t multTrk{tracks.size()};
    if (cfgOpenFullEventQA)
      fillEventQAhistBe(collision, multTrk, centrality);
    if (!eventSelBasic(collision, multTrk, centrality, true))
      return;
    if (cfgOpenFullEventQA)
      fillEventQAhistAf(collision, multTrk, centrality);
    float q2FT0C{std::sqrt(collision.qvecFT0CReVec()[0] * collision.qvecFT0CReVec()[0] + collision.qvecFT0CImVec()[0] * collision.qvecFT0CImVec()[0]) * std::sqrt(collision.sumAmplFT0C())};
    float q2Tarx{0.};
    float q2Tary{0.};
    float q2Refx{0.};
    float q2Refy{0.};
    int multiTar{0};
    int multiRef{0};
    fillHistosQvec(collision);
    fillESECandidates(collision, tracks, q2Tarx, q2Tary, multiTar, q2Refx, q2Refy, multiRef);
    float q2Tar{calculateq2(q2Tarx, q2Tary, multiTar)};
    float q2Ref{calculateq2(q2Refx, q2Refy, multiRef)};
    if (cfgOpenESEQA) {
      ese_parameters::hESEQATar1D[0]->Fill(multiTar);
      ese_parameters::hESEQATar1D[1]->Fill(q2Tar);
      ese_parameters::hESEQATar1D[2]->Fill(std::hypot(q2Tarx, q2Tary));
      ese_parameters::hESEQATar1D[3]->Fill(std::sqrt(static_cast<float>(multiTar)));
      ese_parameters::hESEQATar2D[0]->Fill(q2Tar, centrality);
      ese_parameters::hESEQARef1D[0]->Fill(multiRef);
      ese_parameters::hESEQARef1D[1]->Fill(q2Ref);
      ese_parameters::hESEQARef1D[2]->Fill(std::hypot(q2Refx, q2Refy));
      ese_parameters::hESEQARef1D[3]->Fill(std::sqrt(static_cast<float>(multiRef)));
      ese_parameters::hESEQARef2D[0]->Fill(q2Ref, centrality);
      ese_parameters::hESEQAFT0C1D[0]->Fill(q2FT0C);
      ese_parameters::hESEQAFT0C2D[0]->Fill(q2FT0C, centrality);
    }
    if (cfgOpenv2Ref) {
      for (const auto& c : ese_parameters::eseReferences) {
        if (c.signRef > 0) {
          ese_parameters::hv2Ref[0]->Fill(c.ptRef, centrality, q2Ref, c.v2Ref);
        } else {
          ese_parameters::hv2Ref[1]->Fill(c.ptRef, centrality, q2Ref, c.v2Ref);
        }
      }
    }
    if (multiTar == 0)
      return;
    if (cfgOpenESEQA) {
      ese_parameters::hESEQATar1D[4]->Fill(multiTar);
      ese_parameters::hESEQATar1D[5]->Fill(q2Tar);
      ese_parameters::hESEQATar1D[6]->Fill(std::hypot(q2Tarx, q2Tary));
      ese_parameters::hESEQATar1D[7]->Fill(std::sqrt(static_cast<float>(multiTar)));
      ese_parameters::hESEQATar2D[1]->Fill(q2Tar, centrality);
      ese_parameters::hESEQARef1D[4]->Fill(multiRef);
      ese_parameters::hESEQARef1D[5]->Fill(q2Ref);
      ese_parameters::hESEQARef1D[6]->Fill(std::hypot(q2Refx, q2Refy));
      ese_parameters::hESEQARef1D[7]->Fill(std::sqrt(static_cast<float>(multiRef)));
      ese_parameters::hESEQARef2D[1]->Fill(q2Ref, centrality);
      ese_parameters::hESEQAFT0C1D[1]->Fill(q2FT0C);
      ese_parameters::hESEQAFT0C2D[1]->Fill(q2FT0C, centrality);
    }
    for (const auto& c : ese_parameters::eseCandidates) {
      eseTable(c.vtz, c.centFT0C, c.psi2FT0C, q2Tar, q2Ref, q2FT0C, c.signTar, c.tpcInnerParamTar, c.tpcSignalTar, c.ptTar, c.etaTar, c.phiTar, c.dcaXYTar, c.dcaZTar, c.tpcNclsTar, c.itsNclsTar, c.tpcChi2NDFTar, c.itsChi2NDFTar, c.tpcNSigmaTar, c.tofNSigmaTar, c.itsNSigmaTar, c.itsClusSizeTar);
      if (cfgOpenESE) {
        ese_parameters::hESETar->Fill(c.ptTar, c.centFT0C, q2Ref, std::cos(2 * (c.phiTar - c.psi2FT0C)), c.signTar);
      }
      if (cfgOpenv2Tar) {
        if (c.signTar > 0) {
          ese_parameters::hv2Tar[0]->Fill(c.ptTar, c.centFT0C, q2Tar, std::cos(2 * (c.phiTar - c.psi2FT0C)));
        } else {
          ese_parameters::hv2Tar[1]->Fill(c.ptTar, c.centFT0C, q2Tar, std::cos(2 * (c.phiTar - c.psi2FT0C)));
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEsePHe3>(cfgc),
  };
}
