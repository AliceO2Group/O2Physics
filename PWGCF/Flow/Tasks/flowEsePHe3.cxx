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

#include "DataFormatsTPC/BetheBlochAleph.h"

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

namespace ese_parameters
{
constexpr uint8_t kProton = 0;
constexpr uint8_t kDeuteron = 1;
constexpr uint8_t kTriton = 2;
constexpr uint8_t kHe3 = 3;
constexpr uint8_t kAlpha = 4;
constexpr int kFT0AV0ASigma = 5;
constexpr int kRMSMode = 0;
constexpr int kTPCMode = 1;
constexpr int kTOFOnlyMode = 2;
constexpr float Amplitudelow = 1e-8;
constexpr float Charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};
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
constexpr double TPCnSigmaCutDefault[5][2]{
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.}};
constexpr double ITSnSigmaCutDefault[5][2]{
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.},
  {-3., 3.}};
constexpr double PtPreselection[5][2]{
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.},
  {0.15, 99.}};
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
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
static const std::vector<std::string> pidPtNames{"p_{T} Low", "p_{T} High"};
std::vector<ESECandidate> eseCandidates;
// Tar ptr
std::shared_ptr<TH1> hPIDQATar1D[12];
std::shared_ptr<TH2> hPIDQATar2D[4];
std::shared_ptr<TH3> hPIDQATar3D[3];
std::shared_ptr<TProfile2D> hv2Tar[2];
std::shared_ptr<TH1> hESEQATar1D[2];
std::shared_ptr<TH2> hESEQATar2D;
std::shared_ptr<THnSparse> hESETar;
// Ref ptr
std::shared_ptr<TH1> hPIDQARef1D[12];
std::shared_ptr<TH2> hPIDQARef2D[4];
std::shared_ptr<TH3> hPIDQARef3D[3];
std::shared_ptr<TProfile2D> hv2Ref[2];
std::shared_ptr<TH1> hESEQARef1D[2];
std::shared_ptr<TH2> hESEQARef2D;
} // namespace ese_parameters

using TracksPIDFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>;

struct FlowEsePHe3 {

  EventPlaneHelper helperEP;
  o2::aod::ITSResponse itsResponse;
  HistogramRegistry histsESE{"histsESE", {}, OutputObjHandlingPolicy::AnalysisObject};
  // process POI control
  Configurable<std::string> cfgTarName{"cfgTarName", "kHe3", "Name of the v2 particle: kProton, kDeuteron, kTriton, kHe3, kAlpha"};
  Configurable<std::string> cfgRefName{"cfgRefName", "kProton", "Name of the q2 reference particle: kProton, kDeuteron, kTriton, kHe3, kAlpha"};
  // total control config
  Configurable<bool> cfgOpenAllowCrossTrack{"cfgOpenAllowCrossTrack", false, "Allow one track to be identified as different kind of PID particles"};
  Configurable<bool> cfgOpenFullEventQA{"cfgOpenFullEventQA", true, "Open full QA plots for event QA"};
  Configurable<bool> cfgOpenPIDQA{"cfgOpenPIDQA", true, "Open PID QA plots"};
  Configurable<bool> cfgOpenv2{"cfgOpenv2", true, "Open v2(EP)and q calculation for Proton and He3"};
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
  // PID configs
  Configurable<LabeledArray<double>> cfgPtPreselection{"cfgPtPreselection", {ese_parameters::PtPreselection[0], 5, 2, ese_parameters::names, ese_parameters::pidPtNames}, "Pt preselection for light nuclei"};
  Configurable<std::vector<float>> cfgnSigmaCutTOFProton{"cfgnSigmaCutTOFProton", {-1.5, 1.5}, "TOF nsigma cut limit for Proton"};
  Configurable<std::vector<float>> cfgnSigmaCutRMSProton{"cfgnSigmaCutRMSProton", {-3, 3}, "RMS nsigma cut limit for Proton"};
  Configurable<bool> cfgUseSelfnSigmaTPCProton{"cfgUseSelfnSigmaTPCProton", false, "Use self nSigma TPC for Proton PID"};
  Configurable<int> cfgProtonPIDMode{"cfgProtonPIDMode", 2, "Proton PID mode: 0 for TPC + RMS(TPC,TOF), 1 for TPC only, 2 for TOF only"};
  Configurable<LabeledArray<double>> cfgnSigmaTPC{"cfgnSigmaTPC", {ese_parameters::TPCnSigmaCutDefault[0], 5, 2, ese_parameters::names, ese_parameters::pidTPCnSigmaNames}, "TPC nSigma selection for light nuclei"};
  Configurable<LabeledArray<double>> cfgnSigmaITS{"cfgnSigmaITS", {ese_parameters::ITSnSigmaCutDefault[0], 5, 2, ese_parameters::names, ese_parameters::pidITSnSigmaNames}, "ITS nSigma selection for light nuclei"};
  // PID BBself paras config
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {ese_parameters::BetheBlochDefault[0], 5, 6, ese_parameters::names, ese_parameters::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {ese_parameters::BbMomScalingDefault[0], 5, 2, ese_parameters::names, ese_parameters::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", true, "If true, divide tpcInnerParam by the electric charge"};
  // Axias configs
  ConfigurableAxis cfgrigidityBins{"cfgrigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis cfgdedxBins{"cfgdedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis cfgnSigmaBinsTPC{"cfgnSigmaBinsTPC", {200, -5.f, 5.f}, "Binning for n sigma TPC"};
  ConfigurableAxis cfgnSigmaBinsTOF{"cfgnSigmaBinsTOF", {200, -5.f, 5.f}, "Binning for n sigma TOF"};
  ConfigurableAxis cfgnSigmaBinsITS{"cfgnSigmaBinsITS", {200, -5.f, 5.f}, "Binning for n sigma ITS"};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, "Binning for P_{t}"};
  ConfigurableAxis cfgaxisetaPID{"cfgaxisetaPID", {90, -0.9, 0.9}, "Binning for Pt QA"};
  ConfigurableAxis cfgaxisDCAz{"cfgaxisDCAz", {200, -1, 1}, "Binning for DCAz"};
  ConfigurableAxis cfgaxisDCAxy{"cfgaxisDCAxy", {100, -0.5, 0.5}, "Binning for DCAxy"};
  ConfigurableAxis cfgaxisChi2Ncls{"cfgaxisChi2Ncls", {100, 0, 30}, "Binning for Chi2Ncls TPC/ITS"};
  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};
  ConfigurableAxis cfgaxisNch{"cfgaxisNch", {4000, 0, 4000}, "N_{ch}"};
  ConfigurableAxis cfgaxisT0C{"cfgaxisT0C", {70, 0, 70000}, "N_{ch} (T0C)"};
  ConfigurableAxis cfgaxisT0A{"cfgaxisT0A", {200, 0, 200000}, "N_{ch} (T0A)"};
  ConfigurableAxis cfgaxisNchPV{"cfgaxisNchPV", {4000, 0, 4000}, "N_{ch} (PV)"};
  ConfigurableAxis cfgaxisq2{"cfgaxisq2", {120, 0, 12}, "Binning for P_{t} PID"};
  ConfigurableAxis cfgaxiscos{"cfgaxiscos", {102, -1.02, 1.02}, ""};

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
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(0u, 0u), cfgBetheBlochParams->get(0u, 1u), cfgBetheBlochParams->get(0u, 2u), cfgBetheBlochParams->get(0u, 3u), cfgBetheBlochParams->get(0u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(0u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kDeuteron: {
        const double bgScaling[2]{ese_parameters::Charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / ese_parameters::Masses[1], ese_parameters::Charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / ese_parameters::Masses[1]};
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(1u, 0u), cfgBetheBlochParams->get(1u, 1u), cfgBetheBlochParams->get(1u, 2u), cfgBetheBlochParams->get(1u, 3u), cfgBetheBlochParams->get(1u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(1u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kTriton: {
        const double bgScaling[2]{ese_parameters::Charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / ese_parameters::Masses[2], ese_parameters::Charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / ese_parameters::Masses[2]};
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(2u, 0u), cfgBetheBlochParams->get(2u, 1u), cfgBetheBlochParams->get(2u, 2u), cfgBetheBlochParams->get(2u, 3u), cfgBetheBlochParams->get(2u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(2u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kHe3: {
        const double bgScaling[2]{ese_parameters::Charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / ese_parameters::Masses[3], ese_parameters::Charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / ese_parameters::Masses[3]};
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(3u, 0u), cfgBetheBlochParams->get(3u, 1u), cfgBetheBlochParams->get(3u, 2u), cfgBetheBlochParams->get(3u, 3u), cfgBetheBlochParams->get(3u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(3u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      case ese_parameters::kAlpha: {
        const double bgScaling[2]{ese_parameters::Charges[4] * cfgMomentumScalingBetheBloch->get(4u, 0u) / ese_parameters::Masses[4], ese_parameters::Charges[4] * cfgMomentumScalingBetheBloch->get(4u, 1u) / ese_parameters::Masses[4]};
        double expBethe{tpc::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScaling[iC]), cfgBetheBlochParams->get(4u, 0u), cfgBetheBlochParams->get(4u, 1u), cfgBetheBlochParams->get(4u, 2u), cfgBetheBlochParams->get(4u, 3u), cfgBetheBlochParams->get(4u, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(4u, 5u)};
        double nSigmaTPC{static_cast<float>((track.tpcSignal() - expBethe) / expSigma)};
        return nSigmaTPC;
      }
      default:
        return -99.f;
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
    } else {
      LOGF(warning, "Unknown POI: %s", POI.value.c_str());
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
  }

  template <typename CollType>
  void fillHistosQvec(const CollType& collision)
  {
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (2 - 2);
    int refAInd = refAId * 4 + cfgnTotalSystem * 4 * (2 - 2);
    int refBInd = refBId * 4 + cfgnTotalSystem * 4 * (2 - 2);
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
    if (track.pt() < cfgPtPreselection->get(POI, 0u) || track.pt() > cfgPtPreselection->get(POI, 1u)) {
      return false;
    }
    float nSigmaTPC = 0.f;
    float nSigmaITS = 0.f;
    switch (POI) {
      case ese_parameters::kProton:
        if (cfgProtonPIDMode == ese_parameters::kRMSMode) { // RMS mode
          float nSigmaUse = (track.pt() > cfgPtMaxforTPCOnlyPIDPrton) ? std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr()) : track.tpcNSigmaPr();
          if (nSigmaUse < cfgnSigmaCutRMSProton.value[0] || nSigmaUse > cfgnSigmaCutRMSProton.value[1]) {
            return false;
          }
        } else if (cfgProtonPIDMode == ese_parameters::kTPCMode) { // TPC mode
          nSigmaTPC = (cfgUseSelfnSigmaTPCProton ? getNSigmaTPCSelfBB(track, ese_parameters::kProton) : track.tpcNSigmaPr());
        } else if (cfgProtonPIDMode == ese_parameters::kTOFOnlyMode) { // TOF only mode
          if (!track.hasTOF())
            return false;
          if (track.tofNSigmaPr() < cfgnSigmaCutTOFProton.value[0] || track.tofNSigmaPr() > cfgnSigmaCutTOFProton.value[1]) {
            return false;
          }
        }
        nSigmaITS = itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
        break;

      case ese_parameters::kDeuteron:
        nSigmaTPC = getNSigmaTPCSelfBB(track, ese_parameters::kDeuteron);
        nSigmaITS = itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track);
        break;

      case ese_parameters::kTriton:
        nSigmaTPC = getNSigmaTPCSelfBB(track, ese_parameters::kTriton);
        nSigmaITS = itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
        break;

      case ese_parameters::kHe3:
        nSigmaTPC = getNSigmaTPCSelfBB(track, ese_parameters::kHe3);
        nSigmaITS = itsResponse.nSigmaITS<o2::track::PID::Helium3>(track);
        break;

      case ese_parameters::kAlpha:
        nSigmaTPC = getNSigmaTPCSelfBB(track, ese_parameters::kAlpha);
        nSigmaITS = itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
        break;

      default:
        LOGF(error, "Unknown POI: %d", POI);
        return false;
    }
    if (nSigmaTPC < cfgnSigmaTPC->get(POI, 0u) || nSigmaTPC > cfgnSigmaTPC->get(POI, 1u)) {
      return false;
    }
    if (nSigmaITS < cfgnSigmaITS->get(POI, 0u) || nSigmaITS > cfgnSigmaITS->get(POI, 1u)) {
      return false;
    }
    return true;
  }

  template <typename Tcoll, typename Ttrks>
  void fillESECandidates(Tcoll const& collision, Ttrks const& tracks, float& q2Tarx, float& q2Tary, int& multiTar, float& q2Refx, float& q2Refy, int& multiRef)
  {
    float psi2 = helperEP.GetEventPlane(collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3], 2);
    float q2Tarinit{0.f};
    float q2Refinit{0.f};
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
      if (!cfgOpenAllowCrossTrack && kIsTar && kIsRef) {
        if (getNSigmaTPCSelfBB(track, poiTar) < getNSigmaTPCSelfBB(track, poiRef)) {
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
        float nSigmaTPCTar{(poiTar == ese_parameters::kProton && !cfgUseSelfnSigmaTPCProton) ? track.tpcNSigmaPr() : getNSigmaTPCSelfBB(track, poiTar)};
        float nSigmaTOFTar{getNSigmaTOF(track, poiTar)};
        float nSigmaITSTar{getNSigmaITS(track, poiTar)};
        if (cfgOpenPIDQA) {
          ese_parameters::hPIDQATar1D[0]->Fill(track.pt());
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
          ese_parameters::hPIDQATar2D[1]->Fill(nSigmaTPCTar, track.pt());
          ese_parameters::hPIDQATar2D[2]->Fill(nSigmaTOFTar, track.pt());
          ese_parameters::hPIDQATar2D[3]->Fill(nSigmaITSTar, track.pt());
          if (cfgOpen3DPIDPlots->get(0u)) {
            ese_parameters::hPIDQATar3D[0]->Fill(nSigmaTOFTar, nSigmaITSTar, track.pt());
          }
          if (cfgOpen3DPIDPlots->get(1u)) {
            ese_parameters::hPIDQATar3D[1]->Fill(nSigmaITSTar, nSigmaTPCTar, track.pt());
          }
          if (cfgOpen3DPIDPlots->get(2u)) {
            ese_parameters::hPIDQATar3D[2]->Fill(nSigmaTOFTar, nSigmaTPCTar, track.pt());
          }
        }
        if (cfgOpenv2) {
          if (track.sign() > 0) {
            ese_parameters::hv2Tar[0]->Fill(track.pt(), collision.centFT0C(), std::cos(2 * (track.phi() - psi2)));
          } else {
            ese_parameters::hv2Tar[1]->Fill(track.pt(), collision.centFT0C(), std::cos(2 * (track.phi() - psi2)));
          }
        }
        ese_parameters::eseCandidates.emplace_back(ESECandidate{
          collision.posZ(), collision.centFT0C(), psi2, q2Tarinit, q2Refinit, static_cast<int8_t>(track.sign()), correctedTpcInnerParam, track.tpcSignal(), track.pt(), track.eta(), track.phi(),
          track.dcaXY(), track.dcaZ(), static_cast<uint8_t>(track.tpcNClsFound()), track.itsNCls(), track.tpcChi2NCl(), track.itsChi2NCl(),
          nSigmaTPCTar, nSigmaTOFTar, nSigmaITSTar, track.itsClusterSizes()});
      }
      if (kIsRef) {
        multiRef++;
        q2Refx += std::cos(2 * track.phi());
        q2Refy += std::sin(2 * track.phi());
        bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
        float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
        float nSigmaTPCRef{(poiRef == ese_parameters::kProton && !cfgUseSelfnSigmaTPCProton) ? track.tpcNSigmaPr() : getNSigmaTPCSelfBB(track, poiRef)};
        float nSigmaTOFRef{getNSigmaTOF(track, poiRef)};
        float nSigmaITSRef{getNSigmaITS(track, poiRef)};
        if (cfgOpenPIDQA) {
          ese_parameters::hPIDQARef1D[0]->Fill(track.pt());
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
          ese_parameters::hPIDQARef2D[1]->Fill(nSigmaTPCRef, track.pt());
          ese_parameters::hPIDQARef2D[2]->Fill(nSigmaTOFRef, track.pt());
          ese_parameters::hPIDQARef2D[3]->Fill(nSigmaITSRef, track.pt());
          if (cfgOpen3DPIDPlots->get(0u)) {
            ese_parameters::hPIDQARef3D[0]->Fill(nSigmaTOFRef, nSigmaITSRef, track.pt());
          }
          if (cfgOpen3DPIDPlots->get(1u)) {
            ese_parameters::hPIDQARef3D[1]->Fill(nSigmaITSRef, nSigmaTPCRef, track.pt());
          }
          if (cfgOpen3DPIDPlots->get(2u)) {
            ese_parameters::hPIDQARef3D[2]->Fill(nSigmaTOFRef, nSigmaTPCRef, track.pt());
          }
        }
        if (cfgOpenv2) {
          if (track.sign() > 0) {
            ese_parameters::hv2Ref[0]->Fill(track.pt(), collision.centFT0C(), std::cos(2 * (track.phi() - psi2)));
          } else {
            ese_parameters::hv2Ref[1]->Fill(track.pt(), collision.centFT0C(), std::cos(2 * (track.phi() - psi2)));
          }
        }
      }
    }
  }

  void init(InitContext const&)
  {
    poiTar = getPOI(cfgTarName);
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
    // hists for event level QA
    histsESE.add("EventQA/histEventCount", ";Event Count;Counts", {HistType::kTH1F, {{100, 0, 100}}});
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
      histsESE.add("EventQA/hist_globalTracks_centT0C_before", "before cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, cfgaxisNch}});
      histsESE.add("EventQA/hist_PVTracks_centT0C_before", "before cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, cfgaxisNchPV}});
      histsESE.add("EventQA/hist_globalTracks_PVTracks_before", "before cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histsESE.add("EventQA/hist_globalTracks_multT0A_before", "before cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histsESE.add("EventQA/hist_globalTracks_multV0A_before", "before cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histsESE.add("EventQA/hist_multV0A_multT0A_before", "before cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histsESE.add("EventQA/hist_multT0C_centT0C_before", "before cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, cfgaxisT0C}});
      histsESE.add("EventQA/hist_globalTracks_centT0C_after", "after cut;Centrality T0C;mulplicity global tracks", {HistType::kTH2D, {axisCentForQA, cfgaxisNch}});
      histsESE.add("EventQA/hist_PVTracks_centT0C_after", "after cut;Centrality T0C;mulplicity PV tracks", {HistType::kTH2D, {axisCentForQA, cfgaxisNchPV}});
      histsESE.add("EventQA/hist_globalTracks_PVTracks_after", "after cut;mulplicity PV tracks;mulplicity global tracks", {HistType::kTH2D, {cfgaxisNchPV, cfgaxisNch}});
      histsESE.add("EventQA/hist_globalTracks_multT0A_after", "after cut;mulplicity T0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histsESE.add("EventQA/hist_globalTracks_multV0A_after", "after cut;mulplicity V0A;mulplicity global tracks", {HistType::kTH2D, {cfgaxisT0A, cfgaxisNch}});
      histsESE.add("EventQA/hist_multV0A_multT0A_after", "after cut;mulplicity T0A;mulplicity V0A", {HistType::kTH2D, {cfgaxisT0A, cfgaxisT0A}});
      histsESE.add("EventQA/hist_multT0C_centT0C_after", "after cut;Centrality T0C;mulplicity T0C", {HistType::kTH2D, {axisCentForQA, cfgaxisT0C}});
    }
    histsESE.add("PlanQA/histQvec_CorrL0_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL1_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL2_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histQvec_CorrL3_V2", ";#it{Q_{x}};#it{Q_{y}};Centrality", {HistType::kTH3F, {cfgaxisQvecF, cfgaxisQvecF, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL0_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL1_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL2_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histEvtPl_CorrL3_V2", ";EventPlane angle;Centrality", {HistType::kTH2F, {axisEvtPl, cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_SigRefAV2", ";Centrality;Cos(Sig-RefA)", {HistType::kTProfile, {cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_SigRefBV2", ";Centrality;Cos(Sig-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    histsESE.add("PlanQA/histQvecRes_RefARefBV2", ";Centrality;Cos(RefA-RefB)", {HistType::kTProfile, {cfgaxisCent}});
    // hists for track level QA
    histsESE.add("TrackQA/hist_dEdxTPC_All", ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", {HistType::kTH2F, {cfgrigidityBins, cfgdedxBins}});
    histsESE.add("TrackQA/hist_pt_All", ";#it{p}_{T};counts", {HistType::kTH1F, {cfgaxispt}});
    histsESE.add("TrackQA/hist_eta_All", ";#it{#eta};counts", {HistType::kTH1F, {cfgaxisetaPID}});
    histsESE.add("TrackQA/hist_phi_All", ";#it{#phi};counts", {HistType::kTH1F, {axisPhi}});
    histsESE.add("TrackQA/hist_ITSNcls_All", ";ITSNcls;counts", {HistType::kTH1F, {axisITSNcls}});
    histsESE.add("TrackQA/hist_TPCNcls_All", ";TPCNcls;counts", {HistType::kTH1F, {axisTPCNcls}});
    histsESE.add("TrackQA/hist_ITSChi2NDF_All", ";ITS#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
    histsESE.add("TrackQA/hist_TPCChi2NDF_All", ";TPC#it{#chi^{2}}/NDF;counts", {HistType::kTH1F, {cfgaxisChi2Ncls}});
    histsESE.add("TrackQA/hist_DCAxy_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAxy}});
    histsESE.add("TrackQA/hist_DCAz_All", ";#it{DCA_{xy}};counts", {HistType::kTH1F, {cfgaxisDCAz}});
    // v2 and ESEPlots
    /// QA plots
    if (cfgOpenPIDQA) {
      ese_parameters::hPIDQATar2D[0] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_dEdxTPC_%s", cfgTarName.value.c_str()), ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", HistType::kTH2F, {cfgrigidityBins, cfgdedxBins});
      ese_parameters::hPIDQATar1D[0] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_pt_%s", cfgTarName.value.c_str()), ";#it{p}_{T};counts", HistType::kTH1F, {cfgaxispt});
      ese_parameters::hPIDQATar1D[1] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_eta_%s", cfgTarName.value.c_str()), ";#it{#eta};counts", HistType::kTH1F, {cfgaxisetaPID});
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
      ese_parameters::hPIDQARef2D[0] = histsESE.add<TH2>(Form("ESE/TrackQA/hist_dEdxTPC_%s", cfgRefName.value.c_str()), ";#it{p}^{TPC}/#it{z} (GeV/c);d#it{E}/d#it{x}", HistType::kTH2F, {cfgrigidityBins, cfgdedxBins});
      ese_parameters::hPIDQARef1D[0] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_pt_%s", cfgRefName.value.c_str()), ";#it{p}_{T};counts", HistType::kTH1F, {cfgaxispt});
      ese_parameters::hPIDQARef1D[1] = histsESE.add<TH1>(Form("ESE/TrackQA/hist_eta_%s", cfgRefName.value.c_str()), ";#it{#eta};counts", HistType::kTH1F, {cfgaxisetaPID});
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
      if (cfgOpen3DPIDPlots->get(0u)) {
        ese_parameters::hPIDQATar3D[0] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFITSPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxispt});
        ese_parameters::hPIDQARef3D[0] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFITSPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}ITS;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsITS, cfgaxispt});
      }
      if (cfgOpen3DPIDPlots->get(1u)) {
        ese_parameters::hPIDQATar3D[1] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaITSTPCPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxispt});
        ese_parameters::hPIDQARef3D[1] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaITSTPCPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}ITS;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsITS, cfgnSigmaBinsTPC, cfgaxispt});
      }
      if (cfgOpen3DPIDPlots->get(2u)) {
        ese_parameters::hPIDQATar3D[2] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFTPCPt_%s", cfgTarName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxispt});
        ese_parameters::hPIDQARef3D[2] = histsESE.add<TH3>(Form("ESE/TrackQA/hist_nSigmaTOFTPCPt_%s", cfgRefName.value.c_str()), ";n_{#sigma}TOF;n_{#sigma}TPC;#it{p}_{T}", HistType::kTH3F, {cfgnSigmaBinsTOF, cfgnSigmaBinsTPC, cfgaxispt});
      }
    }
    // v2 plots
    if (cfgOpenv2) {
      ese_parameters::hv2Tar[0] = histsESE.add<TProfile2D>(Form("ESE/V2/hist_%sPosV2", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%)", HistType::kTProfile2D, {cfgaxispt, cfgaxisCent});
      ese_parameters::hv2Tar[1] = histsESE.add<TProfile2D>(Form("ESE/V2/hist_%sNegV2", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%)", HistType::kTProfile2D, {cfgaxispt, cfgaxisCent});
      ese_parameters::hv2Ref[0] = histsESE.add<TProfile2D>(Form("ESE/V2/hist_%sPosV2", cfgRefName.value.c_str()), ";#it{p}_{T};Centrality (%)", HistType::kTProfile2D, {cfgaxispt, cfgaxisCent});
      ese_parameters::hv2Ref[1] = histsESE.add<TProfile2D>(Form("ESE/V2/hist_%sNegV2", cfgRefName.value.c_str()), ";#it{p}_{T};Centrality (%)", HistType::kTProfile2D, {cfgaxispt, cfgaxisCent});
    }
    // ESE plots
    if (cfgOpenESEQA) {
      ese_parameters::hESEQATar1D[0] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum", cfgTarName.value.c_str()), ";Num_{Proton}/Event;counts", HistType::kTH1F, {{100, 0, 100}});
      ese_parameters::hESEQATar1D[1] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2", cfgTarName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2});
      ese_parameters::hESEQATar2D = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent", cfgTarName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2, cfgaxisCent});
      ese_parameters::hESEQARef1D[0] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sNum", cfgRefName.value.c_str()), ";Num_{He3}/Event;counts", HistType::kTH1F, {{10, 0, 10}});
      ese_parameters::hESEQARef1D[1] = histsESE.add<TH1>(Form("ESE/ESEQA/hist_%sq2", cfgRefName.value.c_str()), ";#it{q}_{2};counts", HistType::kTH1F, {cfgaxisq2});
      ese_parameters::hESEQARef2D = histsESE.add<TH2>(Form("ESE/ESEQA/hist_%sq2_Cent", cfgRefName.value.c_str()), ";#it{q}_{2};Centrality (%)", HistType::kTH2F, {cfgaxisq2, cfgaxisCent});
    }
    if (cfgOpenESE) {
      ese_parameters::hESETar = histsESE.add<THnSparse>(Form("ESE/ESE/histESE_%s", cfgTarName.value.c_str()), ";#it{p}_{T};Centrality (%);#it{q}_{2};cos(#phi-#Psi_{2});Charge", HistType::kTHnSparseF, {cfgaxispt, cfgaxisCent, cfgaxisq2, cfgaxiscos, axisCharge});
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::Qvectors>>::iterator const& collision, TracksPIDFull const& tracks)
  {
    ese_parameters::eseCandidates.clear();
    const float centrality{collision.centFT0C()};
    const int64_t multTrk{tracks.size()};
    if (cfgOpenFullEventQA)
      fillEventQAhistBe(collision, multTrk, centrality);
    if (!eventSelBasic(collision, multTrk, centrality, true))
      return;
    if (cfgOpenFullEventQA)
      fillEventQAhistAf(collision, multTrk, centrality);
    float q2Tarx{0.};
    float q2Tary{0.};
    float q2Refx{0.};
    float q2Refy{0.};
    int multiTar{0};
    int multiRef{0};
    fillESECandidates(collision, tracks, q2Tarx, q2Tary, multiTar, q2Refx, q2Refy, multiRef);
    float q2Tar{calculateq2(q2Tarx, q2Tary, multiTar)};
    float q2Ref{calculateq2(q2Refx, q2Refy, multiRef)};
    if (cfgOpenESEQA) {
      ese_parameters::hESEQATar1D[0]->Fill(multiTar);
      ese_parameters::hESEQATar1D[1]->Fill(q2Tar);
      ese_parameters::hESEQATar2D->Fill(q2Tar, centrality);
      ese_parameters::hESEQARef1D[0]->Fill(multiRef);
      ese_parameters::hESEQARef1D[1]->Fill(q2Ref);
      ese_parameters::hESEQARef2D->Fill(q2Ref, centrality);
    }
    if (multiTar == 0)
      return;
    for (const auto& c : ese_parameters::eseCandidates) {
      eseTable(c.vtz, c.centFT0C, c.psi2FT0C, q2Tar, q2Ref, c.signTar, c.tpcInnerParamTar, c.tpcSignalTar, c.ptTar, c.etaTar, c.phiTar, c.dcaXYTar, c.dcaZTar, c.tpcNclsTar, c.itsNclsTar, c.tpcChi2NDFTar, c.itsChi2NDFTar, c.tpcNSigmaTar, c.tofNSigmaTar, c.itsNSigmaTar, c.itsClusSizeTar);
      if (cfgOpenESE) {
        ese_parameters::hESETar->Fill(c.ptTar, c.centFT0C, q2Ref, std::cos(2 * (c.phiTar - c.psi2FT0C)), c.signTar);
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
